######################
###RMark Processing###
######################
	###Processing Data in RMark formatting
	
	agout.proc<-RMark::process.data(data=eh.df,
		model="PoissonMR",                                                       
		counts=list(
			"Unmarked Seen"=eh.aux['unmark',],	
			"Marked Unidentified"=eh.aux['mark.n.id',] ,
			"Known Marks"=kmark),
		begin.time=1,
		time.intervals=interval_info$time.interval.months[-1]
		)
		
	naive.agout.proc<-RMark::process.data(data=eh.df,
		model="PoissonMR",                                                       
		counts=list(
			"Unmarked Seen"=eh.aux['unmark',],	
			"Marked Unidentified"=eh.aux['mark.n.id',],
			"Known Marks"=naive.kmark),
		begin.time=1,
		time.intervals=interval_info$time.interval.months[-1]
		)
	
	#####
	###Generating Design Data information
	#####
	agout.ddl<-RMark::make.design.data(agout.proc)
	
	#####
	###Adding RMark time information to interval summary
	#####
	interval_info$time<-force.numeric(agout.ddl$alpha$time[1:no_int])
	interval_info$survey.mean<-rep(
		tapply(interval_info$time,interval_info$sample_month,mean),
		each=5)
	interval_info$survey.beg<-rep(
		tapply(interval_info$time,interval_info$sample_month,min),
		each=5)
	interval_info$survey.end<-rep(
		tapply(interval_info$time,interval_info$sample_month,max),
		each=5)
	survey.trans<-
		(interval_info$survey.end[1:(length(interval_info$survey.end)-1)] + 
		interval_info$survey.beg[2:length(interval_info$survey.beg)]) / 2
	interval_info$survey.trans<-interval_info$survey.end
	interval_info$survey.trans[1:(length(survey.trans))]<-survey.trans
	
	################################################################
	###Adding individual covariates do Processed Data RMark Object##
	################################################################
	##Born in Captitivy
		agout.proc$data$captive%<>%as.numeric
	##Gender
		agout.proc$data$sex%<>%as.numeric
		agout.proc$data$sex[is.na(agout.proc$data$sex)]<-mean(agout.proc$data$sex,na.rm=TRUE)
	##Marked as Young		
		agout.proc$data$marked.as.young%<>%as.numeric
	##Age at each interval 't'
		##Creating a template matrix to be filled with age values (0 = young; 1 = adult)
		temp<-matrix(ncol=no_int,nrow=nrow(agout.proc$data),dimnames=list(
			rownames=rownames(agout.proc$data),colnames=paste("age.t",interval_info$time,sep="")))
		##Looping over each individual in the analysis
		for (i in 1:nrow(agout.proc$data)){
			##vector of primary intervals before capture (with NA values)
				a<-rep(NA,length(which(interval_info$start<agout.proc$data$resig.cohort[i])))
			##vector of time in days since release for each given individual in a given interval
				tsr<-difftime(interval_info$start,agout.proc$data$date.mark.rel[i],units="days")
				##negative values are replaced as 'NA'
					tsr[tsr<0]<-NA
			##vector of primary intervals over which an individual was young;
			##given that it was captured and mark while young; he was considered young for 60 days
			##after release
				b<-rep((agout.proc$data$marked.as.young[i]==F)*1,length(which(tsr<=60)))
			##vector primary intervals over which an individual was an adult;
			##it is defined as all the intervals that were not neither in 'a' or 'b' vector
				c<-rep(1,no_int-length(c(a,b)))
			##replace each row in the template column with the respective age information
				temp[i,]<-c(a,b,c)
			##replace age values after removal of an individual as NA 
				temp[i,which(interval_info$start>=agout.proc$data$date.rm[i])]<-NA
		}
		##Vector of mean age of marked individuals over each interval
		age.means<-apply(temp,2,mean,na.rm=TRUE)
		##Replacing missing values (NA's) for the mean age of of marked individuals
		for(c in 1:ncol(temp)){temp[is.na(temp[,c]),c]<-age.means[c]}
		##Binding age data to Processed Data RMark object
			agout.proc$data<-cbind(agout.proc$data,temp)

	####################################################################	
	###Adding primary interval covariates to Design Data RMark Object###
	####################################################################	
		agout.ddl$Phi<-RMark::merge_design.covariates(
			ddl=agout.ddl$Phi,df=interval_info,bytime=T)
		agout.ddl$alpha<-RMark::merge_design.covariates(
			ddl=agout.ddl$alpha,df=interval_info,bytime=T)
		agout.ddl$sigma<-RMark::merge_design.covariates(
			ddl=agout.ddl$sigma,df=interval_info,bytime=T)
		agout.ddl$U<-RMark::merge_design.covariates(
			ddl=agout.ddl$U,df=interval_info,bytime=T)
		agout.ddl$GammaDoublePrime<-RMark::merge_design.covariates(
			ddl=agout.ddl$GammaDoublePrime, df=interval_info,bytime=T)
		agout.ddl$GammaPrime<-RMark::merge_design.covariates(
			ddl=agout.ddl$GammaPrime,df=interval_info,bytime=T)
	#######################
	###Fixing Parameters###
	#######################
		agout.ddl.full<-agout.ddl
		##Creating columns to determine fixed values
		agout.ddl$Phi$fix<-NA
		agout.ddl$GammaPrime$fix<-NA
		agout.ddl$GammaDoublePrime$fix<-NA
		
	##Fixing Parameters of non-existing release cohorts (those after the first primary interval of each survey)
		agout.ddl$Phi$fix[
			!agout.ddl$Phi$cohort%in%interval_info$survey.beg]<-1
		agout.ddl$GammaPrime$fix[
			!agout.ddl$GammaPrime$cohort%in%interval_info$survey.beg]<-0	
		agout.ddl$GammaDoublePrime$fix[
			!agout.ddl$GammaDoublePrime$cohort%in%interval_info$survey.beg]<-0
			
	##Saving Design Data object of 'Open Structure'
	op.agout.ddl<-agout.ddl

	##Fixing parameters of closed populations (transitions between primary intervals from the same survey)
		agout.ddl$Phi$fix[
			!agout.ddl$Phi$time%in%interval_info$survey.end]<-1
		agout.ddl$GammaPrime$fix[
			!agout.ddl$GammaPrime$time%in%interval_info$survey.end]<-0
		agout.ddl$GammaDoublePrime$fix[
			!agout.ddl$GammaDoublePrime$time%in%interval_info$survey.end]<-0
			
	##Saving Design Data object of 'Closed Structure'
	cl.agout.ddl<-agout.ddl
	
	##Declaring which covariates are individual, temporal, categoric or continuous
		##Temporal, continous covariates
		time.cont.cov<-c("effort.trapdays","recordloss.agout.total","effort.var","recordloss.agout.var","age.t")
		##Temporal, categoric covariates
		time.catg.cov<-c("time","sample_month")
		##Individual covariates
		indiv.cov<-c("captive","age.t","sex")
	
	####
	###Data.frame containing the conservative estimate of number of marked individuals:
	###The minium number (of marked) known alive
	####	
	survey.mnka<-data.frame(mnka=naive.kmark,interval=names(naive.kmark))%>%
		merge(interval_info[,c('start','sample_month','time','interval')],by='interval')%>%
		group_by(sample_month)%>%
		mutate(surv.mnka=max(mnka))		