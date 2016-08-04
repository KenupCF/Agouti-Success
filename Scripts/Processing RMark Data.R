######################
###RMark Processing###
######################
	###Processing Data in RMark formatting
	
	agout.proc<-RMark::process.data(data=eh.df,
		model="PoissonMR",                                                       
		counts=list(
			"Unmarked Seen"=eh.aux['unmark',],	
			"Marked Unidentified"=eh.aux['mark.n.id',] 
			"Known Marks"=kmark),
		begin.time=1,
		time.intervals=interval_cov$time.interval.months[-1]
		)
		
	naive.agout.proc<-RMark::process.data(data=eh.df,
		model="PoissonMR",                                                       
		counts=list(
			"Unmarked Seen"=eh.aux['unmark',],	
			"Marked Unidentified"=eh.aux['mark.n.id',] 
			"Known Marks"=naive.kmark),
		begin.time=1,
		time.intervals=interval_cov$time.interval.months[-1]
		)
			
	###Generating Design Data information
	agout.ddl<-RMark::make.design.data(agout.proc)
	
	###Adding RMark time information to interval summary
	interval_cov$time<-force.numeric(agout.ddl$alpha$time[1:no_int])
	interval_cov$survey.mean<-rep(
		tapply(interval_cov$time,interval_cov$sample_month,mean),
		each=5)
	interval_cov$survey.beg<-rep(
		tapply(interval_cov$time,interval_cov$sample_month,min),
		each=5)
	interval_cov$survey.end<-rep(
		tapply(interval_cov$time,interval_cov$sample_month,max),
		each=5)
	survey.trans<-(survey.end[1:(length(survey.end)-1)] + survey.beg[2:length(survey.beg)]) / 2
	interval_cov$survey.trans<-interval_cov$survey.end
	interval_cov$survey.trans[1:(length(survey.trans)*5)]<-rep(survey.trans,each=5)	
	
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
			rownames=rownames(agout.proc$data),colnames=paste("age.t",interval_cov$time,sep="")))
		##Looping over each individual in the analysis
		for (i in 1:nrow(agout.proc$data)){
			##vector of primary intervals before capture (with NA values)
				a<-rep(NA,length(which(interval_cov$start<agout.proc$data$resig.cohort[i])))
			##vector of time in days since release for each given individual in a given interval
				tsr<-difftime(interval_cov$start,agout.proc$data$date.mark.rel[i],units="days")
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
				temp[i,which(interval_cov$start>=agout.proc$data$date.rm[i])]<-NA
		}
		##Vector of mean age of marked individuals over each interval
		age.means<-apply(temp,2,mean,na.rm=TRUE)
		##Replacing missing values (NA's) for the mean age of of marked individuals
		for(c in 1:ncol(temp)){temp[is.na(temp[,c]),c]<-age.means[c]}
		##Binding age data to Processed Data RMark object
			agout.proc$data<-cbind(agout.proc$data,temp)
			{# dummy.agout.proc$data<-agout.proc$data
			}
	####################################################################	
	###Adding primary interval covariates to Design Data RMark Object###
	####################################################################	
		agout.ddl$Phi<-RMark::merge_design.covariates(
			ddl=agout.ddl$Phi,df=interval_cov,bytime=T)
		agout.ddl$alpha<-RMark::merge_design.covariates(
			ddl=agout.ddl$alpha,df=interval_cov,bytime=T)
		agout.ddl$sigma<-RMark::merge_design.covariates(
			ddl=agout.ddl$sigma,df=interval_cov,bytime=T)
		agout.ddl$U<-RMark::merge_design.covariates(
			ddl=agout.ddl$U,df=interval_cov,bytime=T)
		agout.ddl$GammaDoublePrime<-RMark::merge_design.covariates(
			ddl=agout.ddl$GammaDoublePrime, df=interval_cov,bytime=T)
		agout.ddl$GammaPrime<-RMark::merge_design.covariates(
			ddl=agout.ddl$GammaPrime,df=interval_cov,bytime=T)
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
			!agout.ddl$Phi$cohort%in%interval_cov$survey.beg]<-1
		agout.ddl$GammaPrime$fix[
			!agout.ddl$GammaPrime$cohort%in%interval_cov$survey.beg]<-0	
		agout.ddl$GammaDoublePrime$fix[
			!agout.ddl$GammaDoublePrime$cohort%in%interval_cov$survey.beg]<-0
			
	##Saving Design Data object of 'Open Structure'
	op.agout.ddl<-agout.ddl

	##Fixing parameters of closed populations (transitions between primary intervals from the same survey)
		agout.ddl$Phi$fix[
			!agout.ddl$Phi$time%in%interval_cov$survey.end]<-1
		agout.ddl$GammaPrime$fix[
			!agout.ddl$GammaPrime$time%in%interval_cov$survey.end]<-0
		agout.ddl$GammaDoublePrime$fix[
			!agout.ddl$GammaDoublePrime$time%in%interval_cov$survey.end]<-0
			
	##Saving Design Data object of 'Closed Structure'
	cl.agout.ddl<-agout.ddl