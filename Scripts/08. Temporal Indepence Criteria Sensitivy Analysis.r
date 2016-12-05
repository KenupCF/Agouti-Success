##Author: Caio Kenup#############
##email: caio.kenup@gmail.com####

#########################
##Temporal Independence##
#########################

###Vector of temporal independence criteria, in minutes
criteria.vec<-c(30,45,60,90,120,60*6,60*12,60*24)

##Creating empty lists
fotos.ind.ls<-list()
mark.results<-list()

##Looping over each independence criteria
for (zz in 1:length(criteria.vec)){
#######################################################
##Temporal Independence Criteria, further Formatting###
#######################################################

	z<-criteria.vec[zz]	#Selecting independence criterion, in minutes
	##Determining temporally independent records
	photos.ind<-ind_record(
		#Vector of species' names
		sp=photos$species,
		#Value of independence criteria
		z=z,
		#Vector of camera trap stations
		st=photos$station,
		#Vector of individuals' ages
		age=photos$age,
		#Vector of individuals' names
		# ind=photos$individuo,
		ind=photos$true.indiv,
		#Vector of time of records
		time=photos$time,
		#data.frame object with appendable variables to each records
		var=data.frame( 
			#Primary interval of records
			interval=photos$interval,
			#Resighting survey of records
			sample_month=photos$sample_month,
			#'True' name of individual
			#(different from 'individuo' since discarded individuals were assigned as 'unmark')
			true.ind=photos$true.indiv))
			
	##Ordering new data.frame by time
	photos.ind%<>%arrange(time)
	photos.ind[photos.ind$age=="unknown","age"]<-NA
	photos.ind$age%<>%force.numeric
	
		
###	
##General formatting
###
		##Include among 'station' levels those without any records
		levels(photos.ind$station)<-union(levels(photos.ind$station),levels(trap_hist$station))
		##Convert sample_month to factor
		interval_info$sample_month%<>%factor
		##Save object with resighting surveys
		months<-levels(interval_info$sample_month)			
		##Keeping only 'Dasyprocta leporina' records##
		agout.ind<-photos.ind%>%
			dplyr::filter(sp=="Dasyprocta.leporina")
			
	####
	##Dealing with 'Unidentifiable' and 'Marked, but not Identified' records;
	##These can be temporally dependent of other, identified records
	##(only to marked individuals' records, in the case of 'Marked, but not Identified')
	####
	
	##Ordering data.frame by: Species, Station, Age and Time
		agout.ind<-agout.ind[order(agout.ind$sp, agout.ind$station, agout.ind$age, agout.ind$time, decreasing=F),]
		
	##Defining independent, 'Unidentifiable' or 'Marked, but not Identified' records' ID
		dubious.records<-rownames(agout.ind[
			agout.ind$independent==T & 
			(agout.ind$indiv=="na" | agout.ind$indiv=="mark.n.id")
			,])

	pb <- winProgressBar(title="Revoking independence from dubious records",
		label="0% done", min=0, max=100, initial=0)
	##Loop over all 'Dubious Records' 
		for (i in dubious.records){
			###If an records' age is UNKNOWN
			if(is.na(agout.ind[i,"age"])){
			
			##Compute the minimum temporal distance to another independent record
			##from the same species at the same station
			min.dist.time<-difftime(
				agout.ind[i,"time"],
				agout.ind[
					rownames(agout.ind)!=i & 
					agout.ind$sp==agout.ind[i,"sp"] & 
					agout.ind$station==agout.ind[i,"station"],"time"],
				units="mins")%>%
			abs%>%min%>%as('numeric')
			
			##Compute the minimum temporal distance to another independent record from a marked individual 
			##of the same species at the same station
			min.dist.time.mark<-difftime(
				agout.ind[i,"time"],
				agout.ind[
					rownames(agout.ind)!=i & 
					agout.ind$sp==agout.ind[i,"sp"] & 
					agout.ind$station==agout.ind[i,"station"]  & 
					str_detect(agout.ind$indiv,"unmark")!=T & 
					agout.ind$indiv!="na","time"],	
				units="min")%>%
			abs%>%min%>%as('numeric')
			
			###If an records' age is KNOWN
			} else{
			
			##Compute the minimum temporal distance to another independent record
			##from the same species at the same station, with the same (or of unknown) age
			min.dist.time<-difftime(
				agout.ind[i,"time"],
				agout.ind[
					rownames(agout.ind)!=i & 
					agout.ind$sp==agout.ind[i,"sp"] & 
					agout.ind$station==agout.ind[i,"station"]  & 
					(agout.ind$age==agout.ind[i,"age"] | is.na(agout.ind$age) ),"time"],
				units="mins")%>%
			abs%>%min%>%as('numeric')
			
			##Compute the minimum temporal distance to another independent record from a marked individual 
			##of the same species at the same station, with the same  (or of unknown) age
			min.dist.time.mark<-difftime(
				agout.ind[i,"time"],
				agout.ind[
					rownames(agout.ind)!=i & 
					agout.ind$sp==agout.ind[i,"sp"] & 
					agout.ind$station==agout.ind[i,"station"]  & 
					(agout.ind$age==agout.ind[i,"age"] | is.na(agout.ind$age)==T) &
					(!str_detect(agout.ind$indiv,"unmark")) & 
					agout.ind$indiv!="na","time"],
				units="min")%>%
			abs%>%min%>%as('numeric')
			
			}
			
			###Revoke independence from 'Unidentifiable' individuals closer than 'z' minutes
			###from any other independent individual
			if (agout.ind[i,"indiv"]=="na" & (min.dist.time<=z | min.dist.time==Inf)){
				agout.ind[i,"independent"]<-F}
				
			###Revoke independence from 'Marked, but Unidentified' individuals closer than 'z' minutes
			###from any other independent marked individual				
			if (agout.ind[i,"indiv"]=="mark.n.id" & (min.dist.time.mark<=z | min.dist.time.mark==Inf)){
				agout.ind[i,"independent"]<-F}
			
			progress<-which(dubious.records==i)/length(dubious.records)
			info <- sprintf("%.2f%% done", progress*100)
			setWinProgressBar(pb, progress*100, label=info)
			}
	close(pb)
	##Keeping just the indepent records
	agout.ind%<>%
		dplyr::filter(independent)

####
##Individuals with the same individual identifcation
##('Unmarked', 'Unidentifiable' or 'Marked, but Unidentified')
##present on the same record received prefixes ('b_', 'c_' and so on)
##so as to prevent being wrongfully discarded by the independence criteria.
##Here we remove these preffixes
####

agout.ind$indiv[agout.ind$indiv%in%c("b_unmark","c_unmark","d_unmark")]<-"unmark"
agout.ind$indiv[agout.ind$indiv%in%c("b_mark.n.id","c_mark.n.id","d_mark.n.id")]<-"mark.n.id"
agout.ind$indiv[agout.ind$indiv%in%c("b_na","c_na","d_na")]<-"na"

##Creating trap_summary object (a list of matrices with 'j' time intervals x 's' stations)
trap_summary<-list()
temp<-data.frame()
temp.eff<-data.frame()
summary.template<-matrix(data=0, ncol=length (levels(agout.ind$station)) ,nrow=no_int,
	dimnames=list(rownames=1:no_int,colnames=levels(trap_hist$station)))
trap_summary$effort.trapdays<-summary.template
trap_summary$total.records<-summary.template
trap_summary$unused.records<-summary.template
trap_summary$recordloss<-summary.template

##Looping over each primary interval
	for (j in 1:no_int){
		#subset of records in the primary interval j
		temp<-agout.ind[agout.ind$interval==j,] 
		#subset of trap effort history in the primary interval j
		temp.eff<-trap_hist[trap_hist$interval==j,]
		##Defining total effort per station in the interval 'j'
		trap_summary$effort.trapdays[j,]<-tapply(temp.eff$effort,temp.eff$station,sum)
		##Defining number of indepedent records per station in the interval 'j'
		trap_summary$total.records[j,]<-tapply(temp$independent,temp$station,sum)
		trap_summary$total.records[which(is.na(trap_summary$total.records))]<-0
		##Defining number of discarded, indepedent records per station in the interval 'j'
		trap_summary$unused.records[j,]<-tapply(temp$independent[temp$indiv=="na"],temp$station[temp$indiv=="na"],sum)
		trap_summary$unused.records[which(is.na(trap_summary$unused.records))]<-0
		##Defining rate of record loss (proportion of discarded records) per station in the interval 'j'
		trap_summary$recordloss[j,]<-trap_summary$unused.records[j,]/trap_summary$total.records[j,]
		trap_summary$recordloss[which(is.na(trap_summary$recordloss))]<-0
	}
	
###Calculating interval and trapping sessions covariates###
	trapping_sessions$duration.days<-(difftime(trapping_sessions$end,trapping_sessions$start,"days")%>%as('numeric') )+ 1
	interval_info$duration.days<-(difftime(interval_info$end,interval_info$start,"days")%>%as('numeric') ) + 1
	interval_info$time_elaps.days[2:no_int]<-difftime(interval_info$start[2:no_int],interval_info$end[1:no_int-1],"days") 		
	interval_info$time.interval.days[2:no_int]<-difftime(interval_info$start[2:no_int],interval_info$start[1:no_int-1],"days")
	interval_info$time.interval.months<-round(interval_info$time.interval.days/30,digits=3)
	interval_info$time.interval.years<-round(interval_info$time.interval.days/365,digits=4)
	interval_info$cum.time_elaps.days[2:no_int]<-cumsum(interval_info$time_elaps.days[-1])
	interval_info$effort.trapdays<-apply(trap_summary$effort.trapdays,1,sum)
	interval_info$failure.rate<-interval_info$effort.trapdays/(6*33)
	interval_info$effort.trapdays.log<-log(interval_info$effort.trapdays)
	interval_info$effort.var<-apply(trap_summary$effort.trapdays,1,var)
	interval_info$recordloss.agout.mean<-apply(trap_summary$recordloss,1,mean)
	interval_info$recordloss.agout.total<-apply(trap_summary$unused.records,1,sum)/apply(trap_summary$total.records,1,sum)
	interval_info$recordloss.agout.var<-apply(trap_summary$recordloss,1,var)
	interval_info$agout.records<-apply(trap_summary$total.records,1,sum)
	rownames(interval_info)<-interval_info$start
	
	
	###Generate RDPNE encounterhistory
		source(".\\Scripts\\03.0.Generating Encounter History formatted Log-Poisson Mark Resight models.R")
	
	###RMark Processing Data
		source(".\\Scripts\\04.0.Processing RMark Data.R")
	
####
##Running simple models, all-constant, to have initial values to start more complex ones##
####
	naive.dot.fix.siman<-RMark::mark(naive.agout.proc,cl.agout.ddl,
		options="SIMANNEAL",
		threads=-1)
	dot.fix.init.siman<-RMark::mark(agout.proc,cl.agout.ddl,
		initial=naive.dot.fix.siman,
		options="SIMANNEAL",
		threads=-1)
		
#Running model with temporal variation on 'alpha'
	mark.model<-RMark::mark(agout.proc,cl.agout.ddl,
		initial=dot.fix.init.siman,
		threads=-1,
		model.parameters=list(
			U=list(formula=~sample_month),
			alpha=list(formula=~1+time)))
			
##Saving results from zzth independence criterion
	mark.results[[zz]]<-mark.model			
}	
