##Author: Caio Kenup#############
##email: caio.kenup@gmail.com####

################################################################################	
####Generating Encounter History formatted "Log-Poisson Mark Resight" models####
################################################################################

	##Character vector of all individuals recorded in the study (even not Identified)
	ind<-union(levels(factor(agout.ind$indiv)),levels(factor(indiv.summ$indiv)))

	##Summarising number of encounter for each individual in each interval
	eh.temp<-agout.ind%>%
		dplyr::filter(!is.na(interval))%>%
		dplyr::group_by(indiv,interval)%>%
		dplyr::summarise(no=n())
	class(eh.temp)<-'data.frame'
	
	##Filling unexisting combinations of individuals and intervals 
	existing.combinations<-apply(eh.temp[,1:2],1,paste,collapse=" ")
	indiv.interval<-expand.grid(indiv=ind,interval=1:no_int) #all combinations of individuals and intervals
	indiv.interval$indiv%<>%as('character')
	unexisting.combinations<-apply(indiv.interval,1,function(x){
		!paste0(x,collapse=" ") %in% existing.combinations})
	eh.temp<-rbind(eh.temp,
		indiv.interval%>%
			dplyr::filter(unexisting.combinations)%>%
			dplyr::mutate(no=0))
	
	##Aditional data.frame formatting (re-factoring intervals)
		eh.temp%<>%
		#refactoring intervals
		dplyr::mutate(interval=factor(interval))%>%	
		#variable to split marked individuals from not marked ou not indentified groups 	
		dplyr::mutate(is.indiv=!indiv%in%c('unmark','mark.n.id','na'))%>%
		#merging in date information of intervals
		merge(interval_info[,c('interval','start')],by='interval')%>%
		#merging in individual information
		merge(indiv.summ[,
				c('indiv','resig.cohort','true.cohort','date.mark.rel','date.rm')],
			by='indiv',all=TRUE)%>%
		#removing records out of intervals (NA)
		dplyr::filter(!is.na(interval))%>%
		#removing discarded individuals (these are present only as 0's, since their records are considered 'unmarked')
		dplyr::filter(!indiv%in%unusable.ind)%>%
		#sorting data.frame
		dplyr::arrange(indiv,interval)
		
#####
#In RDPNE, absence of encounters are registered using three different notations:
#####
	###
	#Intervals in which a given individual was not yet marked, or after it was removed,
	#are written as '..'
	###
	eh.temp[eh.temp$is.indiv &
		( (eh.temp$start<eh.temp$resig.cohort)%>%na2false |
			(eh.temp$start>eh.temp$date.rm)%>%na2false),'no']<-'..'	
	###
	#Intervals in which there was no resighting for an individuals, yet he was known to be in the population
	#(e.g. immediately after marking and release), are written as '+0'
	###
	eh.temp[eh.temp$is.indiv &
		eh.temp$no==0 & 
		eh.temp$start==eh.temp$true.cohort & 
		difftime(eh.temp$true.cohort,eh.temp$date.mark.rel,units='days')<15 &
		eh.temp$indiv%in%captures$indiv,'no']<-'+0'	
	###
	#Intervals in which there was no resighting of an individuals hence it is not known if he is in the population
	#(having exited by either mortality, permanent imigration or temporary imigration)
	#are written as '-0'
	###
	eh.temp[eh.temp$is.indiv &
		eh.temp$no==0,'no']<-'-0'
###
#Here we divided the encounter history, separating marked, individually identified animals, of
#grouping of not individually identifiable classes
###
	eh.temp.split<-split(eh.temp,eh.temp$is.indiv)

##Additional formatting individually identifiable records
	eh.temp.split[['TRUE']]%<>%
		dplyr::group_by(indiv)%>%
		dplyr::filter(sum(no%in%c('-0','..'))!=no_int)%>%
		ungroup%>%
		dplyr::mutate(indiv=factor(indiv),
			value=zero_pad(no,2))
			
##Generating matrix containing records of NOT individually identifiable classes
	eh.aux<-apply(acast(data=eh.temp.split[['FALSE']]%<>%
				mutate(indiv=factor(indiv),value=no),
			formula=indiv~interval),
		2,as.numeric)
	rownames(eh.aux)<-levels(eh.temp.split[['FALSE']]$indiv)	

##Generating matrix containing records of individually identifiable animals
	eh.mat<-acast(data=eh.temp.split[['TRUE']], formula=indiv~interval)
	ch<-apply(eh.mat,1,paste0,collapse="")

##Formatting encounter history as data.frame, with all encounter collapsed as character strings,
##as require by package RMark
	eh.df<-data.frame(indiv=names(ch),ch=ch)%>%merge(indiv.summ,by='indiv',all.y=FALSE)
	rownames(eh.df)<-eh.df$indiv
	eh.df$ch%<>%as.character
#####	
###Vectors cointaning the number of exact number of marks in the population
###('0' means exact number of marks in not know)
#####

	##Naive vector of known marks (considering the minimum number of marks present as the exact number of marks)
	naive.kmark<-apply(eh.mat,2,function(x){
		sum(x%in%c('+0',zero_pad(1:99,2)))
		})		
	##True vector of known marks
	kmark<-apply(eh.mat,2,function(x){
		sum(x%in%c('+0',zero_pad(1:99,2)))*(!any(str_detect(x,"-0")))
		})