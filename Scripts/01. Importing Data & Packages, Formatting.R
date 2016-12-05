##Author: Caio Kenup#############
##email: caio.kenup@gmail.com####

##Packages and External Functions##
##Statistical packages
	require(RMark)
	require(msm)
	require(combinat)
	require(gtools)
##'Data-management' packages
	require(lubridate)
	require(gdata)
	require(stringr)
	require(zoo)
	require(dplyr)
	require(magrittr)
	require(devtools)
	require(reshape2)
###Spatial tools packages
	require(maptools)
	require(rgdal)
	require(rgeos)
###Importing Custom Functions	
	##Transform interval records to point records
	devtools::source_url('https://raw.githubusercontent.com/KenupCF/Agouti-Success/master/Functions/interval2point%20v1.1.R')
	##Define temporally independent records
	devtools::source_url('https://raw.githubusercontent.com/KenupCF/Agouti-Success/master/Functions/Independent_records%20v3.0.R')
	##Quick functions to make life easier
	devtools::source_url('https://raw.githubusercontent.com/KenupCF/Agouti-Success/master/Functions/Aux%20Functions.R')

Sys.setlocale("LC_TIME", "English")

	
##Importing and formatting data##
	##Photographic records
		photos<-read.csv(".\\Data\\Photographic Records Subsetted.csv",T)
	##Capture information
		captures<-read.csv(".\\Data\\Live Capture Records.csv",header=TRUE)
	##Individual information
		indiv.summ<-read.csv(".\\Data\\Individual Summary.csv",T,row.names=1)
		indiv.summ$indiv<-rownames(indiv.summ)
	##Resighting interval information
		interval_info<-read.csv(".\\Data\\Resighting Intervals Information.csv",T)
		no_int<-nrow(interval_info)
	##Trapping sessions information
		trapping_sessions<-read.csv(".\\Data\\Trapping Sessions.csv",T)
	##Camera Trapping Effort
		trap_hist_1314<-read.csv(".\\Data\\Camera Trap Effort 2013-2014.csv",T)
		trap_hist_15<-read.csv(".\\Data\\Camera Trap Effort 2015.csv",T)
		trap_hist<-rbind(trap_hist_1314,trap_hist_15)
	##Spatial data###
		spatial.data.wd<-'.\\Data\\Spatial Files'
		sp.obj<-list()
		sp.obj$camera.traps<-readOGR(dsn=spatial.data.wd,layer="Camera Trap Locations",
			encoding='UTF-8',use_iconv=TRUE)				
		sp.obj$study.area<-readOGR(dsn=spatial.data.wd,layer="TNP Limits",
			encoding='UTF-8',use_iconv=TRUE)
		sp.obj$release.pens<-readOGR(dsn=spatial.data.wd,layer="Release Pen Locations",
			encoding='UTF-8',use_iconv=TRUE)		
		sp.obj$coimbra.sites<-readOGR(dsn=spatial.data.wd,
			layer="Coimbra-Filho Release Sites (Estimated)",encoding='UTF-8',use_iconv=TRUE)
		
##Formatting Time Data as POSIXct##
tz<-'Etc/GMT-3' #Study Area timezone code
photos$time_end<-as.POSIXct(paste(photos$date, substr(photos$time_end,1,5)),format="%d/%m/%Y %H:%M", tz=tz)
photos$time_start<-as.POSIXct(paste(photos$date, substr(photos$time_start,1,5)) ,format="%d/%m/%Y %H:%M", tz=tz)
photos$date%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
indiv.summ$date.rm%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
indiv.summ$date.mark.rel%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$date.capture%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$date.handling%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$date.release%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$month<-myear(captures$date.capture)
trapping_sessions$start%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
trapping_sessions$end%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
trapping_sessions$month<-myear(trapping_sessions[,c('start')])%>%as.character
interval_info$start%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
interval_info$end%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)        
interval_info<-merge(interval_info,
	data.frame(									
		month=yearmon(tapply(myear(interval_info$start),interval_info$sample_month,function(x){
			x<-x[1]
			return(x)})),
	sample_month=unique(interval_info$sample_month)),
	by="sample_month")
trap_hist$date<-as.POSIXct(trap_hist$date,format="%d/%m/%Y", tz=tz) 

##Formating Capture Information##
captures$age<-NA
captures$weight.kg%<>%force.numeric
captures$age[captures$weight.kg < 2]<-"Young"
captures$age[captures$weight.kg >= 2]<-"Adult"

#####
##Formatting photographic records
###Some individuals' marks were not permanentstror unambigous. 
###hese individuals were considered as 'Unmarked' for the purpose of the analysis
#####
unusable.ind<-c("sloth","chico","primo","neguinha","negrinha") #vector of unusable individuals
photos$true.indiv<-photos$indiv
photos$indiv[photos$indiv%in%unusable.ind]<-"unmark"

#####
##Photographic data were formatted with 'start' and 'end' times for each record;
##To facilitade data management, we transformed each record into a series of records
##spaced one minute apart, from its starting to ending time.
#####
photos<-interval2point(photos,photos$time_start,photos$time_end)		

#Defining primary interval of each record
	photos$interval<-sapply(photos$time, function(x){
		if(any(y <- (x>=interval_info$start & x<=interval_info$end+86399))>0) {interval_info$interval[which(y)]} 
		else {NA}})
#Appending interval information to each record
	photos<-merge(photos,interval_info[,c("interval","sample_month","start")],
		by="interval",all.x=TRUE)
		

#Assigning primary interval to trapping effort data
trap_hist$interval<-sapply(trap_hist$date, function(x){
	if(any(y <- (x>=interval_info$start & x<=interval_info$end+86399))>0) {interval_info$interval[which(y)]} 
	else {NA}})

##Removing unused stations
	#Character vector containing discarded stations
	unused.st<-c('AFV25','AFV05','AFV10','AFV21')
	#Removing such stations from effort records
	trap_hist<-trap_hist[!trap_hist$station%in%unused.st,]
	#Removing photographic records from the same stations
	photos<-photos[!photos$station%in%unused.st,]
	#Removing stations' coordinates from spatial object
	sp.obj$camera.traps<-sp.obj$camera.traps[				
		!sp.obj$camera.traps$label%in%unused.st,]
				
##Dealing with redundant stations
	###Merging "AFV05a" and "AFV05b"
		sp.obj$camera.traps@coords<-rbind(sp.obj$camera.traps@coords,
			apply(sp.obj$camera.traps@coords[str_detect(sp.obj$camera.traps$label,"AFV05"),],2,mean))
		sp.obj$camera.traps@data<-rbind(sp.obj$camera.traps@data,
			sp.obj$camera.traps@data[sp.obj$camera.traps@data$label=="AFV05a",])
		sp.obj$camera.traps@data[nrow(sp.obj$camera.traps@data),"label"]<-"AFV05"
		sp.obj$dropped.st<-sp.obj$camera.traps[sp.obj$camera.traps@data$label%in%c("AFV05a","AFV05b"),]
		sp.obj$camera.traps<-sp.obj$camera.traps[sp.obj$camera.traps$label%in%c("AFV05a","AFV05b")==FALSE,]
		sp.obj$camera.traps$label<-factor(sp.obj$camera.traps$label)
		sp.obj$camera.traps<-sp.obj$camera.traps[order(sp.obj$camera.traps$label,decreasing=FALSE),]
		
##Removing subscript information from stations
	trap_hist$station<-factor(substr(trap_hist$station,1,5))
	photos$station<-factor(substr(photos$station,1,5))
	sp.obj$camera.traps@data$label<-factor(substr(sp.obj$camera.traps@data$label,1,5))
				
##Spatial Analysis
	##Defining stations which were deployed earlier in the study	
		sp.obj$camera.traps$grid<-FALSE
		sp.obj$camera.traps$grid[1:21]<-TRUE
		sp.obj$fst.grid$distmatrix<-as.matrix(
			dist(sp.obj$camera.traps@coords[sp.obj$camera.traps$grid==TRUE,],method="euclidean"))
		sp.obj$fst.grid$distmatrix[sp.obj$fst.grid$distmatrix==0]<-NA
		sp.obj$fst.grid$mean_spacing<-mean(apply(sp.obj$fst.grid$distmatrix,1,min,na.rm=TRUE))
		sp.obj$fst.grid$sd_spacing<-sd(apply(sp.obj$fst.grid$distmatrix,1,min,na.rm=TRUE))
		sp.obj$scd.grid$distmatrix<-as.matrix(dist(sp.obj$camera.traps@coords,method="euclidean"))
		sp.obj$scd.grid$distmatrix[sp.obj$scd.grid$distmatrix==0]<-NA
		sp.obj$scd.grid$mean_spacing<-mean(apply(sp.obj$scd.grid$distmatrix,1,min,na.rm=TRUE))
		sp.obj$scd.grid$sd_spacing<-sd(apply(sp.obj$scd.grid$distmatrix,1,min,na.rm=TRUE))
	##Convex Hull from Traps
		sp.obj$fst.grid$MPC<-gConvexHull(sp.obj$camera.traps[sp.obj$camera.traps$grid==TRUE,])
		sp.obj$fst.grid$MPCarea<-gArea(sp.obj$fst.grid$MPC)
		sp.obj$scd.grid$MPC<-gConvexHull(sp.obj$camera.traps)		
		sp.obj$scd.grid$MPCarea<-gArea(sp.obj$scd.grid$MPC)
		
####
##Formatting individual summary 
####

###Defining individuals born in the captivity
	indiv.summ%<>%
		dplyr::mutate(captive=origin!='pnt')
###Defining release cohort (concerning mark resighting analyses)
	indiv.summ$true.cohort<-NA	
	##Loop over each individuals
		for (n in 1:nrow(indiv.summ)){
			##Define 'true.cohort' as the starting date of the first primary interval after marking
			indiv.summ$true.cohort[n]<-interval_info$start[
				min(which(interval_info$start>indiv.summ$date.mark.rel[n]))]
		} 
	##Converting to POSIXct		
	indiv.summ$true.cohort%<>%as.POSIXct(origin="1970-01-01", tz=tz(interval_info$start[1])) 
	
###Delaying inclusion of individuals recently released in the analysis
#Defining time (in days) to delay inclusion in the analysis
	y<-45 
	##Defining which individuals are to be post-poned (logic vector)
	post.poneds<-indiv.summ$captive==T &  							#captives individuals WHICH
		!is.na(indiv.summ$date.mark.rel) &							#have a date of mark and release AND
		indiv.summ$date.mark.rel >= min(interval_info$start) & 		#this date falls between the start and
		indiv.summ$date.mark.rel < max(interval_info$end)			#end of sampling
		
	##Adjusting date of release to recalculate release cohort
		indiv.summ$date.mark.rel[post.poneds==T]<-indiv.summ$date.mark.rel[post.poneds==T] + (y*86400)
		
	##Creating new cohort column
		indiv.summ$resig.cohort<-NA	
	##Loop over each individual
	for (n in 1:nrow(indiv.summ)){ 
		indiv.summ$resig.cohort[n]<-interval_info$start[min(which(interval_info$start>=indiv.summ$date.mark.rel[n]))]
		}
	##Converting to POSIXct		
	indiv.summ$resig.cohort%<>%as.POSIXct(origin="1970-01-01", tz=tz(interval_info$start[1])) 
	##Recovering original release dates
	indiv.summ$date.mark.rel[post.poneds==T]<-indiv.summ$date.mark.rel[post.poneds==T] - (y*86400)
	##Saving post poned information as new column
	indiv.summ$post.poneds<-post.poneds
	indiv.summ%<>%
		dplyr::filter(!is.na(indiv))