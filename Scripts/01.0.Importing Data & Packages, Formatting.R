##Author: Caio Kenup#############
##email: caio.kenup@gmail.com####

##Packages and External Functions##
##Statistical packages
	require(RMark)
	require(msm)
	require(combinat)
##'Data-management' packages
	require(lubridate)
	# require(gtools)
	require(gdata)
	require(stringr)
	require(zoo)
	require(dplyr)
	require(magrittr)
	require(devtools)
###Spatial tools packages
	require(maptools)
	require(rgdal)
	require(rgeos)
###Importing Custom Functions	
	##Transform interval records to point records
	devtools::source_url('https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/point2interval%20v1.1.R')
	##Define temporally independent records
	devtools::source_url('https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Independent_records%20v3.0.R')
	##Quick functions to make life easier (:
	devtools::source_url('https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Aux%20Functions.R')
	
{#@#remove this
##General parameters##
setwd("C:\\Users\\Kenup\\Dropbox\\03-Trabalho\\01-Science\\01-Artigos sendo Escritos\\Sucesso das Cutias\\GitHub Repo")
# setwd("K:\\03-Trabalho\\01-Science\\01-Artigos sendo Escritos\\Sucesso das Cutias")
}		
##Importing and formatting data##
	##Photographic records
		fotos<-read.csv(".\\Data\\Registros Armadilhagem Fotográfica.csv",T)
	##Capture information
		captures<-read.csv(".\\Data\\Capturas Triagem.csv",header=TRUE)
	##Individual information
		dossie<-read.csv(".\\Data\\Dossiê Cutias Resumido.csv",T,row.names=1)
		dossie$indiv<-rownames(dossie)
	##Resighting interval information
		interval_cov<-read.csv(".\\Data\\Informações de Intervalo.csv",T)
		no_int<-nrow(interval_cov)
	##Trapping sessions information
		trapping_sessions<-read.csv(".\\Data\\Sessões de Captura.csv",T)
	##Camera Trapping Effort
		trap_hist_1314<-read.csv(".\\Data\\Esforço 2013-2014.csv",T)
		trap_hist_15<-read.csv(".\\Data\\Esforço 2015.csv",T)
		trap_hist<-rbind(trap_hist_1314,trap_hist_15)
	##Spatial data###
		spatial.data.wd<-'.\\Data\\Spatial Files'
		sp.obj<-list()
		sp.obj$spatial.sampling<-readOGR(dsn=spatial.data.wd,layer="Amostragem v2.0")				
		sp.obj$study.area<-readOGR(dsn=spatial.data.wd,layer="limite_PARNA_TIJUCA")
		sp.obj$landmarks<-readOGR(dsn=spatial.data.wd,layer="Landmarks")
		sp.obj$camera.traps<-sp.obj$spatial.sampling[substr(sp.obj$spatial.sampling$label,1,3)=="AFV",]
		sp.obj$release.pens<-sp.obj$landmarks[str_detect(sp.obj$landmarks$label,"Cercado"),]
		sp.obj$coimbra.sites<-readOGR(dsn=spatial.data.wd,layer="Coimbra Sites Estimated")
			proj.utm<-"+proj=utm +datum=WGS84 +zone=23 +south"
			sp.obj$coimbra.sites%<>%spTransform(CRS(proj.utm))
		
##Formatting Time Data as POSIXct##
tz<-'Etc/GMT-3' #Study Area timezone code
fotos$horario_fim<-as.POSIXct(paste(fotos$data, substr(fotos$horario_fim,1,5)),format="%d/%m/%Y %H:%M", tz=tz)
fotos$horario_inicio<-as.POSIXct(paste(fotos$data, substr(fotos$horario_inicio,1,5)) ,format="%d/%m/%Y %H:%M", tz=tz)
fotos$data%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
dossie$date.rm%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
dossie$date.mark.rel%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$date.capture%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$date.handling%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$date.release%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
captures$month<-myear(captures$date.capture)
trapping_sessions$start%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
trapping_sessions$end%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
trapping_sessions$month<-myear(trapping_sessions[,c('start')])%>%as.character
interval_cov$start%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)
interval_cov$end%<>%as.POSIXct(format="%d/%m/%Y", tz=tz)        
interval_cov<-merge(interval_cov,
	data.frame(									
		month=yearmon(tapply(myear(interval_cov$start),interval_cov$sample_month,function(x){
			x<-x[1]
			return(x)})),
	sample_month=unique(interval_cov$sample_month)),
	by="sample_month")
trap_hist$data<-as.POSIXct(trap_hist$data,format="%d/%m/%Y", tz=tz) 

##Formating Capture Information##
captures$age<-NA
captures$weight.kg%<>%force.numeric
captures$age[captures$weight.kg < 2]<-"Young"
captures$age[captures$weight.kg >= 2]<-"Adult"

#####
##Formatting photographic records
###Some individuals' marks were not permanent or unambigous. 
###hese individuals were considered as 'Unmarked' for the purpose of the analysis
#####
unusable.ind<-c("sloth","chico","primo","neguinha","negrinha") #vector of unusable individuals
fotos$true.indiv<-fotos$individuo
fotos$individuo[fotos$individuo%in%unusable.ind]<-"unmark"

#####
##Photographic data were formatted with 'start' and 'end' times for each record;
##To facilitade data management, we transformed each record into a series of records
##spaced one minute apart, from its starting to ending time.
#####
fotos<-interval2point(fotos,fotos$horario_inicio,fotos$horario_fim)		

#Defining primary interval of each record
	fotos$interval<-sapply(fotos$time, function(x){
		if(any(y <- (x>=interval_cov$start & x<=interval_cov$end+86399))>0) {interval_cov$interval[which(y)]} 
		else {NA}})
#Appending interval information to each record
	fotos<-merge(fotos,interval_cov[,c("interval","sample_month","start")],
		by="interval",all.x=TRUE)
		

#Assigning primary interval to trapping effort data
trap_hist$interval<-sapply(trap_hist$data, function(x){
	if(any(y <- (x>=interval_cov$start & x<=interval_cov$end+86399))>0) {interval_cov$interval[which(y)]} 
	else {NA}})

{#@# rethink this
# age.struct<-data.frame(idade=c("filhote","juvenil","adulto","desconhecida"),age=c(0,1,1,"desconhecida"))
		# fotos<-merge(fotos,age.struct,by="idade",all.x=TRUE)
			#gambiarra, trocando NA's por 1 na idade
			# fotos[is.na(fotos$age),"age"]<-1
}

##Removing unused stations
	#Character vector containing discarded stations
	unused.st<-c("AFV25")
	#Removing such stations from effort records
	trap_hist<-trap_hist[!trap_hist$estacao%in%unused.st,]
	#Removing photographic records from the same stations
	fotos<-fotos[!fotos$estacao%in%unused.st,]
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
	trap_hist$estacao<-factor(substr(trap_hist$estacao,1,5))
	fotos$estacao<-factor(substr(fotos$estacao,1,5))
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
	dossie%<>%
		dplyr::mutate(captive=origin!='pnt')
###Defining release cohort (concerning mark resighting analyses)
	dossie$true.cohort<-NA	
	##Loop over each individuals
		for (n in 1:nrow(dossie)){
			##Define 'true.cohort' as the starting date of the first primary interval after marking
			dossie$true.cohort[n]<-interval_cov$start[
				min(which(interval_cov$start>dossie$date.mark.rel[n]))]
		} 
	##Converting to POSIXct		
	dossie$true.cohort%<>%as.POSIXct(origin="1970-01-01", tz=tz(interval_cov$start[1])) 
	
###Delaying inclusion of individuals recently released in the analysis
#Defining time (in days) to delay inclusion in the analysis
	y<-45 
	##Defining which individuals are to be post-poned (logic vector)
	post.poneds<-dossie$captive==T &  							#captives individuals WHICH
		!is.na(dossie$date.mark.rel) &							#have a date of mark and release AND
		dossie$date.mark.rel >= min(interval_cov$start) & 		#this date falls between the start and
		dossie$date.mark.rel < max(interval_cov$end)			#end of sampling
		
	##Adjusting date of release to recalculate release cohort
		dossie$date.mark.rel[post.poneds==T]<-dossie$date.mark.rel[post.poneds==T] + (y*86400)
		
	##Creating new cohort column
		dossie$resig.cohort<-NA	
	##Loop over each individual
	for (n in 1:nrow(dossie)){ 
		dossie$resig.cohort[n]<-interval_cov$start[min(which(interval_cov$start>=dossie$date.mark.rel[n]))]
		}
	##Converting to POSIXct		
	dossie$resig.cohort%<>%as.POSIXct(origin="1970-01-01", tz=tz(interval_cov$start[1])) 
	##Recovering original release dates
	dossie$date.mark.rel[post.poneds==T]<-dossie$date.mark.rel[post.poneds==T] - (y*86400)
	##Saving post poned information as new column
	dossie$post.poneds<-post.poneds
	dossie%<>%
		dplyr::filter(!is.na(indiv))