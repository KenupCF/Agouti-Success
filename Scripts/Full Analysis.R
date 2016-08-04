##Author: Caio Kenup#############
##email: caio.kenup@gmail.com###

##Packages and External Functions##
##Statistical packages
	require(RMark)
	require(msm)
	require(combinat)
##'Data-management' packages
	require(Matrix)
	require(lubridate)
	require(tidyr)
	require(gtools)
	require(gdata)
	require(stringr)
	require(zoo)
	require(devtools)
	require(dplyr)
	require(magrittr)
###Graphical packages
	require(plotrix)
	require(grDevices)
###Spatial tools packages
	require(maptools)
	require(rgdal)
	require(rgeos)
###Importing Custom Functions	
source_url('https://raw.githubusercontent.com/KenupCF/Agouti-Success/master/Independent_records%20v3.0.R')
source_url('https://raw.githubusercontent.com/KenupCF/Agouti-Success/master/point2interval%20v1.1.R')

####ver quais funções do quick functions tão sendo usadas
source("C:\\Users\\Kenup\\Dropbox\\03-Trabalho\\03-Programa?\\01-R\\Kenup-Repo\\Quick Functions.R")
##General parameters##
# setwd(C:\\Users\\Kenup\\Dropbox\\03-Trabalho\\01-Science\\01-Artigos sendo Escritos\\Sucesso das Cutias)
# discard.mark.n.id<-FALSE #se devemos descartar os registros de animais 'marcados, mas não identificados'
# cut.off<-as.POSIXct("01/05/2015" ,format="%d/%m/%Y", tz="Etc/GMT-3") ##cut-off date to interrupt analysis
		
##Importing and formatting data##
##Photographic records
fotos<-read.csv(".\\Dados Usados\\Armadilhagem Fotogr?ca\\NAO ALTERAR - Registros Armadilhagem Fotogr?ca.csv",T)
##Capture information
captures<-read.table(".\\Dados Usados\\Captura Viva\\Capturas Triagem.txt",header=TRUE,sep="\t")
##Individual information
dossie<-read.table(".\\Dados Usados\\Armadilhagem Fotográfica\\Dossiê Cutias Resumido.txt",T,row.names=1)
dossie$indiv<-rownames(dossie)
# dossie$marked.as.young[is.na(dossie$marked.as.young)==T]<-FALSE	####gambiarrando, trocando NAs por FALSE
##Resighting interval information
interval_cov<-read.table(".\\Dados Usados\\Armadilhagem Fotogr?ca\\Informa?s de Intervalo.txt",T)
no_int<-nrow(interval_cov)
interval_cov$interval<-1:no_int
##Trapping sessions information
trapping_sessions<-read.table(".\\Dados Usados\\Captura Viva\\Sess??de Captura.txt",T)
##Camera Trapping Effort
trap_hist_1314<-read.csv(".\\Dados Usados\\Armadilhagem Fotogr?ca\\Esfor?2013-2014.csv",T)	##planilha de esforço amostral 2013
trap_hist_1516<-read.csv(".\\Dados Usados\\Armadilhagem Fotogr?ca\\Esfor?2015-2016.csv",T)	##planilha de esforço amostral 2015
trap_hist<-rbind(trap_hist_1314,trap_hist_1516)									####juntando esforços amostrais
##Spatial data###
spatial.data.wd<-'.\\Dados Usados\\Shapes'
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
fotos$horario_fim<-as.POSIXct(paste(fotos$data, substr(fotos$horario_fim,1,5)),format="%d/%m/%Y %H:%M", tz="Etc/GMT-3")fotos$horario_inicio<-as.POSIXct(paste(fotos$data, substr(fotos$horario_inicio,1,5)) ,format="%d/%m/%Y %H:%M", tz="Etc/GMT-3")
fotos$data<-as.POSIXct(fotos$data ,format="%d/%m/%Y", tz="Etc/GMT-3")
dossie$date.rm<-as.POSIXct(dossie$date.rm,format="%d/%m/%Y", tz="Etc/GMT-3")
dossie$date.mark.rel<-as.POSIXct(dossie$date.mark.rel,format="%d/%m/%Y", tz="Etc/GMT-3")
captures$date.capture<-as.POSIXct(captures$date.capture,format="%d/%m/%Y", tz="Etc/GMT-3")
captures$date.handling<-as.POSIXct(captures$date.handling,format="%d/%m/%Y", tz="Etc/GMT-3")
captures$date.release<-as.POSIXct(captures$date.release,format="%d/%m/%Y", tz="Etc/GMT-3")
captures$month<-myear(captures$date.capture)
trapping_sessions$start<-as.POSIXct(trapping_sessions$start,format="%d/%m/%Y", tz="Etc/GMT-3")    ##Convertendo datas
trapping_sessions$end<-as.POSIXct(trapping_sessions$end,format="%d/%m/%Y", tz="Etc/GMT-3")        ##para POSIXct
trapping_sessions$month<-myear(trapping_sessions[,c('start')])%>%as.character
interval_cov$start<-as.POSIXct(interval_cov$start,format="%d/%m/%Y", tz="Etc/GMT-3")    ##Convertendo datas
interval_cov$end<-as.POSIXct(interval_cov$end,format="%d/%m/%Y", tz="Etc/GMT-3")        ##para POSIXct
interval_cov<-merge(interval_cov,data.frame(										##adicionando "mes ano" para cada intervalo
		month=yearmon(tapply(myear(interval_cov$start),interval_cov$sample_month,function(x){
			x<-x[1]
			return(x)})),
		sample_month=levels(interval_cov$sample_month)),
	by="sample_month")
trap_hist$data<-as.POSIXct(trap_hist$data,format="%d/%m/%Y", tz="Etc/GMT-3") 
{
# REMOVENDO DADOS E METADADOS POSTERIORES AO CUT-OFF
		# fotos%<>%dplyr::filter(horario_fim < cut.off)
		# dossie%<>%dplyr::filter(date.mark.rel < cut.off)
		# interval_cov<-interval_cov[interval_cov$start < cut.off,]
		# no_int<-nrow(interval_cov)  ##número de intervalos de amostragem (atualizando depois de fazer o cut-off)
		# trap_hist%<>%dplyr::filter(data < cut.off)
		# captures%<>%dplyr::filter(date.capture < cut.off)
		# trapping_sessions%<>%dplyr::filter(start < cut.off)
}	

##Formating Capture Information##
# captures$weight.kg%<>%force.numeric
captures$Peso.kg%<>%force.numeric
captures$age<-NA
captures$age[captures$Peso.kg < 2]<-"Young" 
captures$age[captures$Peso.kg >= 2]<-"Adult"
# captures$age[captures$weight.kg < 2]<-"Young"
# captures$age[captures$weight.kg >= 2]<-"Adult"

#####
##Formatting photographic records
###Some individuals' marks were not permanent or unambigous. 
###hese individuals were considered as 'Unmarked' for the purpose of the analysis
#####
unusable.ind<-c("sloth","chico","primo","neguinha","negrinha") #vector of unusable individuals ()
fotos$true.indiv<-fotos$individuo
fotos$individuo[fotos$individuo%in%unusable.ind]<-"unmark"
{
		# if(discard.mark.n.id==TRUE) {
				# fotos$individuo[fotos$individuo=="mark.n.id"]<-"na"}
}		

##Photographic data were formatted with a 'start' and 'end' time for each record;
##So as to exclude temporally dependent records, we transformed each record into a series of records
##spaced one minute apart, from its starting to ending time.
fotos<-interval2point(fotos,fotos$horario_inicio,fotos$horario_fim)		

#Defining primary interval of each record
fotos$interval<-sapply(fotos$time, function(x){
	if(any(y <- (x>=interval_cov$start & x<=interval_cov$end+86399))>0) {interval_cov$interval[which(y)]} 
	else {NA}})
#Appending interval information to each record
fotos<-merge(fotos,interval_cov[,c("interval","sample_month","start")],by="interval",
	all.x=TRUE)
{
# age.struct<-data.frame(idade=c("filhote","juvenil","adulto","desconhecida"),age=c(0,1,1,"desconhecida"))
		# fotos<-merge(fotos,age.struct,by="idade",all.x=TRUE)
			#gambiarra, trocando NA's por 1 na idade
			# fotos[is.na(fotos$age),"age"]<-1
}

##Assigning primary interval to trapping effort data
trap_hist$interval<-sapply(trap_hist$data, function(x){
	if(any(y <- (x>=interval_cov$start & x<=interval_cov$end+86399))>0) {interval_cov$interval[which(y)]} 
	else {NA}})

##Removing unused stations
unused.st<-c("AFV25","AFV21","AFV10","AFV05") #Character vector containing discarded stations
trap_hist<-trap_hist[trap_hist$estacao%in%unused.st==F,]#Removing such stations from effort records
fotos<-fotos[fotos$estacao%in%unused.st==F,]			#Removing photographic records from the same stations
sp.obj$camera.traps<-sp.obj$camera.traps[				#Removing stations' coordinates from spatial object
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
				

{##CALCULANDO NÍVEL DE CONCORDÂNCIA	(SÓ IMPLEMENTÁVEL QUANDO TODAS AS FOTOS FOREM IDENTIFICADAS POR 2 PESSOAS INDEPENDENTEMENTE)
		# fotos<-fotos[,str_detect(colnames(fotos),"id2")==F] #removendo colunas do segundo identificador (temporário, enquanto elas não estão completas)
		#indiv.agree<-rownames(fotos[fotos$individuo==fotos$individuo.id2,])
		#agree<-fotos$individuo==fotos$individuo.id2
		#indiv.disag<-rownames(fotos[fotos$individuo!=fotos$individuo.id2,])
		#ind.agree.lvl<-length(indiv.disag)/length(c(indiv.disag,indiv.agree))
		#age.agree<-rownames(fotos[fotos$age.id2==fotos$age.id2,])
		#age.disag<-rownames(fotos[fotos$age.id2!=fotos$age.id2,])
		#age.agree.lvl<-length(age.disag)/length(c(age.disag,age.agree))
		#
		#TROCANDO IDENTIDADES PARA A CATEGORIA MENOS INFORMATIVA
			#fotos$individuo[agree==FALSE & (fotos$individuo=="na" | fotos$individuo.id2=="na")]<-"na"
			#fotos$individuo[agree==FALSE & (fotos$individuo=="unmark" | fotos$individuo.id2=="unmark")]<-"na"
			#fotos$individuo[agree==FALSE & (fotos$individuo%in%rownames(dossie) & fotos$individuo.id2%in%rownames(dossie))]<-"mark.n.id"
			#fotos$individuo[agree==FALSE & (fotos$individuo%in%rownames(dossie) & fotos$individuo.id2=="mark.n.id")]<-"mark.n.id"
			#fotos$individuo[agree==FALSE & (fotos$individuo=="mark.n.id" & fotos$individuo.id2%in%rownames(dossie))]<-"mark.n.id"
}

###############
####PART II####
###############

##Temporal Independence##
	criteria.vec<-c(30,45,90,120,60*6,60*12,60*24,60)
		fotos.ind.ls<-list()
		mark.results<-list()
		encounter.histories<-list()
		for (zz in 
			1:length(criteria.vec)){
		z<-criteria.vec[zz]	##definindo critério de independência temporal, em minutos

##Determining temporally independent records
	fotos.ind.ls[[zz]]<-ind_record(
		sp=fotos$especie,
		z=z,
		st=fotos$estacao,
		age=fotos$age,
		ind=fotos$individuo,
		time=fotos$time,
		#data.frame object with appendable variables
		var=data.frame( 
			interval=fotos$interval,
			sample_month=fotos$sample_month,
			true.ind=fotos$true.indiv))
	##Ordering new data.frame by time
	fotos.ind.ls[[zz]]<-fotos.ind.ls[[zz]][order(fotos.ind.ls[[zz]]$time),] ##ordenando planilha nova por data	
	fotos.ind.ls[[zz]][fotos.ind.ls[[zz]]$age=="desconhecida","age"]<-NA
	fotos.ind.ls[[zz]]$age%<>%force.numeric
	##Save data.frame object outside of list
		fotos.ind<-fotos.ind.ls[[zz]]
	##General formatting
		##Include among 'station' levels those without any records
		levels(fotos.ind$station)<-union(levels(fotos.ind$station),levels(trap_hist$estacao))
		##Convert sample_month to factor
		interval_cov$sample_month%<>%factor
		##Save object with resighting surveys
		months<-levels(interval_cov$sample_month)			
	##Keeping only 'Dasyprocta leporina' records##
		agout.ind<-fotos.ind%>%
			dplyr::filter(sp=="Dasyprocta.leporina")
	
	##Dealing with 'Unidentifiable' and 'Marked, but not Identified' records;
	##These can be temporally dependent of other, identified records
	##(only to marked individuals' records, in the case of 'Marked, but not Identified')
	
	##Ordering data.frame by: Species, Station, Age and Time
		agout.ind<-agout.ind[order(agout.ind$sp, agout.ind$station, agout.ind$age, agout.ind$time, decreasing=F),]
	##Loop over all independent, 'Unidentifiable' or 'Marked, but not Identified' records
		for (i in rownames(agout.ind[agout.ind$independent==T & (agout.ind$indiv=="na" | agout.ind$indiv=="mark.n.id"),])){
			if(is.na(agout.ind[i,"age"])){
			
			min.dist.time<-difftime(
				agout.ind[i,"time"],
				agout.ind[
					rownames(agout.ind)!=i & 
					agout.ind$sp==agout.ind[i,"sp"] & 
					agout.ind$station==agout.ind[i,"station"],"time"],
				units="mins")%>%
			abs%>%min%>%as('numeric')
			
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
			
			} else{
			
			min.dist.time<-difftime(
				agout.ind[i,"time"],
				agout.ind[
					rownames(agout.ind)!=i & 
					agout.ind$sp==agout.ind[i,"sp"] & 
					agout.ind$station==agout.ind[i,"station"]  & 
					(agout.ind$age==agout.ind[i,"age"] | is.na(agout.ind$age) ),"time"],
				units="mins")%>%
			abs%>%min%>%as('numeric')
			
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
			
			if (agout.ind[i,"indiv"]=="na" & (min.dist.time<=z | min.dist.time==Inf)){
				agout.ind[i,"independent"]<-F}
			if (agout.ind[i,"indiv"]=="mark.n.id" & (min.dist.time.mark<=z | min.dist.time.mark==Inf)){
				agout.ind[i,"independent"]<-F}
			}

{			
##Creating a separate object for dropped records	
	#dropped.agout.rec<-agout.ind[!agout.ind$independent,] 
}

##Keeping just the indepent records
agout.ind<-agout.ind[agout.ind$independent,]

####
##Individuals with the same individual identifcation
##('Unmarked', 'Unidentifiable' or 'Marked, but Unidentified')
##present on the same record received prefixes ('b_', 'c_' and so on)
##so as to prevent being wrongfully discarded by the independence criteria.
##Here we remove these preffixes
####
agout.ind$indiv[agout.ind$indiv%in%c("b_unmark","c_unmark","d_unmark")]<-"unmark"
agout.ind$indiv[agout.ind$indiv%in%c("b_na","c_na","d_na")]<-"na"

##Creating trap_summary object (a list of 'j' time intervals x 's' stations matrices)
trap_summary<-list()
temp<-data.frame()
temp.eff<-data.frame()
summary.template<-matrix(data=0, ncol=length (levels(agout.ind$station)) ,nrow=no_int,
	dimnames=list(rownames=1:no_int,colnames=levels(trap_hist$estacao)))
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
		trap_summary$effort.trapdays[j,]<-tapply(temp.eff$esforco,temp.eff$estacao,sum)
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
	interval_cov$duration.days<-(difftime(interval_cov$end,interval_cov$start,"days")%>%as('numeric') ) + 1
	interval_cov$time_elaps.days[2:no_int]<-difftime(interval_cov$start[2:no_int],interval_cov$end[1:no_int-1],"days") 		
	interval_cov$time.interval.days[2:no_int]<-difftime(interval_cov$start[2:no_int],interval_cov$start[1:no_int-1],"days")
	interval_cov$time.interval.months<-round(interval_cov$time.interval.days/30,digits=3)
	interval_cov$time.interval.years<-round(interval_cov$time.interval.days/365,digits=4)
	interval_cov$cum.time_elaps.days[2:no_int]<-cumsum(interval_cov$time_elaps.days[-1])
	interval_cov$effort.trapdays<-apply(trap_summary$effort.trapdays,1,sum)
	interval_cov$failure.rate<-interval_cov$effort.trapdays/(6*33)
	interval_cov$effort.trapdays.log<-log(interval_cov$effort.trapdays)
	interval_cov$effort.var<-apply(trap_summary$effort.trapdays,1,var)
	interval_cov$recordloss.agout.mean<-apply(trap_summary$recordloss,1,mean)
	interval_cov$recordloss.agout.total<-apply(trap_summary$unused.records,1,sum)/apply(trap_summary$total.records,1,sum)
	interval_cov$recordloss.agout.var<-apply(trap_summary$recordloss,1,var)
	interval_cov$agout.records<-apply(trap_summary$total.records,1,sum)
	rownames(interval_cov)<-interval_cov$start
	
##Formatting individual summary 
###Defining individuals born in the captivity
	dossie$captive<-FALSE
	dossie$captive[dossie$origin!="pnt"]<-TRUE
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
	post.poneds<-dossie$captive==T & 
		!is.na(dossie$date.mark.rel) &
		dossie$date.mark.rel >= min(interval_cov$start) & 
		dossie$date.mark.rel < max(interval_cov$end)
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
	{
		##GAMBIARRAS
			dossie["x","resig.cohort"]<-dossie["negra","resig.cohort"]
	}
		
################################################################################	
####Generating Encounter History formatted "Log-Poisson Mark Resight" models####
################################################################################

	##Character vector of all individuals recorded in the study (even not Identified)
	ind<-union(levels(factor(agout.ind$indiv)),levels(factor(dossie$indiv)))

	##Summarising number of encounter for each individual in each interval
	eh.temp<-agout.ind%>%
		dplyr::filter(!is.na(interval))%>%
		dplyr::group_by(indiv,interval)%>%
		dplyr::summarise(no=n())
	
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
		merge(interval_cov[,c('interval','start')],by='interval')%>%
		#merging in individual information
		merge(dossie[,
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
	eh.df<-data.frame(indiv=names(ch),ch=ch)%>%merge(dossie,by='indiv',all.y=FALSE)
	rownames(eh.df)<-eh.df$indiv

	{
	####Gambiarras#####
				eh.df["x","captive"]<-TRUE
				eh.df["x","resig.cohort"]<-eh.df["negra","resig.cohort"]
	}
	
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
####################
###RMark Analysis###
####################
	###Processing Data in RMark formatting
	agout.proc<-RMark::process.data(data=eh.df,model="PoissonMR",                                                       
		counts=list(
			"Unmarked Seen"=eh.aux['unmark',],	
			"Marked Unidentified"=eh.aux['mark.n.id',] 
			"Known Marks"=kmark),
		begin.time=1,
		time.intervals=interval_cov$time.interval.months[-1]
		)
	naive.agout.proc<-RMark::process.data(data=eh.df,model="PoissonMR",                                                       
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

####
##Running simple models, all-constant, to have initial values to start more complex ones##
####
	naive.dot.fix.siman<-RMark::mark(naive.agout.proc,op.agout.ddl,options="SIMANNEAL",threads=-1)
	dot.fix.init.siman<-RMark::mark(agout.proc,op.agout.ddl,
		initial=naive.dot.fix.siman,
		options="SIMANNEAL",threads=-1)
		
##Assessing population closure 
	#Running a model with open parameters structure
	open.model<-RMark::mark(agout.proc,op.agout.ddl,initial=dot.fix.init.siman,threads=-1,
			model.parameters=list(
				U=list(formula=~time),
				alpha=list(formula=~1+time)))
				
	#Running a model with closed parameters structure
	clos.model<-RMark::mark(agout.proc,cl.agout.ddl,initial=dot.fix.init.siman,threads=-1,
		model.parameters=list(
			U=list(formula=~sample_month),
			alpha=list(formula=~1+time)))
			
	#Collecting both models for model selection
	open.close.models<-collect.models(lx=c("clos.model","open.model"),type="PoissonMR")
		
		{
		if(clos.model$results$AICc - open.model$results$AICc > 2){closed.months<-FALSE} else {closed.months<-TRUE}
		if (closed.months){true.agout.ddl<-cl.agout.ddl} else{true.agout.ddl<-op.agout.ddl}
				mark.results[[zz]]<-open.close.models
		}		
				
	} #fim do loop de analise de sensibilidade de independendência temporal

###AQUI ACABA A ANALISE DE SENSIBILIDADE E O CLOSURE TEST	
