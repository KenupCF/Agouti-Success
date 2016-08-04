########################################################
###################Interval To Point####################
###################(interval2point)#####################
#takes a data.frame where each record has a defined,####
#temporal begining and end, and converts each record in#
#a sequence of records 1-min apart###################### 
###Version 1.1##########################################
###by Caio Kenup########################################
###contact: caio.kenup@gmail.com########################
#Arguments:#############################################
###df: dataframe to be extended#########################
###start: POSIXct vector containing interval beginings##
###end:POSIXct vector containing interval ends##########
########################################################
require(lubridate)
interval2point<-function(df,start,end){
	run.start<-now()
	pb <- winProgressBar(title="Running interval2point", label="0% done", min=0, max=100, initial=0)
	if(is(start)[1]!="POSIXct" | is(end)[1]!="POSIXct"){
		stop("Time must be inserted as POSIXct class object")} 
	if(tz(start) != tz(end)){
		stop("start and end arguments must be in the same timezone!")}
	resu<-df[0,]
	n<-1
	x<-1
	r<-0
	while(x <= nrow(df)){
		while(r<=as.numeric(difftime(end[x], start[x], units="mins"))){
			resu[n,]<-df[x,]
			resu$time[n]<-start[x]+ 60*r
			n<-n+1
			r<-r+1
			}
		rownames(resu)<-1:nrow(resu)
		r<-0
	progress<-x/nrow(df)
	info <- sprintf("%.2f%% done", progress*100)
	setWinProgressBar(pb, progress*100, label=info)
		x<-x+1
	}
	resu$time<-as.POSIXct(as.numeric(resu$time), origin="1970-01-01", tz=tz(start))
	resu$horario_inicio<-NULL
	resu$horario_fim<-NULL
	close(pb)
	run.end<-now()
	return(resu)
	message("Function ran with success! Total running time:",abs(difftime(run.start,run.end)))
}