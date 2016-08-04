###################################################################
#######################INDEPENDENT_RECORDS V3.0####################
###Arguments#######################################################
######sp = vector of species#######################################
######indiv = vector of individuals ###############################
########(use a dummy vector when not applicable)###################
######station = vector of trapping stations########################
######age = vector of individuals' ages (default is NULL)##########
######time = POSIXct vector of date-time of records################
######var = variables to append to output##########################
###################################################################

ind_record<-function(sp,indiv,station,age=NULL, time, z, var=NULL){
	if(is(time)[1]!="POSIXct"){stop("Time must be inserted as POSIXct class object")}
	if (is.null(age)){age<-rep(1,length(time))} ##Gambiarrinha pra fingir que todas as linhas possuem idades iguais (logo, a idade nao importa)
	data<-data.frame(sp, indiv, station, age, time)		#creates main data.frame
		if(!is.null(var)){data<-cbind(data,var)}  # pool desirable variables together with photo data
	data<-data[order(data$sp, data$indiv, data$age, data$station, data$time, decreasing=F),] #order dataframe in relevant manner
	independent<-T #logical vector, 1st element TRUE
	for (i in 2:nrow(data)){ #for every row but the first:
		li<-max(which(independent==T))  #define index of the last independent record
		if(data$sp[i]==data$sp[li] & data$station[i]==data$station[li] & data$indiv[i]==data$indiv[li] & data$age[i]==data$age[li]){
		 if(abs(as.numeric(difftime(data$time[li],data$time[i],units="mins")))<z){independent[i]<-F}
			else{independent[i]<-T}
		}
		else{independent[i]<-T}
	}
	data<-cbind(data,independent)
	return(data)
}