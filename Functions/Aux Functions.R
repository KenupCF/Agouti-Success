##Author: Caio Kenup#############
##email: caio.kenup@gmail.com####

###General
##Simpler conversion from factor to numeric
force.numeric<-function(x){as.numeric(as.character(x))}

###Strings
##Pad digits with a given numbers of zeroes to force a string of length 'w'
zero_pad<-function(x,w){
	require(stringr)
	y<-str_pad(x, w, side = "left", pad = "0")
	return(y)}
	
###Math
##Calculate the inverse of a number
inv<-function(x){1/x} 
##Calculate the complement of a number	
compl<-function(x) { 	
	if(x<0 | x>1){stop("x must be between 0 and 1")
	}else{
	1-x}
} 

###Time
##Return month and year from a POSIXct object
	myear<-function(x,abbr=TRUE){
		if(is.POSIXct(x)==FALSE){stop("'x' must be a POSIXct object")}
		require(lubridate)
		y<-paste(month(x,label=TRUE,abbr=TRUE),year(x),sep=" ")
		y<-as.yearmon(y, format="%b %Y",LANGUAGE=en)
		return(y)}

###RMark
##Determine impossible combinations of covariates in a model
imposib.mod<-function(x){  #x is a vector of covariate names
	any(c(
		(any(x%in%time.catg.cov) & any(x%in%time.cont.cov)),#models with a categorical and continuous temporal effect
		(any(x=="sample_month") & any(x=="time"))			#models with an effect of both surveys and sampling occasions
	))==FALSE}		
				
				

			
