####################################################
####PARAMETRIC BOOTSTRAPING (DERIVED PARAMETERS)####
####################################################

	##########
	####Estimated parameters and their VCOV matrices
	##########
	{
	###List of estimates and variance-covariance matrices of number of Unmarked Individuals
	U.sizes<-list(estimates=best.model[[1]]$results$real%>%
			mutate(rownam=rownames(best.model[[1]]$results$real))%>%
			filter(str_detect(rownam,'U\\sg'))%>%
			mutate(time=gsub(x=rownam,"^.*t","")))
	U.sizes$estimates<-merge(U.sizes$estimates,interval_info,by="time")%>%
		arrange(force.numeric(time))
	U.index<-which(str_detect(rownames(best.model[[1]]$design.matrix),"U\\s"))
	U.sizes$vcv<-best.model[[1]]$results$real.vcv[U.index,U.index]

	###List of estimates and variance-covariance matrices of conservative estimate of total population size (U + MNKA)
	cons.N.sizes<-list()
	cons.N.sizes$estimates<-U.sizes$estimates%>%merge(survey.mnka%>%mutate(time=force.numeric(time)),by='time',all.y=FALSE)
	cons.N.sizes$estimates%<>%mutate(estimate=estimate+surv.mnka,lcl=lcl+surv.mnka,ucl=ucl+surv.mnka)%>%
		arrange(force.numeric(time))
	cons.N.sizes$vcv<-U.sizes$vcv

	###List of estimates and variance-covariance matrices of survival estimates
	survival<-model.average(best.model,"Phi",vcv=T,drop=F)
	survival$estimates<-merge(survival$estimates,interval_info,by="time")%>%
		arrange(force.numeric(time))
	survival$estimates<-survival$estimates[survival$estimates$fixed!="Fixed",] #removendo valores fixados de Phi
	survival$estimates<-survival$estimates[!duplicated(survival$estimates[,c("time")]),]
	survival$vcv<-survival$vcv[
		as.character(survival$estimates$par.index),
		as.character(survival$estimates$par.index)]
	}
	##########
	####Boostraping preamble (defining empty objects to be filled, number of simulations and intervals' length)
	##########
	{
	N.mc<-list()
	Phi.mc<-list()
	Recrut.mc<-list()
	Fin.rate.inc.mc<-list()
	Trend.mc<-list()
	n.simu<-10000
	intervals<-na.omit(agout.proc$time.intervals[cumsum(c(agout.proc$begin.time,agout.proc$time.intervals))%in%cons.N.sizes$estimates$survey.end])
	}
	##########
	####Sampling N and Phi values from its distributions
	##########
	{
	###Looping over each 6-day interval
	for (n in 1:nrow(cons.N.sizes$estimates)){
	
		###Sample 10 thousand values of N from a log-normal distribution
		N.mc[[n]]<-rlnorm(n.simu, meanlog = log(cons.N.sizes$estimates[n,"estimate"]), 
			###The log-transform Standard Deviation of N is determined through the Delta Method
			sdlog = deltamethod(g=~log(x1),
				mean=cons.N.sizes$estimates[n,"estimate"],
				cov=cons.N.sizes$vcv[n,n],ses=TRUE))
				
		###For every 6-day interval but the last, sample 10 thousand logit-transformed Phi values from a normal distribution
		if(n!=nrow(cons.N.sizes$estimates)){	
		Phi.mc[[n]]<-rnorm(n.simu,
			mean = logit(survival[[1]][n,"estimate"]), 
			###The logit transformed Standard Deviation of Phi is determined through the Delta Method
			sd = deltamethod(g=~log(x1/(1-x1)),
					mean=survival[[1]][n,"estimate"],
					cov=survival[[2]][
						as.character(survival[[1]][n,"par.index"]),
						as.character(survival[[1]][n,"par.index"])],
					ses=TRUE))
			}
		}
	}
	##########
	####Generating derived parameters distributions from sampled N and Phi values
	##########
	{
	###Looping over each 6-day interval but the last
	for(n in 1:(nrow(cons.N.sizes$estimates)-1)){ ###
		###Finite rate of increase, `lambda`
			Fin.rate.inc.mc[[n]]<-N.mc[[n+1]]/N.mc[[n]]
		###Recruitment
			Recrut.mc[[n]]<-N.mc[[n+1]]-((inv.logit(Phi.mc[[n]])^intervals[n]) * N.mc[[n]])
			##Truncate negative recruitment (define as 0)
			Recrut.mc[[n]][Recrut.mc[[n]]<0]<-0			
			##Transform absolute recruitment into monthly rate
			Recrut.mc[[n]]<-Recrut.mc[[n]]/intervals[n] 
			}
	}
	##########
	###Create summaries of derived parameters
	##########
	{
	Recrut.estm2<-data.frame(
		mean=sapply(Recrut.mc,quantile,.5),
		sd=sapply(Recrut.mc,sd),
		lcl=sapply(Recrut.mc,quantile,.025),
		ucl=sapply(Recrut.mc,quantile,.975),
		time=force.numeric(cons.N.sizes$estimates$time[1:n]))%>%
		merge(interval_info,by='time')
	Fin.rate.inc.estm<-data.frame(
		mean=sapply(Fin.rate.inc.mc,mean),
		sd=sapply(Fin.rate.inc.mc,sd),
		lcl=sapply(Fin.rate.inc.mc,quantile,.025),
		ucl=sapply(Fin.rate.inc.mc,quantile,.975),
		time=force.numeric(cons.N.sizes$estimates$time[1:n]))%>%
		merge(interval_info,by='time')		
	#####
	###Overall growth during the study period
	#####
		Fin.rate.overall<-N.mc[[nrow(cons.N.sizes$estimates)]]/N.mc[[1]]
		Fin.rate.ov.mean<-mean(Fin.rate.overall)
		Fin.rate.ov.se<-sd(Fin.rate.overall)
		Fin.rate.ov.lcl<-quantile(Fin.rate.overall,0.025)
		Fin.rate.ov.ucl<-quantile(Fin.rate.overall,0.975)
	}