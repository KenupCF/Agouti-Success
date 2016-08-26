	###############################
	####Stepwise Model Selection###
	###############################
	
	#########################################
	##Running RDPNE Models (step 1, alpha)###
	#########################################
	{
	###Run all possible alpha models inside a single function
	agout.models.alpha<-function(threads,run,initial,data,options,ddl){
		##'U' is always modelled as a function of time
			if(closed.months==TRUE){U.time<-list(formula=~sample_month)} else{
				U.time<-list(formula=~1+time)}
		##Creating formula objects to model 'alpha'
			alpha.covariates<-c(
						"captive","age.t","sex",
						"effort.trapdays",
						"recordloss.agout.total",
						"time")			
			formulas<-"~1" #first element of 'formulas' is the formula for a constant parameter
			mod.list<-list()
			
		##Generating all possible combinations of alpha covariates
		##(while removing unreasonable ones; see 'Aux Functions.R')
			for(m in 1:length(alpha.covariates)){
				combin<-combn(alpha.covariates,m)
					if(m==length(alpha.covariates)){combin%<>%t%>%t}
				combin<-combin[,apply(combin,2,imposib.mod)]
					if(m==1){combin%<>%t}
					if(m!=1 & is.null(dim(combin))){combin%<>%t%>%t}
					if(!any(dim(combin)==0)){	
						mod.list[[m]]<-combin
						formulas<-c(formulas,paste("~1 +",apply(mod.list[[m]],2,paste0,collapse=" + ",sep="")))
			}}					
				
			##Converting strings to formula objects
			alpha<-lapply(formulas,function(x){
				y<-list(formula=as.formula(x))
				environment(y[[1]])<-environment(U.time[[1]])
				return(y)
				})
			names(alpha)<-paste0('alpha.',1:length(alpha))
			##Assigning formulas to function enviroment
			for(x in names(alpha)){assign(x,alpha[[x]])}
			
			##Rest of parameters is kept constant
			sigma.dot<-list(formula=~1)
			Phi.dot=list(formula=~1)
			GammaPrime.dot<-list(formula=~1)
			GammaDoublePrime.dot<-list(formula=~1)
			##Wrapping function up
				cml<-RMark::create.model.list("PoissonMR")	
				results=RMark::mark.wrapper(cml,data=data,
				adjust=TRUE,
				initial=initial, 
				threads=threads, realvcv=TRUE,
				ddl=ddl,run=run) 
			return(results) 
			}
	agout.results.alpha<-agout.models.alpha(
		run=TRUE,
		data=agout.proc,
		ddl=true.agout.ddl,
		initial=dot.fix.init.siman,
		threads=-1)
	###################
	##Model Selection##
	###################
	agout.results.alpha.bkp<-agout.results.alpha #Backup of full results
	##Temp object with selection model table
		temp<-agout.results.alpha$model.table[order(agout.results.alpha$model.table$weight,decreasing=T),]
	##Removing models without convergence
		if(any(temp$AICc==-Inf)){
			agout.results.alpha<-remove.mark(agout.results.alpha,as.numeric(rownames(agout.results.alpha$model.table[agout.results.alpha$model.table$AICc==-Inf,])))}
	##Defining models to be excluded from further model selection (those with Delta AICc larger than 2)			
		cut<-as.numeric(rownames(temp[temp$DeltaAICc > 2,]))
	##Removing such models
		agout.results.alpha<-RMark::remove.mark(agout.results.alpha,cut)
	}		
			
	#########################################
	##Running RDPNE Models (step 2, sigma)###
	#########################################
	{
	###Run all possible sigma models inside a single function
		agout.models.sigma<-function(threads,run,initial,data,options,ddl){
			##'U' is always modelled as a function of time
				if(closed.months==TRUE){U.time<-list(formula=~sample_month)} else{
					U.time<-list(formula=~1+time)}
			##Creating formula objects to model 'sigma'
				sigma.covariates<-c(
					"time","sample_month",
					"effort.var","recordloss.agout.var")	
			formulas<-"~1" #first element of 'formulas' is the formula for a constant parameter
			mod.list<-list()
			
	##Generating all possible combinations of sigma covariates
	##(while removing unreasonable ones; see 'Aux Functions.R')
		for(m in 1:length(sigma.covariates)){
			combin<-combn(sigma.covariates,m)
				if(m==length(sigma.covariates)){combin%<>%t%>%t}
			combin<-combin[,apply(combin,2,imposib.mod)]
				if(m==1){combin%<>%t}
				if(m!=1 & is.null(dim(combin))){combin%<>%t%>%t}
				if(!any(dim(combin)==0)){	
					mod.list[[m]]<-combin
					formulas<-c(formulas,paste("~1 +",apply(mod.list[[m]],2,paste0,collapse=" + ",sep="")))
		}}					
			
	##Converting strings to formula objects
	sigma<-lapply(formulas,function(x){
		y<-list(formula=as.formula(x))
		environment(y[[1]])<-environment(U.time[[1]])
		return(y)
		})
	names(sigma)<-paste0('sigma.',1:length(sigma))
	##Assigning formulas to function enviroment
	for(x in names(sigma)){assign(x,sigma[[x]])}
	##Add an extra model to the enviroment: Sigma Zero (No heterogeneity)
		sigma.zero<-list(formula=~1,fixed=0) ##sem heterogeneidade alguma!
	##Rest of parameters is kept constant
		Phi.dot=list(formula=~1)
		GammaPrime.dot<-list(formula=~1)
		GammaDoublePrime.dot<-list(formula=~1)	
	
	##Extracting models selected in the first step (for alpha)
		alpha<-list()
		for (m in 1:nrow(agout.results.alpha$model.table)){
			alpha[[m]]<-agout.results.alpha[[m]]$model.parameters$alpha}
			names(alpha)<-paste("alpha",1:length(alpha),sep=".")
		for (l in 1:length(alpha)){environment(alpha[[l]]$formula)<-environment(sigma.zero$formula)}
		for (l in 1:length(alpha)){assign(names(alpha)[l], alpha[[l]])}	
		###Define names of effort and rloss alpha and sigma models
		effort.alpha.models<-paste("alpha",as.character(which(sapply(alpha,function(x){
			any(str_detect(as.character(x$formula),"effort"))}))),sep=".")
		recloss.alpha.models<-paste("alpha",as.character(which(sapply(alpha,function(x){
			any(str_detect(as.character(x$formula),"recordloss"))}))),sep=".")
		effort.sigma.models<-paste("sigma",as.character(which(sapply(sigma,function(x){
			any(str_detect(as.character(x$formula),"effort"))}))),sep=".")
		recloss.sigma.models<-paste("sigma",as.character(which(sapply(sigma,function(x){
			any(str_detect(as.character(x$formula),"recordloss"))}))),sep=".")	
		##Wrapping function up
			cml<-create.model.list("PoissonMR")
			##Keeping only sigma 'eff var' or 'rloss var' that have the corresponding alpha covariate
			##('eff' or 'rloss')
			cml<-cml[!(
				( cml$sigma%in%effort.sigma.models & (cml$alpha%in%effort.alpha.models==FALSE) ) |
				( cml$sigma%in%recloss.sigma.models & (cml$alpha%in%recloss.alpha.models==FALSE) ) ),]
					results=mark.wrapper(cml,data=data,
					adjust=TRUE,
					initial=initial, 
					threads=threads, realvcv=TRUE,
					ddl=ddl,run=run) 
				return(results) 
				}
		agout.results.sigma<-agout.models.sigma(
			run=TRUE,
			data=agout.proc,
			ddl=true.agout.ddl,
			initial=dot.fix.init.siman,
			threads=-1)
	###################
	##Model Selection##
	###################
		agout.results.sigma.bkp<-agout.results.sigma #Backup of full results
	##Temp object with selection model table
		temp<-agout.results.sigma$model.table[order(agout.results.sigma$model.table$weight,decreasing=T),]
	##Removing models without convergence
		if(any(temp$AICc==-Inf)){
			agout.results.sigma<-remove.mark(agout.results.sigma,as.numeric(rownames(agout.results.sigma$model.table[agout.results.sigma$model.table$AICc==-Inf,])))}
	##Defining models to be excluded from further model selection (those with Delta AICc larger than 2)			
		cut<-as.numeric(rownames(temp[temp$DeltaAICc > 2,]))
	##Removing such models
		agout.results.sigma<-RMark::remove.mark(agout.results.sigma,cut)
	}
		
	#########################################
	##Running RDPNE Models (step 3, Phi)#####
	#########################################
	{
	###Run all possible Phi models inside a single function
	agout.models.Phi<-function(threads,run,initial,data,options,ddl){
			##'U' is always modelled as a function of time
			if(closed.months==TRUE){U.time<-list(formula=~sample_month)} else{
				U.time<-list(formula=~1+time)}
			##Creating formula objects to model 'Phi'
				Phi.covariates<-c(
					"captive","sample_month",
					#"mo.rfall","mo.rfall.30d_delay","mo.rfall.60d_delay", #climatic variables (not used)
					"age.t","sex")
				formulas<-"~1" #first element of 'formulas' is the formula for a constant parameter
				mod.list<-list()
				
			##Generating all possible combinations of Phi covariates
			##(while removing unreasonable ones; see 'Aux Functions.R')
				for(m in 1:length(Phi.covariates)){
					combin<-combn(Phi.covariates,m)
						if(m==length(Phi.covariates)){combin%<>%t%>%t}
					combin<-combin[,apply(combin,2,imposib.mod)]
						if(m==1){combin%<>%t}
						if(m!=1 & is.null(dim(combin))){combin%<>%t%>%t}
						if(!any(dim(combin)==0)){	
							mod.list[[m]]<-combin
							formulas<-c(formulas,paste("~1 +",apply(mod.list[[m]],2,paste0,collapse=" + ",sep="")))
				}}					
		
			##Converting strings to formula objects
			Phi<-lapply(formulas,function(x){
				y<-list(formula=as.formula(x))
				environment(y[[1]])<-environment(U.time[[1]])
				return(y)
				})
			names(Phi)<-paste0('Phi.',1:length(Phi))
			##Assigning formulas to function enviroment
			for(x in names(Phi)){assign(x,Phi[[x]])}
	
			##Extracting models selected in the second step (for alpha and sigma)
				alpha<-list()
				sigma<-list()
				for (m in 1:nrow(agout.results.sigma$model.table)){
					alpha[[m]]<-agout.results.sigma[[m]]$model.parameters$alpha
					sigma[[m]]<-agout.results.sigma[[m]]$model.parameters$sigma}
				names(alpha)<-paste("alpha",1:length(alpha),sep=".")
				names(sigma)<-paste("sigma",1:length(sigma),sep=".")
				for (l in 1:length(alpha)){
					environment(alpha[[l]]$formula)<-environment(Phi.1$formula)
					environment(sigma[[l]]$formula)<-environment(Phi.1$formula)
					assign(names(alpha)[l], alpha[[l]])
					assign(names(sigma)[l], sigma[[l]])}
		##Rest of parameters is kept constant
			GammaPrime.dot<-list(formula=~1)
			GammaDoublePrime.dot<-list(formula=~1)				
		##Wrapping function up
			cml<-create.model.list("PoissonMR")
		##Ensuring that only previously selected combinations of alpha and sigma are kept for the third step
			cml<-cml[substr(cml$alpha,6,nchar(cml$alpha))==substr(cml$sigma,6,nchar(cml$sigma)),]
			results=mark.wrapper(cml,data=data,
			adjust=TRUE,
			initial=initial, 
			threads=threads, realvcv=TRUE,
			ddl=ddl,run=run) 	
}	
	agout.results.Phi<-agout.models.Phi(run=TRUE,data=agout.proc,ddl=true.agout.ddl,initial=dot.fix.init.siman,threads=-1)
	###################
	##Model Selection##
	###################
	agout.results.Phi.bkp<-agout.results.Phi #Backup of full results
	##Temp object with selection model table
		temp<-agout.results.Phi$model.table[order(agout.results.Phi$model.table$weight,decreasing=T),]
	##Removing models without convergence
		if(any(temp$AICc==-Inf)){
			agout.results.Phi<-remove.mark(agout.results.Phi,as.numeric(rownames(agout.results.Phi$model.table[agout.results.sigma$model.table$AICc==-Inf,])))}
	##Defining models to be excluded from further model selection (those with Delta AICc larger than 2)			
		cut<-as.numeric(rownames(temp[temp$DeltaAICc > 2,]))
	##Removing such models
		agout.results.Phi<-RMark::remove.mark(agout.results.Phi,cut)
	}
	#########################################
	##Running RDPNE Models (step 4, Gamma)###
	#########################################
	{
	###Run all possible Gamma models inside a single function
	agout.models.Gamma<-function(threads,run,initial,data,options,ddl){
			##MODELOS PARA 'U': UM SÓ, EM FUNÇÃO APENAS DO TEMPO (INTERVALO PRIMÁRIO OU SECUNDÁRIO)
			if(closed.months==TRUE){U.time<-list(formula=~sample_month)} else{
				U.time<-list(formula=~1+time)}
					
			##Creating formula objects to model 'Gamma'
				GammaPrime.covariates<-c("sample_month")
				formulas<-"~1" #first element of 'formulas' is the formula for a constant parameter
				mod.list<-list()
				
			##Generating all possible combinations of Gamma covariates
			##(while removing unreasonable ones; see 'Aux Functions.R')
				for(m in 1:length(GammaPrime.covariates)){
					combin<-combn(GammaPrime.covariates,m)
						if(m==length(GammaPrime.covariates)){combin%<>%t%>%t}
					combin<-combin[,apply(combin,2,imposib.mod)]
						if(m==1){combin%<>%t}
						if(m!=1 & is.null(dim(combin))){combin%<>%t%>%t}
						if(!any(dim(combin)==0)){	
							mod.list[[m]]<-combin
							formulas<-c(formulas,paste("~1 +",apply(mod.list[[m]],2,paste0,collapse=" + ",sep="")))
				}}	
			##Converting strings to formula objects
			GammaPrime<-lapply(formulas,function(x){
				y<-list(formula=as.formula(x))
				environment(y[[1]])<-environment(U.time[[1]])
				return(y)
				})
			GammaDoublePrime<-GammaPrime
			names(GammaPrime)<-paste0('GammaPrime.',1:length(GammaPrime))
			names(GammaDoublePrime)<-paste0('GammaDoublePrime.',1:length(GammaPrime))
			##Assigning formulas to function enviroment
			for(x in 1:length(names(GammaPrime))){
				assign(names(GammaPrime)[x],GammaPrime[[x]])
				assign(names(GammaDoublePrime)[x],GammaDoublePrime[[x]])
			}
			##Extracting models selected in the third step (for alpha, sigma and Phi)
				alpha<-list()
				sigma<-list()
				Phi<-list()
				for (m in 1:nrow(agout.results.Phi$model.table)){
					alpha[[m]]<-agout.results.Phi[[m]]$model.parameters$alpha
					sigma[[m]]<-agout.results.Phi[[m]]$model.parameters$sigma
					Phi[[m]]<-agout.results.Phi[[m]]$model.parameters$Phi}
				names(alpha)<-paste("alpha",1:length(alpha),sep=".")
				names(sigma)<-paste("sigma",1:length(sigma),sep=".")
				names(Phi)<-paste("Phi",1:length(Phi),sep=".")
				for (l in 1:length(alpha)){
					environment(alpha[[l]]$formula)<-environment(GammaPrime.1$formula)
					environment(sigma[[l]]$formula)<-environment(GammaPrime.1$formula)
					environment(Phi[[l]]$formula)<-environment(GammaPrime.1$formula)
					assign(names(alpha)[l], alpha[[l]])
					assign(names(sigma)[l], sigma[[l]])
					assign(names(Phi)[l], Phi[[l]])}
			##Add an extra model to the enviroment: GammaPrime and GammaDouble Zero (No temporary migration)
					GammaPrime.zero<-list(formula=~1,fixed=0)
					GammaDoublePrime.zero<-list(formula=~1,fixed=0)	
			##Wrapping function up
				cml<-create.model.list("PoissonMR")
			##Ensuring that only previously selected combinations of alpha, sigma and Phi are kept
			##for the fourth step
				cml<-cml[
					(substr(cml$alpha,6,nchar(cml$alpha))==substr(cml$sigma,6,nchar(cml$sigma))) &
					(substr(cml$sigma,6,nchar(cml$sigma))==substr(cml$Phi,4,nchar(cml$Phi))),]
			##Ensuring that only corresponding models of GammaPrime and GammaDoublePrime are run
				cml<-cml[substr(cml$GammaPrime,11,nchar(cml$GammaPrime))==
					substr(cml$GammaDoublePrime,17,nchar(cml$GammaDoublePrime)),]	
				results=mark.wrapper(cml,data=data,
				adjust=TRUE,
				initial=initial, 
				threads=threads, realvcv=TRUE,
				ddl=ddl,run=run) 	
	}
	agout.results.Gamma<-agout.models.Gamma(run=TRUE,data=agout.proc,ddl=true.agout.ddl,initial=dot.fix.init.siman,threads=-1)
	###################
	##Model Selection##
	###################
	agout.results.Gamma.bkp<-agout.results.Gamma #Backup of full results
	##Temp object with selection model table
		temp<-agout.results.Gamma$model.table[order(agout.results.Gamma$model.table$weight,decreasing=T),]
	##Removing models without convergence
		if(any(temp$AICc==-Inf)){
			agout.results.Gamma<-remove.mark(agout.results.Gamma,as.numeric(rownames(agout.results.Gamma$model.table[agout.results.sigma$model.table$AICc==-Inf,])))}
	##Defining models to be excluded from further model selection (those with Delta AICc larger than 2)			
		cut<-as.numeric(rownames(temp[temp$DeltaAICc > 2,]))
	##Removing such models
		agout.results.subset<-remove.mark(agout.results.Gamma,as.numeric(rownames(temp[temp$DeltaAICc > 2,])))
	##All final models
		agout.results<-agout.results.Gamma
	}
	
	####
	##Keeping the most supported model for subsequent analysis
	####
	best.model<-remove.mark(agout.results,
		as.numeric(rownames(agout.results$model.table[agout.results$model.table$DeltaAICc != 0,])))