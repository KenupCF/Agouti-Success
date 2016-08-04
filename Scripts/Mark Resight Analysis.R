	###Importing Data, Formatting
		source(".\\RDPNE.R")
		
	###Temporal Indepence; Formatting Pt 2
		source(".\\RDPNE.R")
	
	###Generate RDPNE encounterhistory
		source(".\\RDPNE.R")
	
	###RMark Processing Data
		source(".\\RMark Process")
		
	###Defining which covariates are individual, temporal, categoric or continuous
		##Temporal, continous covariates
		time.cont.cov<-c("effort.trapdays","recordloss.agout.total","effort.var","recordloss.agout.var","age.t")
		##Temporal, categoric covariates
		time.catg.cov<-c("time","sample_month")
		##Individual covariates
		indiv.cov<-c("captive","age.t","sex")
	
	
	####Stepwise Model Selection###
	
	###First step: Modelling alpha (resighting rates)
	
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
			names(alpha)<-paste0('alpha.',zero_pad(1:length(alpha),nchar(length(alpha))))
			##Assigning formulas to function enviroment
			for(x in names(alpha)){assign(x,alpha[[x]])}
			
			##Rest of parameters is kept constant
			sigma.dot<-list(formula=~1)
			Phi.dot=list(formula=~1)
			GammaPrime.dot<-list(formula=~1)
			GammaDoublePrime.dot<-list(formula=~1)
			##Wrapping function up
				cml<-RMark::create.model.list("PoissonMR")	
				cml<-cml[
					substr(cml$GammaPrime,11,nchar(cml$GammaPrime)) == 
					substr(cml$GammaDoublePrime,17,nchar(cml$GammaDoublePrime))
					,]
				results=RMark::mark.wrapper(cml,data=data,
				adjust=TRUE,
				initial=initial, 
				threads=threads, realvcv=TRUE,
				ddl=ddl,run=run) 
			return(results) 
			}
			
	##Running RDPNE Models (step 1, alpha)
	agout.results.alpha<-agout.models.alpha(
		run=TRUE,
		data=agout.proc,
		ddl=true.agout.ddl,
		initial=dot.fix.init.siman,
		threads=-1)
	##Model Selection	
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
			
			
	###Second step: Modelling sigma (individual heterogeneity in resighting rates)
	
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
	names(sigma)<-paste0('sigma.',zero_pad(1:length(sigma),nchar(length(sigma))))
	##Assigning formulas to function enviroment
	for(x in names(sigma)){assign(x,sigma[[x]])}
	##Add an extra model to the enviroment: Sigma Zero (No heterogeneity)
		sigma.zero<-list(formula=~1,fixed=0) ##sem heterogeneidade alguma!
	##Rest of parameters is kept constant
		Phi.dot=list(formula=~1)
		GammaPrime.dot<-list(formula=~1)
		GammaDoublePrime.dot<-list(formula=~1)	
	
	####@#@#@#@#@#@##@#@# PAREI AQUI!!!!!#@#@#@#@#@#@#@#@#@#@#@#@
	
	##Extracting models selected in the first step (for alpha)
		alpha<-list()
						for (m in 1:nrow(agout.results.alpha$model.table)){
							alpha[[m]]<-agout.results.alpha[[m]]$model.parameters$alpha}
							names(alpha)<-paste("alpha",1:length(alpha),sep=".")
						for (l in 1:length(alpha)){environment(alpha[[l]]$formula)<-environment(sigma.dot$formula)}
						for (l in 1:length(alpha)){assign(names(alpha)[l], alpha[[l]])}	
						effort.alpha.models<-paste("alpha",as.character(which(sapply(alpha,function(x){
							any(str_detect(as.character(x$formula),"effort"))}
						))),sep=".")
						recloss.alpha.models<-paste("alpha",as.character(which(sapply(alpha,function(x){
							any(str_detect(as.character(x$formula),"recordloss"))}
						))),sep=".")
				##RESTO CONSTANTE
					Phi.dot=list(formula=~1)
					GammaPrime.dot<-list(formula=~1)
					GammaDoublePrime.dot<-list(formula=~1)				
				###RESTO DA FUNÇÃO	
					cml<-create.model.list("PoissonMR")
					cml<-cml[(
						( cml$sigma%in%effort.sigma.models & (cml$alpha%in%effort.alpha.models==FALSE) ) |
						( cml$sigma%in%recloss.sigma.models & (cml$alpha%in%recloss.alpha.models==FALSE) ) )==FALSE,]
					results=mark.wrapper(cml,data=data,
					adjust=TRUE,
					initial=initial, 
					threads=threads, realvcv=TRUE,
					ddl=ddl,run=run) 
				return(results) 
				}
	
	##Running RDPNE Models (step 2, sigma)
	agout.results.sigma<-agout.models.sigma(
		run=TRUE,
		data=agout.proc,
		ddl=true.agout.ddl,
		initial=dot.fix.init.siman,
		threads=-1)
	##Model Selection	
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
		
		
	#################################		
	###SELEÇÃO DE MODELOS PARA PHI###
	#################################
		agout.models.Phi<-function(threads,run,initial,data,options,ddl){
			##MODELOS PARA 'U': UM SÓ, EM FUNÇÃO APENAS DO TEMPO (INTERVALO PRIMÁRIO OU SECUNDÁRIO)
			if(closed.months==TRUE){U.time<-list(formula=~sample_month)} else{
				U.time<-list(formula=~1+time)}
		########################	
		####MODELOS PARA PHI####
		########################
				Phi.dot<-list(formula=~1)
				Phi.covariates<-c(
					"captive","sample_month",
					#"mo.rfall","mo.rfall.30d_delay","mo.rfall.60d_delay", #climatic variables (not used)
					"age.t","sex")
				mod.list<-list()
				formulas<-character()
				f<-1
				for(m in 2:length(Phi.covariates)){
					temp<-as.matrix(combn(Phi.covariates,m))
					temp<-temp[,apply(temp,2,imposib.mod)]
					mod.list[[m]]<-as.matrix(temp)
					if(any(dim(mod.list[[m]])==0)){
						mod.list[[m]]<-NULL
						break} else {
					formulas[f:(f-1+ncol(mod.list[[m]]))]<-paste("~1 +",apply(mod.list[[m]],2,paste0,collapse="",sep="+"))
					f<-length(formulas)+1
					}
				}
				formulas<-gsub('.{1}$', '', formulas)
				formulas[f:(f-1+length(Phi.covariates))]<-paste("~1",Phi.covariates,sep="+")
				Phi<-list()
				for (f in 1:length(formulas)){
					Phi[[f]]<-list(formula=as.formula(formulas[f]))
					names(Phi)[f]<-paste("Phi",f,sep=".")
					environment(Phi[[f]]$formula)<-environment(Phi.dot$formula)
					assign(names(Phi)[f], Phi[[f]])}
		##RECUPERANDO OS MODELOS SELECIONADOS PARA 'alpha' E 'sigma'
				alpha<-list()
				for (m in 1:nrow(agout.results.sigma$model.table)){
					alpha[[m]]<-agout.results.sigma[[m]]$model.parameters$alpha}
					names(alpha)<-paste("alpha",1:length(alpha),sep=".")
				for (l in 1:length(alpha)){environment(alpha[[l]]$formula)<-environment(Phi.dot$formula)}
				for (l in 1:length(alpha)){assign(names(alpha)[l], alpha[[l]])}
				sigma<-list()
				for (m in 1:nrow(agout.results.sigma$model.table)){
					sigma[[m]]<-agout.results.sigma[[m]]$model.parameters$sigma}
					names(sigma)<-paste("sigma",1:length(sigma),sep=".")
				for (l in 1:length(sigma)){environment(sigma[[l]]$formula)<-environment(Phi.dot$formula)}
				for (l in 1:length(sigma)){assign(names(sigma)[l], sigma[[l]])}
		##RESTO CONSTANTE
			GammaPrime.dot<-list(formula=~1)
			GammaDoublePrime.dot<-list(formula=~1)				
		###RESTO DA FUNÇÃO	
			cml<-create.model.list("PoissonMR")
			cml<-cml[substr(cml$alpha,6,nchar(cml$alpha))==substr(cml$sigma,6,nchar(cml$sigma)),]
			results=mark.wrapper(cml,data=data,
			adjust=TRUE,
			initial=initial, 
			threads=threads, realvcv=TRUE,
			ddl=ddl,run=run) 	
}	
	start.computing<-now()
	agout.results.Phi<-agout.models.Phi(run=TRUE,data=agout.proc,ddl=true.agout.ddl,initial=dot.fix.init.siman,threads=-1)
	run.time.Phi<-difftime(now(),start.computing,units="hours")
	save.image(".\\Caio.RData")
	###REALIZANDOO SELEÇÃO DOS MELHORES MODELOS, Phi:	
		agout.results.Phi.bkp<-agout.results.Phi #salva backup dos resultados originais
		temp<-agout.results.Phi$model.table[order(agout.results.Phi$model.table$weight,decreasing=T),] #temp com tabela dos modelos
		###Removendo modelos zoados###
			if(any(temp$AICc==-Inf)){
				agout.results.pt1<-remove.mark(agout.results.pt1,as.numeric(rownames(agout.results.pt1$model.table[agout.results.pt1$model.table$AICc==-Inf,])))}
		cum.weight<-cumsum(temp$weight)
		#if (min(cum.weight)<0.95) {cut<-as.numeric(rownames(temp[cum.weight>0.95,]))} else{cut<-as.numeric(rownames(temp[1,]))}
		cut<-as.numeric(rownames(temp[temp$DeltaAICc > 2,]))
		drop.agout.results.Phi<-remove.mark(agout.results.Phi,as.numeric(setdiff(rownames(temp),cut)))
		agout.results.Phi<-remove.mark(agout.results.Phi,cut)
		temp<-agout.results.Phi$model.table[order(agout.results.Phi$model.table$weight,decreasing=T),]
		if(any(duplicated(temp$AICc)==TRUE)){
			dupli.agout.results.Phi<-remove.mark(agout.results.Phi,as.numeric(rownames(temp[(duplicated(temp$AICc)==F),])))
			agout.results.Phi<-remove.mark(agout.results.Phi,as.numeric(rownames(temp[(duplicated(temp$AICc)==T),])))}
		temp<-agout.results.Phi$model.table[order(agout.results.Phi$model.table$weight,decreasing=T),]
	###################################		
	###SELEÇÃO DE MODELOS PARA GAMMA###
	###################################
		agout.models.Gamma<-function(threads,run,initial,data,options,ddl){
			##MODELOS PARA 'U': UM SÓ, EM FUNÇÃO APENAS DO TEMPO (INTERVALO PRIMÁRIO OU SECUNDÁRIO)
			if(closed.months==TRUE){U.time<-list(formula=~sample_month)} else{
				U.time<-list(formula=~1+time)}
			####MODELOS PARA GAMMA##
				GammaPrime.dot<-list(formula=~1)
					GammaDoublePrime.dot<-list(formula=~1)
				GammaPrime.zero<-list(formula=~1,fixed=0)
					GammaDoublePrime.zero<-list(formula=~1,fixed=0)
				GammaPrime.covariates<-c(
					# "captive","age.t","sex",
					#"mo.rfall","mo.rfall.30d_delay","mo.rfall.60d_delay", #climatic variables (not used)
					"sample_month")
				mod.list<-list()
				formulas<-character()
				f<-1
				if(length(GammaPrime.covariates)>1){
				for(m in 2:length(GammaPrime.covariates)){
					temp<-as.matrix(combn(GammaPrime.covariates,m))
					temp<-temp[,apply(temp,2,imposib.mod)]
					mod.list[[m]]<-as.matrix(temp)
					if(any(dim(mod.list[[m]])==0)){
						mod.list[[m]]<-NULL
						break} else {
					formulas[f:(f-1+ncol(mod.list[[m]]))]<-paste("~1 +",apply(mod.list[[m]],2,paste0,collapse="",sep="+"))
					f<-length(formulas)+1
					}
				} }
				formulas<-gsub('.{1}$', '', formulas)
				formulas<-c(formulas,paste("~1",GammaPrime.covariates,sep="+"))
				GammaPrime<-list()
				GammaDoublePrime<-list()
				for (f in 1:length(formulas)){
					GammaPrime[[f]]<-list(formula=as.formula(formulas[f]))
						GammaDoublePrime[[f]]<-GammaPrime[[f]]
					names(GammaPrime)[f]<-paste("GammaPrime",f,sep=".")
						names(GammaDoublePrime)[f]<-paste("GammaDoublePrime",f,sep=".")
					environment(GammaPrime[[f]]$formula)<-environment(GammaPrime.dot$formula)
						environment(GammaDoublePrime[[f]]$formula)<-environment(GammaPrime.dot$formula)
					assign(names(GammaPrime)[f], GammaPrime[[f]])
						assign(names(GammaDoublePrime)[f], GammaDoublePrime[[f]])}
			##RECUPERANDO OS MODELOS SELECIONADOS PARA 'alpha' E 'sigma' E 'Phi'
					alpha<-list()
					for (m in 1:nrow(agout.results.Phi$model.table)){
						alpha[[m]]<-agout.results.Phi[[m]]$model.parameters$alpha}
						names(alpha)<-paste("alpha",1:length(alpha),sep=".")
					for (l in 1:length(alpha)){environment(alpha[[l]]$formula)<-environment(GammaPrime.dot$formula)}
					for (l in 1:length(alpha)){assign(names(alpha)[l], alpha[[l]])}
					sigma<-list()
					for (m in 1:nrow(agout.results.Phi$model.table)){
						sigma[[m]]<-agout.results.Phi[[m]]$model.parameters$sigma}
						names(sigma)<-paste("sigma",1:length(sigma),sep=".")
					for (l in 1:length(sigma)){environment(sigma[[l]]$formula)<-environment(GammaPrime.dot$formula)}
					for (l in 1:length(sigma)){assign(names(sigma)[l], sigma[[l]])}				
					Phi<-list()
					for (m in 1:nrow(agout.results.Phi$model.table)){
						Phi[[m]]<-agout.results.Phi[[m]]$model.parameters$Phi}
						names(Phi)<-paste("Phi",1:length(Phi),sep=".")
					for (l in 1:length(Phi)){environment(Phi[[l]]$formula)<-environment(GammaPrime.dot$formula)}
					for (l in 1:length(Phi)){assign(names(Phi)[l], Phi[[l]])}
			###RESTO DA FUNÇÃO	
				cml<-create.model.list("PoissonMR")
				cml<-cml[
					(substr(cml$alpha,6,nchar(cml$alpha))==substr(cml$sigma,6,nchar(cml$sigma))) &
					(substr(cml$sigma,6,nchar(cml$sigma))==substr(cml$Phi,4,nchar(cml$Phi))),]
				cml<-cml[substr(cml$GammaPrime,11,nchar(cml$GammaPrime))==substr(cml$GammaDoublePrime,17,nchar(cml$GammaDoublePrime)),]	
				results=mark.wrapper(cml,data=data,
				adjust=TRUE,
				initial=initial, 
				threads=threads, realvcv=TRUE,
				ddl=ddl,run=run) 	
	}
	start.computing<-now()
	agout.results.Gamma<-agout.models.Gamma(run=TRUE,data=agout.proc,ddl=true.agout.ddl,initial=dot.fix.init.siman,threads=-1)
	run.time.Gamma<-difftime(now(),start.computing,units="hours")
	run.time.total<-sum(run.time.alpha+run.time.sigma+run.time.Phi+run.time.Gamma)
	save.image(".\\Caio.RData")
	###REALIZANDOO SELEÇÃO DOS MELHORES MODELOS, Gamma:	
		agout.results.Gamma.bkp<-agout.results.Gamma #salva backup dos resultados originais
		temp<-agout.results.Gamma$model.table[order(agout.results.Gamma$model.table$weight,decreasing=T),] #temp com tabela dos modelos
		###Removendo modelos zoados###
			if(any(temp$AICc==-Inf)){
				agout.results.Gamma<-remove.mark(agout.results.Gamma,
					as.numeric(rownames(agout.results.Gamma$model.table[agout.results.Gamma$model.table$AICc==-Inf,])))}
		cum.weight<-cumsum(temp$weight)
		#if (min(cum.weight)<0.95) {cut<-as.numeric(rownames(temp[cum.weight>0.95,]))} else{cut<-as.numeric(rownames(temp[1,]))}
		# cut<-as.numeric(rownames(temp[temp$weight == 0,]))
		# drop.agout.results.Gamma<-remove.mark(agout.results.Gamma,as.numeric(setdiff(rownames(temp),cut)))
		# agout.results.Gamma<-remove.mark(agout.results.Gamma,cut)
		temp<-agout.results.Gamma$model.table[order(agout.results.Gamma$model.table$weight,decreasing=T),]
		if(any(duplicated(temp$AICc)==TRUE)){
			dupli.agout.results.Gamma<-remove.mark(agout.results.Gamma,as.numeric(rownames(temp[(duplicated(temp$AICc)==F),])))
			agout.results.Gamma<-remove.mark(agout.results.Gamma,as.numeric(rownames(temp[(duplicated(temp$AICc)==T),])))}
		temp<-agout.results.Gamma$model.table[order(agout.results.Gamma$model.table$weight,decreasing=T),]
		agout.results<-agout.results.Gamma
		agout.results.subset<-remove.mark(agout.results.Gamma,as.numeric(rownames(temp[temp$DeltaAICc > 2,])))


###Criando Objeto com Número Minimo de Marcas Vivas
		
survey.mnka<-data.frame(mnka=naive.kmark,interval=names(naive.kmark))%>%
	merge(interval_cov[,c('start','sample_month','time','interval')],by='interval')%>%
	group_by(sample_month)%>%
	mutate(surv.mnka=max(mnka))
		
######################################################################		
#####FORCANDO A DECISAO DE REALIZAR OU NAO A PONDERACAO DE MODELOS####		
# if(sum(agout.results$model.table$DeltaAICc <= 2)==1) {
	model.av<-FALSE
# } else {
	# model.average<-TRUE
	# }
##################################################
######################################################################	
###VALORES ESTIMADOS COM UM MODELO SÓ		
	
	if(!model.av){
		best.model<-remove.mark(agout.results,as.numeric(rownames(agout.results$model.table[agout.results$model.table$DeltaAICc != 0,])))
		pop.sizes<-list()
			pop.sizes$estimates<-best.model[[1]]$results$derived[[1]][1:no_int,]
			pop.sizes$vcv<-best.model[[1]]$results$derived.vcv[[1]]
			pop.sizes$estimates$interval<-1:no_int
				pop.sizes$estimates<-merge(pop.sizes$estimates,interval_cov,by="interval")
				
		U.sizes <- list()
			# temp.U<-model.average(best.model,"U",vcv=T,drop=F)
			temp.U<-list(estimates=best.model[[1]]$results$real%>%
				mutate(rownam=rownames(best.model[[1]]$results$real))%>%
				filter(str_detect(rownam,'U\\sg'))%>%
				mutate(time=gsub(x=rownam,"^.*t","")))
				temp.U$estimates<-merge(temp.U$estimates,interval_cov,by="time")
				# temp.U$estimates<-temp.U$estimates[temp.U$estimates$fixed!="Fixed",] #removendo valores fixados de Phi
				# temp.U$estimates<-temp.U$estimates[duplicated(temp.U$estimates[,c("time")])==FALSE,]
			U.sizes $estimates<-temp.U$estimates
				U.index<-which(str_detect(rownames(best.model[[1]]$design.matrix),"U\\s"))
				U.sizes$vcv<-best.model[[1]]$results$real.vcv[U.index,U.index]
				# temp<-as.character(U.sizes$estimates$par.index)
				# U.sizes$vcv<-temp.U$vcv[temp,temp]
			
		survival<-list()
			temp.surv<-model.average(best.model,"Phi",vcv=T,drop=F)
				temp.surv$estimates<-merge(temp.surv$estimates,interval_cov,by="time")
				temp.surv$estimates<-temp.surv$estimates[temp.surv$estimates$fixed!="Fixed",] #removendo valores fixados de Phi
				temp.surv$estimates<-temp.surv$estimates[duplicated(temp.surv$estimates[,c("time")])==FALSE,]
			survival$estimates<-temp.surv$estimates
				temp<-as.character(survival$estimates$par.index)
				survival$vcv<-temp.surv$vcv[temp,temp]
	}		
####COM DOIS OU MAIS MODELOS - MODEL-AVERAGING
	if(model.av){ #Se houver mais de 1 modelo com Delta menor que 2
		##PARÂMETROS REAIS
			agout.Phi.mod.avg<-model.average(agout.results,"Phi",vcv=T,drop=F)
				agout.Phi.mod.avg$estimates<-merge(agout.Phi.mod.avg$estimates,interval_cov,by="time")
			agout.alpha.mod.avg<-model.average(agout.results,"alpha",vcv=T,drop=F)
				agout.alpha.mod.avg$estimates<-merge(agout.alpha.mod.avg$estimates,interval_cov,by="time")
			agout.sigma.mod.avg<-model.average(agout.results,"sigma",vcv=T,drop=F)
				agout.sigma.mod.avg$estimates<-merge(agout.sigma.mod.avg$estimates,interval_cov,by="time")
			agout.U.mod.avg<-model.average(agout.results,"U",vcv=T,drop=F)
				agout.U.mod.avg$estimates<-merge(agout.U.mod.avg$estimates,interval_cov,by="time")
			agout.GammaP.mod.avg<-model.average(agout.results,"GammaPrime",vcv=T,drop=F)
				agout.GammaP.mod.avg$estimates<-merge(agout.GammaP.mod.avg$estimates,interval_cov,by="time")
			agout.GammaDP.mod.avg<-model.average(agout.results,"GammaDoublePrime",vcv=T,drop=F)
				agout.GammaDP.mod.avg$estimates<-merge(agout.GammaDP.mod.avg$estimates,interval_cov,by="time")
		##PARÂMETROS DERIVADOS
			###Criando objetos vazios para preencher com os valores estimados de cada modelo
				estm.N<-matrix(																			#matriz para valores de ^N				
					0,ncol=length(agout.results[[1]]$results$derived$`N Population Size`$estimate),
					nrow=nrow(agout.results$model.table))
				estm.Lmb<-matrix(																		#matriz para valores de Lambda
					0,ncol=length(agout.results[[1]]$results$derived$`Lambda Mean Resightings`$estimate),
					nrow=nrow(agout.results$model.table))	
				vcv.N<-vector("list", length=nrow(agout.results$model.table))							#lista com matriz de VCOV de N
				vcv.Lmb<-vector("list", length=nrow(agout.results$model.table))							#lista com matriz de VCOV de Lambda
				for(i in 1:nrow(agout.results$model.table)){						#'for' passando por todos os modelos (dos selecionados)					
					mod.num<-as.numeric(row.names(agout.results$model.table))		#índices dos modelos
					x<-agout.results[[mod.num[i]]]$results							#objeto com os resultados do modelo 'i'
					estm.N[i,]<-x$derived$`N Population Size`$estimate				#cada linha da matriz são as estimativas de um modelo 'i'
					estm.Lmb[i,]<-x$derived$`Lambda Mean Resightings`$estimate		#cada linha da matriz são as estimativas de um modelo 'i'
					vcv.N[[i]]<-x$derived.vcv$`N Population Size`					#cada objeto da lista é a VCOV de um modelo 'i'
					vcv.Lmb[[i]]<-x$derived.vcv$`Lambda Mean Resightings`			#cada objeto da lista é a VCOV de um modelo 'i'
					}																#fim do 'for'
			temp<-model.average(list(estimate=estm.N, weight=agout.results$model.table$weight, vcv=vcv.N))
			agout.Nhat.mod.avg<-list(estimates=data.frame(estimate=temp$estimate, se=temp$se),
				vcv.derived=temp$vcv)
			temp<-model.average(list(estimate=estm.Lmb, weight=agout.results$model.table$weight, vcv=vcv.Lmb))
			agout.Lambdahat.mod.avg<-list(estimates=data.frame(estimate=temp$estimate, se=temp$se),
				vcv.derived=temp$vcv)	
			agout.Nhat.mod.avg$estimates<-merge(agout.Nhat.mod.avg$estimates,agout.ddl.fix$U[,c("time","start","sample_month","group","captive")],by=0)
				agout.Nhat.mod.avg$estimates<-agout.Nhat.mod.avg$estimates[order(agout.Nhat.mod.avg$estimates$group,agout.Nhat.mod.avg$estimates$time),]
					agout.Nhat.mod.avg$estimates$sample_month<-factor(agout.Nhat.mod.avg$estimates$sample_month)
					rownames(agout.Nhat.mod.avg$estimates)<-agout.Nhat.mod.avg$estimates$Row.names
					agout.Nhat.mod.avg$estimates$Row.names<-NULL
			agout.Nhat.mod.avg$estimates<-merge(agout.Nhat.mod.avg$estimates,interval_cov,by="time")		
				agout.Nhat.mod.avg$estimates$dummy.index<-rownames(agout.Nhat.mod.avg$estimates)
					rownames(agout.Nhat.mod.avg$vcv.derived)<-agout.Nhat.mod.avg$estimates$dummy.index
					colnames(agout.Nhat.mod.avg$vcv.derived)<-agout.Nhat.mod.avg$estimates$dummy.index
			agout.Lambdahat.mod.avg$estimates<-merge(agout.Lambdahat.mod.avg$estimates,agout.ddl.fix$U[,c("time","start","sample_month","group","captive")],by=0)
				agout.Lambdahat.mod.avg$estimates<-agout.Lambdahat.mod.avg$estimates[order(agout.Lambdahat.mod.avg$estimates$group,agout.Lambdahat.mod.avg$estimates$time),]
					agout.Lambdahat.mod.avg$estimates$sample_month<-factor(agout.Lambdahat.mod.avg$estimates$sample_month)
					rownames(agout.Lambdahat.mod.avg$estimates)<-agout.Lambdahat.mod.avg$estimates$Row.names
					agout.Lambdahat.mod.avg$estimates$Row.names<-NULL
				agout.Lambdahat.mod.avg$estimates$dummy.index<-rownames(agout.Lambdahat.mod.avg$estimates)
					rownames(agout.Lambdahat.mod.avg$vcv.derived)<-agout.Lambdahat.mod.avg$estimates$dummy.index
					colnames(agout.Lambdahat.mod.avg$vcv.derived)<-agout.Lambdahat.mod.avg$estimates$dummy.index
				###FIXANDO VALORES DE N-HAT COMO IGUAIS PARA O MESMO MES, CASO CADA 30 DIAS SEJA FECHADO.
					if (closed.months==TRUE){
						rownames(agout.Nhat.mod.avg$vcv.derived)<-agout.Nhat.mod.avg$estimates$sample_month
						colnames(agout.Nhat.mod.avg$vcv.derived)<-agout.Nhat.mod.avg$estimates$sample_month
						agout.Nhat.mod.avg$vcv.wild<-agout.Nhat.mod.avg$vcv.derived[1:no_int,1:no_int]
						agout.Nhat.mod.avg$vcv.captive<-agout.Nhat.mod.avg$vcv.derived[(no_int+1):(2*no_int),(no_int+1):(2*no_int)]
						agout.clNhat.mod.avg<-list()
						temp.df<-agout.Nhat.mod.avg$estimates[0,]
						for (l in c(FALSE,FALSE)){
							if(l==TRUE){temp.vcv<-agout.Nhat.mod.avg$vcv.captive} else{temp.vcv<-agout.Nhat.mod.avg$vcv.wild}
							temp<-agout.Nhat.mod.avg$estimates[agout.Nhat.mod.avg$estimates$captive==l,]
							closed.estim<-tapply(temp$estimate,temp$sample_month,mean)
								closed.estim<-data.frame(closed.estim=closed.estim,sample_month=names(closed.estim))
								temp<-merge(temp,closed.estim,by="sample_month")
							closed.se<-numeric()	
							for (m in 1:length(months)){
								temp.vcv.sub<-temp.vcv[rownames(temp.vcv)==months[m],colnames(temp.vcv)==months[m]]
								closed.se[m]<-deltamethod.special("sum",temp$estimate[temp$sample_month==months[m]],temp.vcv.sub,ses=T) * 1/nrow(temp.vcv.sub)
								names(closed.se)[m]<-months[m]}
							closed.se<-data.frame(closed.se=closed.se,sample_month=names(closed.se))
							temp<-merge(temp,closed.se,by="sample_month")
							temp.df<-rbind(temp.df,temp)}
						agout.clNhat.mod.avg$estimates<-temp.df}
				##CALCULANDO INTERVALOS DE CONFIANCA PARA ESTIMATIVAS
					##N-HAT
						temp<-agout.Nhat.mod.avg$estimates
							Mt1<-c(dummy.kmark.g1,dummy.kmark.g2)
							f0<-temp$estimate-Mt1 
							C<-exp(1.96*sqrt(log(1+(temp$se/f0)^2)))
						temp$lcl<-Mt1+f0/C
						temp$ucl<-Mt1+f0*C
							temp$lcl[is.na(temp$lcl)]<-temp$estimate[is.na(temp$lcl)]
							temp$ucl[is.na(temp$ucl)]<-temp$estimate[is.na(temp$ucl)]
						if (closed.months==TRUE){	
							c.f0<-temp$closed.estim-Mt1
							c.C<-exp(1.96*sqrt(log(1+(temp$closed.se/c.f0)^2)))	
							temp$closed.lcl<-Mt1+c.f0/c.C
							temp$closed.ucl<-Mt1+c.f0*c.C
								temp$closed.lcl[is.na(temp$closed.lcl)]<-temp$closed.estim[is.na(temp$closed.lcl)]
								temp$closed.ucl[is.na(temp$closed.ucl)]<-temp$closed.estim[is.na(temp$closed.ucl)]}
						agout.Nhat.mod.avg$estimates<-temp
					##LAMBDA-HAT
						temp<-agout.Lambdahat.mod.avg$estimates
							C<-exp(1.96*sqrt(log(1+(temp$se/temp$estimate)^2)))
						temp$lcl<-temp$estimate/C
						temp$ucl<-temp$estimate*C
						agout.Lambdahat.mod.avg$estimates<-temp
			
			pop.sizes<-list()
				pop.sizes$estimates<-agout.Nhat.mod.avg$estimates[agout.Nhat.mod.avg$estimates$captive==FALSE,]
					pop.sizes$vcv<-agout.Nhat.mod.avg$vcv.derived[
						as.character(pop.sizes$estimates$dummy.index),as.character(pop.sizes$estimates$dummy.index)]
			survival<-list(wild=list(),captive=list())
			temp.surv<-agout.Phi.mod.avg$estimates[agout.Phi.mod.avg$estimates$fixed!="Fixed",] #removendo valores fixados de Phi e salvando em um objeto temporário
			temp.surv<-temp.surv[duplicated(temp.surv[,c("time")])==FALSE,] #removendo valores duplicados em relação a tempo (logo, com cohortes diferentes)
			survival$estimates<-temp.surv[temp.surv$captive==FALSE,]
				temp<-as.character(survival$estimates$par.index)
				survival$vcv<-agout.Phi.mod.avg$vcv[temp,temp]
		}
		
	pop.sizes$estimates$dummy.index<-rownames(pop.sizes$estimates)
		rownames(pop.sizes$vcv)<-pop.sizes$estimates$dummy.index
		colnames(pop.sizes$vcv)<-pop.sizes$estimates$dummy.index		
	pop.sizes$singular.estm<-pop.sizes$estimates[pop.sizes$estimates$se==0 | pop.sizes$estimates$se > 1000 ,]
		pop.sizes$singular.vcv<-pop.sizes$vcv[
			pop.sizes$estimates$se==0 | pop.sizes$estimates$se > 1000,
			pop.sizes$estimates$se==0 | pop.sizes$estimates$se > 1000]
	pop.sizes$vcv<-pop.sizes$vcv[
		pop.sizes$estimates$se!=0 & pop.sizes$estimates$se <= 1000,
		pop.sizes$estimates$se!=0 & pop.sizes$estimates$se <= 1000]
	pop.sizes$estimates<-pop.sizes$estimates[
		pop.sizes$estimates$se!=0 & pop.sizes$estimates$se <= 1000,]
	pop.sizes$estimates<-pop.sizes$estimates[order(pop.sizes$estimates$time),]
		no.est_int<-nrow(pop.sizes$estimates)
	

	###Criando estimativa conservadora de N!!!
	cons.N.sizes<-list()
	cons.N.sizes$estimates<-U.sizes$estimates%>%merge(survey.mnka%>%mutate(time=force.numeric(time)),by='time',all.y=FALSE)
	cons.N.sizes$estimates%<>%mutate(estimate=estimate+surv.mnka,lcl=lcl+surv.mnka,ucl=ucl+surv.mnka)%>%
		arrange(force.numeric(time))
	cons.N.sizes$vcv<-U.sizes$vcv
	
####################################################
####PARAMETRIC BOOTSTRAPING (DERIVED PARAMETERS)####
####################################################
###Preamble###
	N.mc<-list()
	Phi.mc<-list()
	Recrut.mc<-list()
	Fin.rate.inc.mc<-list()
	Trend.mc<-list()
	Trend.stats<-list()
	N.mc.stats<-list()
	Phi.mc.stats<-list()
	Recrut.stats<-list()
	Fin.rate.inc.stats<-list()
	N.mean<-numeric()
	N.se<-numeric()
	N.log.se<-numeric()
	Phi.mean<-numeric()
	Phi.se<-numeric()
	x<-numeric()
	n.simu<-10000
	Phi.temp<-survival$estimates
	Phi.temp<-Phi.temp[order(Phi.temp$interval,decreasing=FALSE),]
	Phi.temp.vcv<-survival$vcv[as.character(Phi.temp$par.index),as.character(Phi.temp$par.index)]
	end.N<-which(cons.N.sizes$estimates$survey.end==cons.N.sizes$estimates$time)
	beg.N<-which(cons.N.sizes$estimates$survey.beg==cons.N.sizes$estimates$time)
	end.Phi<-which(Phi.temp$survey.end==Phi.temp$time)
	beg.Phi<-which(Phi.temp$survey.beg==Phi.temp$time)
	intervals<-na.omit(agout.proc$time.intervals[cumsum(c(agout.proc$begin.time,agout.proc$time.intervals))%in%cons.N.sizes$estimates$survey.end])
	####Generating random values 
		for (n in 1:nrow(cons.N.sizes$estimates)){		
			N.mc[[n]]<-rlnorm(n.simu, meanlog = log(cons.N.sizes$estimates[n,"estimate"]), 
				sdlog = deltamethod(g=~log(x1),
					mean=cons.N.sizes$estimates[n,"estimate"],
					cov=cons.N.sizes$vcv[n,n],ses=TRUE))
				N.mc.stats$mean[n]<-mean(N.mc[[n]])
				N.mc.stats$se[n]<-sd(N.mc[[n]])				
			if(n!=nrow(cons.N.sizes$estimates)){	
			Phi.mc[[n]]<-rnorm(n.simu,
				mean = logit(Phi.temp[n,"estimate"]), 
				sd = deltamethod(g=~log(x1/(1-x1)),
							mean=Phi.temp[n,"estimate"],
							cov=Phi.temp.vcv[
								as.character(Phi.temp[n,"par.index"]),as.character(Phi.temp[n,"par.index"])],
							ses=TRUE))		
					Phi.mc.stats$mean[n]<-mean(Phi.mc[[n]])
					Phi.mc.stats$se[n]<-sd(Phi.mc[[n]])
			}
		}
	####GERANDO VALORES DERIVADOS A PARITR DAS DISTRIBUIÇÕES DOS PARAMETROS ESTIMADOS	
		for(n in 1:(nrow(cons.N.sizes$estimates)-1)){ ###
		###Taxa finita de crescimento, `lambda`
			Fin.rate.inc.mc[[n]]<-N.mc[[n+1]]/N.mc[[n]]
				Fin.rate.inc.stats$mean[n]<-mean(Fin.rate.inc.mc[[n]])
				Fin.rate.inc.stats$se[n]<-sd(Fin.rate.inc.mc[[n]])
				Fin.rate.inc.stats$lcl[n]<-quantile(Fin.rate.inc.mc[[n]],0.025)
				Fin.rate.inc.stats$ucl[n]<-quantile(Fin.rate.inc.mc[[n]],0.975)
				Fin.rate.inc.stats$time[n]<-force.numeric(cons.N.sizes$estimates$time[n])
		###Recrutamento
			Recrut.mc[[n]]<-N.mc[[n+1]]-((inv.logit(Phi.mc[[n]])^intervals[n]) * N.mc[[n]])
			Recrut.mc[[n]][Recrut.mc[[n]]<0]<-0
			
			###[[ATENÇÃO! Essa linha transforma o recrutamento em recrutamento mensal
			Recrut.mc[[n]]<-Recrut.mc[[n]]/intervals[n] 
			
				Recrut.stats$mean[n]<-quantile(Recrut.mc[[n]],0.5)
				Recrut.stats$se[n]<-sd(Recrut.mc[[n]])
				Recrut.stats$lcl[n]<-quantile(Recrut.mc[[n]],0.025)
				Recrut.stats$ucl[n]<-quantile(Recrut.mc[[n]],0.975)
				Recrut.stats$time[n]<-force.numeric(cons.N.sizes$estimates$time[n])
			}
			Recrut.estm<-merge(as.data.frame(Recrut.stats),interval_cov,by="time")
			Fin.rate.inc.estm<-merge(as.data.frame(Fin.rate.inc.stats),interval_cov,by="time")
		###Tendência anual
			Fin.rate.overall<-N.mc[[nrow(cons.N.sizes$estimates)]]/N.mc[[1]]
			Fin.rate.ov.mean<-mean(Fin.rate.overall)
			Fin.rate.ov.se<-sd(Fin.rate.overall)
			Fin.rate.ov.lcl<-quantile(Fin.rate.overall,0.025)
			Fin.rate.ov.ucl<-quantile(Fin.rate.overall,0.975)
			Trend.mc<-apply(simplify2array(Fin.rate.inc.mc),1,prod)^(1/(max(interval_cov$cum.time_elaps.days,na.rm=TRUE)/365))
			Trend.stats$mean<-mean(Trend.mc)
			Trend.stats$se<-sd(Trend.mc)
			Trend.stats$lcl<-quantile(Trend.mc,c(0.025))
			Trend.stats$ucl<-quantile(Trend.mc,c(0.975))
			Trend.estm<-as.data.frame(Trend.stats)		