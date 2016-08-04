	###Importing Data, Formatting
		source(".\\RDPNE.R")
		
	###Temporal Indepence; Formatting Pt 2
		source(".\\RDPNE.R")
	
	###Generate RDPNE encounterhistory
		source(".\\RDPNE.R")
	
	###RMark Processing Data
		source(".\\RMark Process")
		
####
##Running simple models, all-constant, to have initial values to start more complex ones##
####
	naive.dot.fix.siman<-RMark::mark(naive.agout.proc,cl.agout.ddl,
		options="SIMANNEAL",
		threads=-1)
	dot.fix.init.siman<-RMark::mark(agout.proc,cl.agout.ddl,
		initial=naive.dot.fix.siman,
		options="SIMANNEAL",
		threads=-1)
		
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
	population.closure.test<-collect.models(lx=c("clos.model","open.model"),type="PoissonMR")
		
	####
	##Defining parameter structure to keep: 
	##If the closed strucutre has delta AICc greater than two, 
	##keep the open parameters strucutre.
	##If not, keep the closed parameters strucuture.
	####
	
		if(clos.model$results$AICc - open.model$results$AICc > 2){
			closed.months<-FALSE
		} else {
			closed.months<-TRUE}
			
		if (closed.months){
			true.agout.ddl<-cl.agout.ddl
		} else{
			true.agout.ddl<-op.agout.ddl}