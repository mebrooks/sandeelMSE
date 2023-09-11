########################################################################
stock.name="SAN_area_4_C"
big_output=TRUE

mse.stock.dir = paste0("~/Sandeel_MSE/MSE_code/", stock.name)
code.path = "~/Sandeel_MSE/smsR"
#mse.stock.dir = paste0("~/Desktop/DTU_postdoc/Sandeel_MSE/MSE_code/", stock.name)
#code.path = "~/Desktop/DTU_postdoc/Sandeel_MSE/smsR"

load(file.path(mse.stock.dir, "organized data.Rdata"))
source(file.path(mse.stock.dir, "global_stock_settings.R"))

library(smsR)
require(MASS)
library(parallel)
library(reshape)
library(plyr)
source("utilities.R")
library(TMB)
library(abind)
library(ggplot2); theme_set(theme_bw())
# Load the TMB model
TMB::compile(file.path(code.path, "src/smsR.cpp"))
dyn.load(file.path(code.path, dynlib("src/smsR")))

########################################################################
# INITIALIZATION SECTION
########################################################################
#size of MSE
nyears=20 #years to simulate forward
ntrials=1000
Fcap=seq(0.05, 0.15, by=0.01)
TACmax= 1000000
TACmin= c(0, 5000)

no_obs_error = FALSE
no_init_error= FALSE
if(no_obs_error) Sigma_survey[]=0
if(no_init_error) Sigma[]=0

grid=expand.grid(trial=1:ntrials, Fcap=Fcap, TACmax=TACmax, TACmin=TACmin, 
								 Fhist = Fhist)

#simulate M, PM, weca, and west for all trials to save time
set.seed(111)
sampyr=replicate(sample(minyear:startyear, ntrials, replace=TRUE), n = nyears) #ntrials x nyears

#set up constants
last.years.states = distmat$value
last.years.state.cov=Sigma

#scale old effort by Fbar to get seasonal pattern indept of Fbar
oldeff = df.tmb$effort/apply(apply(cond$F0[2:3,,],MARGIN=2:3, FUN=sum), MARGIN=1, FUN=mean)
oldeff=oldeff*is.finite(oldeff)
eff = matrix(colMeans(oldeff[df.tmb$bidx==max(df.tmb$bidx),], na.rm = TRUE), #assume seasonal effort is mean of past years
						 byrow = TRUE, ncol=nseason, nrow = nages)

#prep data for stock assessment
df.assess=df.tmb

#Allocate storage
Fscaling = array(dim=c(nyears), dimnames=list(year=startyear+(1:nyears)))
TAC = array(dim=c(nyears), dimnames=list(year=startyear+(1:nyears))) #Allowable Catch
Yield = array(dim=c(nyears), dimnames=list(year=startyear+(1:nyears))) #Taken Catch
bb = rep(-0.1, nyears); bb[!(1:nyears %% 2)] = bb[!(1:nyears %% 2)]* -1; bb=c(0, bb) #banking and borrowing years
Ntrue = array(dim=c(length(ages), nseason, nyears), dimnames=list(age=ages, q=1:nseason, year=startyear+(1:nyears)))
Ftrue = array(dim=c(length(ages), nseason, nyears), dimnames=list(age=ages, q=1:nseason, year=startyear+(1:nyears)))
Nest = array(dim=c(length(ages), 1, nyears), dimnames=list(age=ages, q=1, year=startyear+(1:nyears))) 
Fest = array(dim=c(length(ages), nseason, nyears), dimnames=list(age=ages, q=1:nseason, year=startyear+(1:nyears)))
SSBtrue = array(dim=c(nseason,nyears), dimnames=list(q=1:nseason, year=startyear+(1:nyears)))
BioParstrue = array(dim=c(4, nyears, length(ages), nseason), dimnames=list(par=c("M", "PM", "weca","west"), year=startyear+(1:nyears), age=ages, q=1:nseason))

########################################################################
#CALCULATION SECTION
########################################################################
cl=makeCluster(min(c(20, nrow(grid))), type="FORK")
triallist = parApply(cl, grid, 1, function(grid){
	trial=grid["trial"]
	Fcap=grid["Fcap"]
	Fhist =grid["Fhist"]
	TACmax=grid["TACmax"]
	TACmin=grid["TACmin"]

	set.seed(20+trial)
	
	#simulate multivariate log-normal N, E, R, Q for this trial (only use N in first year)
	true.states = mvrnorm(1, mu=last.years.states, Sigma=Sigma)
	#pick out starting stock numbers
	Ntrue.next = exp(matrix(true.states[Nidx], nrow= nages))
	
	#pick out exploitation pattern and resacale it
	E.true = exp(true.states[Fidx])
  #rescale exploitation to make Fscaling = Fbar
  Ftrue0=matrix(E.true, nrow= nages)*eff #only temporary to get Etrue, overwritten later with Fscaling included
  Ftemp=sum(Ftrue0[2:3,])/2 #ages 1 and 2
  Etrue=Ftrue0/Ftemp
  
  #pick out past recruitment
  R.true = exp(true.states[Ridx]) #keep as vector
  
  #pick out catchability and reshape it
  
  Q.mle = getCatchability(df.tmb, sas)
  Q.true = Q.mle #for structure
  Q.true$Q = exp(true.states[Qidx])
  Q.true = cast(Q.true, age~survey, value="Q", drop=FALSE)[,-1, drop=FALSE]
  
	#simulate M, PM, weca, and west	for this trial
	BioParstrue["M",,,]=M_array[as.character(sampyr[trial,]), , ]
	BioParstrue["PM",,,]=PM_array[as.character(sampyr[trial,]), , ]
	BioParstrue["weca",,,]=weca_array[as.character(sampyr[trial,]), , ]
	BioParstrue["west",,,]=west_array[as.character(sampyr[trial,]), , ]

	for(yr in 1:nyears) {
		
		#MANAGEMENT PROCEDURE
		#do a stock assessment
		x=try({
			# Set structure and initial values of parameters
			inits = getParms(df.assess)
			# Get non-estimated parameters, based on info in df.assess
			mps = getMPS(df.assess, inits)
			smsmod = runAssessment(df.assess, parms = inits, mps = mps)
			sdr = smsmod$reps
		})
		if(inherits(x, "try-error")) break()
			
		#save the estimated state on the real scale
		Nest[,,yr] = exp(matrix(sdr$value[names(sdr$value)=="term_logN_next"],nrow=nages))
		Fest[,,yr] = exp(matrix(sdr$value[names(sdr$value)=="log_exp_pattern"],nrow=nages))*eff

		#the forecast uses geom mean recruitment in the period from 
		# f.rec years before the assessment year to l.rec years before the assessment year
		Rec = data.frame(estR = getR(df.assess, smsmod)$R[1:df.assess$nyears], 
										 year=df.assess$years,
										 estSSB = getSSB(df.assess, smsmod)$SSB[1:df.assess$nyears])
		ayear=max(df.assess$years)#assessment year
		R.for = exp(mean(log(subset(Rec, year >= (ayear-f.rec()+1) & year <=(ayear-l.rec))$estR)))
		
		#organize bio pars for forecast according to stock annex vals in global_stock_settings.R 
		nyia=df.assess$nyears #(n)umber (y)ears (i)n (a)ssessment
		M.for=apply(df.assess$M[,(nyia-f.m):(nyia-l.m),, drop=FALSE], MARGIN = c(1,3), mean)
		PM.for=apply(df.assess$Mat[,(nyia-f.pm):(nyia-l.pm),, drop=FALSE], MARGIN = c(1,3), mean)
		west.for=apply(df.assess$west[,(nyia-f.west):(nyia-l.west),, drop=FALSE], MARGIN = c(1,3), mean)
		weca.for=apply(df.assess$weca[,(nyia-f.weca):(nyia-l.weca),, drop=FALSE], MARGIN = c(1,3), mean)

		#rescale estimated exploitation to make Fscaling = Fbar
		Ftemp = sum(Fest[2:3,,yr])/2 #F values in matrix are half-yearly rates, Fbar is annual.
		Fest[,,yr]= Fest[,,yr]/Ftemp

		#do a forecast to suggest a TAC based on the estimated state and Fcap
		Nestmat = cbind(Nest[,,yr], rep(0, nages))#assumes nseasons=2
		Forecast = do.forecast(N=Nestmat, E=Fest[,,yr], Fcap=Fcap, M=M.for, PM=PM.for, west=west.for, weca=weca.for, Bmsy=Bmsy, R=R.for)

		#HCR
		TACsuggested = Forecast$TAC
		if(TACsuggested>TACmax) TACsuggested=TACmax
		if(TACsuggested<TACmin) TACsuggested=TACmin
		
		#store the allowable catch		
		TAC[yr] = TACsuggested

		#Implement Banking and borrowing
		if(yr==1){ Yield[yr] = TAC[yr] } # bb[yr]=0 
		if(yr>1){	Yield[yr] = TAC[yr]*(1 + bb[yr]) - bb[yr-1]*TAC[yr-1]}
		if(Yield[yr]<0){ Yield[yr]=0}
		
		#check what Fbar corresponds to the yield given the estimated state
		# and seasonal effort
		Forecast$Fbar = calc.fterm(target.TAC= Yield[yr], E=Fest[,,yr], N=Nestmat, 
															 M=M.for, PM=PM.for, west=west.for, weca=weca.for, Fcap= Fhist, R=R.for)$minimum
		
		
		#save what the MP expected
		Fest[,,yr]=Fest[,,yr]*Forecast$Fbar

		
		# OPERATING MODEL
		#set up bio pars to use for this year's dynamics
		M = BioParstrue["M",yr,,]
		PM = BioParstrue["PM",yr,,]
		weca = BioParstrue["weca",yr,,]
		west = BioParstrue["west",yr,,]

		#simulate true recruitment from season 1 SSB
		Ntrue1=cbind(Ntrue.next, #N in q1 for this TAC year
				rep(0, nages))#assumes only 2 seasons
		SSB1=sum((Ntrue1*west*PM)[,1])
		Rtrue = sim.recruitment(SSB=SSB1, R=R.true[years%in%years_Rec_above_Blim], hockey.break=Blim)

		#check what Fbar corresponds to the suggested TAC given the true state
		# and seasonal effort
		Fbar = calc.fterm(target.TAC= Yield[yr], E=Etrue, N=Ntrue1, 
											weca=weca, west=west, PM=PM, M=M, Fcap= Fhist, R=Rtrue)$minimum
		
		Fscaling[yr]=Fbar

		#update the true stock numbers
		tmp = update.stock(N=Ntrue1, M=M, f=Fbar*Etrue, PM=PM, west=west, R=Rtrue)
		Ntrue[,,yr]=tmp$N #this year's dynamics
		Ntrue.next=tmp$N_next[,1] #N in q1 of next year
		Ctrue=tmp$C
		
		#save the true state on the real scale
		Ftrue[,,yr]=Fbar*Etrue

		#save SSB for performance statistics 
		SSBtrue[,yr] =apply(PM * west * Ntrue[,,yr], 2, sum)

		# OBSERVATION SIMULATOR
		#simulate catch and survey data
		Csim=sim.catch(C=Ctrue, N=Ntrue[,,yr], Q = Q.true, Sigma_survey=Sigma_survey, Z=M+Ftrue[,,yr])
		if(stock.name=="SAN_area_4_C") Csim$Sobs[,1]=-1 #stopped doing that dredge survey
		
		Up = update.in.files(df=df.assess, pars=inits, Csim=Csim, weca=weca, west=west, 
												 PM=PM, M=M, effort=Fbar*eff[1,])
		
		df.assess = Up$df

	} #end year for-loop
	########################################################################
	# OUTPUT SECTION
	########################################################################

	#change structures to be easier to use later 
	BioParstrue = recast(melt(BioParstrue), age+ q+ year ~par, id.var=1:4)
	
	Fscaling = melt(Fscaling) #by TACyr
	names(Fscaling)[grep("value", names(Fscaling))]="Fbar"

	TAC = melt(TAC) #by TACyr
	names(TAC)[grep("value", names(TAC))]="TAC"
	
	Yield = melt(Yield) #by TACyr
	names(Yield)[grep("value", names(Yield))]="Yield"
	
	Ntrue = melt(Ntrue) #by age, q, TACyr
	names(Ntrue)[grep("value", names(Ntrue))]="Ntrue"
	
	Ftrue = melt(Ftrue) #by age, q, TACyr
	names(Ftrue)[grep("value", names(Ftrue))]="Ftrue"

	Nest = melt(Nest) #by age, q (1), TACyr
	names(Nest)[grep("value", names(Nest))]="Nest"

	Fest = melt(Fest) #by age, q, TACyr
	names(Fest)[grep("value", names(Fest))]="Fest"

	SSBtrue = melt(SSBtrue) #by q, TACyr
	names(SSBtrue)[grep("value", names(SSBtrue))]="SSBtrue"

	
	if(big_output){
		dat=join(join(join(join(join(join(join(join(Ntrue, Ftrue), Nest), Fest), BioParstrue), Fscaling), TAC), Yield), SSBtrue)
	}else{
		dat=join(join(join(join(join(Ntrue,  BioParstrue), Fscaling), TAC), Yield), SSBtrue)
	}	
	return(dat)
})#end parApply on grid
stopCluster(cl)

dat= combine.trials(triallist, grid)

if(!big_output) save(dat, grid, sampyr, file=paste0("MSE_varNERQ_BnB_", stock.name, "_small.Rdata"))
if(big_output) save(dat, grid, sampyr, file=paste0("MSE_varNERQ_BnB_", stock.name, ".Rdata"))
