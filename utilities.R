
#R is a vector of past recruitment estimates on the real scale when MLE SSB was >hockey.break
sim.recruitment =function(SSB, R, hockey.break=110000, seed=NULL){
	
	if(! is.null(seed)){set.seed(seed)}
	
		b = mean(R)
		var = var(R)
		n = length(SSB)
		
		#Calculate expected recruitment on real scale
		m = ifelse(SSB>hockey.break, b, SSB*b/hockey.break)
		
		#variance is specified on the natural scale
		mu = log(m/sqrt(1+var/(m^2))) #log-scale mean
		sig = sqrt(log(1+var/(m^2))) #log-scale std dev
		z = rnorm(n)

		x = mu +z *sig
		Rec = exp(x)
		
		return(Rec)
}

##' Get the correlation matrix for "term_N_next" and "exploi_pattern"
##' @export
getcovNERQ=function(sdr, years, ages, nseason) {
	#Covariance matrix
	cov0 = sdr$cov
	
	ind = c(which(names(sdr$value) == "term_logN_next"),
					which(names(sdr$value) == "log_exp_pattern"),
					which(names(sdr$value) == "logRec"),
					which(names(sdr$value) == "logQsurv"))

	covmat = cov0[ind,ind]
	
	colnames(covmat) = names(sdr$value)[ind]
	rownames(covmat) = names(sdr$value)[ind]
	
	return(covmat)
}

##' Get the coefficients of variation for "term_logN_next", "log_exp_pattern", "logRec", "logQsurv"
##' @export
getdistNERQ=function(sdr, years, ages, nseason) {
	par0 = sdr$value
	
	ind = c(which(names(par0) == "term_logN_next"),
					which(names(par0) == "log_exp_pattern"),
					which(names(par0) == "logRec"),
					which(names(par0) == "logQsurv"))
	
	
	par1 = par0[ind]
	mat = data.frame("name"=names(par1), "value"= par1, "std.dev"= sdr$sd[ind])
	mat$value[mat$value==0] = -1000
	mat = transform(mat, cv = abs(std.dev/value))
	
	#	info = info0[ind, c("par", "year", "quarter", "age")]
	#	mat2=cbind(mat, info[,-1])
	mat$year=c(rep(max(years)+1, length(ages)), #term_logN_next
						 rep(max(years), length(ages)*nseason),
						 years,
						 rep(max(years), nrow(getCatchability(df.tmb, sas)))
						 )
	
	mat[is.na(mat)]=0
	mat[mat==-1]=NA
	return(mat)
}

#Read refence points from file reference_points.in 
Read.reference.points=function(data.in.path){
    a=scan(file.path(data.in.path, "reference_points.in"),comment.char = "#",quiet = TRUE)
    b=matrix(a,nrow=1,ncol=4,byrow=TRUE)
    colnames(b)=c("Flim","Fpa","Blim","Bpa")   
    b
}


# move the N forward to the end of the year while taking the catch
move.on<-function(N,M,f,R){
    Z=M+f
    C=Z   #copy structure
    #put recruitment in the right place
    N[f.age, rec.season]=R
    for (q in (f.season:l.season)) {
      for (a in (f.age:l.age)){
        if (Z[a,q]>0) C[a,q]=f[a,q]*N[a,q]*(1-exp(-Z[a,q]))/Z[a,q]
        if (q<l.season) { 
        	if (!(q<rec.season & a==f.age)) N[a,q+1]=N[a,q]*exp(-Z[a,q]) 
        }
      }
    }
   list("N"=N,"C"=C)
}

#move the N forward from q4 this year to q1 of next year 
birthday=function(N, Z, PM, west) {
  N[,1]=0
  for (a in (l.age:(f.age+1))){
   if (a==l.age && plus.group==1) N[a,f.season]=N[a-1,l.season]*exp(-Z[a-1,l.season])+N[a,l.season]*exp(-Z[a,l.season])
   if (a<l.age || plus.group==0) N[a,f.season]=N[a-1,l.season]*exp(-Z[a-1,l.season])
  }
  N[,(f.season+1):l.season]=0  
  SSB0=sum(N*PM*west)
  list("N"=N,"SSB"=SSB0)
}

#predict the escaped SSB and catch from a given Fbar
do.prediction=function(Fbar, E, N, weca, PM, M, west, R){   
    f=E*Fbar
    res=move.on(N=N, M=M, f=f, R=R)
    N=res$N
    C=res$C
    TAC=sum(C*weca, rm.na=T)
    # longterm data 
    res=birthday(N=N, Z=M+f, PM=PM, west=west)
    list("SSB"=res$SSB, "TAC"=TAC)
}

#
calc.fterm=function(target.SSB=NULL, target.TAC=NULL, E, Fcap, N, weca, west, PM, M, R) {
	if(is.null(target.TAC)) {
		mini = function(x, ...) {
		    a=do.prediction(Fbar=x, E=E, N=N, weca=weca, west=west, PM=PM, M=M, R=R) 
		    (a$SSB-target.SSB)^2
		  }
	} else {
		mini = function(x, ...) {
		    a=do.prediction(Fbar=x, E=E, N=N, weca=weca, west=west, PM=PM, M=M, R=R) 
		    (a$TAC-target.TAC)^2
		  }
	}
  res = optimize(f=mini, interval=c(0, Fcap), E=E, N=N, weca=weca, west=west, PM=PM, M=M, R=R)
  res
}

do.forecast =function(N, R, E, M, PM, west, weca, Bmsy, Fcap) {
	if(round(Fcap, 2)==0){
		Fbar=0
	} else{ 
		# MSY, find F 
		res=calc.fterm(target.SSB=Bmsy, E=E, Fcap= Fcap, N=N, weca=weca, west=west, PM=PM, M=M, R=R)
		Fbar=res$minimum
		obj=res$objective
	}
	tmp=do.prediction(Fbar, E=E, N=N, weca=weca, west=west, PM=PM, M=M, R=R)
	SSB1=tmp$SSB
	TAC=tmp$TAC
	
	return(list("TAC"=TAC, "Fbar"=Fbar, "SSB"=SSB1))
}

#Input: true stock in q1, f required to take the TAC, PM, 
#Output true stock in all q of this year and q1 of next year
update.stock=function(N, M, f, PM, west, R) {

	#move from q1 to q2
	tmp1= move.on(N, M, f, R)
	N1=tmp1$N
	C=tmp1$C
	
	#birthday from q4 to q1
	tmp2 = birthday(N=N1, Z=M+f, PM=PM, west=west)
	N2=tmp2$N

	return(list("N"=N1, "N_next"=N2, "C"=C))
}


assessment.emulator = function(state.true, state.cov) {
	mvrnorm(1, mu= state.true, Sigma= state.cov, 1e-5)
}


## generate a list with names equal to values
named.list  =  function (...) {
    L  =  list(...)
    snm  =  sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm  =  names(L)))
        nm  =  snm
    if (any(nonames  =  nm == ""))
        nm[nonames]  =  snm[nonames]
    setNames(L, nm)
}

combine.trials=function(triallist, grid) {
	do.call(rbind, lapply(1:length(triallist), 
	function(i) { 
		data.frame(triallist[[i]], grid[i,]) 
	}))
}

#N0 is the stock number at the beginning of the season.
#t is how far into the season to extend the average (<=1)
#Z is total mortality rate over the season
calc.Nbar=function(N0, t, Z)
{
	if(t==0) return(N0)
	#else
	newN = N0/Z*(1-exp(-Z*t))/t
	indx = which(Z==0)
	newN[indx]=N0[indx]
	return(newN)
}
	
sim.catch=function(C, N, Q, Sigma_survey, Z, est)
{
	
	S0 = df.tmb$Surveyobs[,1,, drop=FALSE] #copy shape, fill in pieces, then abind
	Snew = S0
	
	# Expected surveys
	for(i in 1:df.tmb$nsurvey)
	{
		s = df.tmb$surveySeason[i]
		
		#calc N at start of survey
		if(df.tmb$surveyStart[i]==0) N0=N[, s]
		if(df.tmb$surveyStart[i]>0) N0=N[, s]*exp(-Z[, s]*df.tmb$surveyStart[i])
		
		duration = df.tmb$surveyEnd[i]-df.tmb$surveyStart[i]
		
		Nbar = calc.Nbar(N0, duration, Z[, s])

		catchability =	Q[,i] #cond$Qsurv[,i]
		Snew[,1,i]=Nbar*catchability
	}

	#Add error to expected value
	sim.resids = mvrnorm(1, mu=rep(0, ncol(Sigma_survey)), Sigma=Sigma_survey, tol=1e-1)
	sim.resids = data.frame(agesurveys=names(sim.resids), resid=sim.resids)
	sim.resids = join(data.frame(sim.resids), Sigma_survey_struct)	
	
	Cdf=melt(C)
	Cdf$survey="Catch"
	names(Cdf)[2]="s"
	
	Sdf=melt(Snew)
	
	sim.resids.survey=subset(sim.resids, survey!="Catch")
	sim.resids.catch=subset(sim.resids, survey=="Catch")
	
	Cdf=join(Cdf, sim.resids.catch)
	Sdf=join(Sdf, sim.resids.survey)
	
	both=rbind(Cdf, Sdf[,names(Cdf)])
	both=transform(both, obs = exp(log(value)+resid))
	
	Sobs= cast(subset(both, survey!="Catch")[,c("age", "survey", "obs")], age~survey, value="obs")
	Cobs= cast(subset(both, survey=="Catch")[,c("age", "s", "obs")], age~s, value="obs")[,-1]
	
	#put surveys in the same order as before
	Sobs=Sobs[,dimnames(Snew)[[3]], drop=FALSE]

	Cobs[is.na(Cobs)]=0
	Sobs[is.na(Sobs)]=-1
	
	return(list(Cobs=Cobs, Sobs=Sobs))
}

update.in.files=function(df, pars, Csim, weca, west, PM, M, effort) {
		
	df$weca = abind(abind(df$weca[,1:df$nyears,], weca, along = 2), weca, along = 2)
	df$west = abind(abind(df$west[,1:df$nyears,], west, along = 2), west, along = 2)
	df$M = abind(abind(df$M[,1:df$nyears,], M, along = 2), M, along = 2)
	df$Mat = abind(abind(df$Mat[,1:df$nyears,], PM, along = 2), PM, along = 2)
	
	df$Catchobs = abind(df$Catchobs, Csim$Cobs, along=2)
#	df$Surveyobs = abind(df$Surveyobs, Csim$Sobs, along=2)
	df$Surveyobs = abind(df$Surveyobs, array(unname(Csim$Sobs), dim=c(dim(df$Surveyobs)[1],1,dim(df$Surveyobs)[3])), along = 2)

	df$nocatch = rbind(df$nocatch, colSums(Csim$Cobs)>100000)
#	df$nocatch = apply(df$Catchobs, MARGIN=2:3, FUN=function(x){as.numeric(sum(x)>100000)})
	df$effort = rbind(df$effort, effort)
	df$bidx = c(df$bidx, tail(df$bidx, 1))
	df$scv = abind(df$scv, df$scv[,df$nyears-1,, drop=FALSE], along=2)
	df$propM = abind(df$propM, df$propM[,df$nyears-1,, drop=FALSE], along=2)
	df$propF = abind(df$propF, df$propF[,df$nyears-1,, drop=FALSE], along=2)
	df$no = df$no + matrix(data=c(0,0, 2,2,  2,2), byrow=TRUE, nrow=3, ncol=2) #hack that only works here
	
	df$nyears = df$nyears+1
	df$years = c(df$years, max(df$years)+1)
	
  return(list(df=df, pars=pars))
}

#Seasonal data on biological parameters for resampling in future years
bioparfun = function(bparray){
	dimnames(bparray)=list(age=ages, year=c(years, max(years)+1), s=1:nseason)
	cast(subset(melt(bparray), year>=minyear &year<=startyear), year~age~s)
}
