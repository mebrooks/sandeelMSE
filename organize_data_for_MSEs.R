library(smsR)
library(plyr)
library(reshape)
library(MASS)
library(ggplot2); theme_set(theme_bw())

source("utilities.R")

stocks = c("SAN_area_1r_C", 
					 "SAN_area_2r_C", 
					 "SAN_area_3r_C", 
					 "SAN_area_4_C")


for(stock.name in stocks){
	wd = file.path("assessments", stock.name)
	output.path=wd
	
	load(file.path(wd, "assessment_objects.Rdata"))#sas and df.tmb
	
	
	# fix how years vary by stock
	first_year_rec = min(df.tmb$years) #first year to include in stock-recruit info
	startyear  = max(df.tmb$years)
	minyear = startyear-10+1 #first year of biological parameters to use
	
	###################################################################
	# Reorganize pieces from the assessment for conditioning
	
	#dimensions of some objects
	ages=df.tmb$age
	years=df.tmb$years
	nseason = df.tmb$nseason
	nages=length(ages)
	
	#reorganize bio pars into arrays with dim= year, age, season
	M_array=bioparfun(df.tmb$M)
	west_array=bioparfun(df.tmb$west)
	weca_array=bioparfun(df.tmb$weca)
	PM_array=bioparfun(df.tmb$Mat)
	M_mean=cast(melt(M_array), age~s, fun.aggregate=mean)
	
	
	#####################
	distmat0=getdistNERQ(sas$reps, years, ages, nseason)
	covmat0=getcovNERQ(sas$reps, years, ages, nseason)
	
	keep=distmat0$year>= first_year_rec #drop est rec before first_year_rec
	distmat=distmat0[keep,]
	Sigma=covmat0[keep,keep]
	
	Fidx=which(distmat$name == "log_exp_pattern")
	Nidx=which(distmat$name == "term_logN_next")
	Ridx=which(distmat$name == "logRec")
	Qidx=which(distmat$name == "logQsurv")
	
	tmp=mvrnorm(1, mu=distmat[,"value"], Sigma= Sigma, tol=1e-5) #just a test
	
	
	tmp=Read.reference.points(wd)
	Bmsy=tmp[1,'Bpa']
	Blim=tmp[1,'Blim']
	
	#Plot Recruitment
	SSB <- exp(sas$reps$value[names(sas$reps$value) == 'logSSB'])
	rec <- sas$reps$value[names(sas$reps$value) == 'Rsave']
	FF <- exp(array(sas$reps$value[names(sas$reps$value) == 'logF0'],
									dim = c(df.tmb$nage, df.tmb$nyears,df.tmb$nseason),
									dimnames=list(age=ages, year=years, s=1:nseason)))
	
	FF[FF==1]=0
	FFm=subset(melt(FF), age %in% 1:2)
	Fbar=ddply(FFm, ~year, summarize, F=sum(value)/2)
	
	ydat=data.frame(year=years, SSB=SSB[-length(SSB)], Rec=rec,  Fbar=Fbar$F)
	
	#Recruitment residuals to be used for simulating recruitment in OM
	Rec_above_Blim=subset(ydat, SSB>Blim & year %in% first_year_rec:min(years)-1)$Rec
	years_Rec_above_Blim = subset(ydat, SSB>Blim & year %in% first_year_rec:min(years)-1)$year
	minRec = min(ydat$Rec, na.rm=TRUE)
	
	ggplot(subset(ydat, year %in% 2000:startyear), aes(SSB, Rec))+
		geom_text(aes(label=year)) +
		geom_vline(xintercept=Blim, color="red")+
		ylab("Recruitment")+
		ggtitle(stock.name)
	
	ggsave(file.path(output.path, paste0("recruitment vs SSB 2000 to ",startyear, ".png")), height=4, width=4)
	
	ggplot(subset(ydat, year %in% first_year_rec:startyear), aes(SSB, Rec))+
		geom_text(aes(label=year)) +
		geom_vline(xintercept=Blim, color="red")+
		ylab("Recruitment")+
		ggtitle(stock.name)
	
	ggsave(file.path(output.path,paste0("recruitment vs SSB ", first_year_rec," to ", startyear, ".png")), height=4, width=4)
	
	save(ydat, file=file.path(output.path,"yearly.Rdata"))
	
	
	#Historical max F
	Fhist=max(ydat$Fbar, na.rm=TRUE)
	
	
	#########TMB resids
	x = sas$obj$report(sas$obj$env$last.par)
	cond = x #values to use for CONDitioning the operating model
	
	#survey residuals
	resids1=x$resid_survey #(dimensions: age, year, survey)
	dimnames(resids1)=dimnames(df.tmb$Surveyobs)
	resids1=melt(resids1)
	resids1$s = df.tmb$surveySeason[resids1$survey] #survey1 in s1, survey2 in s2
	
	#catch residuals
	resids2=x$resid_catch
	dimnames(resids2)=dimnames(FF)
	resids2=melt(resids2)
	resids2$survey="Catch"
	
	tmp0=subset(rbind(resids1[,c("age","year","s","survey","value")], resids2[,c("age","year","s","survey","value")]), value>-99)
	names(tmp0)[5]="residual"
	tmp=transform(tmp0,
								agesurveys=paste0("age", age, "_", survey, "_s", s))
	tmpwide=reshape(tmp[, c("agesurveys", "residual","year")] , direction="wide", timevar ="agesurveys", v.names="residual", idvar=c("year"))
	#colnames(dwide)[1]="year"
	#tmpwide=tmpwide[,colnames(dwide)]
	#Sigma_survey2 =cov(tmpwide[, -1], use="pairwise.complete.obs")
	Sigma_survey =cov(tmpwide[, -1], use="na.or.complete")#"pairwise.complete.obs")#use=
	
	#Save structure to help reorganize simulated residuals later
	Sigma_survey_struct=unique(tmp[,c("age", "survey", "s", "agesurveys")])
	
	if(stock.name %in%c("SAN_area_3r_C", "SAN_area_4_C"))
	{
		Sigma_survey =cov(tmpwide[, -1], use="pairwise.complete.obs")#
		
		Sigma_survey[!outer(Sigma_survey_struct$survey, Sigma_survey_struct$survey, "==")	]=0
	}
	test=mvrnorm(1, mu=rep(0, ncol(Sigma_survey)), Sigma= Sigma_survey, tol=1e-1) #just a test
	
	rownames(Sigma_survey)= gsub("residual.", "", rownames(Sigma_survey))
	colnames(Sigma_survey)= gsub("residual.", "", colnames(Sigma_survey))
	
	if(!all(rownames(Sigma_survey) %in% Sigma_survey_struct$agesurveys)) warning("Problem with Sigma_survey!")
	
	logFavg=data.frame(m=sas$reps$value[names(sas$reps$value)=="logFavg"],
										 sd=sas$reps$sd[names(sas$reps$value)=="logFavg"],
										 year=c(df.tmb$years))
	
	logRec=data.frame(m=sas$reps$value[names(sas$reps$value)=="logRec"],
										sd=sas$reps$sd[names(sas$reps$value)=="logRec"],
										year=df.tmb$years)
	
	logCatch=data.frame(m=sas$reps$value[names(sas$reps$value)=="logCatchtot"],
											sd=sas$reps$sd[names(sas$reps$value)=="logCatchtot"],
											year=df.tmb$years)
	
	
	logSSB=data.frame(m=sas$reps$value[names(sas$reps$value)=="logSSB"],
										sd=sas$reps$sd[names(sas$reps$value)=="logSSB"],
										year=c(df.tmb$years, max(df.tmb$years)+1))
	
	
	#Save results
	save(distmat, Sigma, Bmsy, Blim, Rec_above_Blim, years_Rec_above_Blim, 
			 # qdat, hdat, nages,
			 #	N_array , E_array, dist_FSSB,
			 M_array, M_mean, west_array, weca_array, PM_array,
			 Fidx, Nidx, Ridx, Qidx, minyear, startyear,
			 Sigma_survey, Sigma_survey_struct, Fhist,
			 df.tmb, cond, sas,
			 ages, years, nseason, nages,
			 logCatch, logRec, logSSB, logFavg,
			 file=file.path(output.path, "organized_data.Rdata"))
	
}
