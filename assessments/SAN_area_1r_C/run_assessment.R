
#devtools::install_github('https://github.com/nissandjac/smsR', dependencies = F)  #or dependencies = T

library(smsR)

maxage <- 4
years = 1983:2021
nyear <- length(years)
seasons <- 1:2

  dat <- getDataSMS(wd,
                        maxage = maxage,
                        survey.age = list(0:1, 1:3), # Ages in the two surveys
                        survey.years = list(2004:2021, 2011:2020),
                        survey.names = c('Dredge','RTM'),
                        survey.quarter = c(2,1),
                        years = years,
                        seasons = seasons

  )


  Qminage = c(0,1) # Qminage = c(0,1) minimum age in surveys
  Qmaxage = c(1,3) #Qmaxage = c(1,3)
  surveyStart = c(0.75,0) #c(0.75,0)
  surveyEnd =  c(1,0) # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = c(2,1) # c(2,1)Which seasons do the surveys occur in
  surveyCV =  list(c(0,1),
                   c(1,2)) #c(1,2)),
  powers = list(NA,NA)  #c(0)
# Load packages and files #


ages <- 0:maxage
nseason <- 2 # Number of seasons
beta <- 110000

# Fishing effort
# Get the effort data
Feffort <- matrix(scan(file.path(wd,'effort.in'), comment.char = '#'),
                  ncol = 2, nrow = nyear)

# Normalize effort to 1
# Format input data matrices into TMB style
# mtrx <- sumtable_to_matrix(sandeel.age)
Surveyobs <- survey_to_matrix(dat[['survey']], years)
Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)

nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#', skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files

# Save some data for package example
dat$effort <- Feffort
dat$nocatch <- nocatch#*0+1
#dat$nocatch[dat$effort == 0] <- 0



df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  nseason = nseason, # Number of seasons
  useEffort = 1,
  ages = ages, # Ages of the species
  recseason = 2, # Season where recruitment occurs
  CminageSeason = c(1,1),
  Fmaxage = 3, # Fully selected fishing mortality age
  Qminage = Qminage, # Qminage = c(0,1) minimum age in surveys
  Qmaxage = Qmaxage, #Qmaxage = c(1,3)
  Fbarage = c(1,2),
  isFseason = c(1,0), # Seasons to calculate fishing in
  effort = Feffort,
  powers = powers,
  blocks = c(1983,1999), # Blocks with unique selectivity c(1983,1989,1999,2005,2010)
  endFseason = 2, # which season does fishing stop in the final year of data
  nocatch = as.matrix(nocatch),
  surveyStart = surveyStart, #c(0.75,0)
  surveyEnd =  surveyEnd, # c(1,0) Does the survey last throughout the season it's conducted?
  surveySeason = surveySeason, # c(2,1)Which seasons do the surveys occur in
  surveyCV =  surveyCV, #c(1,2)),
  catchCV = list(c(0,1,3),
                 c(0,1,3)),
  recmodel = 2, # Chose recruitment model (2 = estimated)
  estCV = c(0,2,0), # Estimate
  beta = 110000, # Hockey stick plateau
  nllfactor = c(1,1,0.05) # Factor for relative strength of log-likelihood

)



# Get initial parameter structure
parms <- getParms(df.tmb)

# Get non-estimated parameters, based on info in df.tmb
mps <- getMPS(df.tmb, parms)

# Set boundaries
# This model works best if SDsurvey is mapped
#parms[names(parms) == 'SDsurvey'][[1]] <- c(0.5, 0.3, 0.3,.3)

# Try initial parameters in weird spots
sas <- runAssessment(df.tmb, parms = parms, mps = mps)
# save(df.tmb, sas, file="assessment_objects.Rdata")
#
save(df.tmb, sas, file=file.path(wd, "assessment_objects.Rdata"))



