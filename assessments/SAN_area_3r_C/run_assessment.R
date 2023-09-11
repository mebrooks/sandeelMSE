# Run sprat from R
#devtools::install_github('https://github.com/nissandjac/smsR', dependencies = TRUE)
library(smsR)

maxage <- 4
years = 1986:2021
nyear <- length(years)
seasons <- 1:2
#
dat <- getDataSMS(wd,
                      maxage = maxage,
                      survey.age = list(0:1, 1:4), # Ages in the two surveys
                      survey.years = list(2006:2021, 2009:2021),
                      survey.names = c('Dredge','Acoustic'),
                      survey.quarter = c(2,1), # ? not sure if this is correct
                      years = years,
                      seasons = seasons

)
#
# #
# dat <- getDataSandeel(wd,
#                       maxage = maxage,
#                       survey.age = list(0:1), # Ages in the two surveys
#                       survey.years = list(2006:2021),
#                       survey.names = c('Dredge'),
#                       survey.quarter = c(2), # ? not sure if this is correct
#                       years = years,
#                       seasons = seasons
#
# )


ages <- 0:maxage
nseason <- 2 # Number of seasons
beta <- 80000

# Fishing effort
# Get the effort data
Feffort <- matrix(scan(file.path(wd,'effort.in'), comment.char = '#'),
                  ncol = 2, nrow = nyear)

# Normalize effort to 1
# Format input data matrices into TMB style
# mtrx <- sumtable_to_matrix(sandeel.age)
Surveyobs <- survey_to_matrix(dat[['survey']], years)
Catchobs <- df_to_matrix(dat[['canum']], season =  1:2)
#Catchobs[,]

#Catchobs[,years == 2013,] <- 0
nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#', skip = 3)

#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
#nocatch[29,2] <- 0 # Test  this

df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  startYear = 1987,
  nseason = nseason, # Number of seasons
  useEffort = 1,
  ages = ages, # Ages of the species
  recseason = 2, # Season where recruitment occurs
  Fmaxage = 3, # Fully selected fishing mortality age
  Fbarage = c(1,2),
  CminageSeason = c(1,1),
  Qminage = c(0,1), # minimum age in surveys
  Qmaxage = c(1,4),
  Qlastage = c(1,3),
  isFseason = c(1,0), # Seasons to calculate fishing in
  effort = Feffort,
  blocks = FALSE, #c(1986,1999),
  powers = list(NA,NA),
  #leavesurveyout=c(1,1),
  endFseason = 2, # which season does fishing stop in the final year of data
  nocatch = as.matrix(nocatch),
  surveyStart = c(0.75,.5),
  surveyEnd =  c(1,.7), # Does the survey last throughout the http://127.0.0.1:35111/graphics/plot_zoom_png?width=1200&height=900season it's conducted?
  surveySeason = c(2,1), # Which seasons do the surveys occur in
  surveyCV =  list(c(0,1),
                   c(1)),
  catchCV = list(c(0,1,3),
                 c(0,1,3)),
  estCV = c(0,2,0),
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,0.05) # Factor for relative strength of log-likelihood

)

# Get initial parameter structure
parms <- getParms(df.tmb)

# Get non-estimated parameters, based on info in df.tmb

#mps <- getMPS(df.tmb, parms, mapExtra = 'SDsurvey')

mps <- getMPS(df.tmb, parms)
sas <- runAssessment(df.tmb, parms = parms, mps = mps)
# save(df.tmb, sas, file="assessment_objects.Rdata")

save(df.tmb, sas, file=file.path(wd, "assessment_objects.Rdata"))
