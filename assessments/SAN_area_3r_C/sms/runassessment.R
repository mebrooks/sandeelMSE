# Run sprat from R
#devtools::install_github('https://github.com/nissandjac/smsR', dependencies = TRUE)
library(smsR)
# Working directory with SMS files
wd <- 'M:/Tobis/Tobis_assessment/SMS_2022_Benchmark/SAN_area_3r_Current/sms'
## Load packages and files #

maxage <- 4
years = 1986:2021
nyear <- length(years)
seasons <- 1:2
#
dat <- getDataSandeel(wd,
                      maxage = maxage,
                      survey.age = list(0:1, 1:4), # Ages in the two surveys
                      survey.years = list(2004:2021, 2009:2021),
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

nocatch <- read.table(file.path(wd,'zero_catch_year_season.in'), comment = '#', skip = 3)
#Feffort[which(df.tmb$years == 2013),2] <- 0 # This is different two places in the input files
#nocatch[29,2] <- 0 # Test  this



df.tmb <- get_TMB_parameters(
  mtrx = dat[['mtrx']], # List that contains M, mat, west, weca
  Surveyobs = Surveyobs, # Survey observations (dimensions age, year, quarter, number of surveys)
  Catchobs = Catchobs, # Catch observations  (dimensions age, year, quarter)
  years = years, # Years to run
  nseason = nseason, # Number of seasons
  useEffort = 1,
  ages = ages, # Ages of the species
  recseason = 2, # Season where recruitment occurs
  Fmaxage = 3, # Fully selected fishing mortality age
  Fbarage = c(1,2),
  CminageSeason = c(1,0),
  Qminage = c(0,1), # minimum age in surveys
  Qmaxage = c(1,4),
  Qlastage = c(0,3),
  isFseason = c(1,0), # Seasons to calculate fishing in
  effort = Feffort,
  blocks = c(1986, 1999),
  powers = list(NA,NA),
  endFseason = 2, # which season does fishing stop in the final year of data
  nocatch = as.matrix(nocatch),
  surveyStart = c(0.75,.5),
  surveyEnd =  c(1,.7), # Does the survey last throughout the season it's conducted?
  surveySeason = c(2,1), # Which seasons do the surveys occur in
  surveyCV =  list(c(0,1),
                   c(1,3)),
  catchCV = list(c(0,1,3),
                 c(0,1,3)),
  recmodel = 2, # Chose recruitment model (2 = estimated)
  estCV = c(0,2,0),
  beta = beta, # Hockey stick plateau
  nllfactor = c(1,1,0.01) # Factor for relative strength of log-likelihood

)

# Get initial parameter structure
parms <- getParms(df.tmb)

# Get non-estimated parameters, based on info in df.tmb
mps <- getMPS(df.tmb, parms)

#parms$SDsurvey <- rep(0.3, length(parms$SDsurvey))
# Set boundaries
# This model works best if SDsurvey is mapped
#parms[names(parms) == 'SDsurvey'][[1]] <- c(0.5, 0.3, 0.3,.3)
#parms$logQ <- c(1,2,3,4)
# Try initial parameters in weird spots
sas <- runAssessment(df.tmb, parms = parms, mps = mps)
#
sas$reps

p1 <- smsPlots(df.tmb = df.tmb,sas)
# # Plot numbers at age. Is it F or N that goes wrong
 mr <- mohns_rho(df.tmb, peels = 5, parms, mps, plotfigure = TRUE)
# saveOutput(df.tmb, sas , MR = mr, save = TRUE, wd = wd)
# Compare with sms
reps <- sas$reps

wd.sms <- "M:/Tobis/Tobis_assessment/SMS_2022/SAN-area-3r/"

sms.admb <- read.table(file.path(wd.sms,'summary_table_raw.out'), header = TRUE)


df.plot.tmb <- data.frame(SSB = getSSB(df.tmb,sas)$SSB,
                     R = getR(df.tmb,sas)$R,
                     Catchtot = c(getCatch(df.tmb,sas)$Catch,NA),
                     model = 'TMB',
                     year = c(df.tmb$years,max(df.tmb$years)+1)

)
sms.plot <- sms.admb %>% select(SSB, Rec, Yield.hat, Year) %>% rename('Catchtot' = Yield.hat,
                                                                      'year' = Year,
                                                                      'R' = Rec) %>%
  mutate(model = 'admb')


df.plot <- rbind(df.plot.tmb, sms.plot) %>% pivot_longer(1:3) %>%
  ggplot(aes(x=  year, y = value, color = model))+geom_line()+facet_wrap(~name, scales = 'free_y')+
  theme_classic()
df.plot

ggsave(filename = 'Area4.png', df.plot, width = 16, height = 10, units ='cm')

