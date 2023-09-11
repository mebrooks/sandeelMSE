## Add the new data
library(smsR)
library(tidyverse)
wd <- 'M:/Tobis/Tobis_assessment/SMS_2022_Benchmark/SAN_area_3r_Current/'
wddat <- 'M:/Tobis/Tobis_assessment/SMS_2022_Benchmark/data/'
#setwd(wd)
# Load packages and files #
df <- read.csv(file.path(wddat,'Total_catch_in_numbers_and_mean_weight_WKSAND16.csv'))
df <- df[df$area == '3r',]
# Remove the 1982 year
df <- df[df$aar > 1985,]

# Note that SSB pre 1983 is not correct then

maxage <- 4
years = 1986:2021
nyear <- length(years)
seasons <- 1:2



# Weca first

weca <- df %>% select(paste('mw',0:4, sep=''), aar, hy) %>% arrange(aar) %>%
  mutate('hide' = '#') %>% select(paste('mw',0:4, sep = ''), hide, aar, hy)

# Weca and west are the same (excluded 5 year avg)
# Weca and west are the same (excluded 5 year avg)
write.table(weca, file = file.path(wd,'weca.in'), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(weca, file = file.path(wd,'west.in'), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Now do canum
canum <- df %>% select(paste('n',0:4, sep=''), aar, hy) %>% arrange(aar) %>%
  mutate('hide' = '#') %>% select(paste('n',0:4, sep = ''), hide, aar, hy)

write.table(canum, file = file.path(wd,'canum.in'), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Test the data


dat <- getDataSandeel(wd,
                      maxage = maxage,
                      survey.age = list(0:1,1:4), # Ages in the two surveys
                      survey.years = list(2006:2021, 2009:2021),
                      survey.names = c('Dredge','Acoustic'),
                      survey.quarter = c(2,1),
                      years = years,
                      seasons = seasons

)


