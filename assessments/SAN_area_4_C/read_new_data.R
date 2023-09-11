## Add the new data
library(smsR)
library(tidyverse)
wd <- 'M:/Tobis/Tobis_assessment/SMS_2022_Benchmark/SAN_area_4_Current/'
wddat <- 'M:/Tobis/Tobis_assessment/SMS_2022_Benchmark/data/'
# Load packages and files #
df <- read.csv(file.path(wddat,'Total_catch_in_numbers_and_mean_weight_WKSAND16.csv'))
df <- df[df$area == '4',]

# Remove the 1982 year
# Remove the 1982 year
df <- df[df$aar > 1992,]

# Note that SSB pre 1983 is not correct then

maxage <- 4
years = 1993:2021
nyear <- length(years)
seasons <- 1:2



# Weca first

weca <- df %>% select(paste('mw',0:4, sep=''), aar, hy) %>% arrange(aar) %>%
  mutate('hide' = '#') %>% select(paste('mw',0:4, sep = ''), hide, aar, hy)

myear <- (max(weca$aar)-4):max(weca$aar)

s1end <- as.numeric(colMeans(weca[weca$hy == 1 & weca$aar %in% myear,] %>% select(paste('mw',0:4, sep = ''))))
s2end <- as.numeric(colMeans(weca[weca$hy == 2 & weca$aar %in% myear,] %>% select(paste('mw',0:4, sep = ''))))
# Add two years of 5 year average
wmean1 <- data.frame(mw0 = s1end[1],
                     mw1 = s1end[2],
                     mw2 = s1end[3],
                     mw3= s1end[4],
                     mw4= s1end[5],
                     hide = '#', aar = 5, hy = 1)

wmean2 <- data.frame(mw0 = s2end[1],
                     mw1 = s2end[2],
                     mw2 = s2end[3],
                     mw3= s2end[4],
                     mw4= s2end[5],
                     hide = '#', aar = 5, hy = 1)

weca <- bind_rows(weca, wmean1,wmean2)


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
                      survey.age = list(0:1), # Ages in the two surveys
                      survey.years = list(2008:2021),
                      survey.names = c('Dredge'),
                      survey.quarter = c(2),
                      years = years,
                      seasons = seasons

)



df.test <-data.frame(year = df[df$Area == 1 & df$hy == 2 & df$aar>1982,]$aar,
                     catchage = df[df$Area == 1 & df$hy == 2 & df$aar>1982,]$n2)

canum <- dat[['canum']]

ggplot(canum ,aes(x = year, y= catchage, color = factor(Quarter)))+geom_line()+
  geom_point(data = df.test, col = 'red')+theme_classic()+
  facet_wrap(~Age, scales = 'free_y')

