#forecast recipe from stock annex
#average over this many years ago to get pars for forecast
f.west <<- 5 # weight in stock (in past 5 years)
l.west <<- 0 #
f.weca <<- 5 # weight in catch (in past 5 years)
l.weca <<- 0 #
f.rec <<- function(){df.assess$nyears} #recruitment (all past years)
l.rec <<- 0
f.pm <<- 0 #prop mature (most recent)
l.pm <<- 0
f.m <<- 0 #natural mortality (most recent)
l.m <<- 0


#General info on stock to be used in forecast code	
f.age <<- 1 #indicies in R, not ages
l.age <<- 5 #indicies in R, not ages
f.season <<- 1
l.season <<- df.tmb$nseason
plus.group <<- TRUE
rec.season <<- df.tmb$recseason
qq_power <<- FALSE
