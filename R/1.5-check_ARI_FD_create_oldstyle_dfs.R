#note some parts of tidyverse don't install on
#cl servers
#can just use dplyr for this script instead
#if rquired
#library(tidyverse)
library(dplyr)
library(tidyr)
source("R/0-functions.R")

#Read data in from JSON files using Maff Gretton's script
#(requires JSON lite)
#saves to csvs which are loaded below
library(jsonlite)
source("R/1-ConvertResultsToCSV.R")

ARdats <- read.csv("data/AR_ARI_results.csv")

#remove practice trials and first block of 
#stop trials (which we consider stop practice)
ARdats$practice <- ARdats$global_trial_num<285

ARdats <- ARdats %>% group_by(id) %>%  mutate(SSD = lag(ss_delay,1))

#t test to check whether 
#RT 'significantly' slower on ARI trials than AR trials
# (we want it to be faster but decided to settle for just
#'not substantially slower')
ts <- ARdats %>% group_by(id) %>% 
  filter(response_type!="NO_RESPONSE" & !practice) %>% 
  summarise(t= t.test(response_time[trial_type=="AR"] ,
                      response_time[trial_type=="ARI"])[["statistic"]],
            p=t.test(response_time[trial_type=="AR"] ,
                     response_time[trial_type=="ARI"])[["p.value"]])

#participants to exclude for context independence issue
ts %>% filter(t<0 & p<.05)

#check stop rates for reasonability
ARI_stoprates <- ARdats %>% group_by(id, trial_type) %>% filter (!practice) %>%  
  summarise(stoprate=mean(response_type=="NO_RESPONSE"))

#all stop rates are within 10% of 50%
ARI_stoprates %>% filter (abs(stoprate-0.5) >0.1 & trial_type=="ARI")

#Only a couple of participants who stop on go trials and not that often
#note p23 will be removed anyway due to FD reasons below
ARI_stoprates %>% filter (stoprate>0.05 & trial_type=="AR")

# 

FDdats <- read.csv("data/FD_results.csv")
FDdats$practice <- FDdats$global_trial_num<285
FDdats$ss <- factor(FDdats$ss, levels=c("FALSE", "TRUE"),
                               labels = c("go", "stop"))


FDdats <- FDdats %>% group_by(id) %>%  mutate(SSD = lag(ss_delay,1))

#check context independence

FDRTs <- FDdats %>% group_by(id, bias, ss) %>% 
  filter(response_type!="NO_RESPONSE" &!practice) %>% 
  summarise(meanRT=mean(response_time))

FDRTs_spread <- spread(FDRTs, ss, meanRT)

FDRTs_spread %>% filter(go<stop)

ts <- FDdats %>% group_by(id, bias) %>% 
  filter(response_type!="NO_RESPONSE" & !practice) %>% 
  summarise(t= t.test(response_time[ss=="go"] ,
                      response_time[ss=="stop"])[["statistic"]],
            p=t.test(response_time[ss=="go"] ,
                     response_time[ss=="stop"])[["p.value"]])

#just 8- blue violates context independence
ts %>% filter (t <0 & p <.05)


FDaccs <- FDdats %>% group_by(id, bias, ss) %>% 
  filter(response_type!="NO_RESPONSE" &!practice) %>% 
  summarise(meanacc=mean(substr(response_type,10,20)==bias))

#Note p47 is low acc - might want to exclude
FDaccs %>% filter(meanacc <0.8)


#check stop rates for reasonability
FD_stoprates <- FDdats %>% group_by(id, ss) %>% filter (!practice) %>%  
  summarise(stoprate=mean(response_type=="NO_RESPONSE"))

FD_stoprates %>% filter (abs(stoprate-0.5) >0.1 & ss=="stop")

#quite a few severely bad stop rates- p 8, p23, p46, p54
#currently these have been excluded from standard analyses
FD_stoprates %>% filter (abs(stoprate-0.5) >0.2 & ss=="stop")

#plots below reveal various problems with the stopping
plot_ssdelay(FDdats, "8")
plot_ssdelay(FDdats, "23")
plot_ssdelay(FDdats, "46")
plot_ssdelay(FDdats, "54")

#NOTE: Also some pathologies with SSD for p6
#and 47. Maybe should also be excluded?
plot_ssdelay(FDdats, "6")
plot_ssdelay(FDdats, "47")

###CREATE FINAL ARI AND FD DATA FRAMES

bad <- c("8", "23", "46", "54", "4", "48")

#total exclusions (counting replaced files):
length(bad) + 3

#proportion excluded - quite high already
(length(bad) + 3) / (length(unique(ARdats$id))+3)

#remove bad participants and convert back to data frame
#to avoid potential issues with DMC
AR_final <- as.data.frame(ARdats %>% filter(!(id %in% bad)))
FD_final <- as.data.frame(FDdats %>% filter(!(id %in% bad)))

table(AR_final$id, AR_final$practice)
table(FD_final$id, FD_final$practice)

save(AR_final, FD_final, file="data/final_dfs.RData")


is.stop <- FDdats$ss=="stop"
godat <- FDdats[!is.stop & !FDdats$practice,]
stopdat <- FDdats[is.stop& !FDdats$practice,] 
stopdat$id <- factor(stopdat$id)
stopdat[stopdat$response_type=="NO_RESPONSE","response_time"] <- NA


SSRTs <- get_ssrt_FD(godat, stopdat)

SSRT_integrate <- do.call("rbind", lapply(SSRTs, function(x) x$SSRTintegrate))

write.csv(SSRT_integrate, file= "img/FDSSRTs.csv")



SSRTs <- get_ssrt(ARdats[ARdats$trial_type=="ARI"& !ARdats$practice,],
                  ARdats[ARdats$trial_type=="AR" & ARdats$response_type=="RESPONSE" & !ARdats$practice,])

SSRTs_integrate <-
  as.data.frame(do.call("rbind", lapply(SSRTs, function(x)
    x$SSRTintegrate)))

SSRTs_integrate$pid <-names(SSRTs)

write.csv(SSRTs_integrate, file="img/ari_SSRTs_integrate.csv", row.names=F)


