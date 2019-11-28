#note some parts of tidyverse don't install on
#cl servers
#can just use dplyr for this script instead
library(tidyverse)
source("R/0-functions.R")

#Justifying replacements of certain data files

#Below justifies replacing participant 42
# - violation of context independence
# greater AR response time than ARI
# significant according to t-test
ARdats <- read.csv("data/AR_ARI_Replaced_results.csv")

#remove practice trials and first block of 
#stop trials (which we consider stop practice)
ARdats$practice <- ARdats$global_trial_num<285

ari_RTs <- ARdats %>% group_by(id, trial_type) %>% 
  filter(response_type!="NO_RESPONSE" & !practice) %>% 
  summarise(meanRT=mean(response_time))

ari_RTs_SPSS <- spread(ari_RTs, trial_type, meanRT)

#check for context independence
round(ari_RTs_SPSS$AR - ari_RTs_SPSS$ARI, 4)

#p42 shows a potential violation of context independence (10ms)
#Luke note had p41 written here before - typo?

#t test check

ts <- ARdats %>% group_by(id) %>% 
  filter(response_type!="NO_RESPONSE" & !practice) %>% 
  summarise(t= t.test(response_time[trial_type=="AR"] ,
                               response_time[trial_type=="ARI"])[["statistic"]],
            p=t.test(response_time[trial_type=="AR"] ,
                   response_time[trial_type=="ARI"])[["p.value"]])



ts

#Below justifies replacing p4 and p15 based on 
#FD stop rates

FDdats <- read.csv("data/FD_Replaced_results.csv")

FDdats$practice <- FDdats$global_trial_num<285

stop_rates <- FDdats %>% group_by(id) %>% filter (ss==TRUE & !practice) %>%  
  summarise(stoprate=mean(response_type=="NO_RESPONSE"))

#4 giving up and not stopping at all in later blocks
plot_ssdelay(FDdats, "4")

#showing p15 drift up to a strong proactive slowing strategy
plot_ssdelay(FDdats, "15")

p15_FDRTs <- FDdats %>% mutate (trialnums = cut(global_trial_num,
                                 breaks=8)) %>% 
  filter(response_type!="NO_RESPONSE" & !ss & id=="15" & !practice)  %>%
  group_by(trialnums) %>% 
  summarise(meanRT=mean(response_time))

plot(p15_FDRTs$trialnums, p15_FDRTs$meanRT, xlab="Trial Range (trials 0-285 excluded)",
     ylab="FD RT")

