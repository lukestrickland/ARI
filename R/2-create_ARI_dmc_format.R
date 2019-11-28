#note some parts of tidyverse don't install on
#cl servers
#can just use dplyr for this script instead
#if rquired
#library(tidyverse)
library(dplyr)
library(tidyr)
source("R/0-functions.R")

#Read data in from JSON files using Maff Gretton's script
#and save in csvs(requires JSON lite)
#source("R/1-ConvertResultsToCSV.R")

ARdatsnew <- read.csv("data/AR_ARI_results.csv") %>% 
                    mutate(id=as.character(id))

#Extra data which was replaced on the basis of FD results
#(FD results now irrelevant)
ARreplaced <- read.csv("data/AR_ARI_Replaced_results.csv") %>% 
                        mutate(id=paste(id, "_old", sep=""))

#Only p42 from ARreplaced was excluded based on ARI so add other 
#two into full data frame

ARdats <- full_join(ARdatsnew, ARreplaced %>% filter(id %in% c("4_old", "15_old")))

#remove practice trials and first block of 
#stop trials (which we consider stop practice)
ARdats$practice <- ARdats$global_trial_num<285

ARdats <- ARdats %>% group_by(id) %>%  mutate(SSD = lag(ss_delay,1))
ARdats$SSD[ARdats$trial_type=="AR"] <- Inf

##Dora's re-formatting for dmc

dat <- subset(ARdats,!ARdats$practice)
dat <- data.frame(s = dat$id,
                  level=dat$level,
                  TRIAL=dat$block_trial_num,
                  S=dat$trial_type,
                  SS=dat$trial_type,
                  R=dat$response_type,
                  RT=dat$response_time,
                  SSD=dat$SSD,
                  diff=dat$diff)

### make S column
dat$S <- as.vector(as.character(dat$S))
dat$S<-"ar"
dat$S <- factor(dat$S,levels=c("ar","ar2"))

### make R column
dat$R <- as.vector(as.character(dat$R))
dat$R[dat$R=="NO_RESPONSE"]<-"NR"
dat$R[dat$R=="RESPONSE"]<-"AR"
dat$R <- factor(dat$R,levels=c("NR","AR","AR2"))

### make SS column
levels(dat$SS) <- c("go","stop")

###Fix RT column: replace ~1 for NO_RESPONSE with NR
dat$RT[dat$R=="NR"]<-NA

### Fix s column
dat$s <- factor(dat$s)

### Fix SSD column
dat$SSD <- round(dat$SSD,3)

#Check my refactoring of parse didn't affect anything
#Old subjects who dora didn't look at filtered and
#df rearranged appropriately
load("data/ARIDat_Dora.RData")
all.equal(dat %>% filter(!grepl("old", s))%>% 
            arrange(as.numeric(as.character(s))) ,
          ARIdat,
          check.attributes=FALSE)

ARIdat <- dat
save(ARIdat,file="data/ARIDat.RData")
