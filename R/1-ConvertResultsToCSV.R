#Script by Maff Gretton to read JSON files
#Run top to bottom
rm(list = ls())
# install.packages("jsonlite")

library(jsonlite)

results_dir <- "Results"

fnams <- dir(path=results_dir, pattern="^Pa.*.json$")
raw.dats <-  vector(mode="list",length=length(fnams))
names(raw.dats) <- lapply(strsplit(fnams,".",fixed=TRUE),function(x){x[[1]]})

# Load
for (i in 1:length(fnams)) {
  f <- file(paste0(results_dir, "\\", fnams[i]))
  text <- suppressWarnings(readLines (f))
  raw.dats[[i]] <- fromJSON (text)
  
  close(f)
}

source ("R/0-ExperimentConstants.R")
source ("R/0-functions.R")

# Split into FD files and ARI files
session_type <- vector()
for (i in 1:length(fnams)) {
    if (raw.dats[[i]]$session_type == TASKTYPE_ARI)      { session_type[i] <- TASKTYPE_ARI }
    if (raw.dats[[i]]$session_type == TASKTYPE_FLASHING) { session_type[i] <- TASKTYPE_FLASHING }
}
    

ari.dats <- raw.dats[session_type == TASKTYPE_ARI]
fd.dats  <- raw.dats[session_type == TASKTYPE_FLASHING]

BuildCSVData <- function (dats, factoriser)
{
    csvdat <- data.frame() 
    
    for (i in 1:length(dats)) 
    {
        dats[[i]]$trials <- factoriser (dats[[i]]$trials) 
        
        row <- cbind(dats[[i]]$participant_id, dats[[i]]$task_version, dats[[i]]$session, dats[[i]]$session_type, dats[[i]]$trials)
        
        csvdat <- rbind(csvdat, row, make.row.names=F)  
    }
    
    csvdat
}

# Factorise and build the AR/ARI data for CSV export
ari.csvdat <- BuildCSVData (ari.dats, FactoriseARARI)
fd.csvdat <- BuildCSVData (fd.dats, FactoriseFD)

# Reorganise columns into a nicer order
ari.csvdat <- ari.csvdat[c(5,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16)]
fd.csvdat  <- fd.csvdat[c(5,1,2,3,4, 17, 6,7,8,9,10, 18, 19, 11,12,13,14, 16)]

    
names (ari.csvdat) <- c("timestamp", "id", "version", "session", "session_type", 
                        "level", "level_type", "global_trial_num", "block_trial_num", "trial_type", 
                        "response_type", "response_correct", "response_time", "ss_delay", "diff", 
                        "score")

names (fd.csvdat)  <- c("timestamp", "id", "version", "session", "session_type", "keys_lr",
                        "level", "level_type", "global_trial_num", "block_trial_num", 
                        "trial_type", "bias", "ss",
                        "response_type", "response_correct", "response_time", "ss_delay", "score")

write.csv (ari.csvdat, "data/AR_ARI_Results.csv")
write.csv (fd.csvdat, "data/FD_Results.csv")


##LUKE DO SAME THING FOR REPLACED DATA
#all just copy+pasted but different results_dir


results_dir <- "Results/REPLACED"

fnams <- dir(path=results_dir, pattern="^Pa.*.json$")
raw.dats <-  vector(mode="list",length=length(fnams))
names(raw.dats) <- lapply(strsplit(fnams,".",fixed=TRUE),function(x){x[[1]]})

# Load
for (i in 1:length(fnams)) {
  f <- file(paste0(results_dir, "\\", fnams[i]))
  text <- suppressWarnings(readLines (f))
  raw.dats[[i]] <- fromJSON (text)
  
  close(f)
}

# Split into FD files and ARI files
session_type <- vector()
for (i in 1:length(fnams)) {
  if (raw.dats[[i]]$session_type == TASKTYPE_ARI)      { session_type[i] <- TASKTYPE_ARI }
  if (raw.dats[[i]]$session_type == TASKTYPE_FLASHING) { session_type[i] <- TASKTYPE_FLASHING }
}


ari.dats <- raw.dats[session_type == TASKTYPE_ARI]
fd.dats  <- raw.dats[session_type == TASKTYPE_FLASHING]

BuildCSVData <- function (dats, factoriser)
{
  csvdat <- data.frame() 
  
  for (i in 1:length(dats)) 
  {
    dats[[i]]$trials <- factoriser (dats[[i]]$trials) 
    
    row <- cbind(dats[[i]]$participant_id, dats[[i]]$task_version, dats[[i]]$session, dats[[i]]$session_type, dats[[i]]$trials)
    
    csvdat <- rbind(csvdat, row, make.row.names=F)  
  }
  
  csvdat
}

# Factorise and build the AR/ARI data for CSV export
ari.csvdat <- BuildCSVData (ari.dats, FactoriseARARI)
fd.csvdat <- BuildCSVData (fd.dats, FactoriseFD)

# Reorganise columns into a nicer order
ari.csvdat <- ari.csvdat[c(5,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16)]
fd.csvdat  <- fd.csvdat[c(5,1,2,3,4, 17, 6,7,8,9,10, 18, 19, 11,12,13,14, 16)]


names (ari.csvdat) <- c("timestamp", "id", "version", "session", "session_type", 
                        "level", "level_type", "global_trial_num", "block_trial_num", "trial_type", 
                        "response_type", "response_correct", "response_time", "ss_delay", "diff", 
                        "score")

names (fd.csvdat)  <- c("timestamp", "id", "version", "session", "session_type", "keys_lr",
                        "level", "level_type", "global_trial_num", "block_trial_num", 
                        "trial_type", "bias", "ss",
                        "response_type", "response_correct", "response_time", "ss_delay", "score")

write.csv (ari.csvdat, "data/AR_ARI_Replaced_Results.csv")
write.csv (fd.csvdat, "data/FD_Replaced_Results.csv")

