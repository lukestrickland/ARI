# Included by ConvertResultsToCSV.R

FactoriseARARI <- function (trials)
{ 
  trials$trial_type    <- factor(trials$trial_type, levels=c(TASKCODE_ARI_AR_TRIAL,TASKCODE_ARI_ARI_TRIAL), labels=c("AR", "ARI"))
  
  block_type <- bitwAnd (trials$block_type, 63)
  trials$block_type    <- factor(block_type, levels=c(BLOCKTYPE_AR_ONLY,BLOCKTYPE_AR_ARI), labels=c("AR_ONLY", "AR_ARI")) 
  
  trials$response_type <- factor(trials$response_type, levels=c(1,2), labels=c("NO_RESPONSE", "RESPONSE"))
  trials$response_correct <- factor(trials$response_correct, levels=c(1,2), labels=c("CORRECT", "INCORRECT"))
  
  trials
}

FactoriseFD <- function (trials, participant_num)
{
  trials$bias          <- factor(bitwAnd (trials$trial_type, 3), levels=c(1,2), labels=c('ORANGE', 'BLUE')) 
  trials$SS            <- factor(bitwAnd (trials$trial_type, 8),   levels=c(0,8), labels=c('FALSE', 'TRUE'))  
  
  # TODO: trial_type factoring assumes one difficulty-level only
  trials$trial_type    <- factor(trials$trial_type, levels=c(1,2,9,10), labels=c("ORANGE", "BLUE", "SS_ORANGE", "SS_BLUE"))
  
  
  # TODO: block_type factoring NYI for Elise's conditions
  block_type <- bitwAnd (trials$block_type, 63)
  trials$block_type    <- factor(block_type, levels=c(3,4), labels=c("FLASHING_ONLY", "FLASHING_SS")) 
  
  trials$response_type <- factor(trials$response_type, levels=c(1,2,3), labels=c("NO_RESPONSE", "RESPONSE_ORANGE", "RESPONSE_BLUE"))
  trials$response_correct <- factor(trials$response_correct, levels=c(1,2), labels=c("CORRECT", "INCORRECT"))
  
  #   Not needed    
  #    trials$accuracy      <- factor(trials$accuracy, levels=c(-1000,1000), labels=c("INCORRECT", "CORRECT"))
  
  keys_lr <- bitwAnd (trials$counter_bits, 1)
  trials$counter_bits <- factor(keys_lr, levels=c(0,1), labels=c("TRUE", "FALSE")) 
  
  trials
}

plot_ssdelay <- function(dats, id) {
  plot(rownames(dats[dats$id==id,]), dats$ss_delay[dats$id==id])
  lines(x = rownames(dats[dats$id==id,]), dats$ss_delay[dats$id==id])
}

get_ssrt <- function(aridat,ardat) {
  aridat$id <- factor(aridat$id)
  out <- vector(mode="list",length=length(levels(aridat$id)))
  names(out) <- levels(aridat$id)
  for (i in levels(aridat$id)) {
    ari <- aridat[aridat$id==i,]
    ar <- ardat[ardat$id==i,]
    nri <- tapply(ari$response_type=="RESPONSE",ari$SSD,sum)
    ni <- table(ari$SSD)
    out[[i]] <- ssrt(nri,ni,goi=ar$response_time,ssd=ari$SSD)
  } 
  out
}

ssrt <- function(nri,ni,goi,ssd,lowp=0,hip=1,lown=1,do.weight=TRUE) 
  # SSRT estimates, nri=counts of signal responds for each ssd
  #                 ni=corresponding trial counts, goi=vector of go rts
  # removes cases where n<=lown, and pr(signal respond) outside of lop-hip
  # Calcualtes SSRTobs = SSRT for each SSD by integration method
  #            SSRTav = average of SSRTobs (possibly weighted by n)
  #            SSRTmean = SSRT by mean method
  #            SSRTmedian = SSRT by median method (possibly weighted by n)
{  
  SSRTintegrate <- quantile(goi,probs=sum(nri)/sum(ni)) - mean(ssd)
  nri <- nri[!is.na(ni)]
  ni <- ni[!is.na(ni)]
  pr <- nri/ni
  ok <- ni>lown & pr >= lowp & pr <= hip
  SSRTobs <- quantile(goi,probs=pr[ok])-as.numeric(names(pr[ok]))
  if (!do.weight) SSRTav <- mean(SSRTobs) else
    SSRTav <- sum(SSRTobs*ni[ok]/sum(ni[ok]))
  dpr <- diff(pr)
  SSRTmean <- mean(goi) - sum(dpr*as.numeric(names(dpr)))/
    (pr[length(pr)]-pr[1])
  data <- as.data.frame(cbind(y=as.numeric(names(pr)),x=pr))
  if (!do.weight)
    lmj <- lm(y~x,data) else
      lmj <- lm(y~x,data,weights=ni)
  SSRTmed <- median(goi) - predict(lmj,data.frame(x=.5))
  list(n=ni,pr=pr,ok=ok,SSRTintegrate=SSRTintegrate,SSRTobs=SSRTobs,SSRTav=SSRTav,
       SSRTmean=SSRTmean,SSRTmed=SSRTmed)  
}

get_ssrt_FD <- function(godat,stopdat) {
  colnames(godat)[15] <- "response"
  colnames(stopdat)[15] <- "response"
  out <- vector(mode="list",length=length(levels(stopdat$id)))
  names(out) <- levels(stopdat$id)
  for (i in levels(stopdat$id)) {
    ari <- stopdat[stopdat$id==i,]
    ar <- godat[godat$id==i,]
    nri <- tapply(ari$response!="NO_RESPONSE",ari$SSD,sum)
    ni <- table(ari$SSD)
    out[[i]] <- ssrt(nri,ni,goi=ar$response_time,ssd=ari$SSD)
  } 
  out
}

