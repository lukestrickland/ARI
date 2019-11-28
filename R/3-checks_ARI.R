source("dmc/dmc.R")

dtmp <- load("data/ARIDat.RData")
#[1] "ARIdat"

### ------------  Some checks:
table(ARIdat$S)
# ar   ar2 
# 56120     0 
table(ARIdat$R)
# NR    AR   AR2 
# 7573 48547     0 
table(ARIdat$SS)
# go  stop 
# 42090 14030 

# Number of participants:
unique(ARIdat$s) # 61

tmp <- ARIdat$R; levels(tmp) <- c(NA,"ar","ar2")
crct <- ARIdat$S==tmp; mean(crct,na.rm=TRUE)
is.go <- ARIdat$SS=="go"

# Number of go and stop trials/participant
tapply(ARIdat$RT[is.go],ARIdat$s[is.go],length)   # 690 Go
tapply(ARIdat$RT[!is.go],ARIdat$s[!is.go],length) # 230 Stop (0.25%)
 
# Accuracy on go and failed stop trials
go_error <- round(1-tapply(crct[is.go],ARIdat$s[is.go],mean,na.rm=TRUE),3) # Go trials
failed_stop_error <- round(1-tapply(crct[!is.go],ARIdat$s[!is.go],mean,na.rm=TRUE),3) # Stop trials

# Go omissions
go_omissions <- round(tapply(ARIdat$RT[is.go],ARIdat$s[is.go],function(x) mean(is.na(x))),3)

# Inhibition rate on stop trials
inhibition_rate <- round(tapply(ARIdat$RT[!is.go],ARIdat$s[!is.go],function(x) mean(is.na(x))),3)

layout(matrix(1:4,2,2))
hist(go_error) # No choice errors as expected; sanity check
hist(failed_stop_error) # No choice errors as expected; sanity check
hist(go_omissions) # All under 7%
hist(inhibition_rate) #  All between 0.485 and 0.515 
# which(inhibition_rate < 0.25 | inhibition_rate > 0.75)

### -------- Independence: mean signal respond RT (failed stop) < mean go RT:
mean_RTs <- tapply(ARIdat$RT,cbind.data.frame(Subj=ARIdat$s,Go=is.go),mean,na.rm=TRUE)
cbind(mean_RTs,mean_RTs[,1]-mean_RTs[,2],mean_RTs[,1]>mean_RTs[,2])[which(mean_RTs[,1]>mean_RTs[,2]),]

# FALSE      TRUE               
# 3  0.8106137 0.8095657 1.048083e-03 1
# 4  0.8099535 0.8022698 7.683703e-03 1 # big
# 5  0.8085611 0.8065323 2.028788e-03 1
# 16 0.8136104 0.8123804 1.230044e-03 1
# 19 0.8120674 0.8114658 6.016362e-04 1
# 22 0.8052780 0.8052155 6.258334e-05 1
# 27 0.8020593 0.8011687 8.906416e-04 1
# 35 0.8036024 0.8021893 1.413120e-03 1
# 48 0.8265579 0.8124326 1.412535e-02 1 # big
# 50 0.8251921 0.8239024 1.289766e-03 1
# 53 0.8054167 0.8037904 1.626312e-03 1
# 62 0.8030417 0.8013033 1.738390e-03 1

### Go RTs Collapsed across subjects: lower cut-off? outliers?
pdf("img/goRT_ARI_collapsed.pdf",width=10,height=6)
layout(matrix(1:2,1))
hist(ARIdat$RT[is.go],breaks=40,main="GO RTs")
hist(ARIdat$RT[!is.go],breaks=40,main="Signal-respond RTs")
dev.off()

#head(sort(ARIdat$RT[is.go]))

# Plot dsitribution of SSDs
hist(ARIdat$SSD,breaks=40)
range(ARIdat[ARIdat$SS=="stop","SSD"])
#[1] 0.32 0.74 # Quite narrow range

### Check funny-looking median signal-respond RTs:
tmp_funny <- ARIdat[ARIdat$s=="4",]
layout(1)
plot_SS_srrt.dmc(tmp_funny)
tapply(!is.na(tmp_funny$RT),tmp_funny[,c("SS","SSD")],sum)[2,]  # How many RTs /SSD?
#tapply(tmp_funny$RT,tmp_funny[,c("SS","SSD")],mean,na.rm=TRUE)[2,] # Means /SSD

# Conclusions:
# 2: Probably just noisy / can't judge
# 4: : Probably just noisy / can't judge
# 5: Bad
# 6: Bad
# 10: Bad
# 13: Probably just noisy / can't judge
# 16: Probably just noisy / can't judge
# 19: Bad
# 22: Probably just noisy / can't judge / TF
# 23: Probably just noisy / can't judge
# 27: Bad
# 34: Bad
# 38: Bad
# 40: Bad
# 42: Bad
# 48: Bad
# 52: Probably just noisy / can't judge
# 56: Probably just noisy / can't judge
# 58: Probably just noisy / can't judge
# 62: Bad
# 72: Probably just noisy / can't judge

bad = c(5,6,10,19,27,34,38,40,42,48,62)
pdf("img/all_ARI.pdf",width=26, height=20)

par(cex.main=3,cex.lab=3)
layout(matrix(1:24,4,6,byrow=TRUE))
for(i in levels(ARIdat$s)){
  
  if(i %in% bad) 
    hist(ARIdat$RT[is.go & ARIdat$s==i],main=i,xlab="go RT",breaks=40,col="red")
  else
    hist(ARIdat$RT[is.go & ARIdat$s==i],main=i,xlab="go RT",breaks=40)
  
  hist(ARIdat[ARIdat$s==i,"SSD"],breaks=40, main=i,xlab="SSD distribution")
  
  ssd <- ARIdat[ARIdat$s==i,"SSD"]
  plot(1:length(ssd),ssd,main=i,xlab="SSD per trial")
  
  plot_SS_if.dmc(ARIdat[ARIdat$s==i,],main=i)
  
  plot_SS_srrt.dmc(ARIdat[ARIdat$s==i,],main=i)
  
  tmp_go <- ARIdat$RT[ARIdat$s==i & ARIdat$SS=="go"]
  tmp_go <- tmp_go[!is.na(tmp_go)]
  trials <- 1:length(tmp_go)
  plot(trials,tmp_go,ylim=c(0,2),xlim=c(0,700),ylab="go RT by trial",xlab="trial",main=i)
  abline(lm(tmp_go~trials),col="red",lwd=2)
  }
dev.off()

### Bisset plot to check independence:
source("R/0-Plot.R")

levels(ARIdat$SS) <- c("GO","SS")

# Extract summary statistics
viol<- summary.violation(data=ARIdat[,c("s","SS","RT","SSD")])

pdf("img/independence_ARI.pdf",width=12, height=8)
layout(matrix(1:2,2))
bissett.plot(violation.summary=viol,main="individuals")
aggregate.plot(violation.summary=viol,main="average")
dev.off()

#Luke plot RTs

pdf("img/RTdists.pdf",width=12, height=20)
ggplot(data=ARIdat, aes(x=RT)) + 
  geom_histogram() +facet_wrap(~s, ncol=3)
dev.off()

ARIdat %>% group_by(s) %>% filter(SS=="go") %>% summarise(count_long=mean(R=="NR"))

test <- ARIdat %>% group_by(s) %>% filter(SS=="go") %>% summarise(nrs=mean(R=="NR"))

test[test$nrs>0.03,]
