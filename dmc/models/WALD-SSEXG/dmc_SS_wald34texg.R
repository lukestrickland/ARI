##################  DMC Stop Signal Model Lesson
#
rm(list=ls())
source ("dmc/dmc.R")
#
# This lesson illustrates the use of the lower truncated  
# 1) waldSSexg:  
# 2) waldSSesN: The same as (1) but with a variable number of accumulators
# 3) waldSSexg_probit: The same as (1) but with probabilites on probit
#
# The trunaction minEXG is hard coded to zero in the model files, but can be 
# changed to another value there (in two places, random and likelihood)

# Use a pure Wald model (no start point variability, so A=0, AS=0).


# is.tf <- TRUE
# is.gf <- TRUE
# is.ts <- FALSE
# use.staircase <- TRUE

########################################################## 1)

# load_model ("WALD-SSEXG","waldSSexg.R")
load_model ("WALD-SSEXG","waldSStexg.R") 

# Basic Model
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(A=0,ts=0,tf=0,gf=0)) # No start point noise or TRIAL effects

# # Choose small mu and large sigma so substantial truncation e.g.
# pEXG(0,mu=0.1,sigma=0.2,tau=0.1)
# [1] 0.1838131
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.1,sigma=0.2,tau=0.1,t0=.2)
# pEXG(0,mu=0.2,sigma=0.2,tau=0.1)
# [1] 0.08495332
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.2,sigma=0.2,tau=0.1,t0=.2)
# pEXG(0,mu=0.3,sigma=0.2,tau=0.1)
# [1] 0.03228198
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.3,sigma=0.2,tau=0.1,t0=.2)
# pEXG(0,mu=0.15,sigma=0.1,tau=0.2)
# [1] 0.01223247
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.15,sigma=0.1,tau=0.2,t0=.2)
# pEXG(0,mu=0.5,sigma=0.1,tau=0.1)
# [1] 4.524153e-08
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.5,sigma=0.1,tau=0.1,t0=.2)


# With gf and tf
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(A=0,ts=0)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "v.true"  "v.false" "B"       "mu"      "sigma"   "tau"     "t0"
# [8] "tf"      "gf"
#
# Constants are (see attr(,"constants") ):
#  A ts
#  0  0
#
# Model type = waldss

#
# # Choose small mu and large sigma so substantial truncation e.g.
# pEXG(0,mu=0.1,sigma=0.2,tau=0.1)
# [1] 0.1838131
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.1,sigma=0.2,tau=0.1,t0=.2,tf=.1,gf=.1)
# pEXG(0,mu=0.5,sigma=0.1,tau=0.1)
# [1] 4.524153e-08
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.5,sigma=0.1,tau=0.1,t0=.2,tf=.1,gf=.1)

              

# # check.p.vector(p.vector,model)
# # print.cell.p(p.vector,model)
# n=6
# data <- data.model.dmc(simulate.dmc(p.vector,model,
#                        n=n,SSD=c(Inf,Inf,.25,.25)),model)
# likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25)),model) 

# mean(data$RT[data$SS=="GO"])
# mean(data$RT[data$SS!="GO"],na.rm=TRUE)

# # SSDs
# sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# # Broken down by SSD
# tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("mu",.01,.8,p.vector,data)
profile.dmc("sigma",.01,.5,p.vector,data)
profile.dmc("tau",.01,.5,p.vector,data)
profile.dmc("v.true",2,3,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)

profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)

profile.dmc("A",.25,.75,p.vector,data)

# recovery all good expect sigma consistently a little under-estiamted.

##################################  Mixed 2 and 3 accumulator example

load_model ("WALD-SSEXG","waldSStexgN.R") 

# Just test full model
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2","s3"),SS=c("GO","SS"), # Go stimulus, go and stop trials
               NC=c("a2","a3")),               # 2 choice vs. 3 choice case
  responses <- c("NR","r1","r2","r3"),         # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             N="NC",                           # Number of accumulators
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s3="r3",s1="NR")), # NR mapping ignored
  constants = c(ts=0,N.a2=3,N.a3=4)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "v.true"  "v.false" "B"       "A"       "mu"      "sigma"   "tau"    
#  [8] "t0"      "tf"      "gf"     
# 
# Constants are (see attr(,"constants") ):
#   ts N.a2 N.a3 
#    0    3    4 
# 
# Model type = waldss 


# # Choose small mu and large sigma so substantial truncation e.g.
# pEXG(0,mu=0.1,sigma=0.2,tau=0.1)
# [1] 0.1838131
p.vector  <- c(v.true=2.5,v.false=.5,B=1,A=.5,mu=0.1,sigma=0.2,tau=0.1,t0=.2,tf=.1,gf=.1)
# pEXG(0,mu=0.5,sigma=0.1,tau=0.1)
# [1] 4.524153e-08
p.vector  <- c(v.true=2.5,v.false=.5,B=1,A=.5,mu=0.5,sigma=0.1,tau=0.1,t0=.2,tf=.1,gf=.1)


# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
# Full factorial 12 cells, but remove s3 in a2 for both GO and SS
SSD <- c(rep(Inf,3),rep(.25,3),rep(Inf,3),rep(.25,3)) # Per cell, two ignored
# SSD <- c(rep(Inf,2*6),rep(.25,2*6),rep(Inf,3*6),rep(.25,3*6)) # Per value

# Check likelihood
n <- c(6,6,0,6,6,0,6,6,6,6,6,6)  
data <- data.model.dmc(simulate.dmc(p.vector,model,n=n,SSD=SSD),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,0,2.5e3,2.5e3,0,7.5e3,7.5e3,7.5e3,2.5e3,2.5e3,2.5e3)/10
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=SSD),model) 

# # SSDs
# sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# # Broken down by SSD
# tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("mu",.01,.8,p.vector,data)
profile.dmc("sigma",.01,.5,p.vector,data)
profile.dmc("tau",.01,.5,p.vector,data)
profile.dmc("v.true",2,3,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)
profile.dmc("A",.25,.75,p.vector,data)

# All checks out!


############### Probit

load_model ("WALD-SSEXG","waldSStexg_probit.R") 

# Full model (didnt test ts = trial slope as slow!)
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(ts=0)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "v.true"  "v.false" "B"       "A"       "mu"      "sigma"  
#  [7] "tau"     "t0"      "tf"      "gf"     
# 
# Constants are (see attr(,"constants") ):
# ts 
#  0 
# 
# Model type = waldss 

#

# # Choose small mu and large sigma so substantial truncation e.g.
# pEXG(0,mu=0.1,sigma=0.2,tau=0.1)
# [1] 0.1838131
p.vector  <- c(v.true=2.5,v.false=.5,B=1,A=.5,mu=0.1,sigma=0.2,tau=0.1,t0=.2,
  tf=qnorm(.1),gf=qnorm(.1))
# pEXG(0,mu=0.5,sigma=0.1,tau=0.1)
# [1] 4.524153e-08
p.vector  <- c(v.true=2.5,v.false=.5,B=1,A=.5,mu=0.5,sigma=0.1,tau=0.1,t0=.2,
  tf=qnorm(.1),gf=qnorm(.1))

# # check.p.vector(p.vector,model)
# # print.cell.p(p.vector,model)
# n=6
# data <- data.model.dmc(simulate.dmc(p.vector,model,
#                        n=n,SSD=c(Inf,Inf,.25,.25)),model)
# likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25)),model) 

# mean(data$RT[data$SS=="GO"])
# mean(data$RT[data$SS!="GO"],na.rm=TRUE)

# # SSDs
# sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# # Broken down by SSD
# tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("mu",.01,.8,p.vector,data)
profile.dmc("sigma",.01,.5,p.vector,data)
profile.dmc("tau",.01,.5,p.vector,data)
profile.dmc("v.true",2,3,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("gf",qnorm(.05),qnorm(.15),p.vector,data)
profile.dmc("tf",qnorm(.05),qnorm(.15),p.vector,data)
profile.dmc("A",.25,.75,p.vector,data)

### Do some sampling, use a simpler model with n A 

model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(ts=0,A=0)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "v.true"  "v.false" "B"       "mu"      "sigma"   "tau"    
# [7] "t0"      "tf"      "gf"       
# 
# Constants are (see attr(,"constants") ):
# ts 
#  0 
# 
# Model type = waldss 

# # Choose small mu and large sigma so substantial truncation e.g.
# pEXG(0,mu=0.1,sigma=0.2,tau=0.1)
# [1] 0.1838131
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.1,sigma=0.2,tau=0.1,t0=.2,
  tf=qnorm(.1),gf=qnorm(.1))


n <- c(7.5e3,7.5e3,2.5e3,2.5e3)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25)),model) 

# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS)
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Truncated normal priors:
p1 <- p.vector; p1[1:length(p1)] <- c(2,1,1,.3,.1,.1,.3,-2,-2)
p2=c(2,1,1,1,.5,.5,.5,2,2)
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p.vector,p2=p2, 
  lower=c(rep(0,length(p1)-2),-8,-8),upper=c(rep(NA,3),rep(1,4),rep(6,2)) 
)
# par(mfcol=c(2,5)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Start sampling
samples <- samples.dmc(nmc=100,p.prior,data)
samples <- run.dmc(samples,p.migrate=.05,report=1,cores=9)
plot.dmc(samples,pll.chain = TRUE)

samples <- run.unstuck.dmc(samples.dmc(samples=samples,nmc=100),
  report=10,cores=9,p.migrate=.05,verbose=TRUE)
plot.dmc(samples,pll.chain=TRUE,start=80)
samples1 <- run.converge.dmc(samples.dmc(samples=samples,nmc=120,thin=5),
      report=10,cores=9,verbose=TRUE,nmc=40,minN=500,max.try=20)
samples2 <- run.dmc(samples.dmc(samples=samples1,nmc=250),cores=27)

# Parameter chains look like fat hairy catapillars:
plot.dmc(samples2,layout=c(2,6))
plot.dmc(samples2,pll.chain=TRUE)

# Rhat not great needs thinning.
gelman.diag.dmc(samples2)

# Plot the posteriors:
plot.dmc(samples2,layout=c(2,6),p.prior=p.prior)

# Good parameter recovery
check.recovery.dmc(samples2,p.vector)
#                v.true v.false    B    mu sigma  tau   t0    tf    gf
# True             2.50    0.50 1.00  0.10  0.20 0.10 0.20 -1.28 -1.28
# 2.5% Estimate    2.46    0.47 0.99  0.00  0.01 0.00 0.19 -1.71 -1.33
# 50% Estimate     2.52    0.54 1.02  0.05  0.21 0.10 0.20 -1.25 -1.30
# 97.5% Estimate   2.57    0.61 1.04  0.16  0.32 0.27 0.20 -1.01 -1.27
# Median-True      0.02    0.04 0.02 -0.05  0.01 0.00 0.00  0.03 -0.02

# Good fits, as expected for simulated data: 
pp <- post.predict.dmc(samples2,n.post=1000)
pp2 <- post.predict.dmc(samples2,n.post=1000,save.simulation=TRUE)

# Fit all look good.
plot.pp.dmc(pp,style="cdf",layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="")
layout(1)
plot_SS_if.dmc(data=samples2$data,sim=pp2)
# Inhibition function 
#     SSD    n     p
# 1  0.00   21 0.220
# 2  0.05   55 0.363
# 3  0.10  166 0.879
# 4  0.15  410 0.847
# 5  0.20  730 0.541
# 6  0.25  977 0.475
# 7  0.30 1022 0.664
# 8  0.35  818 0.347
# 9  0.40  478 0.113
# 10 0.45  217 0.664
# 11 0.50   82 0.161
# 12 0.55   20 0.576
# 13 0.60    4 0.454
layout(1)
plot_SS_srrt.dmc(data=samples2$data,sim=pp2)
# Median signal-respond RTs 
#     SSD nrt     p n.sim
# 1  0.00   6 0.822   992
# 2  0.05  15 0.223  1000
# 3  0.10  40 0.235  1000
# 4  0.15 126 0.730  1000
# 5  0.20 284 0.616  1000
# 6  0.25 446 0.254  1000
# 7  0.30 529 0.626  1000
# 8  0.35 491 0.745  1000
# 9  0.40 327 0.749  1000
# 10 0.45 151 0.277  1000
# 11 0.50  66 0.974  1000
# 12 0.55  16 0.158  1000
# 13 0.60   4 0.523   995



#########  Finally try a verison in which the trucation is 50ms ----

# This is done by changing likelihood.dmc and random.dmc in the model file
load_model ("WALD-SSEXG","waldSSt50exg_probit.R") 

model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(ts=0)) # No start point noise or TRIAL effects

# # Choose small mu and large sigma so substantial truncation e.g.
# pEXG(0,mu=0.1,sigma=0.2,tau=0.1)
# [1] 0.1838131
p.vector  <- c(v.true=2.5,v.false=.5,B=1,A=.5,mu=0.1,sigma=0.2,tau=0.1,t0=.2,
  tf=qnorm(.1),gf=qnorm(.1))

# Now try a version in which t0 < truncation
p.vector  <- c(v.true=2.5,v.false=.5,B=1,A=.5,mu=0.1,sigma=0.2,tau=0.1,t0=0,
  tf=qnorm(.1),gf=qnorm(.1))

# # check.p.vector(p.vector,model)
# # print.cell.p(p.vector,model)
# n=6
# data <- data.model.dmc(simulate.dmc(p.vector,model,
#                        n=n,SSD=c(Inf,Inf,.25,.25)),model)
# likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25)),model) 

# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("mu",.01,.8,p.vector,data)
profile.dmc("sigma",.01,.5,p.vector,data)
profile.dmc("tau",.01,.5,p.vector,data)
profile.dmc("v.true",2,3,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("gf",qnorm(.05),qnorm(.15),p.vector,data)
profile.dmc("tf",qnorm(.05),qnorm(.15),p.vector,data)
profile.dmc("A",.25,.75,p.vector,data)
