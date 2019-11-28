rm(list=ls())
source ("dmc/dmc.R")


# QUICK SIMULATION EXAMPLE TO GET A FEEL ...

# Wald vs. ExG race.
load_model ("WALD-SSEXG","waldSSexg_probit.R") 

# Need to declare two go accumulators, but make sure false accumulator never
# wins by setting its rate = 0
model <- model.dmc(type="waldss",              
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),
  responses = c("NR","r1","r2"),
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")),  # NR mapping ignored
  p.map=list(A="1",B="1",v="M", t0="1", gf="1",
             mu="1",sigma="1",tau="1",tf="1",ts="1"),                                   
  constants = c(A=0,ts=0,v.false=0)) # False never wins + no start point noise or TRIAL effects

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "B"      "v.true" "t0"     "gf"     "mu"     "sigma"  "tau"    "tf"    
# 
# Constants are (see attr(,"constants") ):
#       A      ts v.false 
#       0       0       0 
# 
# Model type = waldss 

# Note tf and gf on probit scale
# Big vlaues of B and v promote symmetry and low variance like AR
p.vector <- c(B=9,v.true=15,t0=.2,mu=.3,sigma=.1,tau=.1,tf=qnorm(.1),gf=qnorm(.1))

# Runs a staircase, starting at 0.25
dat <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,n=c(1e3,0,1e3,0),
                  SSD=c(Inf,Inf,.25,.25)),model)
# Only s1 trials
any(dat$S=="s2")
# False never wins
any((dat$S=="s1" & dat$R=="r2"))
# Good distribution of SSDs
table(dat$SSD)
# Staircase gets 50%
mean(dat$R[dat$SS=="SS"]=="NR")

# Distribution around right speed and low variance
mean(dat$RT,na.rm=TRUE)
hist(dat$RT,breaks="fd")

# All look OK
par(mfrow=c(2,4))
profile.dmc("B",8,10,p.vector,dat)
profile.dmc("v.true",14,16,p.vector,dat)
profile.dmc("t0",.1,.3,p.vector,dat)
profile.dmc("mu",.2,.4,p.vector,dat)
profile.dmc("sigma",.05,.25,p.vector,dat)
profile.dmc("tau",.05,.15,p.vector,dat)
profile.dmc("gf",qnorm(.05),qnorm(.15),p.vector,dat)
profile.dmc("tf",qnorm(.05),qnorm(.15),p.vector,dat)


####  THE REAL THING! 

rm(list=ls())
source ("dmc/dmc.R")

# Truncated Wald vs. ExG race.
load_model ("WALD-SSEXG","waldSSexg_probit.R") 

# Check and clean data

tmp=load("ARIDat.RData")
# Remove irrelevant columns
dat <- ARIdat[,c("s","S","SS","R","RT","SSD")]
# lapply(dat,levels)
# # $s
# #  [1] "1"      "10"     "11"     "12"     "13"     "14"     "15"    
# #  [8] "15_old" "16"     "17"     "18"     "19"     "2"      "20"    
# # [15] "21"     "22"     "23"     "24"     "25"     "26"     "27"    
# # [22] "28"     "29"     "3"      "30"     "31"     "32"     "33"    
# # [29] "34"     "35"     "36"     "37"     "38"     "39"     "4"     
# # [36] "4_old"  "40"     "41"     "42"     "43"     "44"     "45"    
# # [43] "46"     "47"     "48"     "49"     "5"      "50"     "51"    
# # [50] "52"     "53"     "54"     "55"     "56"     "58"     "59"    
# # [57] "6"      "62"     "66"     "7"      "72"     "8"      "9"     
# # 
# # $S
# # [1] "ar"  "ar2"
# # 
# # $SS
# # [1] "go"   "stop"
# # 
# # $R
# # [1] "NR"  "AR"  "AR2"

# # Note that ar2 and AR2 are dummies, dont exist in data
# table(dat$S)
# #    ar   ar2 
# # 57879     0 
# table(dat$R)
#  #   NR    AR   AR2 
#  # 7820 50059     0 
 
# Levels I like (anything works, uses SSD=Inf to identify go trials)
levels(dat$SS) <- c("GO","SS")
# Fix levels, above has one level a substring of another, that bugs
levels(dat$S) <- c("s1","s2")
# Must call first level NR, otherwise ok if not substring
levels(dat$R) <- c("NR","r1","r2")


# Remove fast responses
bad <- dat$RT < .6
bad[is.na(bad)] <- FALSE
# round(sort(100*tapply(bad,dat$s,mean)),2)
#   #    1     12     13 15_old     16     17     18     19      2     21 
#   # 0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00 
#   #   22     25     27     28     29     31     32     33     35     36 
#   # 0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00 
#   #   37     38      4  4_old     40     41     42     44     46     48 
#   # 0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00 
#   #   50     51     52     53     54     55     56     59     72      8 
#   # 0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00 
#   #   14     20     24      3     30     34     43     45     47     58 
#   # 0.11   0.11   0.11   0.11   0.11   0.11   0.11   0.11   0.11   0.11 
#   #   62      9     23     39      5     49     66      7     10     11 
#   # 0.11   0.11   0.22   0.22   0.22   0.33   0.33   0.33   0.43   0.54 
#   #    6     26     15 
#   # 0.65   1.52   2.72 
# mean(bad)*100
# # [1] 0.1397516 
dat <- dat[!bad,]

load_model ("WALD-SSEXG","waldSStexg_probit.R") 

# # Later could try these if it looks like EXG estimates produce mass < 0
# # Zero Truncated Wald vs. ExG race.
# load_model ("WALD-SSEXG","waldSStexg_probit.R") 
# load_model ("WALD-SSEXG","waldSSt50exg_probit.R") 



# Need to declare two go accumulators, but make sure false accumulator never
# wins by setting its rate = 0
model <- model.dmc(type="waldss",              
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),
  responses = c("NR","r1","r2"),
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")),  # NR mapping ignored
  p.map=list(A="1",B="1",v="M", t0="1", gf="1",
             mu="1",sigma="1",tau="1",tf="1",ts="1"),                                   
  constants = c(A=0,ts=0,v.false=0)) # False never wins + no start point noise or TRIAL effects

dWEXG <- data.model.dmc(dat,model)

# PRIOR
# Using the exploration from above to set prior, 
# COULD BE COMPLETELY WRONG! FINE TUNE after RUN.dmc
p.vector <- c(B=9,v.true=15,t0=.2,mu=.3,sigma=.1,tau=.1,tf=qnorm(.1),gf=qnorm(.1))

p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p.vector)),
  p1=p.vector,
  p2=c(5,8,1,.5,.25,.25,3,3),
  lower=c(0,0,.05,0,0,0,NA,NA),
  upper=c(NA,NA,1,NA,NA,NA,NA,NA))
# par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior)


sWEXG <- h.samples.dmc(nmc=100,p.prior,dWEXG)
save(dWEXG,sWEXG,file="WEXG.RData")

  
### EXGAUSSIAN VERSION
  
load_model ("EXG-SS","exgSS.R")

model <- model.dmc(type="exgss",              
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),
  responses = c("NR","r1","r2"),
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")),  # NR mapping ignored
  p.map=list(mu="M",sigma="M",tau="M",gf="1",
             muS="1",sigmaS="1",tauS="1",tf="1",ts="1"),             
  # False never wins + no TRIAL effects
  constants = c(ts=0,mu.false=1e6,sigma.false=.001,tau.false=.001)) 
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "mu.true"    "sigma.true" "tau.true"   "gf"         "muS"       
# [6] "sigmaS"     "tauS"       "tf"        
# 
# Constants are (see attr(,"constants") ):
#          ts    mu.false sigma.false   tau.false 
#       0e+00       1e+06       1e-03       1e-03 
# 
# Model type = exgss

dEXG2 <- data.model.dmc(dat,model)

# PRIOR
p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.03,tauS=.03,tf=.1,gf=.1) 

p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p.vector)),
  p1=p.vector,
  p2=c(1,.5,.2,.2,.2,.2,3,3),
  lower=c(0,0,0,0,0,0,NA,NA),
  upper=c(NA,NA,NA,NA,NA,NA,NA,NA))
# par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior)


sEXG2 <- h.samples.dmc(nmc=100,p.prior,dEXG2)
save(dEXG2,sEXG2,file="EXG2.RData")














