# Running on 
rm(list=ls()) 
source ("dmc/dmc.R")
load_model ("WALD-SSEXG","waldSSexg_probit.R") 
tmp=load("WEXG.RData")
cores=63 # Ideally set at number of subjects if cores available

# STAGE 1: indivdiual fits
sWEXG <- h.RUN.dmc(sWEXG,cores=cores)
save(WEXG,sWEXG,file="donesWEXG")
load("donesWEXG")
# Make hierarchical
p.prior <- sWEXG[[1]]$p.prior
p1 <- get.p.vector(sWEXG[[1]])[names(p.prior)]
s.prior <- prior.p.dmc(p1=p1,p2=p1,
  dists=rep("gamma",length(p1)))
pp.prior=list(p.prior,s.prior)
hstart <- make.hstart(sWEXG)
theta1 <- make.theta1(sWEXG)
hWEXG <- h.samples.dmc(nmc=100,p.prior,WEXG,pp.prior,
  hstart.prior=hstart,theta1=theta1,thin=10)
rm(WEXG,sWEXG); gc()

# STAGE 2: Hierarachical
pp.prior <- attr(hWEXG,"hyper")$pp.prior
hWEXG  <- h.run.unstuck.dmc(hWEXG, p.migrate = .05, cores = cores)
save(hWEXG,file="WEXG.RData")
hWEXG <- h.run.converge.dmc(h.samples.dmc(nmc=120, samples=hWEXG),
  thorough=TRUE,nmc=40,cores=cores)
save(hWEXG,file="WEXG.RData")

# # Run these lines if not yet converged (see diagnostics below), repeat as 
# # necessary, adjust nmc, up thin if necessary
# hWEXG <- h.run.dmc(h.samples.dmc(samples=hWEXG,nmc=200,thin=10),cores=cores,report=10)
# save(hWEXG,file="WEXG.RData")

# # Following useful if some subjects have a few bad chains and you dont want 
# # to go back to unstuck, first line can be run to give report of what chains
# # will be removed. Second line does a quick run to settle things down, final
# # line takes another a final sample
# hWEXG <- h.samples.dmc(samples=hWEXG,nmc=40,replace.bad.chains=TRUE)
# hWEXG <- h.run.dmc(hWEXG,cores=cores,report=10)
# hWEXG <- h.run.dmc(h.samples.dmc(samples=hWEXG,nmc=200,thin=20),cores=cores,report=10)
# save(hWEXG,file="WEXG.RData")


# POST PROCESS
pdf("WEXG_chains.pdf",height=6,width = 8)
plot.dmc(hWEXG,hyper=TRUE,pll.chain=TRUE) 
plot.dmc(hWEXG,hyper=TRUE,layout=c(3,3))
plot.dmc(hWEXG,hyper=TRUE,layout=c(3,3),p.prior=pp.prior)
dev.off()

ppWEXG <- h.post.predict.dmc(hWEXG,cores=20,gglist=TRUE)
save(hWEXG,ppWEXG,file="WEXG.RData")
pdf("WEXG_fit.pdf",height=6, width = 8)
plot.pp.dmc(ppWEXG,model.legend = FALSE,layout=c(2,2))
dev.off()

h.IC.dmc(hWEXG,DIC=TRUE)

gelman.diag.dmc(hWEXG,hyper=TRUE)
gelman.diag.dmc(hWEXG)

