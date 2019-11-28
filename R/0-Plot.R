
summary.violation <- function(data,TRIALS=NULL,only.first.stop=TRUE) 
# Data is a data frame with ALL trials (none censored) and columns:
#   s: subject indicator (factor)
#   SS: factor with levels GO and STOP for corresponding trials
#   RT: RT, with NA for non-reponses
#   SSD: SSD, with Inf for GO trials
# TRIALS: trial within block index, used to discard pairs that cross a block
#         break, if NULL they are retained.
#
# only.first.stop=TRUE, default Bissett method, mean stop failure RT from 
#   Trial N minus mean no-stop-signal RT from Trial N-1. Otherwise will match
#   EVERY signal respond with the last go (so same go used in multiple pairs).
{
  
  if (!all((data$SS=="GO")!=is.finite(data$SSD)))
    stop("Something is wrong with data, not all SSD=Inf are coded GO in SS")
  ignore.block.break <- !is.null(TRIALS)
  independence=vector(mode="list",length=length(levels(data$s)))
  names(independence) <- levels(data$s)
  SSDS.all <- sort(unique(data$SSD)[is.finite(unique(data$SSD))])
  SSDSall <- vector(mode="list",length=length(SSDS.all))
  SSDSall <- lapply(SSDSall,function(x){numeric(0)})
  names(SSDSall) <- SSDS.all
  SSDSallc <- SSDSall
  cat("Processing subject\n")
  for(s in levels(data$s)) {
    cat(paste(s," "))
    if (ignore.block.break) trials <- TRIALS[data$s==s] 
    temp <- data[data$s==s,]
    temp <- temp[!is.na(temp$RT),]               #remove all NAs
    rows<-nrow(temp)
    # make sure that dataset starts with at least one go before first stop
    start<-pmax(min(which(temp$SS=="SS")),2) 
    temp<-temp[start:rows,]
    pairs = data.frame("SRRT"=NA,"GORT"=NA,"SSD.SSRT"=NA,"SSD.GORT"=NA)
    for(i in 1:nrow(temp)){
      signal.respond <- temp[i,"SS"]=="SS" & !is.na(temp[i,"RT"])
      back <- NA
      # Always get go for stop even in runs
      if( signal.respond & !only.first.stop ){
        back = 1
        while(temp[i-back,"SS"]=="SS") back=back+1
      } else if ((i>1) && 
          (signal.respond & only.first.stop & (temp[i-1,"SS"]!="SS"))) back = 1
      if ( !is.na(back) ) {
        if (ignore.block.break) {
          ok <- TRIALS[i]>TRIALS[i-back] 
        } else ok=TRUE     
      if (ok & !is.na(back))  
        pairs[i,] = cbind(temp[i,"RT"],temp[i-back,"RT"],
                          temp[i,"SSD"],temp[i-back,"SSD"])
      }
    }
    pairs<-pairs[!is.na(pairs[,1]),]
    pairs<-cbind(pairs$SRRT-pairs$GORT,pairs[,3:4])
    colnames(pairs)<-c("Diff","SSD_SSRT","SSD_GORT")
    for (i in as.character(unique(pairs[,"SSD_SSRT"]))) {
      SSDSall[[i]] <- c(SSDSall[[i]],pairs[,"Diff"])
      SSDSallc[[i]] <- c(SSDSallc[[i]],c(pairs[,"Diff"]-mean(pairs[,"Diff"])))
    }
    violations<-tapply(pairs$Diff,c(SSD=pairs$SSD_SSRT),mean)
    n<-tapply(pairs$Diff,c(SSD=pairs$SSD_SSRT),length)
    SD<-tapply(pairs$Diff,c(SSD=pairs$SSD_SSRT),sd)
    independence[[s]]$SS_G <- violations
    independence[[s]]$n <- n
    independence[[s]]$SD <- SD
  }
  cat("\n")
  SSDSall <- SSDSall[lapply(SSDSall,length)!=0]
  SSDSallc <- SSDSallc[lapply(SSDSallc,length)!=0]
  list(independence=independence,SSDSall=SSDSall,SSDSallc=SSDSallc)
}


bissett.plot <- function(violation.summary,ylim=NA,main="",nMIN=0) 
  # Plots each subject on same graph as lines with no indication of uncertianty
  # Violation summary: list produce by summary.violation
  # ylim: filled in automatically if not NA
  # nMIN: minimum number of observations above which to plot
  # main: Title for plot
{
  independence <- violation.summary$independence
  SSDS <- lapply(independence,function(x){as.numeric(names(x$SS_G))})
  minSSD <- min(unlist(lapply(SSDS,min)))
  maxSSD <- max(unlist(lapply(SSDS,max)))
  if (any(is.na(ylim))) ylim <- 
      c(min(unlist(lapply(independence,function(x){min(x$SS_G[x$n>nMIN])}))),
        max(unlist(lapply(independence,function(x){max(x$SS_G[x$n>nMIN])}))))
  plot(as.numeric(names(independence[[1]][[1]])),independence[[1]][[1]],xlab="SSD",
     ylab="Signal-respond RT -GORT",type="l",lwd=1,ylim=ylim,xlim=c(minSSD,maxSSD),
     col="white",xaxt="n",main=main)
  axis(1,round(seq(minSSD,maxSSD,by=.05),3))
  abline(h=0,lty=2,lwd=4)
  for(i in 1:length(independence)){
    ind<-which(independence[[i]][[2]]>nMIN)
    lines(as.numeric(names(independence[[i]][[1]][ind])),independence[[i]][[1]][ind],
          lwd=2,col=i)
  }  
}


aggregate.plot <- function(violation.summary,ci=95,partial.subject=TRUE,
                           main="",len=.05)
  # Plots mean of all trials aggregated (so each subjects influence is 
  # properly weighted by the number of trials they contribute)
  # Violation summary: list produce by summary.violation
  # ci confidnce interval width as a percent
  # partial.subject = TRUE uses deviations from subject mean (i.e., partials
  #   out subject effect) to get error bars, otherwise just uses SD of same 
  #   data used to get mean.
  # len: length of error bars
  # main: graph title
{
  SSDSall <- violation.summary$SSDSall  
  SSDSallc <- violation.summary$SSDSallc  
  y <- unlist(lapply(SSDSall,mean))
  n <- unlist(lapply(SSDSall,length))
  if (partial.subject) SD <- unlist(lapply(SSDSall,sd)) else
    SD <- unlist(lapply(SSDSallc,sd))
  ssds <- as.numeric(names(y))
  half <- -qt((1-ci/100)/2,n-1)*SD/sqrt(n)
  hi <- y + half
  lo <- y - half
  ylim <- c(min(lo),max(hi))
  plot(ssds,y,xlab="SSD",ylab="Signal-respond RT -GO RT",ylim=ylim,
       xaxt="n",main=main)
  axis(1,as.numeric(names(y)))
  abline(h=0,lty=2)
  for (k in names(lo)) {
    arrows(as.numeric(k),y[k],as.numeric(k),hi[k],angle=90,len=len)
    arrows(as.numeric(k),y[k],as.numeric(k),lo[k],angle=90,len=len)
  }
}



individual.plots <- function(violation.summary,ci=95,len=.05,nMIN=2,
                             layout=c(2,4),main="",plot.it=TRUE) 
  # Individual plots with error bars
  # Violation summary: list produce by summary.violation
  # ci confidnce interval width as a percent
  # len: length of error bars
  # main: graph title
  # plot.it: plot graphs
  #
  # Invisible return: count per subject of nSSD (number of SSDs)
{
  if (main!="") main <- paste(main,": ",sep="")
  independence <- violation.summary$independence
  if (nMIN<2) {
    warning("nMIN must be at least 2")
    nMIN <- 2
  }  
  par(mfrow=layout,mar=c(5, 4, 4, 2) - 1)
  sig <- vector(mode="list",length=length(names(independence)))
  names(sig) <- names(independence)
  SSDS <- sort(unique(unlist(lapply(independence,function(x){
    as.numeric(names(x$SS_G))}))))
  sig.neg <- matrix(NA,nrow=length(names(independence)),ncol=length(SSDS))
  dimnames(sig.neg) <- list(s=names(independence),SSDS=SSDS)
  sig.pos <- sig.neg
  for (i in names(independence)) {
    ssds <- as.numeric(names(independence[[i]]$SS_G))
    half <- numeric(length(ssds))
    for (j in 1:length(ssds)) {
      if (independence[[i]]$n[j] < nMIN)
        half[j] <- NA else
        half[j] <- -qt((1-ci/100)/2,independence[[i]]$n[j]-1)
    }
    y <- independence[[i]]$SS_G
    lo <- y - half*independence[[i]]$SD/sqrt(independence[[i]]$n)
    hi <- y + half*independence[[i]]$SD/sqrt(independence[[i]]$n)
    okssd <- !is.na(lo)
    ylim <- c(min(c(lo[okssd],y)),max(c(hi[okssd],y)))
    if (plot.it) {
      plot(ssds,y,main=paste(main,i),xlab="SSD",ylab="Signal-respond RT -GORT",
         ylim=ylim,xaxt="n")
      axis(1,as.numeric(names(y)))
      abline(h=0,lty=2)
      for (k in names(lo[okssd])) {
        arrows(as.numeric(k),y[k],as.numeric(k),hi[k],angle=90,len=len)
        arrows(as.numeric(k),y[k],as.numeric(k),lo[k],angle=90,len=len)
      }
    }
    sig[[i]] <- c(nSSD=length(y),nOK=sum(okssd),
                  sig.neg=sum(hi[okssd]<0),sig.pos=sum(lo[okssd]>0))
    is.ssd <- names(y[okssd])
    sig.neg[i,is.ssd] <- hi[is.ssd]<0
    sig.pos[i,is.ssd] <- lo[is.ssd]>0
  }
  bad <- apply(sig.pos,2,function(x){all(is.na(x))})
  sig <- do.call(rbind,sig)
  row.names(sig) <- names(independence) 
  par(mar=c(5, 4, 4, 2) + 0.1)
  invisible(list(overall=data.frame(sig),
                 sig.neg=sig.neg[,!bad],sig.pos=sig.pos[,!bad]))
}

sig.summary <- function(sig,digits=2)
  # Summary of output from individual plots
{
  n=apply(sig$sig.neg,2,function(x){sum(!is.na(x))})
  n=c(n,Overall=sum(n))
  neg=apply(sig$sig.neg,2,sum,na.rm=TRUE)
  neg=c(neg,Overall=sum(neg))
  pos=apply(sig$sig.pos,2,sum,na.rm=TRUE)
  pos=c(pos,Overall=sum(pos))
  cbind(n=n,neg=neg,pcnt.neg=round(100*neg/n,digits),
            pos=pos,pcnt.pos=round(100*pos/n,digits))
}

