library(DiceKriging)
library(lhs)
library(mvtnorm)

lhd_mult_norm=function(seed,n,mu,c){
    d=length(mu)
    set.seed(seed)
    lhs_unif=maximinLHS(n,d,eps=1e-8)
    lhs_norm=qnorm(lhs_unif)
    cholc=chol(c)
    lhs_mvnorm=t(t(cholc)%*%t(lhs_norm)+mu)
    lhs_mvnorm
}

### flag forMichigan or Brazil (brs) location
brsflag=T
lhdfilepath=ifelse(brsflag,"latin_hypercube_brs.csv","latin_hypercube.csv")
lhdruns=read.csv(lhdfilepath, row.names = 1)

load("PFT_mean_cov.Rdata")
pftnum=ifelse(brsflag,4,7)
truemean=subPFT_summ[[pftnum]]$mean
truecov=subPFT_summ[[pftnum]]$cov

coords=log(lhdruns[-1,1:3])
GPPmat=lhdruns[-1,-(1:3)]

n.yr=ncol(GPPmat)
start.yr=as.numeric(substr(colnames(GPPmat)[1],2,5))-1
yrs=start.yr+1:n.yr

### 95% CI for each trait 
par.iqr=qnorm(0.975,truemean, sqrt(diag(truecov)))-
    qnorm(0.025,truemean, sqrt(diag(truecov)))
     
### spdep -- twwo values correspond to low and high spatial values respectively
### 0.05 means theta=iqr*0.05, 0.25 means theta=iqr*0.25 i.e. effective range=iqr*0.75
spdep=c(0.05,0.25)

settings=t(t(expand.grid(spdep,spdep,spdep))*(par.iqr))
settings=rbind(settings,rep(1e-10,3))   ### adding the actual estimates from the data analysis

set.seed(1)
fullsamples=rmvnorm(100000,truemean,truecov)
fullsamples=rbind(truemean,log(colMeans(exp(fullsamples))),fullsamples)
row.names(fullsamples)=NULL

#### first see how well LHD does in approximating trait density
traits=colnames(lhdruns)[1:3]
nvec=c(20,50,100,250,500)

svec=c(1,3,8,9) ## we pick settings corresponding to all low spdep, 2low 1 high, all high and all zero
setting.names=c("low","mixed","high","zero")

meanGPP=GPPmean=GPPlogmean=matrix(0,length(svec),n.yr)
for(s in 1:length(svec)){
  theta=settings[svec[s],]
  for(i in 1:n.yr){
    print(paste(s,i))
    y=GPPmat[,i]
    mu=mean(y)
    sd=sd(y)
    ind=which(y<=0) ### leaving out the zeros (for Brazil)
    newind=setdiff(1:length(y),ind)
    model=km(design=coords[newind,], response=y[newind],
      covtype="exp", coef.trend=mu, coef.cov=theta,
      coef.var=sd^2)
    pred=predict(model,fullsamples,type="UK",se.compute=F)
    meanGPP[s,i]=mean(pred$mean[-(1:2)])
    GPPlogmean[s,i]=pred$mean[1]
    GPPmean[s,i]=pred$mean[2]
    }
  }

# write.csv(cbind(settings,meanGPP),"meanGPP_emulations.csv")

# meanGPP=as.matrix(read.csv("meanGPP_emulations.csv"))

nseed=100

lhdGPPsamples=ranGPPsamples=array(NA,c(length(svec),n.yr,length(nvec),nseed))
    
for(j in 1:length(svec)){
    theta=settings[svec[j],]
      for(i in 1:n.yr){
        print(paste(svec[j],i))
        y=GPPmat[,i]
        mu=mean(y)
        sd=sd(y)
        ind=which(y<=0) ### leaving out the zeros (for Brazil)
        newind=setdiff(1:length(y),ind)
        model=km(design=coords[newind,], response=y[newind],
          covtype="exp", coef.trend=mu, coef.cov=theta,
          coef.var=sd^2)
        for(k in 1:length(nvec)) for(seed in 1:nseed){
            n=nvec[k]
            lhdsamples=lhd_mult_norm(seed,n,truemean,truecov)
            set.seed(seed)
            ransamples=rmvnorm(n,truemean,truecov)
            predlhd=predict(model,lhdsamples,type="UK",se.compute=F)
            lhdGPPsamples[j,i,k,seed]=mean(predlhd$mean)
            predran=predict(model,ransamples,type="UK",se.compute=F)
            ranGPPsamples[j,i,k,seed]=mean(predran$mean)
            }
          }
    }

savefile=ifelse(brsflag,"emulations_brs.Rdata","emulations.Rdata")
save.image(savefile)

for(j in 1:length(svec)){
    filename=ifelse(brsflag,paste0("figures/emulations_brs_",setting.names[j],".pdf"),
                    paste0("figures/emulations_",setting.names[j],".pdf"))
    pdf(filename,height=8,width=14)
    par(mfrow=c(ifelse(brsflag,3,4),4),
        oma = c(3,4,0,0) + 0.1,
        mar = c(0,2,2,2) + 0.1)
    for(i in 1:n.yr){
        nlist=rep(1:length(nvec),nseed)
        lhd=as.vector(lhdGPPsamples[j,i,,])
        ran=as.vector(ranGPPsamples[j,i,,])
        plot(nlist,lhd,col='dodgerblue',xaxt='n',
            xlab=ifelse(i > 12, "n", ""),ylab=ifelse(i %in% c(1,5,9,13),
            "GPP mean",""),main=start.yr+i,
            ylim=range(c(lhd,ran,meanGPP[j,i],GPPmean[j,i],GPPlogmean[j,i])))  
        points(nlist+0.1,ran,col='cyan3') 
        abline(h=meanGPP[j,i],col="orange",lwd=2)
        abline(h=0.975*meanGPP[j,i],col="orange",lwd=2,lty=2)
        abline(h=1.025*meanGPP[j,i],col="orange",lwd=2,lty=2)
        #abline(h=GPPlogmean[j,i],col="purple",lwd=2)
        abline(h=GPPmean[j,i],col="black",lwd=2)
        if(i > 12) axis(side=1, at=unique(nlist)+0.05,labels=nvec)
    }
    plot(1:10,1:10,col="white",xaxt='n',yaxt='n',xlab="",ylab="",main="",bty='n')
    # legend("center",c("true GPP mean","+- 2.5% of true GPP mean","GPP_at_log_trait_mean","GPP_at_trait_mean",
    #                   "estimated_GPP_mean_LHD","estimated_GPP_mean_Rand"),
    #                 col=c("orange","orange","purple","black","dodgerblue","cyan3"),
    #                 lwd=c(2,2,2,2,NA,NA),pch=c(NA,NA,NA,NA,1,1),bty='n',cex=1.6,lty=c(1,2,1,1,NA,NA))
    legend("center",c("true GPP mean","+- 2.5% of true GPP mean","GPP_at_trait_mean",
                      "estimated_GPP_mean_LHD","estimated_GPP_mean_Rand"),
           col=c("orange","orange","black","dodgerblue","cyan3"),
           lwd=c(2,2,2,NA,NA),pch=c(NA,NA,NA,1,1),bty='n',cex=1.6,lty=c(1,2,1,NA,NA))
    plot(1:10,1:10,col="white",xaxt='n',yaxt='n',xlab="",ylab="",main="",bty='n')
    dev.off()
}
graphics.off()

#####################################################################################
## repeating emulations for Brazil leaving the samples which produced zero GPP out ##
# brsflag=T

svec=c(1,3,8,9) ## we pick settings corresponding to all low spdep, 2low 1 high, all high and all zero
setting.names=c("low","mixed","high","zero")

#### using log mean trait vs mean log trait makes a huge difference
meanGPP=GPPmean=GPPlogmean=matrix(0,length(svec),n.yr)
for(s in 1:length(svec)){
  theta=settings[svec[s],]
  for(i in 1:n.yr){
    print(paste(s,i))
    y=GPPmat[,i]
    ind=which(y<=0) ### leaving out the zeros
    mu=mean(y)
    sd=sd(y)
    model=km(design=coords[-ind,], response=y[-ind],
             covtype="exp", coef.trend=mu, coef.cov=theta,
             coef.var=sd^2)
    pred=predict(model,fullsamples,type="UK",se.compute=F)
    meanGPP[s,i]=mean(pred$mean[-(1:2)])
    GPPlogmean[s,i]=pred$mean[1]
    GPPmean[s,i]=pred$mean[2]
  }
}

# write.csv(cbind(settings,meanGPP),"meanGPP_emulations.csv")

# meanGPP=as.matrix(read.csv("meanGPP_emulations.csv"))

nseed=100

lhdGPPsamples=ranGPPsamples=array(NA,c(length(svec),n.yr,length(nvec),nseed))

for(j in 1:length(svec)){
  theta=settings[svec[j],]
  for(i in 1:n.yr){
    print(paste(svec[j],i))
    y=GPPmat[,i]
    ind=which(y<=0) ### leaving out the zeros
    mu=mean(y)
    sd=sd(y)
    model=km(design=coords[-ind,], response=y[-ind],
             covtype="exp", coef.trend=mu, coef.cov=theta,
             coef.var=sd^2)
    for(k in 1:length(nvec)) for(seed in 1:nseed){
      n=nvec[k]
      lhdsamples=lhd_mult_norm(seed,n,truemean,truecov)
      set.seed(seed)
      ransamples=rmvnorm(n,truemean,truecov)
      predlhd=predict(model,lhdsamples,type="UK",se.compute=F)
      lhdGPPsamples[j,i,k,seed]=mean(predlhd$mean)
      predran=predict(model,ransamples,type="UK",se.compute=F)
      ranGPPsamples[j,i,k,seed]=mean(predran$mean)
    }
  }
}

savefile=ifelse(brsflag,"emulations_brs_wo_zeros.Rdata","emulations.Rdata")
save.image(savefile)

for(j in 1:length(svec)){
  filename=ifelse(brsflag,paste0("figures/emulations_brs_wo_zeros_",setting.names[j],".pdf"),
                  paste0("figures/emulations_",setting.names[j],".pdf"))
  pdf(filename,height=8,width=14)
  par(mfrow=c(ifelse(brsflag,3,4),4),
      oma = c(3,4,0,0) + 0.1,
      mar = c(0,2,2,2) + 0.1)
  for(i in 1:n.yr){
    nlist=rep(1:length(nvec),nseed)
    lhd=as.vector(lhdGPPsamples[j,i,,])
    ran=as.vector(ranGPPsamples[j,i,,])
    plot(nlist,lhd,col='dodgerblue',xaxt='n',
         xlab=ifelse(i > 12, "n", ""),ylab=ifelse(i %in% c(1,5,9,13),
         "GPP mean",""),main=start.yr+i,
         ylim=range(c(lhd,ran,meanGPP[j,i],GPPmean[j,i],GPPlogmean[j,i])))
    points(nlist+0.1,ran,col='cyan3') 
    abline(h=meanGPP[j,i],col="orange",lwd=2)
    abline(h=0.975*meanGPP[j,i],col="orange",lwd=2,lty=2)
    abline(h=1.025*meanGPP[j,i],col="orange",lwd=2,lty=2)
    #abline(h=GPPlogmean[j,i],col="black",lwd=2)
    abline(h=GPPmean[j,i],col="black",lwd=2)
    if(i > 12) axis(side=1, at=unique(nlist)+0.05,labels=nvec)
  }
  plot(1:10,1:10,col="white",xaxt='n',yaxt='n',xlab="",ylab="",main="",bty='n')
  legend("center",c(expression("GPP"[dist]),expression('+/- 2.5% of '*'GPP'[dist]),
		    "input trait mean","emulated output mean [LHD]","emulated output mean [MVN]"),
         col=c("orange","orange","black","dodgerblue","cyan3"),
         lwd=c(2,2,2,NA,NA),pch=c(NA,NA,NA,1,1),bty='n',cex=1.6,lty=c(1,2,1,NA,NA))
  plot(1:10,1:10,col="white",xaxt='n',yaxt='n',xlab="",ylab="",main="",bty='n')
  dev.off()
}
graphics.off()

