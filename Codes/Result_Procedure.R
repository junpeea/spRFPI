
library(xtable)
library(ggplot2)
library(reshape2)
library(dplyr)


table1.outbut.mat = matrix(NA,nrow=3*8,ncol=2*3)
table1.output.mat = matrix(NA,nrow=2*8,ncol=3*3)

load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f1).RData")
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220405_nspRFPI(f1).RData")
output.cpr.mat = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)
output.len.mat = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)
output.ais.mat = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)

table1.outbut.mat[1:7,1:3] <-
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2),"(",round(apply(output.cpr.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2),"(",round(apply(output.len.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2),"(",round(apply(output.ais.mat,2,function(w){sd(w,na.rm=T)}),2),")")))

table1.output.mat[1:7,1:3] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2))))

png(filename=paste0(Sys.Date(),"LINEAR_Cond_Visual.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(1,3,1,4))
# sel.ind = seq(2,50,by=1)
sel.ind = seq(2,50,by=2)
boxplot(lendmat[,sel.ind],ylim=c(2,2*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25),cex.axis=1.6)
par(new=T);boxplot(lendmat2[,sel.ind],ylim=c(2,2*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=T);boxplot(lendmat3[,sel.ind],ylim=c(2,2*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=TRUE)
plot(colMeans(condmat[,sel.ind],na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
lines(colMeans(condmat2[,sel.ind],na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0,1),pch=1)
lines(colMeans(condmat3[,sel.ind],na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0,1),pch=1)
axis(4,at=c(0.8,0.9,1.0),las=2,cex.axis=1.6)
# par(mfrow=c(1,1),mar=c(1,3,3,3))
# boxplot(output.lendmat[,c(1,6,11,16,2,7,12,16,3,8,13,16,4,9,14,16,5,10,15)],ylim=c(0,2.5*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=rep(col.ind,5))
# par(new=TRUE)
# plot(colMeans(condmat,na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0,1),pch=1,xlab="",ylab="");abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
# lines(colMeans(condmat2,na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0,1),pch=1)
# lines(colMeans(condmat3,na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0,1),pch=1)
# abline(v=1,col=1,lty=2,xlab="",ylab="");abline(v=50,col=1,lty=2,xlab="",ylab="")
# abline(v=10.8,col=1,lty=2,xlab="",ylab="");abline(v=20.6,col=1,lty=2,xlab="",ylab="")
# abline(v=30.4,col=1,lty=2,xlab="",ylab="");abline(v=40.2,col=1,lty=2,xlab="",ylab="")
# axis(3,at=seq(1,50,length.out=6),label=seq(-3,3,length.out=6))
# axis(4,at=c(0.8,0.9,1.0),las=2)
# text(x=6,y=0,labels=c("X1=-3"));text(x=16,y=0,labels=c("X1=-1.5"))
# text(x=26,y=0,labels=c("X1=0"));text(x=36,y=0,labels=c("X1=1.5"));text(x=46,y=0,labels=c("X1=3"))
dev.off()

load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f2)_homo.RData")
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220405_nspRFPI(f2)_homo.RData")
output.cpr.mat = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)
output.len.mat = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)
output.ais.mat = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)

table1.outbut.mat[1:7,4:6] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2))))

table1.output.mat[9:15,1:3] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2))))


png(filename=paste0(Sys.Date(),"SINUSOIDAL(Homo)_Cond_Visual.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(1,3,1,4))
# sel.ind = seq(2,50,by=1)
sel.ind = seq(2,50,by=2)
boxplot(lendmat[,sel.ind],ylim=c(2,1.5*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25),cex.axis=1.6)
par(new=T);boxplot(lendmat2[,sel.ind],ylim=c(2,1.5*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=T);boxplot(lendmat3[,sel.ind],ylim=c(2,1.5*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=TRUE)
plot(colMeans(condmat[,sel.ind],na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
lines(colMeans(condmat2[,sel.ind],na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0,1),pch=1)
lines(colMeans(condmat3[,sel.ind],na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0,1),pch=1)
axis(4,at=c(0.8,0.9,1.0),las=2,cex.axis=1.6)
dev.off()


load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f2)_heavy.RData")
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220405_nspRFPI(f2)_heavy.RData")
output.cpr.mat = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)
output.len.mat = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)
output.ais.mat = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)

table1.outbut.mat[9:15,1:3] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2),"(",round(apply(output.cpr.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2),"(",round(apply(output.len.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2),"(",round(apply(output.ais.mat,2,function(w){sd(w,na.rm=T)}),2),")")))

table1.output.mat[9:15,4:6] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2))))


png(filename=paste0(Sys.Date(),"SINUSOIDAL(HEAVY)_Cond_Visual.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(1,3,1,4))
# sel.ind = seq(2,50,by=1)
sel.ind = seq(2,50,by=2)
boxplot(lendmat[,sel.ind],ylim=c(2,1.5*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25),cex.axis=1.6)
par(new=T);boxplot(lendmat2[,sel.ind],ylim=c(2,1.5*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=T);boxplot(lendmat3[,sel.ind],ylim=c(2,1.5*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=TRUE)
plot(colMeans(condmat[,sel.ind],na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
lines(colMeans(condmat2[,sel.ind],na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0,1),pch=1)
lines(colMeans(condmat3[,sel.ind],na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0,1),pch=1)
axis(4,at=c(0.8,0.9,1.0),las=2,cex.axis=1.6)
dev.off()


load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f2)_hetero.RData")
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220405_nspRFPI(f2)_hetero.RData")
output.cpr.mat = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)
output.len.mat = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)
output.ais.mat = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)

table1.outbut.mat[9:15,4:6] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2),"(",round(apply(output.cpr.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2),"(",round(apply(output.len.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2),"(",round(apply(output.ais.mat,2,function(w){sd(w,na.rm=T)}),2),")")))

table1.output.mat[9:15,7:9] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2))))


png(filename=paste0(Sys.Date(),"SINUSOIDAL(HETERO)_Cond_Visual.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(1,3,1,4))
# sel.ind = seq(2,50,by=1)
sel.ind = seq(2,50,by=2)
boxplot(lendmat[,sel.ind],ylim=c(2,1.5*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25),cex.axis=1.6)
par(new=T);boxplot(lendmat2[,sel.ind],ylim=c(2,1.5*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=T);boxplot(lendmat3[,sel.ind],ylim=c(2,1.5*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=TRUE)
plot(colMeans(condmat[,sel.ind],na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
lines(colMeans(condmat2[,sel.ind],na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0,1),pch=1)
lines(colMeans(condmat3[,sel.ind],na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0,1),pch=1)
axis(4,at=c(0.8,0.9,1.0),las=2,cex.axis=1.6)
dev.off()


load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f3).RData")
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220405_nspRFPI(f3).RData")
output.cpr.mat = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)
output.len.mat = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)
output.ais.mat = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)

table1.outbut.mat[17:23,1:3] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2),"(",round(apply(output.cpr.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2),"(",round(apply(output.len.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2),"(",round(apply(output.ais.mat,2,function(w){sd(w,na.rm=T)}),2),")")))

table1.output.mat[1:7,4:6] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2))))

png(filename=paste0(Sys.Date(),"STEP_Cond_Visual.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(1,3,1,4))
# sel.ind = seq(2,50,by=1)
sel.ind = seq(2,50,by=2)
boxplot(lendmat[,sel.ind],ylim=c(2,1.2*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25),cex.axis=1.6)
par(new=T);boxplot(lendmat2[,sel.ind],ylim=c(2,1.2*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=T);boxplot(lendmat3[,sel.ind],ylim=c(2,1.2*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=TRUE)
plot(colMeans(condmat[,sel.ind],na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
lines(colMeans(condmat2[,sel.ind],na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0,1),pch=1)
lines(colMeans(condmat3[,sel.ind],na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0,1),pch=1)
axis(4,at=c(0.8,0.9,1.0),las=2,cex.axis=1.6)
dev.off()


load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f4).RData")
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220405_nspRFPI(f4).RData")
output.cpr.mat = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)[1:96,]
output.len.mat = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)[1:96,]
output.ais.mat = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)[1:96,]

table1.outbut.mat[17:23,4:6] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2),"(",round(apply(output.cpr.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2),"(",round(apply(output.len.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2),"(",round(apply(output.ais.mat,2,function(w){sd(w,na.rm=T)}),2),")")))

table1.output.mat[1:7,7:9] <- 
  cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.len.mat,na.rm=T),2))),
        rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2))))

png(filename=paste0(Sys.Date(),"FRIEDMAN_Cond_Visual.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(1,3,1,4))
# sel.ind = seq(2,50,by=1)
sel.ind = seq(2,50,by=2)
boxplot(lendmat[,sel.ind],ylim=c(2,0.75*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25),cex.axis=1.6)
par(new=T);boxplot(lendmat2[,sel.ind],ylim=c(2,0.75*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=T);boxplot(lendmat3[,sel.ind],ylim=c(2,0.75*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="",cex.axis=1.6)
par(new=TRUE)
plot(colMeans(condmat[,sel.ind],na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
lines(colMeans(condmat2[,sel.ind],na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0,1),pch=1)
lines(colMeans(condmat3[,sel.ind],na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0,1),pch=1)
axis(4,at=c(0.8,0.9,1.0),las=2,cex.axis=1.6)
dev.off()

# xtable(table1.outbut.mat)
xtable(table1.output.mat)










# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220330_spRFPI(f1).RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f1).RData")
output.alpha.cpr.mar = matrix(NA,length(alphavec),7)
for(alphaind in 1:length(alphavec)){
  output.alpha.cpr.mar[alphaind,] <- output.cpr.mar.list[[alphaind]] %>% colMeans
}
png(filename=paste0(Sys.Date(),"LINEAR_Alpha_Mar.png"),width=300,height=300)
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(alphavec,1-output.alpha.cpr.mar[,1],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=4)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,3],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=3)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,4],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=2)
abline(a=0,b=1,lty=2)
legend("bottomright",col=c(4,3,2),lty=c(1,1,1),legend=c("OOBGK","OOBW","LSCP"))
dev.off()

# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220330_spRFPI(f2)_homo.RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f2)_homo.RData")
output.alpha.cpr.mar = matrix(NA,length(alphavec),7)
for(alphaind in 1:length(alphavec)){
  output.alpha.cpr.mar[alphaind,] <- output.cpr.mar.list[[alphaind]] %>% colMeans
}
png(filename=paste0(Sys.Date(),"SINUSOIDAL(HOMO)_Alpha_Mar.png"),width=300,height=300)
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(alphavec,1-output.alpha.cpr.mar[,1],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=4)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,3],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=3)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,4],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=2)
abline(a=0,b=1,lty=2)
legend("bottomright",col=c(4,3,2),lty=c(1,1,1),legend=c("OOBGK","OOBW","LSCP"))
dev.off()

# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220330_spRFPI(f2)_heavy.RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f2)_heavy.RData")
output.alpha.cpr.mar = matrix(NA,length(alphavec),7)
for(alphaind in 1:length(alphavec)){
  output.alpha.cpr.mar[alphaind,] <- output.cpr.mar.list[[alphaind]] %>% colMeans
}
png(filename=paste0(Sys.Date(),"SINUSOIDAL(HEAVY)_Alpha_Mar.png"),width=300,height=300)
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(alphavec,1-output.alpha.cpr.mar[,1],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=4)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,3],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=3)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,4],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=2)
abline(a=0,b=1,lty=2)
legend("bottomright",col=c(4,3,2),lty=c(1,1,1),legend=c("OOBGK","OOBW","LSCP"))
dev.off()

# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220330_spRFPI(f2)_hetero.RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f2)_hetero.RData")
output.alpha.cpr.mar = matrix(NA,length(alphavec),7)
for(alphaind in 1:length(alphavec)){
  output.alpha.cpr.mar[alphaind,] <- output.cpr.mar.list[[alphaind]] %>% colMeans
}
png(filename=paste0(Sys.Date(),"SINUSOIDAL(HETERO)_Alpha_Mar.png"),width=300,height=300)
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(alphavec,1-output.alpha.cpr.mar[,1],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=4)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,3],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=3)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,4],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=2)
abline(a=0,b=1,lty=2)
legend("bottomright",col=c(4,3,2),lty=c(1,1,1),legend=c("OOBGK","OOBW","LSCP"))
dev.off()

# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220330_spRFPI(f3).RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f3).RData")
output.alpha.cpr.mar = matrix(NA,length(alphavec),7)
for(alphaind in 1:length(alphavec)){
  output.alpha.cpr.mar[alphaind,] <- output.cpr.mar.list[[alphaind]] %>% colMeans
}
png(filename=paste0(Sys.Date(),"STEP_Alpha_Mar.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(alphavec,1-output.alpha.cpr.mar[,1],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=4)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,3],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=3)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,4],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=2)
abline(a=0,b=1,lty=2)
legend("bottomright",col=c(4,3,2),lty=c(1,1,1),legend=c("OOBGK","OOBW","LSCP"))
dev.off()

# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220330_spRFPI(f4).RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f4).RData")
output.alpha.cpr.mar = matrix(NA,length(alphavec),7)
for(alphaind in 1:length(alphavec)){
  output.alpha.cpr.mar[alphaind,] <- output.cpr.mar.list[[alphaind]] %>% colMeans
}
png(filename=paste0(Sys.Date(),"FRIEDMAN_Alpha_Mar.png"),width=360,height=520)
par(mfrow=c(1,1),mar=c(3,3,1,1))
plot(alphavec,1-output.alpha.cpr.mar[,1],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=4)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,3],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=3)
par(new=T);plot(alphavec,1-output.alpha.cpr.mar[,4],type='o',xlim=c(0,0.21),ylim=c(0,0.21),xlab="nominal",ylab="estimated",col=2)
abline(a=0,b=1,lty=2)
legend("bottomright",col=c(4,3,2),lty=c(1,1,1),legend=c("OOBGK","OOBW","LSCP"))
dev.off()


# colnames(output.alpha.cpr.mar) <- names(output.cpr.mar.list[[1]])
# mydat = output.alpha.cpr.mar %>% melt(id.vars =c("id")) %>% filter(Var2 %in% c("oobgk.cpr.mar","oobw.cpr.mar","lscp.cpr.mar"))
# mydat$Var1 = rep(rev(alphavec),3)
# colnames(mydat) <- c("id","Type","cpr")
# ggplot(data=mydat, aes(x=id, y=cpr, group=Type)) +
#   geom_line()+
#   geom_point()

cpr.by.epsilon.mat = matrix(NA,nrow = 10, ncol=11)
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f2)_homo.RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f2)_homo.RData")
cpr.by.epsilon.mat[,1:3] <- lapply(c(1:length(alphavec)),function(w){output.cpr.mar.list[[w]] %>% colMeans}) %>% unlist %>% matrix(ncol=7,byrow=TRUE) %>% .[,c(1,3,4)]
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f2)_heavy.RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f2)_heavy.RData")
cpr.by.epsilon.mat[,5:7] <- lapply(c(1:length(alphavec)),function(w){output.cpr.mar.list[[w]] %>% colMeans}) %>% unlist %>% matrix(ncol=7,byrow=TRUE) %>% .[,c(1,3,4)]
# load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220329_spRFPI(f2)_hetero.RData")
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/20220401_spRFPI(f2)_hetero.RData")
cpr.by.epsilon.mat[,9:11] <- lapply(c(1:length(alphavec)),function(w){output.cpr.mar.list[[w]] %>% colMeans}) %>% unlist %>% matrix(ncol=7,byrow=TRUE) %>% .[,c(1,3,4)]
png(filename=paste0(Sys.Date(),"SINUSOIDAL_Box_Mar.png"),width=450,height=450)
par(mfrow=c(1,1),mar=c(1,3,1,1))
boxplot(cpr.by.epsilon.mat,xaxt="n",ylim=c(0.5,1),col=c(4,3,2,NA),cex.axis=1.6);abline(h=1-0.10,lty=2);abline(v=4,lty=2);abline(v=8,lty=2)
text(x=2,y=0.64,labels=c("HOMO"),cex=2);text(x=6,y=0.64,labels=c("HEAVY"),cex=2);text(x=10,y=0.64,labels=c("HETERO"),cex=2)
legend("bottomright",col=c(4,3,2),lty=c(1,1,1),legend=c("OOBGK","OOBW","LSCP"),cex=1.5)
dev.off()


