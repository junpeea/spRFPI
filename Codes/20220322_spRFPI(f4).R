
rm(list = ls())
library(randomForest);
library(RandomFields)
library(rfinterval)
library(forestError)
library(plot3D)
library(plotrix)
library(dplyr)
library(expm)
library(spatstat)
library(RandomForestsGLS)
# devtools::install_github("mhuiying/lscp")
library(scp)
library(scales)
library(xtable)
# source('RFGLS_support_YBJ.R')
source('RFGLS_support2_YBJ.R')
source('RFGLS_support3_YBJ.R')
################################

nsim = 100;
ntrain = 200;
ntest = 50;
nsample = ntrain + ntest;
sd2 = sqrt(0.6); # marginal variance for spatial error
sd3 = sqrt(3.0 - sd2^2); # marginal variance (nugget) for non-spatial error
nu0 = 0.5; phi = 1
ntree.gls <- 0.25*ntrain; ntree <- 0.25*ntrain
nx = 20 ;           # number of dimensions of covariates
ns = 20 # nodesize
eta = 25 # smooth parameter

#########

# Linear
f1 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  EY = (x[1]+x[2])
  return(EY)
}
# Nonlinear (default)
f2 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  # EY = 20*exp(-abs(x[1])-abs(x[2]))
  EY = 5*sin(pi*x[1]/3)
  return(EY)
}
# Step
f3 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  # EY = ((x[1] >= 0)-(x[1] < 0))+x[2]
  EY = 3*(((x[1] >= 0)-(x[1] < 0))+x[2])
  return(EY)
}
#Friedman
f4 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  EY = 0.005*(10*(sin(pi*x[1]*x[2]) + 20*(x[3]-0.5)^2 + 10*x[4] + 5*x[5]))
  return(EY)
}
#########

mfunc = f4

#########

simulate_data = function(seed,ntrain,ntest,nx,nu0,sd2,sd3,mfunc){
  
  ### Simulate dataset
  # set.seed(seed+1); X = matrix(runif((ntrain+ntest) * nx,min=-3,max=3), ncol = nx)
  # set.seed(seed+2); Loc = cbind(runif((ntrain+ntest),0,1),runif((ntrain+ntest),0,1))
  set.seed(seed+10); X.train = matrix(runif(ntrain * nx,min=-3.5,max=3.5), ncol = nx)
  set.seed(seed+20); X.test  = matrix(runif(ntest  * nx,min=-3.5,max=3.5), ncol = nx)
  X = rbind(X.train, X.test)
  set.seed(seed+30); Loc = cbind(runif((ntrain+ntest),0,1),runif((ntrain+ntest),0,1))
  ### Endow zero/one vector for conditional analysis
  X[(ntrain+1),] <- rep(0,nx); Loc[(ntrain+1),] <- rep(0,2)
  X[(ntrain+2):(ntrain+ntest),1] <- seq(-3,3,length.out=(ntest-1))
  X[(ntrain+2):(ntrain+ntest),2] <- rep(1,(ntest-1))
  # combine spatial information into the covariates
  X.train = cbind(X[1:ntrain,],Loc[1:ntrain,]); colnames(X.train) = c(paste("X-", 1:nx, sep = ""),paste0("s-",1:2))
  X.test  = cbind(X[(ntrain+1):(ntrain+ntest),],Loc[(ntrain+1):(ntrain+ntest),]); colnames(X.test) = c(paste("X-", 1:nx, sep = ""),paste0("s-",1:2))
  if(sd2 ==0){
    set.seed(seed+3); Y = apply(X, 1, mfunc) + rnorm((ntrain+ntest), 0, sd3)
  }else{
    ### Spatial Error generation
    model.option = RMmatern(nu=nu0,var=sd2^2,scale=1)
    epsi= RandomFields::RFsimulate(model=model.option,x=Loc[,1],y=Loc[,2])$variable1
    #homoscadesticity
    set.seed(seed+3); Y = apply(X, 1, mfunc) + epsi + rnorm((ntrain+ntest), 0, sd3)
    #heavy-tail
    # set.seed(seed+3); Y = apply(X, 1, mfunc) + epsi + rt(ntrain, 3)/sqrt(3)
    #heteroscadesticity
    # set.seed(seed+3); Y = apply(X, 1, mfunc) + epsi + rnorm((ntrain+ntest), 0, (sd3+abs(apply(X, 1, mfunc))/4))
  }
  Y.train = Y[1:ntrain]; Y.test = Y[(ntrain+1):(ntrain+ntest)]
  s.train = Loc[1:ntrain,]; s.test = Loc[(ntrain+1):(ntrain+ntest),]
  train_data = data.frame(Y.train,X.train); colnames(train_data)[1] <- c("y")
  test_data  = data.frame(Y.test,X.test);   colnames(test_data)[1] <- c("y")
  
  output = list(train_data,test_data,X.train,Y.train,s.train,X.test,Y.test,s.test)
  return(output)
  
}

### Determine bandwidth parameter : mtry, adj
# We chose the set of bandwidth parameters to minimize the average interval score (Gneiting and Raftery 2007)

#########

seed0 = 123
alpha = 0.10

oobg.cpr.mar   = rep(NA,nsim);oobg.len.mar   = rep(NA,nsim);oobg.ais.mar   = rep(NA,nsim)
oobgk.cpr.mar   = rep(NA,nsim);oobgk.len.mar   = rep(NA,nsim);oobgk.ais.mar   = rep(NA,nsim)
oobw.cpr.mar   = rep(NA,nsim);oobw.len.mar   = rep(NA,nsim);oobw.ais.mar   = rep(NA,nsim);
oob.cpr.mar    = rep(NA,nsim);oob.len.mar    = rep(NA,nsim);oob.ais.mar    = rep(NA,nsim);
sc.cpr.mar     = rep(NA,nsim);sc.len.mar     = rep(NA,nsim);sc.ais.mar     = rep(NA,nsim);
qrf.cpr.mar    = rep(NA,nsim);qrf.len.mar    = rep(NA,nsim);qrf.ais.mar    = rep(NA,nsim);
lscp.cpr.mar    = rep(NA,nsim);lscp.len.mar  = rep(NA,nsim);lscp.ais.mar     = rep(NA,nsim);

OutTable = matrix(NA,nrow=4,ncol=7)
condmat = matrix(NA,nrow=nsim,ncol=ntest)
condmat2= matrix(NA,nrow=nsim,ncol=ntest)
condmat3= matrix(NA,nrow=nsim,ncol=ntest)
lendmat = matrix(NA,nrow=nsim,ncol=ntest)
lendmat2= matrix(NA,nrow=nsim,ncol=ntest)
lendmat3= matrix(NA,nrow=nsim,ncol=ntest)

m.cand = seq(12,20,by = 1)
adj.vec = seq(1,5,by = 0.25)
alphavec = rev(seq(0.02,0.2,by=0.02))
output.cpr.mar.list = list()
output.len.mar.list = list()
output.ais.mar.list = list()

# for(alphaind in 1:length(alphavec)){
  
  # alpha = alphavec[alphaind]
  
  for (k in 1:nsim){
    
    seed <- seed0 + k
    simdata = simulate_data(seed,ntrain,ntest,nx,nu0,sd2,sd3,mfunc)
    train_data = simdata[[1]]; test_data = simdata[[2]]
    X.train    = simdata[[3]]; X.test    = simdata[[6]]
    Y.train    = simdata[[4]]; Y.test    = simdata[[7]]
    s.train    = simdata[[5]]; s.test    = simdata[[8]]
    set.seed(seed);sample.0 = sample(c(1:ntrain),(ntrain-50));sample.v = setdiff(c(1:ntrain),sample.0)
    X.train.0 = X.train[sample.0,];X.train.v = X.train[sample.v,]
    Y.train.0 = Y.train[sample.0]; Y.train.v = Y.train[sample.v]
    s.train.0 = s.train[sample.0,];s.train.v = s.train[sample.v,]
    
    ###########################
    AIS.mat = matrix(NA,nrow=length(m.cand), ncol=length(adj.vec))
    CpabsBias.mat = matrix(NA,nrow=length(m.cand), ncol=length(adj.vec))
    PMSPE.vec  = rep(NA,length(m.cand))
    PMSPE2.vec = rep(NA,length(m.cand))
    for(l in 1:length(m.cand)){
      rf.gls <- RFGLS_estimate_spatial(s.train.0, y = Y.train.0, X = X.train.0, ntree = ntree.gls, param_estimate = TRUE, nthsize = ns, h=16, mtry = m.cand[l])
      oobgk.output.list <- quantGLS(rf.gls, X.train.0, X.train.v, s.train.v, alpha = alpha, adj.vec)
      PMSPE.vec[l] = mean((Y.train.v - oobgk.output.list[[1]]$pred)^2)
      AIS.mat[l,] = lapply(c(1:length(oobgk.output.list)), function(t){
        mean( (oobgk.output.list[[t]][,5] - oobgk.output.list[[t]][,4]) + (2/alpha) * (sapply(oobgk.output.list[[t]][,4] - Y.train.v, function(w){max(w,0)}) + sapply(Y.train.v - oobgk.output.list[[t]][,5], function(w){max(w,0)})), na.rm=T )
      }) %>% unlist
      CpabsBias.mat[l,] = lapply(c(1:length(oobgk.output.list)), function(t){
        abs((1-alpha) - mean(oobgk.output.list[[t]][,4] < Y.train.v & oobgk.output.list[[t]][,5] > Y.train.v))
      }) %>% unlist
      rf <- randomForest::randomForest(X.train.0, Y.train.0, keep.inbag = TRUE, proximity = TRUE, mtry = m.cand[l])
      train_nodes <- findOOBErrors(rf, X.train.0)
      oobw.output <- quantForestError(rf, X.train.0, X.train.v, train_nodes = train_nodes, alpha = alpha)$estimate
      PMSPE2.vec[l] = mean((Y.train.v - oobw.output$pred)^2)
      
    }
    mtry1 = m.cand[order(PMSPE.vec)[1]];mtry1
    mtry2 = m.cand[order(PMSPE2.vec)[1]];mtry2
    adjk   = adj.vec[order(AIS.mat[order(PMSPE.vec)[1],])[1]];adjk
    
    rf.gls <- RFGLS_estimate_spatial(s.train, y = Y.train, X = X.train, ntree = ntree.gls, param_estimate = TRUE, nthsize = ns, h=16, mtry = mtry1)
    rf.oracle.gls <- RFGLS_estimate_spatial(s.train, y = Y.train, X = X.train, ntree = ntree.gls, param_estimate = FALSE, nthsize = ns, h=16, mtry = mtry1,
                                            sigma.sq = sd2^2, tau.sq = sd3^2,  phi = 1, nu = 0.5)
    oobgk.output <- quantGLSForestError(rf.gls, X.train, X.test, alpha = alpha, kernel=TRUE, adjust = adjk)
    oobg.output <- quantGLSForestError(rf.oracle.gls, X.train, X.test, alpha = alpha, kernel=TRUE, adjust = adjk)
    
    # fit random forest to the training data
    rf <- randomForest::randomForest(X.train, Y.train, keep.inbag = TRUE, proximity = TRUE, mtry = mtry2)
    train_nodes <- findOOBErrors(rf, X.train)
    oobw.output <- quantForestError(rf, X.train, X.test,
                                    train_nodes = train_nodes, alpha = alpha)
    
    s0 = s.test; s = s.train; Y = matrix(Y.train - rf$predicted, ncol=1)
    rf.pred = predict(rf,X.test)
    bins     = seq(0.01,0.2,0.01)
    thetaHat = get_theta(s,Y,dists=bins)
    thetaHat[4] <- min(thetaHat[4],10)
    lscp.out = scp(s0,s,Y,alpha = alpha,thetaHat = thetaHat, global = FALSE, eta=eta)
    lscp.out = data.frame(lower=lscp.out[,1]+rf.pred,upper=lscp.out[,2]+rf.pred)
    
    alter.output <- rfinterval(y~., train_data = train_data, test_data = test_data,
                               method = c("oob", "split-conformal", "quantreg"),
                               symmetry = TRUE ,alpha = alpha, params_ranger=list(mtry=mtry2))
    
    oobg.cpr.mar[k]   = mean(oobg.output[,4] < Y.test &  oobg.output[,5] > Y.test)
    oobgk.cpr.mar[k]   = mean(oobgk.output[,4] < Y.test & oobgk.output[,5] > Y.test)
    oobw.cpr.mar[k]   = mean(oobw.output$estimate[,4] < Y.test &  oobw.output$estimate[,5] > Y.test)
    oob.cpr.mar[k]    = mean(alter.output$oob_interval$lo < Y.test & alter.output$oob_interval$up > Y.test)
    sc.cpr.mar[k]     = mean(alter.output$sc_interval$lo < Y.test & alter.output$sc_interval$up > Y.test)
    qrf.cpr.mar[k]    = mean(alter.output$quantreg_interval$lo < Y.test & alter.output$quantreg_interval$up > Y.test)
    lscp.cpr.mar[k]    = mean(lscp.out$lower < Y.test & lscp.out$upper > Y.test)
    
    oobg.len.mar[k]   = mean(oobg.output[,5] - oobg.output[,4])
    oobgk.len.mar[k]  = mean(oobgk.output[,5] - oobgk.output[,4])
    oobw.len.mar[k]   = mean(oobw.output$estimate[,5] - oobw.output$estimate[,4])
    oob.len.mar[k]    = mean(alter.output$oob_interval$up-alter.output$oob_interval$lo)
    sc.len.mar[k]     = mean(alter.output$sc_interval$up-alter.output$sc_interval$lo)
    qrf.len.mar[k]    = mean(alter.output$quantreg_interval$up-alter.output$quantreg_interval$lo)
    lscp.len.mar[k]    = mean(lscp.out$upper-lscp.out$lower)
    
    condmat[k,]  <- oobgk.output[,4] < Y.test & oobgk.output[,5] > Y.test
    condmat2[k,] <- lscp.out$lower < Y.test & lscp.out$upper > Y.test
    condmat3[k,] <- oobw.output$estimate[,4] < Y.test &  oobw.output$estimate[,5] > Y.test
    lendmat[k,]  <- lapply(c(1:ntest),function(w){oobgk.output[w,5] - oobgk.output[w,4]}) %>% unlist
    lendmat2[k,] <- lapply(c(1:ntest),function(w){lscp.out$upper[w]-lscp.out$lower[w]}) %>% unlist
    lendmat3[k,] <- lapply(c(1:ntest),function(w){oobw.output$estimate[w,5] - oobw.output$estimate[w,4]}) %>% unlist
    
    oobg.ais.mar[k]   = mean( (oobg.output[,5] - oobg.output[,4]) + (2/alpha) * (sapply(oobg.output[,4] - Y.test, function(w){max(w,0)}) + sapply(Y.test - oobg.output[,5], function(w){max(w,0)})), na.rm=T )
    oobgk.ais.mar[k]  = mean( (oobgk.output[,5] - oobgk.output[,4]) + (2/alpha) * (sapply(oobgk.output[,4] - Y.test, function(w){max(w,0)}) + sapply(Y.test - oobgk.output[,5], function(w){max(w,0)})), na.rm=T )
    oobw.ais.mar[k]   = mean( (oobw.output$estimate[,5] - oobw.output$estimate[,4]) + (2/alpha) * (sapply(oobw.output$estimate[,4] - Y.test, function(w){max(w,0)}) + sapply(Y.test - oobw.output$estimate[,5], function(w){max(w,0)})), na.rm=T )
    oob.ais.mar[k]    = mean( (alter.output$oob_interval$up-alter.output$oob_interval$lo) + (2/alpha) * (sapply(alter.output$oob_interval$lo - Y.test, function(w){max(w,0)}) + sapply(Y.test - alter.output$oob_interval$up, function(w){max(w,0)})), na.rm=T )
    sc.ais.mar[k]     = mean( (alter.output$sc_interval$up-alter.output$sc_interval$lo)+ (2/alpha) * (sapply(alter.output$sc_interval$lo - Y.test, function(w){max(w,0)}) + sapply(Y.test - alter.output$sc_interval$up, function(w){max(w,0)})), na.rm=T )
    qrf.ais.mar[k]    = mean( (alter.output$quantreg_interval$up-alter.output$quantreg_interval$lo)+ (2/alpha) * (sapply(alter.output$quantreg_interval$lo - Y.test, function(w){max(w,0)}) + sapply(Y.test - alter.output$quantreg_interval$up, function(w){max(w,0)})), na.rm=T )
    lscp.ais.mar[k]   = mean( (lscp.out[,2] - lscp.out[,1]) + (2/alpha) * (sapply(lscp.out[,1] - Y.test, function(w){max(w,0)}) + sapply(Y.test - lscp.out[,2], function(w){max(w,0)})), na.rm=T )
    
    
    par(mfrow=c(1,3))
    yl = c(min(Y.test)-10,max(Y.test)+10); ind = c(2:ntest)
    #############################################################################
    #############################################################################
    plot(X.test[ind,1],Y.test[ind],xlim=c(-3,3),ylim=yl,pch=18,xlab="X1",ylab="OOBW",cex=1.5)
    points(X.test[ind,1],oobgk.output$pred[ind],col=4,pch=18,cex=1.5)
    polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(oobgk.output[ind,4],rev(oobgk.output[ind,5])),
            col=scales::alpha(4,0.25),border=FALSE)
    lines(X.test[ind,1],oobgk.output[ind,5],col=4)
    lines(X.test[ind,1],oobgk.output[ind,4],col=4)
    lines(X.test[ind,1],apply(X.test[ind,], 1, mfunc))
    
    points(X.test[ind,1],oobw.output$estimate[ind,1],col=3,pch=18,cex=1.5)
    polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(oobw.output$estimate[ind,4],rev(oobw.output$estimate[ind,5])),
            col=scales::alpha(3,0.25),border=FALSE)
    lines(X.test[ind,1],oobw.output$estimate[ind,5],col=3)
    lines(X.test[ind,1],oobw.output$estimate[ind,4],col=3)
    
    #############################################################################
    
    #############################################################################
    plot(X.test[ind,1],Y.test[ind],xlim=c(-3,3),ylim=yl,pch=18,xlab="X1",ylab="LSCP",cex=1.5)
    points(X.test[ind,1],oobgk.output$pred[ind],col=4,pch=18,cex=1.5)
    polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(oobgk.output[ind,4],rev(oobgk.output[ind,5])),
            col=scales::alpha(4,0.25),border=FALSE)
    lines(X.test[ind,1],oobgk.output[ind,5],col=4)
    lines(X.test[ind,1],oobgk.output[ind,4],col=4)
    lines(X.test[ind,1],apply(X.test[ind,], 1, mfunc))
    
    points(X.test[ind,1],rowMeans(lscp.out)[ind],col=2,pch=18,cex=1.5)
    polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(lscp.out[ind,1],rev(lscp.out[ind,2])),
            col=scales::alpha(2,0.25),border=FALSE)
    lines(X.test[ind,1],lscp.out[ind,1],col=2)
    lines(X.test[ind,1],lscp.out[ind,2],col=2)
    #############################################################################
    
    boxplot(lendmat[,ind],ylim=c(0,2.5*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25))
    par(new=T);boxplot(lendmat2[,ind],ylim=c(0,2.5*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="")
    par(new=T);boxplot(lendmat3[,ind],ylim=c(0,2.5*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="")
    par(new=TRUE)
    plot(colMeans(condmat,na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0.25,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
    lines(colMeans(condmat2,na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0.25,1),pch=1)
    lines(colMeans(condmat3,na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0.25,1),pch=1)
    axis(4,at=c(0.5,0.6,0.7,0.8,0.9,1.0))
    
    cat(paste(Sys.time(),"alpha",alpha,"simulation", k, "completed"), mtry1, "(",adjk, ")", mtry2, "|",
        paste0(100*round(mean(oobgk.cpr.mar,na.rm=T),3),"(",round(mean(oobgk.len.mar,na.rm=T),3),")","[",round(mean(oobgk.ais.mar,na.rm=T),2),"]"),
        paste0(100*round(mean(oobw.cpr.mar,na.rm=T),3),"(",round(mean(oobw.len.mar,na.rm=T),3),")","[",round(mean(oobw.ais.mar,na.rm=T),2),"]"),
        paste0(100*round(mean(lscp.cpr.mar,na.rm=T),3),"(",round(mean(lscp.len.mar,na.rm=T),3),")","[",round(mean(lscp.ais.mar,na.rm=T),2),"]"),"\n")
  } # for k ends

  # output.cpr.mar.list[[alphaind]]  = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)
  # output.len.mar.list[[alphaind]]  = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)
  output.ais.mar.list[[alphaind]]  = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)
  
# } # for alphaind ends

# output.cpr.mat = data.frame(oobgk.cpr.mar,oobg.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,oob.cpr.mar,sc.cpr.mar,qrf.cpr.mar)
# output.len.mat = data.frame(oobgk.len.mar,oobg.len.mar,oobw.len.mar,lscp.len.mar,oob.len.mar,sc.len.mar,qrf.len.mar)
# output.ais.mat = data.frame(oobgk.ais.mar,oobg.ais.mar,oobw.ais.mar,lscp.ais.mar,oob.ais.mar,sc.ais.mar,qrf.ais.mar)
# 
# outbut.table1.mat = cbind(rev(paste0(round(colMeans(output.cpr.mat,na.rm=T),2),"(",round(apply(output.cpr.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
#                           rev(paste0(round(colMeans(output.len.mat,na.rm=T),2),"(",round(apply(output.len.mat,2,function(w){sd(w,na.rm=T)}),2),")")),
#                           rev(paste0(round(colMeans(output.ais.mat,na.rm=T),2),"(",round(apply(output.ais.mat,2,function(w){sd(w,na.rm=T)}),2),")")))
# 
# xtable(outbut.table1.mat)



# png(filename=paste0(Sys.Date(),"LINEAR_Visualization1.png"),width=720,height=720)
# par(mfrow=c(1,1),mar=c(3,3,1,3))
# 
# yl = c(min(Y.test)-3,max(Y.test)+3); ind = c(2:ntest)
# 
# plot(X.test[ind,1],Y.test[ind],xlim=c(-3,3),ylim=yl,pch=18,xlab="X1",ylab="OOBW",cex=1.5)
# points(X.test[ind,1],oobgk.output$pred[ind],col=4,pch=18,cex=1.5)
# polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(oobgk.output$lower_0.1[ind],rev(oobgk.output$upper_0.1[ind])),
#         col=scales::alpha(4,0.25),border=FALSE)
# lines(X.test[ind,1],oobgk.output$upper_0.1[ind],col=4)
# lines(X.test[ind,1],oobgk.output$lower_0.1[ind],col=4)
# lines(X.test[ind,1],apply(X.test[ind,], 1, mfunc))
# 
# points(X.test[ind,1],oobw.output$estimate[ind,1],col=3,pch=18,cex=1.5)
# polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(oobw.output$estimate[ind,4],rev(oobw.output$estimate[ind,5])),
#         col=scales::alpha(3,0.25),border=FALSE)
# lines(X.test[ind,1],oobw.output$estimate[ind,5],col=3)
# lines(X.test[ind,1],oobw.output$estimate[ind,4],col=3)
# dev.off()
# 
# png(filename=paste0(Sys.Date(),"LINEAR_Visualization2.png"),width=720,height=720)
# par(mfrow=c(1,1),mar=c(3,3,1,3))
# 
# yl = c(min(Y.test)-3,max(Y.test)+3); ind = c(2:ntest)
# 
# plot(X.test[ind,1],Y.test[ind],xlim=c(-3,3),ylim=yl,pch=18,xlab="X1",ylab="LSCP",cex=1.5)
# points(X.test[ind,1],oobgk.output$pred[ind],col=4,pch=18,cex=1.5)
# polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(oobgk.output$lower_0.1[ind],rev(oobgk.output$upper_0.1[ind])),
#         col=scales::alpha(4,0.25),border=FALSE)
# lines(X.test[ind,1],oobgk.output$upper_0.1[ind],col=4)
# lines(X.test[ind,1],oobgk.output$lower_0.1[ind],col=4)
# lines(X.test[ind,1],apply(X.test[ind,], 1, mfunc))
# 
# points(X.test[ind,1],rowMeans(lscp.out)[ind],col=2,pch=18,cex=1.5)
# polygon(c(X.test[ind,1],rev(X.test[ind,1])),c(lscp.out[ind,1],rev(lscp.out[ind,2])),
#         col=scales::alpha(2,0.25),border=FALSE)
# lines(X.test[ind,1],lscp.out[ind,1],col=2)
# lines(X.test[ind,1],lscp.out[ind,2],col=2)
# dev.off()
# 
# png(filename=paste0(Sys.Date(),"LINEAR_Cond_Visual.png"),width=720,height=720)
# par(mfrow=c(1,1),mar=c(3,3,1,3))
# 
# yl = c(min(Y.test)-10,max(Y.test)+10); ind = c(2:ntest)
# 
# boxplot(lendmat[,ind],ylim=c(0,2.0*max(Y.test)),xlab="Test Index",ylab="Width",xaxt="n",col=scales::alpha(4,0.25))
# par(new=T);boxplot(lendmat2[,ind],ylim=c(0,2.0*max(Y.test)),xaxt="n",col=scales::alpha(2,0.25),xlab="",ylab="")
# par(new=T);boxplot(lendmat3[,ind],ylim=c(0,2.0*max(Y.test)),xaxt="n",col=scales::alpha(3,0.25),xlab="",ylab="")
# par(new=TRUE)
# plot(colMeans(condmat,na.rm=T),axes=FALSE, type="o", col=4,ylim=c(0.25,1),pch=1);abline(h=(1-alpha),col=1,lty=2,xlab="",ylab="")
# lines(colMeans(condmat2,na.rm=T),axes=FALSE, type="o", col=2,ylim=c(0.25,1),pch=1)
# lines(colMeans(condmat3,na.rm=T),axes=FALSE, type="o", col=3,ylim=c(0.25,1),pch=1)
# axis(4,at=c(0.5,0.6,0.7,0.8,0.9,1.0))
# dev.off()