
rm(list = ls())
library(randomForest);
library(RandomFields)
library(rfinterval)
library(forestError)
library(plot3D)
library(plotrix)
library(dplyr)
library(expm)
library(scales)
library(spatstat)
library(RandomForestsGLS)
# devtools::install_github("mhuiying/lscp")
library(scp)
# source('RFGLS_support_YBJ.R')
source('RFGLS_support2_YBJ.R')
################################

nsim = 1;
ntrain = 200;
ntest = 50;
nsample = ntrain + ntest;
sd2 = sqrt(2.4); # marginal variance for spatial error
sd3 = sqrt(3.0 - sd2^2); # marginal variance (nugget) for non-spatial error
nu0 = 0.5; phi = 1
ntree.gls <- 0.5*ntrain; ntree <- ntrain
nx = 20 ;           # number of dimensions of covariates
MTRY = max(1,nx/3);
ns = 12 # nodesize
eta = 25 # smooth parameter
#########

# Linear
f1 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  EY = x[1]+x[2]
  return(EY)
}
# Nonlinear (default)
f2 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  # EY = 20*exp(-abs(x[1])-abs(x[2]))
  EY = 10*sin(pi*x[1]/3)
  return(EY)
}
# Step
f3 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  # EY = ((x[1] >= 0)-(x[1] < 0))+x[2]
  EY = 1.5*((x[1] >= 0)-(x[1] < 0))+x[2]
  return(EY)
}
#Friedman
f4 = function(x){
  # 'x' is a vector, covariates corresponding to a case
  EY = (10*(sin(pi*x[1]*x[2]) + 20*(x[3]-0.5)^2 + 10*x[4] + 5*x[5]))
  return(EY)
}
#########

mfunc = f2

#########

simulate_data = function(seed,ntrain,ntest,nx,nu0,sd2,sd3,mfunc){
  
  ### Simulate dataset
  # set.seed(seed+1); X = matrix(runif((ntrain+ntest) * nx,min=-3,max=3), ncol = nx)
  # set.seed(seed+2); Loc = cbind(runif((ntrain+ntest),0,1),runif((ntrain+ntest),0,1))
  set.seed(seed+10); X.train = matrix(runif(ntrain * nx,min=-3,max=3), ncol = nx)
  set.seed(seed+20); X.test  = matrix(runif(ntest  * nx,min=-3,max=3), ncol = nx)
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


#########

seed0 = 123
alpha = 0.1
alphavec = rev(seq(0.1,0.9,by=0.1))
sd2vec   = sqrt(rev(seq(0.2,2.8,by=0.2)))
mvec = seq(1,nx,by=2)
ntree.gls.vec = seq(100,3*ntrain,by=50)

output.cpr.mar.list = list()
output.len.mar.list = list()

k <- 1
rf.gls.list = list()
rf.list = list()
for (adjust in seq(1,9,by=1)){
  
  seed <- seed0 + k
  simdata = simulate_data(seed,ntrain,ntest,nx,nu0,sd2,sd3,mfunc)
  train_data = simdata[[1]]; test_data = simdata[[2]]
  X.train    = simdata[[3]]; X.test    = simdata[[6]]
  Y.train    = simdata[[4]]; Y.test    = simdata[[7]]
  s.train    = simdata[[5]]; s.test    = simdata[[8]]
  
  ##############################
  rf.gls.list[[k]] <- RFGLS_estimate_spatial(s.train, y = Y.train, X = X.train, ntree = ntree.gls, param_estimate = TRUE, nthsize = ns, h=8, mtry = MTRY, phi=1)
  
  # fit random forest to the training data
  rf.list[[k]] <- randomForest::randomForest(X.train, Y.train, nodesize = 5, ntree = ntree, keep.inbag = TRUE, proximity = TRUE)
  
  cat(paste(Sys.time(),"alpha",alpha,"simulation", k, "completed","\n"))
  k <- k + 1
} # for k ends
  
  # output.cpr.mar.list[[alphaind]]  = data.frame(oobgk.cpr.mar,oobw.cpr.mar,lscp.cpr.mar,qrf.cpr.mar)
  # output.len.mar.list[[alphaind]]  = data.frame(oobgk.len.mar,oobw.len.mar,lscp.len.mar,qrf.len.mar)
  
# } # for alphaind ends 
par(mfrow=c(3,3))
for(adj in 1:9){
# for(k in 1:9){
  limind = c(floor(min(rf.list[[k]]$y - rf.list[[k]]$predicted)),ceiling(max(rf.list[[k]]$y - rf.list[[k]]$predicted)))
  oobgk.output <- quantGLSForestError(rf.gls, rf.gls.list[[k]]$y, rf.gls.list[[k]]$X, alpha = 0.1, kernel=TRUE, adjust = adj)
  # hist(rf.gls.list[[k]]$y - rf.gls.list[[k]]$predicted, xlim = limind, ylim = c(0,100), col = alpha(2,0.5), main="",xlab="");
  # par(new=T);hist(rf.list[[k]]$y - rf.list[[k]]$predicted, xlim = limind, ylim = c(0,100), col = alpha(4,0.5), main="",xlab="");
  
}

par(mfrow=c(2,2))
temp.table = output.table[1:50,]
plot(temp.table[1:10,3],col=temp.table[1:10,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[11:20,3],col=temp.table[11:20,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[21:30,3],col=temp.table[21:30,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[31:40,3],col=temp.table[31:40,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[41:50,3],col=temp.table[41:50,1],type='o',ylim=c(5,45),ylab="")

temp.table = output.table[51:100,]
plot(temp.table[1:10,3],col=temp.table[1:10,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[11:20,3],col=temp.table[11:20,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[21:30,3],col=temp.table[21:30,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[31:40,3],col=temp.table[31:40,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[41:50,3],col=temp.table[41:50,1],type='o',ylim=c(5,45),ylab="")

temp.table = output.table[101:150,]
plot(temp.table[1:10,3],col=temp.table[1:10,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[11:20,3],col=temp.table[11:20,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[21:30,3],col=temp.table[21:30,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[31:40,3],col=temp.table[31:40,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[41:50,3],col=temp.table[41:50,1],type='o',ylim=c(5,45),ylab="")

temp.table = output.table[151:200,]
plot(temp.table[1:10,3],col=temp.table[1:10,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[11:20,3],col=temp.table[11:20,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[21:30,3],col=temp.table[21:30,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[31:40,3],col=temp.table[31:40,1],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[41:50,3],col=temp.table[41:50,1],type='o',ylim=c(5,45),ylab="")

par(mfrow=c(1,2))
temp.table = output.table
plot(temp.table[seq(1,200,by=10),3],col=temp.table[seq(1,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(2,200,by=10),3],col=temp.table[seq(2,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(3,200,by=10),3],col=temp.table[seq(3,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(4,200,by=10),3],col=temp.table[seq(4,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(5,200,by=10),3],col=temp.table[seq(5,200,by=10),2],type='o',ylim=c(5,45),ylab="")

temp.table = output.table
plot(temp.table[seq(6,200,by=10),3],col=temp.table[seq(1,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(7,200,by=10),3],col=temp.table[seq(2,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(8,200,by=10),3],col=temp.table[seq(3,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(9,200,by=10),3],col=temp.table[seq(4,200,by=10),2],type='o',ylim=c(5,45),ylab="")
par(new=T);plot(temp.table[seq(10,200,by=10),3],col=temp.table[seq(5,200,by=10),2],type='o',ylim=c(5,45),ylab="")

