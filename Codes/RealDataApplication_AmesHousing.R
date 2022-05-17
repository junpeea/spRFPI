load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/Backups/AmesHousing/AmesHousing220310.RData")

# rm(list=ls())
library(dplyr)
library(RandomForestsGLS)
source('RFGLS_support2_YBJ.R')
source('RFGLS_support3_YBJ.R')
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
library(AmesHousing)

mydat = Comp.Ames_xy_Housing[!is.na(Comp.Ames_xy_Housing$SalePrice),]
mydat$origAddress.lon <- (mydat$origAddress.lon - mean(mydat$origAddress.lon))/sd(mydat$origAddress.lon)
mydat$origAddress.lat <- (mydat$origAddress.lat - mean(mydat$origAddress.lat))/sd(mydat$origAddress.lat)
mydat$SalePrice <- (log(mydat$SalePrice) - mean(log(mydat$SalePrice)))/sd(log(mydat$SalePrice))

alpha = 0.1
ntrain = ceiling(0.9*nrow(mydat))
ntree.gls <- floor(0.25*ntrain); ntree <- floor(0.25*ntrain)
ns = 100
adj = 3.0
eta = 25
nsim  = 30
Final.output = matrix(NA,nsim,6*3)

for(sim in 1:nsim){
  
  set.seed(2205041+sim);train.ind = sample(c(1:nrow(mydat)),ntrain);length(train.ind)
  test.ind  = setdiff(c(1:nrow(mydat)),train.ind);length(test.ind)
  
  train_data = mydat[train.ind,]; test_data = mydat[test.ind,]
  X.train    = train_data[,c(-1,-2,-3)] %>% as.matrix; X.test    = test_data[,c(-1,-2,-3)] %>% as.matrix
  Y.train    = train_data[,c(1)]; Y.test    = test_data[,c(1)]
  s.train    = train_data[,c(2,3)] %>% as.matrix; s.test    = test_data[,c(2,3)] %>% as.matrix
  mtry = floor((16/20) * ncol(X.train))
  
  rf.gls <- RFGLS_estimate_spatial(s.train, y = Y.train, X = X.train, ntree = ntree.gls, param_estimate = TRUE, nthsize = ns, h=16, mtry = mtry)
  oobgk.output <- quantGLSForestError(rf.gls, X.train, X.test, alpha = alpha, kernel=TRUE, adjust = adj)
  
  rf <- randomForest::randomForest(X.train, Y.train, keep.inbag = TRUE, proximity = TRUE, mtry = mtry, nodesize=ns)
  train_nodes <- findOOBErrors(rf, X.train)
  oobw.output <- quantForestError(rf, X.train, X.test,
                                  train_nodes = train_nodes, alpha = alpha)
  
  s0 = s.test; s = s.train; Y = matrix(Y.train - rf$predicted, ncol=1)
  rf.pred = predict(rf,X.test)
  # s0 = s.test; s = s.train; Y = matrix(Y.train, ncol=1)
  bins     = seq(0.01,0.2,0.01)
  thetaHat = get_theta(s,Y,dists=bins)
  thetaHat[4] <- min(thetaHat[4],10)
  sjit  = s+cbind(rnorm(nrow(s.train),0,1e-03),rnorm(nrow(s.train),0,1e-03))
  lscp.out = scp(s0,sjit,Y,alpha = 0.1,thetaHat = thetaHat, global = FALSE, eta=eta)
  lscp.out = data.frame(lower=lscp.out[,1]+rf.pred,upper=lscp.out[,2]+rf.pred)
  # lscp.out = data.frame(lower=lscp.out[,1],upper=lscp.out[,2])
  
  alter.output <- rfinterval(SalePrice~., train_data = train_data, test_data = test_data,
                             method = c("oob", "split-conformal", "quantreg"),
                             symmetry = TRUE ,alpha = alpha, params_ranger = list(mtry=mtry,min.node.size=ns))
  
  oobgk.cpr.mar   = mean(oobgk.output[,4] < Y.test & oobgk.output[,5] > Y.test)
  oobw.cpr.mar    = mean(oobw.output$estimate[,4] < Y.test &  oobw.output$estimate[,5] > Y.test)
  oob.cpr.mar    = mean(alter.output$oob_interval$lo < Y.test & alter.output$oob_interval$up > Y.test)
  sc.cpr.mar     = mean(alter.output$sc_interval$lo < Y.test & alter.output$sc_interval$up > Y.test)
  qrf.cpr.mar    = mean(alter.output$quantreg_interval$lo < Y.test & alter.output$quantreg_interval$up > Y.test)
  lscp.cpr.mar    = mean(lscp.out$lower < Y.test & lscp.out$upper > Y.test)
  
  oobgk.len.mar  = mean(oobgk.output[,5] - oobgk.output[,4])
  oobw.len.mar   = mean(oobw.output$estimate[,5] - oobw.output$estimate[,4])
  oob.len.mar    = mean(alter.output$oob_interval$up-alter.output$oob_interval$lo)
  sc.len.mar     = mean(alter.output$sc_interval$up-alter.output$sc_interval$lo)
  qrf.len.mar    = mean(alter.output$quantreg_interval$up-alter.output$quantreg_interval$lo)
  lscp.len.mar   = mean(lscp.out$upper-lscp.out$lower)
  
  condmat  <- oobgk.output[,4] < Y.test & oobgk.output[,5] > Y.test; mean(condmat)
  condmat2 <- lscp.out$lower < Y.test & lscp.out$upper > Y.test; mean(condmat2)
  condmat3 <- oobw.output$estimate[,4] < Y.test &  oobw.output$estimate[,5] > Y.test; mean(condmat3)
  lendmat  <- lapply(c(1:length(test.ind)),function(w){oobgk.output[w,5] - oobgk.output[w,4]}) %>% unlist
  lendmat2 <- lapply(c(1:length(test.ind)),function(w){lscp.out$upper[w]-lscp.out$lower[w]}) %>% unlist
  lendmat3 <- lapply(c(1:length(test.ind)),function(w){oobw.output$estimate[w,5] - oobw.output$estimate[w,4]}) %>% unlist
  
  oobgk.ais.mar  = mean( (oobgk.output[,5] - oobgk.output[,4]) + (2/alpha) * (sapply(oobgk.output[,4] - Y.test, function(w){max(w,0)}) + sapply(Y.test - oobgk.output[,5], function(w){max(w,0)})), na.rm=T )
  oobw.ais.mar   = mean( (oobw.output$estimate[,5] - oobw.output$estimate[,4]) + (2/alpha) * (sapply(oobw.output$estimate[,4] - Y.test, function(w){max(w,0)}) + sapply(Y.test - oobw.output$estimate[,5], function(w){max(w,0)})), na.rm=T )
  oob.ais.mar   = mean( (alter.output$oob_interval$up-alter.output$oob_interval$lo) + (2/alpha) * (sapply(alter.output$oob_interval$lo - Y.test, function(w){max(w,0)}) + sapply(Y.test - alter.output$oob_interval$up, function(w){max(w,0)})), na.rm=T )
  sc.ais.mar   = mean( (alter.output$sc_interval$up-alter.output$oob_interval$lo) + (2/alpha) * (sapply(alter.output$sc_interval$lo - Y.test, function(w){max(w,0)}) + sapply(Y.test - alter.output$sc_interval$up, function(w){max(w,0)})), na.rm=T )
  qrf.ais.mar   = mean( (alter.output$quantreg_interval$up-alter.output$oob_interval$lo) + (2/alpha) * (sapply(alter.output$quantreg_interval$lo - Y.test, function(w){max(w,0)}) + sapply(Y.test - alter.output$quantreg_interval$up, function(w){max(w,0)})), na.rm=T )
  lscp.ais.mar   = mean( (lscp.out[,2] - lscp.out[,1]) + (2/alpha) * (sapply(lscp.out[,1] - Y.test, function(w){max(w,0)}) + sapply(Y.test - lscp.out[,2], function(w){max(w,0)})), na.rm=T )
  
  
  Final.output[sim,] = c(qrf.cpr.mar,sc.cpr.mar,oob.cpr.mar,lscp.cpr.mar,oobw.cpr.mar,oobgk.cpr.mar,
                         qrf.len.mar,sc.len.mar,oob.len.mar,lscp.len.mar,oobw.len.mar,oobgk.len.mar,
                         qrf.ais.mar,sc.ais.mar,oob.ais.mar,lscp.ais.mar,oobw.ais.mar,oobgk.ais.mar)
  
  cat(paste0(Sys.time(),"_Experiment"),sim,"/",nsim,"completed","\n")
  
}
colMeans(Final.output) %>% round(2)
# write.csv(Final.output,file="database_Jun.csv")
# save.image(file=paste0(Sys.Date(),"_AmesHousing.Rdata"))





