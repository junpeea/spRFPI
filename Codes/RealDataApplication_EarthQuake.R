# rm(list=ls())
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/AmesHousingProject(New)/RealDataAnalysis/220406_EarthQuake.RData")
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

alpha = 0.1
# ntrain = 1000
YYYY = substr(Date,7,10) %>% as.integer
train.ind = which(YYYY < 2016 & YYYY > 2004);
ntree.gls <- 250; ntree <- 1000
ns = 2500
mtry = floor((16/20) * ncol(X.train))
adj = 2.5
eta = 25

# set.seed(220406);train.ind = sample(which(mydat$YYYY2016==0),ntrain);length(train.ind)
test.ind  = which(mydat$YYYY2016==1);length(test.ind)

train_data = mydat[train.ind,]; test_data = mydat[test.ind,]
X.train    = train_data[,c(-1,-2,-3)] %>% as.matrix; X.test    = test_data[,c(-1,-2,-3)] %>% as.matrix
Y.train    = train_data[,c(1)]; Y.test    = test_data[,c(1)]
s.train    = train_data[,c(2,3)] %>% as.matrix; s.test    = test_data[,c(2,3)] %>% as.matrix
mtry = floor((16/20) * ncol(X.train))

start = Sys.time()
seedkey=2205101
set.seed(seedkey)
rf.gls <- RFGLS_estimate_spatial(s.train, y = Y.train, X = X.train, ntree = ntree.gls, param_estimate = TRUE, nthsize = ns, mtry = mtry)
oobgk.output <- quantGLSForestError(rf.gls, X.train, X.test, alpha = alpha, kernel=TRUE, adjust = adj)
save(oobgk.output,file=paste0(Sys.Date(),"_oobgk.output",seedkey,".Rdata"))
time = Sys.time() - start
head(rf.gls$P_matrix)
time

load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/spRFPI_FINAL/2022-05-10_oobgk.output2205101.Rdata")
oobgk.output1 <- oobgk.output
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/spRFPI_FINAL/2022-05-10_oobgk.output2205102.Rdata")
oobgk.output2 <- oobgk.output
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/spRFPI_FINAL/2022-05-10_oobgk.output2205103.Rdata")
oobgk.output3 <- oobgk.output
load("C:/Users/yjun/Desktop/WORK2022/SpRFPI/spRFPI_FINAL/2022-05-10_oobgk.output2205104.Rdata")
oobgk.output4 <- oobgk.output
oobgk.output = (oobgk.output1+oobgk.output2+oobgk.output3+oobgk.output4)/4

rf <- randomForest::randomForest(X.train, Y.train, keep.inbag = TRUE, proximity = TRUE, mtry = mtry)
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
lscp.out

oobgk.cpr.mar   = mean(oobgk.output[,4] < Y.test & oobgk.output[,5] > Y.test)
oobw.cpr.mar    = mean(oobw.output$estimate[,4] < Y.test &  oobw.output$estimate[,5] > Y.test)
lscp.cpr.mar    = mean(lscp.out$lower < Y.test & lscp.out$upper > Y.test)

oobgk.len.mar  = mean(oobgk.output[,5] - oobgk.output[,4])
oobw.len.mar   = mean(oobw.output$estimate[,5] - oobw.output$estimate[,4])
lscp.len.mar   = mean(lscp.out$upper-lscp.out$lower)

condmat  <- oobgk.output[,4] < Y.test & oobgk.output[,5] > Y.test; mean(condmat)
condmat2 <- lscp.out$lower < Y.test & lscp.out$upper > Y.test; mean(condmat2)
condmat3 <- oobw.output$estimate[,4] < Y.test &  oobw.output$estimate[,5] > Y.test; mean(condmat3)
lendmat  <- lapply(c(1:length(test.ind)),function(w){oobgk.output[w,5] - oobgk.output[w,4]}) %>% unlist
lendmat2 <- lapply(c(1:length(test.ind)),function(w){lscp.out$upper[w]-lscp.out$lower[w]}) %>% unlist
lendmat3 <- lapply(c(1:length(test.ind)),function(w){oobw.output$estimate[w,5] - oobw.output$estimate[w,4]}) %>% unlist

oobgk.ais.mar  = mean( (oobgk.output[,5] - oobgk.output[,4]) + (2/alpha) * (sapply(oobgk.output[,4] - Y.test, function(w){max(w,0)}) + sapply(Y.test - oobgk.output[,5], function(w){max(w,0)})), na.rm=T )
oobw.ais.mar   = mean( (oobw.output$estimate[,5] - oobw.output$estimate[,4]) + (2/alpha) * (sapply(oobw.output$estimate[,4] - Y.test, function(w){max(w,0)}) + sapply(Y.test - oobw.output$estimate[,5], function(w){max(w,0)})), na.rm=T )
lscp.ais.mar   = mean( (lscp.out[,2] - lscp.out[,1]) + (2/alpha) * (sapply(lscp.out[,1] - Y.test, function(w){max(w,0)}) + sapply(Y.test - lscp.out[,2], function(w){max(w,0)})), na.rm=T )


Final.output = data.frame(s.test,oobgk.output[,c(1,4,5)],oobw.output$estimates[,c(1,4,5)],rf.pred,lscp.out)
colnames(Final.output) <- c("Latitude","Longitude","OOBGK_pred","OOBGK_lower","OOBGK_upper","OOBW_pred","OOBW_lower","OOBW_upper","LSCP_pred","LSCP_lower","LSCP_upper")
head(Final.output)


write.csv(Final.output,file=paste0(Sys.Date(),"_database_Jun.csv"))

head(Final.output)

mean(Final.output$OOBGK_lower > 5.5)
mean(Final.output$OOBW_lower > 5.5)
mean(Final.output$LSCP_lower > 5.5)






# library(ape)
# earthquake.dists <- as.matrix(dist(cbind(s.train)))
# earthquake.dists.inv <- 1/earthquake.dists
# diag(earthquake.dists.inv) <- 0
# Moran.I(Y.train, earthquake.dists.inv)
# # $observed
# # [1] 0.005592958
# # 
# # $expected
# # [1] -0.001001001
# # 
# # $sd
# # [1] 0.007782364
# # 
# # $p.value
# # [1] 0.3968307

