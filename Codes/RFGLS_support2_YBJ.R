
quantGLSForestError = function(rf.gls, X.train, X.test, alpha, kernel = TRUE, adjust = 1){
  
  ntrain = nrow(X.train); ntest = nrow(X.test); ncx = ncol(X.train)
  
  col_names <- c("pred", "mspe", "bias")
  
  train.terminal.nodes  <- getGLS_TerminalNodes(rf.gls,X.train)
  colnames(train.terminal.nodes) <- c(1:ntree.gls);rownames(train.terminal.nodes) <- c(1:ntrain)
  
  inbag = matrix(0,nrow=ntrain,ncol=ntree.gls)
  for(s in 1:ntree.gls){
    tab.df <- table(rf.gls$P_matrix[,s])
    inbag[as.integer(names(tab.df)[which(names(tab.df) > 0)]),s] <- matrix(tab.df[which(names(tab.df) > 0)],ncol=1)
  }
  rf.gls$inbag <- inbag
  bag.count <- rf.gls$inbag
  
  # oob.errors <- rf.gls$y - rf.gls$predicted
  # tmpdist = dist(rf.gls$coords) %>% as.matrix
  # Gamma = sd2^2*geoR::matern(tmpdist,phi=1,kappa=nu0) + sd3^2*diag(1,nrow(tmpdist))
  
  oob_err_list <- lapply(c(1:ntree.gls),function(k){
    oob_pred_mat <-  rf.gls$predicted_matrix
    a = rf.gls$y[bag.count[,k]==0]
    # b = oob_pred_mat[bag.count[,k]==0,k] + Gamma[bag.count[,k]==0,bag.count[,k]!=0]%*%solve(Gamma[bag.count[,k]!=0,bag.count[,k]!=0])%*%(rf.gls$y[bag.count[,k]!=0]-oob_pred_mat[bag.count[,k]!=0,k])
    b = oob_pred_mat[bag.count[,k]==0,k]
    return(a-b)
  })
  
  oob_err_matrix = matrix(NA,ntrain,ntree.gls)
  # oob_err_matrix = apply(rf.gls$predicted_matrix,2,function(w){rf.gls$y-w})
  for(k in 1:ntree.gls) oob_err_matrix[bag.count[,k]==0,k] <- oob_err_list[[k]]
  oob.errors = rowMeans(oob_err_matrix,na.rm=T)
  
  train.terminal.nodes[bag.count != 0] <- NA
  train_nodes <- data.table::as.data.table(train.terminal.nodes)
  train_nodes[, `:=`(oob_error = oob.errors)]
  train_nodes[, `:=`(rowid_train = c(1:ntrain))]
  train_nodes <- data.table::melt(
    train_nodes,
    id.vars = c("oob_error","rowid_train"),
    # id.vars = c("oob_error"),
    measure.vars = 1:ncol(train.terminal.nodes),
    variable.name = "tree",
    value.name = "terminal_node",
    variable.factor = FALSE,
    na.rm = TRUE)
  
  # collapse the long data.table by unique tree/node
  train_nodes <- train_nodes[,
                             .(node_errs = list(oob_error), id_train = list(rowid_train)),
                             keyby = c("tree", "terminal_node")]
  # get test predictions
  test.preds <- RFGLS_predict_spatial(rf.gls, s.test, X.test, h=12, verbose = FALSE)$prediction
  names(test.preds) <- c(1:ntest)
  
  # get terminal nodes of test observations
  test.terminal.nodes <- getGLS_TerminalNodes(rf.gls,X.test)
  
  # reshape test.terminal.nodes to be a long data.table and
  # add unique IDs and predicted values
  test_nodes <- data.table::melt(
    data.table::as.data.table(test.terminal.nodes)[, `:=`(rowid_test = .I, pred = test.preds)],
    id.vars = c("rowid_test", "pred"),
    measure.vars = 1:ncol(test.terminal.nodes),
    variable.name = "tree",
    variable.factor = FALSE,
    value.name = "terminal_node")
  
  # set key columns for faster indexing
  data.table::setkey(test_nodes, tree, terminal_node)
  
  oob_error_stats <-
    train_nodes[test_nodes,
                .(tree, terminal_node, rowid_test, pred, node_errs,id_train)][,
                                                                              .(mspe = mean(unlist(node_errs) ^ 2),
                                                                                bias = -mean(unlist(node_errs)),
                                                                                all_errs = list(sort(unlist(node_errs))),
                                                                                all_id_train = list(unlist(id_train)[order(unlist(node_errs))])),
                                                                              keyby = c("rowid_test", "pred")]
  
  # for(k in 1:ntest){
  #   temp = c(oob_err_matrix[oob_error_stats$all_id_train[[k]],])
  #   oob_error_stats$all_errs[[k]] <- sort(temp[!is.na(temp)])
  # }
  # 
  
  # format pred2iction interval output
  percentiles <- sort(c(alpha / 2, 1 - (alpha / 2)))
  interval_col_names <- paste0(rep(c("lower_", "upper_"), each = length(alpha)),
                               c(alpha, rev(alpha)))
  
  col_names <- c(col_names, interval_col_names)
  
  if(kernel){
    oob_error_stats[, (interval_col_names) :=  lapply(percentiles, FUN = function(p){
      pred + lapply(c(1:ntest),function(k){
        bins = seq(min(oob_err_matrix,na.rm=T),max(oob_err_matrix,na.rm=T),length.out=1000)
        temp = density(oob_error_stats$all_errs[[k]], adjust=adjust); temp2 = CDF(temp)
        bins[max(which(temp2(bins) <= p))]
      }) %>% unlist 
    })]
  }else{
    # compute prediction intervals
    oob_error_stats[, (interval_col_names) := lapply(percentiles, FUN = function(p) {
      pred + purrr::map_dbl(all_errs, ~.x[ceiling(length(.x) * p)])})]
    
  }
  
  output_df <- data.table::setDF(oob_error_stats[, ..col_names])
  
  return(output_df)
}


# getTREE
getGLS_Tree = function(rf.gls,k=1){
  outlist = rf.gls$RFGLS_object
  output = cbind(outlist$ldaughter[,k],outlist$rdaughter[,k],outlist$mbest[,k],outlist$upper[,k],outlist$nodestatus[,k])
  colnames(output) <- c("left daughter","right daughter","split var","split point","status")
  rownames(output) <- c(1:nrow(output))
  return(output)
}
# test -> FinalNode
# temp = getGLS_Tree(rf.gls,1)

getGLS_terminalnode = function(rf.gls,X.test,ix,it){
  ref_tree <- getGLS_Tree(rf.gls,it)
  x <- X.test[ix,]
  t <- 1; k <- -3; exit <- 1000
  while(k < -2 & exit > 0){
    if(x[as.integer(ref_tree[t,3])] > ref_tree[t,4]){
      t <- as.integer(ref_tree[t,2])
    }else{
      t <- as.integer(ref_tree[t,1])
    }
    k <- ref_tree[t,5]; exit <- exit - 1
  }
  return(t)
}
getGLS_TerminalNodes = function(rf.gls,X.test){
  m = nrow(X.test); n = ntree.gls
  out = matrix(NA,m,n)
  for(i in 1:m){
    for(j in 1:n){
      out[i,j] = getGLS_terminalnode(rf.gls,X.test,i,j)
    }
  }
  colnames(out) <- c(1:n); rownames(out) <- c(1:m)
  return(out)
}
# getGLS_TerminalNodes(rf.gls,X.test)
