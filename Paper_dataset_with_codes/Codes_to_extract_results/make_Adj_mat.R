
make.adj.matrix <-
  function(theta, trunc)
  {
    adj1 = theta; adj1[abs(adj1)>trunc] = 1;  adj1[abs(adj1)<trunc] = 0;
    return(adj1)
  }


TPR_FPR<-function(A1, B1){
  
  K = length(A1)
  true_det = false_det = net_true = net_false = edge_selected = NULL
  for(k in 1:K){
  A = A1[[k]]
  B = B1[[k]]
  true_det = c(true_det,length(which(B[col(B)>row(B)]==1&A[col(A)>row(A)]==1)))
  false_det = c(false_det,length(which(B[col(B)>row(B)]==1&A[col(A)>row(A)]==0)))
  net_true =  c(net_true, length(which(A[col(A)>row(A)]==1)))
  net_false = c(net_false, length(which(A[col(A)>row(A)]==0)))
  edge_selected  = c(edge_selected, length(which(B[col(B)>row(B)]==1)))
  }
  
  TPR = sum(true_det)/sum(net_true)
  FPR = sum(false_det)/sum(net_false)
  return(list(TPR, FPR, P = sum(edge_selected)))
}





TPR_FPR_df<-function(A1, B1){
  
  K = length(A1)
  true_det = false_det = net_true = net_false = edge_selected = NULL
  net_true = net_false  = list()
  for(k in 1:K){
    A = A1[[k]]
    net_true[[k]]  = which(A[col(A)>row(A)]==1)
    net_false[[k]] = which(A[col(A)>row(A)]==0)
  }
  true_diff = setdiff(net_true[[1]],net_true[[2]])
  false_diff = setdiff(net_false[[2]],net_false[[1]])

  det_true = det_false  = list()
  for(k in 1:K){
    B = B1[[k]]
    det_true[[k]]  = which(B[col(B)>row(B)]==1)
    det_false[[k]] = which(B[col(B)>row(B)]==0)
  }
  true_det= setdiff(det_true[[1]],det_true[[2]])
  false_det= setdiff(det_false[[1]],det_false[[2]])
  
 TP =  intersect(union(true_diff,false_diff),  union(true_det,false_det))
 FP =  setdiff(union(true_det,false_det),TP)
 
 TP =  intersect(true_diff,  true_det)
 FP =  setdiff(union(true_det,false_det),TP)
 
  #trueTP = length(TP) / length(true_diff)
  #falseFP = length(FP) / length(false_diff)
  
  return(list(length(TP),length(FP)))
}

PPV<-function(A1, B1){
  
  K = length(A1)
  true_det = false_det = net_true = net_false = edge_selected = NULL
  for(k in 1:K){
    A = A1[[k]]
    B = B1[[k]]
    true_det = c(true_det,length(which(B[col(B)>row(B)]==1&A[col(A)>row(A)]==1)))
    false_det = c(false_det,length(which(B[col(B)>row(B)]==1&A[col(A)>row(A)]==0)))
    net_true =  c(net_true, length(which(A[col(A)>row(A)]==1)))
    net_false = c(net_false, length(which(A[col(A)>row(A)]==0)))
    edge_selected  = c(edge_selected, length(which(B[col(B)>row(B)]==1)))
  }
  
  Recall = sum(true_det)/sum(net_true)
  Precision = ifelse(sum(true_det)+sum(false_det) == 0, 1, sum(true_det)/(sum(true_det)+sum(false_det)))
  return(list(Precision, Recall, P = sum(edge_selected)))
}



