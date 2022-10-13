
#-----Function for making adjacency matrix------
make.adj.matrix <- function(theta, trunc)
  {
    adj1 = theta; adj1[abs(adj1)>trunc] = 1;  adj1[abs(adj1)<trunc] = 0;
    return(adj1)
  }


#-----Function for computing precision/recall from true and estimated adjacency matrices------
PPV <- function(A1, B1)
  {
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
  Precision = ifelse(sum(true_det)+sum(false_det) == 0, 1, 
                     sum(true_det)/(sum(true_det)+sum(false_det)))
  return(list(Precision, Recall, P = sum(edge_selected)))
}



