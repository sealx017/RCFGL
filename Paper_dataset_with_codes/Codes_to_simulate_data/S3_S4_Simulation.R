###------------------Code to simulate datasets under simulation scenarios S3 and S4---------------
##-------------------Simulation scenario S3-------------------

#------------------Defining a couple of data generating functions--------------------------------
unif_dis<-function(){
  y = sample(1:2,1)
  if(y==1){x = runif(1,0.6,0.9)}
  else x = runif(1,-0.9,-0.6)
  return(x)
}

scaler = function(x){
  m = rowSums(x)
  nonz = which(m>0)
  x[nonz,] = x[nonz,]/m[nonz]/1.5
  return(x)
}

#-----------------Generating five blocks of sub-networks for four conditions using Barabasi-Albert model------------
library(igraph)
library(raster)
folder_path  = "/Users/seals/Desktop/CSPH/CFGL/31Aug/Paper_dataset_with_codes/" #the location of the folder on the hard drive
set.seed(5)
size = 100 #size of every block
class1 = class2 = class3 = class4 =  NULL
for(i in 1:5){
  g1 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  g2 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  g3 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  g4 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  
  class1[[i]] = g1
  class2[[i]] = g2
  class3[[i]] = g3
  class4[[i]] = g4
}
p = 5*size # total number of genes 

#------------------Using the generated sub-networks to simulate true covariance/precision matrices---------
#------------------The algorithm below is fully explained in the manuscript---------------------------------

real_net1 = real_net2 = real_net3 = real_net4 = first_cov = second_cov = third_cov = fourth_cov = matrix(0,p,p)

for(i in c(1:5))
{
  adjm1 = as.matrix(as_adjacency_matrix(class1[[i]]))
  real_net1[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm1
  
  adjm2 = as.matrix(as_adjacency_matrix(class2[[i]]))
  real_net2[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm2
  
  adjm3 = as.matrix(as_adjacency_matrix(class3[[i]]))
  real_net3[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm3
  
  adjm4 = as.matrix(as_adjacency_matrix(class4[[i]]))
  real_net4[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm4
  
  adjm1_off = adjm1; adjm1_off[row(adjm1_off)>col(adjm1_off)]=0
  adjm2_off = adjm2; adjm2_off[row(adjm2_off)>col(adjm2_off)]=0
  adjm3_off = adjm3; adjm3_off[row(adjm3_off)>col(adjm3_off)]=0
  adjm4_off = adjm4; adjm4_off[row(adjm4_off)>col(adjm4_off)]=0
  
  set1 = which(adjm1_off==1, arr.ind=T)
  set2 = which(adjm2_off==1, arr.ind=T)
  set3 = which(adjm3_off==1, arr.ind=T)
  set4 = which(adjm4_off==1, arr.ind=T)
  
  for(k in 1:dim(set1)[1]){
    adjm1_off[set1[k,1],set1[k,2]] = unif_dis()
  }
  
  for(k in 1:dim(set2)[1]){
    adjm2_off[set2[k,1],set2[k,2]] = unif_dis()
  }
  
  for(k in 1:dim(set3)[1]){
    adjm3_off[set3[k,1],set3[k,2]] = unif_dis()
  }
  
  for(k in 1:dim(set4)[1]){
    adjm4_off[set4[k,1],set4[k,2]] = unif_dis()
  }
  adjm1_off = adjm1_off + t(adjm1_off);  adjm2_off = adjm2_off + t(adjm2_off);  adjm3_off = adjm3_off + t(adjm3_off); adjm4_off = adjm4_off + t(adjm4_off)
  adjm1 = (adjm1_off)
  adjm2 = (adjm2_off)
  adjm3 = (adjm3_off)
  adjm4 = (adjm4_off)
  
  diag(adjm1)= diag(adjm2) = diag(adjm3) = diag(adjm4) = 1
  eig1 = eigen(adjm1);   eig2 = eigen(adjm2);  eig3 = eigen(adjm3); eig4 = eigen(adjm4)
  adjm1 = adjm1 + diag((abs(min(eig1$values))+0.1),size,size)
  adjm2 = adjm2 + diag((abs(min(eig2$values))+0.1),size,size)
  adjm3 = adjm3 + diag((abs(min(eig3$values))+0.1),size,size)
  adjm4 = adjm4 + diag((abs(min(eig4$values))+0.1),size,size)
  
  #print(matrixcalc::is.positive.semi.definite(adjm1))
  #print(matrixcalc::is.positive.semi.definite(adjm2))
  #print(matrixcalc::is.positive.semi.definite(adjm3))
  #print(matrixcalc::is.positive.semi.definite(adjm4))
  
  adjm1_inv = cov2cor(solve(adjm1))
  adjm2_inv = cov2cor(solve(adjm2))
  adjm3_inv = cov2cor(solve(adjm3))
  adjm4_inv = cov2cor(solve(adjm4))
  
  diag(adjm1_inv) = diag(adjm2_inv) = diag(adjm3_inv) = diag(adjm4_inv) = 1
  
  first_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm1_inv
  #second_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm2_inv
  third_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm3_inv
  #fourth_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm4_inv

}


#------------------Keeping the covariance matrices of the first two conditions exactly equal---------
#------------------and the covariance matrices of the last two conditions exactly equal--------------

second_cov = first_cov
fourth_cov = third_cov


#-----------------Computing the corresponding precision matrices and saving them-------------------

n = 100 #number of subjects
true_theta1 = solve(first_cov)
true_theta2 = solve(second_cov)
true_theta3 = solve(third_cov)
true_theta4 = solve(fourth_cov)

write.csv(true_theta1,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",1,".csv"))
write.csv(true_theta2,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",2,".csv"))
write.csv(true_theta3,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",3,".csv"))
write.csv(true_theta4,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",4,".csv"))

#----------------Generating the expression data using the true covariance matrices under each condition---

for(sim in 1:10){
  
  y1 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = first_cov)
  y2 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = second_cov)
  y3 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = third_cov)
  y4 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = fourth_cov)
  
  write.csv(y1,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",1,"_sim_",sim,".csv"),row.names = F)
  
  write.csv(y2,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",2,"_sim_",sim,".csv"),row.names = F)
  
  write.csv(y3,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",3,"_sim_",sim,".csv"),row.names = F)
  
  write.csv(y4,paste0(folder_path , "Simulated_datasets/4condition/two_twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",4,"_sim_",sim,".csv"),row.names = F)
}



##-------------------Simulation scenario S4-------------------

#-----------------Generating five blocks of sub-networks for four conditions using Barabasi-Albert model------------

set.seed(4)
size = 100 #size of every block
class1 = class2 = class3 = class4 =  NULL
for(i in 1:5){
  g1 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  g2 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  g3 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  g4 = barabasi.game(size, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE,directed=F)
  
  class1[[i]] = g1
  class2[[i]] = g2
  class3[[i]] = g3
  class4[[i]] = g4
}
p = 5*size # total number of genes 

#------------------Using the generated sub-networks to simulate true covariance/precision matrices---------

real_net1 = real_net2 = real_net3 = real_net4 = first_cov = second_cov = third_cov = fourth_cov = matrix(0,p,p)

for(i in c(1:5))
{
  adjm1 = as.matrix(as_adjacency_matrix(class1[[i]]))
  real_net1[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm1
  
  adjm2 = as.matrix(as_adjacency_matrix(class2[[i]]))
  real_net2[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm2
  
  adjm3 = as.matrix(as_adjacency_matrix(class3[[i]]))
  real_net3[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm3
  
  adjm4 = as.matrix(as_adjacency_matrix(class4[[i]]))
  real_net4[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm4
  
  adjm1_off = adjm1; adjm1_off[row(adjm1_off)>col(adjm1_off)]=0
  adjm2_off = adjm2; adjm2_off[row(adjm2_off)>col(adjm2_off)]=0
  adjm3_off = adjm3; adjm3_off[row(adjm3_off)>col(adjm3_off)]=0
  adjm4_off = adjm4; adjm4_off[row(adjm4_off)>col(adjm4_off)]=0
  
  set1 = which(adjm1_off==1, arr.ind=T)
  set2 = which(adjm2_off==1, arr.ind=T)
  set3 = which(adjm3_off==1, arr.ind=T)
  set4 = which(adjm4_off==1, arr.ind=T)
  
  for(k in 1:dim(set1)[1]){
    adjm1_off[set1[k,1],set1[k,2]] = unif_dis()
  }
  
  for(k in 1:dim(set2)[1]){
    adjm2_off[set2[k,1],set2[k,2]] = unif_dis()
  }
  
  for(k in 1:dim(set3)[1]){
    adjm3_off[set3[k,1],set3[k,2]] = unif_dis()
  }
  
  for(k in 1:dim(set4)[1]){
    adjm4_off[set4[k,1],set4[k,2]] = unif_dis()
  }
  adjm1_off = adjm1_off + t(adjm1_off);  adjm2_off = adjm2_off + t(adjm2_off);  adjm3_off = adjm3_off + t(adjm3_off); adjm4_off = adjm4_off + t(adjm4_off)
  adjm1 = (adjm1_off)
  adjm2 = (adjm2_off)
  adjm3 = (adjm3_off)
  adjm4 = (adjm4_off)
  
  diag(adjm1)= diag(adjm2) = diag(adjm3) = diag(adjm4) = 1
  eig1 = eigen(adjm1);   eig2 = eigen(adjm2);  eig3 = eigen(adjm3); eig4 = eigen(adjm4)
  adjm1 = adjm1 + diag((abs(min(eig1$values))+0.1),size,size)
  adjm2 = adjm2 + diag((abs(min(eig2$values))+0.1),size,size)
  adjm3 = adjm3 + diag((abs(min(eig3$values))+0.1),size,size)
  adjm4 = adjm4 + diag((abs(min(eig4$values))+0.1),size,size)
  
  #print(matrixcalc::is.positive.semi.definite(adjm1))
  #print(matrixcalc::is.positive.semi.definite(adjm2))
  #print(matrixcalc::is.positive.semi.definite(adjm3))
  #print(matrixcalc::is.positive.semi.definite(adjm4))
  
  adjm1_inv = cov2cor(solve(adjm1))
  adjm2_inv = cov2cor(solve(adjm2))
  adjm3_inv = cov2cor(solve(adjm3))
  adjm4_inv = cov2cor(solve(adjm4))
  
  diag(adjm1_inv) = diag(adjm2_inv) = diag(adjm3_inv) = diag(adjm4_inv) = 1
  
  first_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm1_inv
  #second_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm2_inv
  third_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm3_inv
  fourth_cov[(((i-1)*size+1):(i*size)),(((i-1)*size+1):(i*size))] = adjm4_inv
  
}


#------------------Keeping the covariance matrices of the first two conditions exactly equal---------
#------------------whereas the covariance matrices of the last two conditions are kept different-----

second_cov = first_cov


#-----------------Computing the corresponding precision matrices and saving them-------------------

n = 100 #number of subjects
true_theta1 = solve(first_cov)
true_theta2 = solve(second_cov)
true_theta3 = solve(third_cov)
true_theta4 = solve(fourth_cov)

write.csv(true_theta1,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",1,".csv"))
write.csv(true_theta2,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",2,".csv"))
write.csv(true_theta3,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",3,".csv"))
write.csv(true_theta4,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_condition_",4,".csv"))

#----------------Generating the expression data using the true covariance matrices under each condition---

for(sim in 1:10){
  
  y1 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = first_cov)
  y2 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = second_cov)
  y3 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = third_cov)
  y4 = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = fourth_cov)
  
  write.csv(y1,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",1,"_sim_",sim,".csv"),row.names = F)
  
  write.csv(y2,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",2,"_sim_",sim,".csv"),row.names = F)
  
  write.csv(y3,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",3,"_sim_",sim,".csv"),row.names = F)
  
  write.csv(y4,paste0(folder_path , "Simulated_datasets/4condition/twosame_morediff",
                      "p_",p,"_n_",n,"_condition_",4,"_sim_",sim,".csv"),row.names = F)
}
