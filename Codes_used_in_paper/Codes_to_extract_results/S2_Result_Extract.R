##------This R code extracts the results of simulation  scenario (S2)--------
#----------------------------------------------------------------------------

library(Matrix)
library(reticulate) #package to read Python result files in R
np = import("numpy")
sp = import("scipy")
folder_path = "/Users/seals/Desktop/CSPH/CFGL/31Aug/Codes_used_in_paper/" #the location of the folder on the hard drive
source(paste0(folder_path, "Codes_to_extract_results/make_Adj_mat.R")) #loading some functions for downstream analysis 



p = 500; n = 100; sim = 10 #number of genes, samples and simulations
lambdas = matrix(c(0.01,0.05,0.01,0.1,0.01,0.2,0.03,0.05,
                   0.03,0.1,0.03,0.2,0.04,0.05,0.04,0.1,
                   0.04,0.2,0.05,0.05,0.05,0.1,0.05,0.2,0.1,0.05,
                   0.1,0.1,0.1,0.2,0.15,0.05,
                   0.15,0.1,0.15,0.2,0.2,0.05,0.2,0.1,
                   0.2,0.2,0.25,0.05,0.25,0.1,
                   0.25,0.2), nrow = 24, ncol = 2, byrow = T) #different combinations of lambda1 and lambda2
#for which the methods were run


true_theta1 = read.csv(paste0(folder_path, 
"Simulated_datasets/3condition/twosame_morediff_theta_",
                              "p_",p,"_n_",n,"_condition_",1,".csv"),header = T)[,-1]

true_theta2 = read.csv(paste0(folder_path, 
"Simulated_datasets/3condition/twosame_morediff_theta_",
                              "p_",p,"_n_",n,"_condition_",2,".csv"),header = T)[,-1]

true_theta3 = read.csv(paste0(folder_path, 
"Simulated_datasets/3condition/twosame_morediff_theta_",
                              "p_",p,"_n_",n,"_condition_",3,".csv"),header = T)[,-1] #loading true precision matrices 

trunc = 1e-2 #threshold to declare an off-diagonal element of a precision matrix as an edge
true_adj1 = make.adj.matrix(true_theta1, trunc);
true_adj2 = make.adj.matrix(true_theta2, trunc);
true_adj3 = make.adj.matrix(true_theta3, trunc); #constructing true adjacency matrices from true precision matrices

lambda1_unique = unique(lambdas[,1]) #unique values of lambda1 considered
lambda2_unique = unique(lambdas[,2]) #unique values of lambda1 considered

##--------Extracting the results of RFGL-----------
#--------------------------------------------------

SSE_lambda2 =  PPV_lambda2 = NULL
for(i in lambda2_unique){
  lambda2 = i;
  all_SSE = all_PPV = NULL
  
  for(j in lambda1_unique){
    lambda1 = j;
    SSE_mat  = matrix(0,sim,3)
    PPV_mat = matrix(0,sim,3)
    for(simn in c(1:sim)){
      #----loading the three estimated precision matrices under three conditions---
      c1 = sp$sparse$load_npz(paste0(folder_path, 
      'Simulated_datasets/Results/Python_res/3condition/RFGL_twosame_new_morediff_theta',
              'p_',as.character(p),'_n_',as.character(n),'_condition_',
              as.character(0),'_sim_',
              as.character(simn),'_condition_',as.character(0),
              '_lambda1_', as.character(format(round(lambda1, 2), nsmall = 1)),
              '_lambda2_', as.character(lambda2),'.npz'))
      
      c2 = sp$sparse$load_npz(paste0(folder_path, 
      'Simulated_datasets/Results/Python_res/3condition/RFGL_twosame_new_morediff_theta',
              'p_',as.character(p),'_n_',as.character(n),'_condition_',
               as.character(1),'_sim_',
               as.character(simn),'_condition_',as.character(1),
               '_lambda1_', as.character(format(round(lambda1, 2), nsmall = 1)),
               '_lambda2_', as.character(lambda2),'.npz'))
      
      c3 = sp$sparse$load_npz(paste0(folder_path, 
      'Simulated_datasets/Results/Python_res/3condition/RFGL_twosame_new_morediff_theta',
              'p_',as.character(p),'_n_',as.character(n),'_condition_',
              as.character(2),'_sim_',
              as.character(simn),'_condition_',as.character(2),
              '_lambda1_', as.character(format(round(lambda1, 2), nsmall = 1)),
              '_lambda2_', as.character(lambda2),'.npz'))
      
      #--------Computing SSE by comparing true and estimated precision matrices------------
      diff_1 = as.matrix(c1)-true_theta1;   diff_2 = as.matrix(c2)-true_theta2; diff_3 = as.matrix(c3)-true_theta3
      SSE = sum(diff_1[col(diff_1)>=row(diff_1)]^2)+sum(diff_2[col(diff_2)>=row(diff_2)]^2)+sum(diff_3[col(diff_3)>=row(diff_3)]^2)
      
      #--------Computing estimated adjacency matrices------------------------------
      adj1 = make.adj.matrix(as.matrix(c1),trunc)
      adj2 = make.adj.matrix(as.matrix(c2),trunc)
      adj3 = make.adj.matrix(as.matrix(c3),trunc)
      
      f1 = TPR_FPR(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      f2 = PPV(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      SSE_mat[simn,] = c(SSE, f1$P,lambda1)
      PPV_mat[simn,] = c(f2[[1]],f2[[2]],lambda1)
    }
    all_SSE = rbind(all_SSE, colMeans(SSE_mat))
    all_PPV = rbind(all_PPV, colMeans(PPV_mat))
  }
  SSE_lambda2 = rbind(SSE_lambda2,data.frame(all_SSE, lambda2))
  PPV_lambda2 = rbind(PPV_lambda2,data.frame(all_PPV, lambda2))
}

#-----------------------SSE extraction-----------------
names(SSE_lambda2) = c("RFGL", "RFGL","lambda1","lambda2")
SSE_long_SSE = tidyr::gather(SSE_lambda2[,c(1,3,4)],Method, SSE,
                             1, factor_key=TRUE)
SSE_long_edges= tidyr::gather(SSE_lambda2[,c(2,3,4)],Method, Edges,
                              1, factor_key=TRUE)
SSE_long2 = merge(SSE_long_edges, SSE_long_SSE,by=c("lambda1","lambda2","Method"))

#-----------------------PPV extraction------------------

names(PPV_lambda2) = c("RFGL", "RFGL","lambda1","lambda2")
Precision = tidyr::gather(PPV_lambda2[,c(1,3,4)],Method, Precision,
                          1, factor_key=TRUE)
Recall = tidyr::gather(PPV_lambda2[,c(2,3,4)],Method, Recall,
                       1, factor_key=TRUE)
PPV_long2 = merge(Precision, Recall,by=c("lambda1","lambda2","Method"))



##--------Extracting the results of RCFGL-----------
#---------------------------------------------------

SSE_lambda2 =  PPV_lambda2 = NULL
for(i in lambda2_unique){
  lambda2 = i;
  all_SSE = all_PPV = NULL
  
  for(j in lambda1_unique){
    lambda1 = j;
    SSE_mat  = matrix(0,sim,3)
    PPV_mat = matrix(0,sim,3)
    for(simn in c(1:sim)){
      #----loading the three estimated precision matrices under three conditions---
      c1 = sp$sparse$load_npz(paste0(folder_path, 
      'Simulated_datasets/Results/Python_res/3condition/RCFGL_twosame_new_morediff_theta',
              'p_',as.character(p),'_n_',as.character(n),'_condition_',
               as.character(0),'_sim_',
               as.character(simn),'_condition_',as.character(0),
               '_lambda1_', as.character(format(round(lambda1, 2), nsmall = 1)),
               '_lambda2_', as.character(lambda2),'.npz'))
      
      c2 = sp$sparse$load_npz(paste0(folder_path, 
      'Simulated_datasets/Results/Python_res/3condition/RCFGL_twosame_new_morediff_theta',
              'p_',as.character(p),'_n_',as.character(n),'_condition_',
               as.character(1),'_sim_',
              as.character(simn),'_condition_',as.character(1),
              '_lambda1_', as.character(format(round(lambda1, 2), nsmall = 1)),
               '_lambda2_', as.character(lambda2),'.npz'))
      
      c3 = sp$sparse$load_npz(paste0(folder_path, 
      'Simulated_datasets/Results/Python_res/3condition/RCFGL_twosame_new_morediff_theta',
              'p_',as.character(p),'_n_',as.character(n),'_condition_',
              as.character(2),'_sim_',
              as.character(simn),'_condition_',as.character(2),
              '_lambda1_', as.character(format(round(lambda1, 2), nsmall = 1)),
              '_lambda2_', as.character(lambda2),'.npz'))
      
      #--------Computing SSE by comparing true and estimated precision matrices------------
      diff_1 = as.matrix(c1)-true_theta1;   diff_2 = as.matrix(c2)-true_theta2; diff_3 = as.matrix(c3)-true_theta3
      SSE = sum(diff_1[col(diff_1)>=row(diff_1)]^2)+sum(diff_2[col(diff_2)>=row(diff_2)]^2)+sum(diff_3[col(diff_3)>=row(diff_3)]^2)
      
      #--------Computing estimated adjacency matrices------------------------------
      adj1 = make.adj.matrix(as.matrix(c1),trunc)
      adj2 = make.adj.matrix(as.matrix(c2),trunc)
      adj3 = make.adj.matrix(as.matrix(c3),trunc)
      
      
      f1 = TPR_FPR(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      f2 = PPV(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      SSE_mat[simn,] = c(SSE, f1$P,lambda1)
      PPV_mat[simn,] = c(f2[[1]],f2[[2]],lambda1)
    }
    all_SSE = rbind(all_SSE, colMeans(SSE_mat))
    all_PPV = rbind(all_PPV, colMeans(PPV_mat))
  }
  SSE_lambda2 = rbind(SSE_lambda2,data.frame(all_SSE, lambda2))
  PPV_lambda2 = rbind(PPV_lambda2,data.frame(all_PPV, lambda2))
}

#-----------------------SSE extraction------------------

names(SSE_lambda2) = c("RCFGL", "RCFGL","lambda1","lambda2")
SSE_long_SSE = tidyr::gather(SSE_lambda2[,c(1,3,4)],Method, SSE,
                             1, factor_key=TRUE)
SSE_long_edges= tidyr::gather(SSE_lambda2[,c(2,3,4)],Method, Edges,
                              1, factor_key=TRUE)
SSE_long2 = merge(SSE_long_edges, SSE_long_SSE,by=c("lambda1","lambda2","Method"))

#-----------------------PPV extraction------------------

names(PPV_lambda2) = c("RCFGL", "RFGL","lambda1","lambda2")
Precision = tidyr::gather(PPV_lambda2[,c(1,3,4)],Method, Precision,
                          1, factor_key=TRUE)
Recall = tidyr::gather(PPV_lambda2[,c(2,3,4)],Method, Recall,
                       1, factor_key=TRUE)
PPV_long2 = merge(Precision, Recall,by=c("lambda1","lambda2","Method"))



##--------Extracting the results of FGL-----------
#--------------------------------------------------

SSE_lambda2 =  PPV_lambda2 = NULL

for(i in lambda2_unique){
  lambda2 = i;
  all_SSE = NULL
  all_PPV = NULL
  
  for(j in lambda1_unique){
    lambda1 = j;
    SSE_mat  = matrix(0,sim,3)
    PPV_mat = matrix(0,sim,3)
    for(simn in c(1:sim)){
      FGL = readRDS(paste0(folder_path, 
      "Simulated_datasets/Results/R_res/3condition/FGL_twosame_morediff_theta_",
                           "p_",p,"_n_",n,"_sim_",simn,
                           "_lamda1_",as.character(lambda1),
                           "_lamda2_", as.character(lambda2),".rdata"))
      
      c1 = FGL[[1]]$theta[[1]]; c2 = FGL[[1]]$theta[[2]]; c3 = FGL[[1]]$theta[[3]]; t1 = FGL[[2]]
      diff_1 = as.matrix(c1)-true_theta1;   diff_2 = as.matrix(c2)-true_theta2; diff_3 = as.matrix(c3)-true_theta3
      
      SSE = sum(diff_1[col(diff_1)>=row(diff_1)]^2)+sum(diff_2[col(diff_2)>=row(diff_2)]^2)+sum(diff_3[col(diff_3)>=row(diff_3)]^2)
      
      adj1 = make.adj.matrix(c1,trunc)
      adj2 = make.adj.matrix(c2,trunc)
      adj3 = make.adj.matrix(c3,trunc)
      
      f1 = TPR_FPR(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      f2 = PPV(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      SSE_mat[simn,] = c(SSE, f1$P,lambda1)
      PPV_mat[simn,] = c(f2[[1]],f2[[2]],lambda1)
    }
    all_SSE = rbind(all_SSE, colMeans(SSE_mat))
    all_TPR_FPR = rbind(all_TPR_FPR, colMeans(TPR_FPR_mat))
    all_PPV = rbind(all_PPV, colMeans(PPV_mat))
  }
  SSE_lambda2 = rbind(SSE_lambda2,data.frame(all_SSE, lambda2))
  TPR_FPR_lambda2 = rbind(TPR_FPR_lambda2,data.frame(all_TPR_FPR, lambda2))
  PPV_lambda2 = rbind(PPV_lambda2,data.frame(all_PPV, lambda2))
}

#-----------------------SSE extraction------------------

names(SSE_lambda2) = c("FGL", "FGL","lambda1","lambda2")
SSE_long_SSE = tidyr::gather(SSE_lambda2[,c(1,3,4)],Method, SSE,
                             1, factor_key=TRUE)
SSE_long_edges= tidyr::gather(SSE_lambda2[,c(2,3,4)],Method, Edges,
                              1, factor_key=TRUE)
SSE_long4 = merge(SSE_long_edges, SSE_long_SSE,by=c("lambda1","lambda2","Method"))

#-----------------------PPV extraction------------------

names(PPV_lambda2) = c("FGL", "FGL","lambda1","lambda2")
Precision = tidyr::gather(PPV_lambda2[,c(1,3,4)],Method, Precision,
                          1, factor_key=TRUE)
Recall = tidyr::gather(PPV_lambda2[,c(2,3,4)],Method, Recall,
                       1, factor_key=TRUE)
PPV_long4 = merge(Precision, Recall,by=c("lambda1","lambda2","Method"))

##--------Extracting the results of CFGL-----------
#--------------------------------------------------

SSE_lambda2 =  PPV_lambda2 = NULL

for(i in lambda2_unique){
  lambda2 = i;
  all_SSE = NULL
  all_PPV = NULL
  
  for(j in lambda1_unique){
    lambda1 = j;
    SSE_mat  = matrix(0,sim,3)
    PPV_mat = matrix(0,sim,3)
    for(simn in c(1:sim)){
      FGL = readRDS(paste0(folder_path, 
      "Simulated_datasets/Results/R_res/3condition/CFGL_twosame_morediff_theta_",
                           "p_",p,"_n_",n,"_sim_",simn,
                           "_lamda1_",as.character(lambda1),
                           "_lamda2_", as.character(lambda2),".rdata"))
      
      c1 = FGL[[1]]$theta[[1]]; c2 = FGL[[1]]$theta[[2]]; c3 = FGL[[1]]$theta[[3]]; t1 = FGL[[2]]
      diff_1 = as.matrix(c1)-true_theta1;   diff_2 = as.matrix(c2)-true_theta2; diff_3 = as.matrix(c3)-true_theta3
      
      SSE = sum(diff_1[col(diff_1)>=row(diff_1)]^2)+sum(diff_2[col(diff_2)>=row(diff_2)]^2)+sum(diff_3[col(diff_3)>=row(diff_3)]^2)
      
      adj1 = make.adj.matrix(c1,trunc)
      adj2 = make.adj.matrix(c2,trunc)
      adj3 = make.adj.matrix(c3,trunc)
      
      f1 = TPR_FPR(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      f2 = PPV(list(true_adj1,true_adj2, true_adj3),list(adj1,adj2,adj3))
      SSE_mat[simn,] = c(SSE, f1$P,lambda1)
      PPV_mat[simn,] = c(f2[[1]],f2[[2]],lambda1)
    }
    all_SSE = rbind(all_SSE, colMeans(SSE_mat))
    all_TPR_FPR = rbind(all_TPR_FPR, colMeans(TPR_FPR_mat))
    all_PPV = rbind(all_PPV, colMeans(PPV_mat))
  }
  SSE_lambda2 = rbind(SSE_lambda2,data.frame(all_SSE, lambda2))
  TPR_FPR_lambda2 = rbind(TPR_FPR_lambda2,data.frame(all_TPR_FPR, lambda2))
  PPV_lambda2 = rbind(PPV_lambda2,data.frame(all_PPV, lambda2))
}

#-----------------------SSE extraction------------------

names(SSE_lambda2) = c("CFGL", "CFGL","lambda1","lambda2")
SSE_long_SSE = tidyr::gather(SSE_lambda2[,c(1,3,4)],Method, SSE,
                             1, factor_key=TRUE)
SSE_long_edges= tidyr::gather(SSE_lambda2[,c(2,3,4)],Method, Edges,
                              1, factor_key=TRUE)
SSE_long5 = merge(SSE_long_edges, SSE_long_SSE,by=c("lambda1","lambda2","Method"))

#-----------------------PPV extraction------------------

names(PPV_lambda2) = c("CFGL", "CFGL","lambda1","lambda2")
Precision = tidyr::gather(PPV_lambda2[,c(1,3,4)],Method, Precision,
                          1, factor_key=TRUE)
Recall = tidyr::gather(PPV_lambda2[,c(2,3,4)],Method, Recall,
                       1, factor_key=TRUE)
PPV_long5 = merge(Precision, Recall,by=c("lambda1","lambda2","Method"))



##-------Coupling all the results together and save them neatly as rdata files-------
#--------------------------------------------------

SSE = rbind(SSE_long2, SSE_long3, SSE_long4, SSE_long5); SSE$lambda2 = as.factor(SSE$lambda2)
PPV = rbind(PPV_long2, PPV_long3, PPV_long4, PPV_long5); PPV$lambda2 = as.factor(PPV$lambda2)

saveRDS(SSE, paste0(folder_path, "Simulated_datasets/Extracted_Results/S2_SSE_results.rdata"))
saveRDS(PPV, paste0(folder_path, "Simulated_datasets/Extracted_Results/S2_PPV_results.rdata"))