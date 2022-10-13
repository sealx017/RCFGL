#---------Code to run FGL and CFGL on simulated datasets under scenario S1-------------

require(JGL)
require(CFGL)

folder_path = "/Users/seals/Desktop/CSPH/CFGL/31Aug/Paper_dataset_with_codes/" #address of the paper data/codes folder

p = 500; n = 100; # number of genes and subjects
lambdas = matrix(c(0.01,0.05,0.01,0.1,0.01,0.2,0.03,0.05,
                   0.03,0.1,0.03,0.2,0.04,0.05,0.04,0.1,
                   0.04,0.2,0.05,0.05,0.05,0.1,0.05,0.2,0.1,0.05,
                   0.1,0.1,0.1,0.2,0.15,0.05,
                   0.15,0.1,0.15,0.2,0.2,0.05,0.2,0.1,
                   0.2,0.2,0.25,0.05,0.25,0.1,
                   0.25,0.2), nrow = 24, ncol = 2, byrow = T) #different combinations of lambda1 and lambda2 for which the methods were run

for(sim in 1:10){
  y1 = read.csv(paste0("/home/seals/CFGL/Codes/April29/simulation_data/3condition/twosame_third_partially_same_morediff",
                       "p_",p,"_n_",n,"_condition_",1,"_sim_",sim,".csv"))
  
  y2 = read.csv(paste0("/home/seals/CFGL/Codes/April29/simulation_data/3condition/twosame_third_partially_same_morediff",
                       "p_",p,"_n_",n,"_condition_",2,"_sim_",sim,".csv"))
  
  y3 = read.csv(paste0("/home/seals/CFGL/Codes/April29/simulation_data/3condition/twosame_third_partially_same_morediff",
                       "p_",p,"_n_",n,"_condition_",3,"_sim_",sim,".csv"))
  
  maxiter = 200
  Y = NULL
  Y[[1]] = y1; Y[[2]] = y2; Y[[3]] = y3;
  
  #-----Run FGL------
  
  for(i in 1:dim(lambdas)[1]){
    lam1 = lambdas[i,1]; lam2 = lambdas[i,2]
    t3<-as.numeric(Sys.time())
    a1 = JGL(Y,penalty="fused",lambda1=lam1 ,lambda2=lam2 ,rho=1,
             weights="sample.size",penalize.diagonal=TRUE,
             maxiter=maxiter,warm=NULL,return.whole.theta=TRUE, screening="fast")
    t3<-as.numeric(Sys.time())-t3
    
    saveRDS(list(a1,t3),paste0(folder_path, "Simulated_datasets/Results/R_res/3condition/
                FGL_twosame_third_partially_same_morediff_theta_",
                "p_",p,"_n_",n,"_sim_",sim,
                "_lamda1_",as.character(lam1),
                "_lamda2_", as.character(lam2),".rdata"))
  }
  #-----Run CFGL------
  
  t2<-as.numeric(Sys.time())
  weight_12 = get_scr_mat(Y[[1]],Y[[2]])
  weight_13 = get_scr_mat(Y[[1]],Y[[3]])
  weight_23 = get_scr_mat(Y[[2]],Y[[3]])
  t2<-as.numeric(Sys.time())-t2
  
  for(i in 1:dim(lambdas)[1]){
    lam1 = lambdas[i,1]; lam2 = lambdas[i,2]
    t3<-as.numeric(Sys.time())
    a1 = CFGL(Y,lambda1=lam1,lambda2=lam2,maxiter = maxiter,loglik.trace=TRUE,
              btc.screening = list(weight_12$scr.mat,weight_13$scr.mat,weight_23$scr.mat))
    t3<-as.numeric(Sys.time())-t3
    saveRDS(list(a1,t2+t3),file = paste0(folder_path, "Simulated_datasets/Results/R_res/3condition/
              CFGL_twosame_third_partially_same_morediff_theta_", "p_",p,
              "_n_",n,"_sim_",sim,"_lamda1_",as.character(lam1), "_lamda2_", as.character(lam2),".rdata"))
  }
print(sim)
}

