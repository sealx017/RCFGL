#---------Code to run FGL on simulated datasets under scenario S3-------------

require(JGL)

folder_path = "/Users/seals/Desktop/CSPH/CFGL/31Aug/Codes_used_in_paper/" #address of the paper data/codes folder

p = 500; n = 100;
lambdas = matrix(c(0.01,0.05,0.01,0.1,0.01,0.2,0.03,0.05,
                   0.03,0.1,0.03,0.2,0.04,0.05,0.04,0.1,
                   0.04,0.2,0.05,0.05,0.05,0.1,0.05,0.2,0.1,0.05,
                   0.1,0.1,0.1,0.2,0.15,0.05,
                   0.15,0.1,0.15,0.2,0.2,0.05,0.2,0.1,
                   0.2,0.2,0.25,0.05,0.25,0.1,
                   0.25,0.2),nrow = 24, ncol = 2, byrow = T)  #different combinations of lambda1 and lambda2 for which the methods were run
for(sim in 1:10){
y1 = read.csv(paste0(folder_path, "Simulated_datasets/4condition/
                     two_twosame_morediff",
                     "p_",p,"_n_",n,"_condition_",1,"_sim_",sim,".csv"))

y2 = read.csv(paste0(folder_path, "Simulated_datasets/4condition/
                     two_twosame_morediff",
                     "p_",p,"_n_",n,"_condition_",2,"_sim_",sim,".csv"))

y3 = read.csv(paste0(folder_path, "Simulated_datasets/4condition/
                     two_twosame_morediff",
                     "p_",p,"_n_",n,"_condition_",3,"_sim_",sim,".csv"))

y4 = read.csv(paste0(folder_path, "Simulated_datasets/4condition/
                     two_twosame_morediff",
                     "p_",p,"_n_",n,"_condition_",4,"_sim_",sim,".csv"))
maxiter = 200
Y = NULL
Y[[1]] = y1; Y[[2]] = y2; Y[[3]] = y3; Y[[4]] = y4

for(i in 1:dim(lambdas)[1]){
  lam1 = lambdas[i,1]; lam2 = lambdas[i,2]
  t3<-as.numeric(Sys.time())
  a1 = JGL(Y,penalty="fused",lambda1=lam1 ,lambda2=lam2 ,rho=1,
           weights="sample.size",penalize.diagonal=TRUE,
           maxiter=maxiter,warm=NULL,return.whole.theta=TRUE, screening="fast")
  t3<-as.numeric(Sys.time())-t3
  
  saveRDS(list(a1,t3),paste0(folder_path, "Simulated_datasets/Results/R_res/
                             4condition/FGL_two_twosame_morediff_theta_",
                             "p_",p,"_n_",n,"_sim_",sim,
                             "_lamda1_", as.character(lam1),
                             "_lamda2_", as.character(lam2),".rdata"))
}
print(sim)
}

