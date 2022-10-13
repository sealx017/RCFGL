##------This R code creates the SSE and PPV plots of the main paper for scenarios S1 and S2--------
#----------------------------------------------------------------------------------------------

library(Matrix)
folder_path = "/Users/seals/Desktop/CSPH/CFGL/31Aug/Paper_dataset_with_codes/" #the location of the folder on the hard drive

#-------Loading the SSE and PPV files of scenario S1----------
SSE = readRDS(paste0(folder_path,"Simulated_datasets/Extracted_Results/S1_SSE_results.rdata"))
PPV = readRDS(paste0(folder_path,"Simulated_datasets/Extracted_Results/S1_PPV_results.rdata"))
SSE = cbind(SSE, "(S1)"); PPV = cbind(PPV, "(S1)")
#-------Loading the SSE and PPV files of scenario S4----------
SSE2 = readRDS(paste0(folder_path,"Simulated_datasets/Extracted_Results/S2_SSE_results.rdata"))
PPV2 = readRDS(paste0(folder_path,"Simulated_datasets/Extracted_Results/S2_PPV_results.rdata"))
SSE2 = cbind(SSE2, "(S2)"); PPV2 = cbind(PPV2, "(S2)")

colnames(SSE)[6] = colnames(PPV)[6] = colnames(SSE2)[6] = colnames(PPV2)[6]  = "case"

#-------Combining the results of S3 and S4 to create the final dataframes-------------
SSE3 = rbind(SSE, SSE2)
PPV3 = rbind(PPV, PPV2)

lam2.labs <- c("lambda2 = 0.05", "lambda2 = 0.1", "lambda2 = 0.2")
names(lam2.labs)<-levels(PPV3$lambda2)
color_codes = c("dodgerblue1","orangered1", "mediumorchid1","seagreen4")

#------Plotting SSE---------------
require(ggplot2)
p = 500; n = 100
png(paste0(folder_path, "Plots/S1_S2_SSE_",
           "p_",p,"_n_",n,".png"),
    height = 1200, width = 1400, res = 220)
ggplot(SSE3, aes(x=Edges, y=SSE,color = Method)) +
  geom_point(size=0, shape=0) + geom_line(aes(linetype=Method)) +  
  theme(legend.position="bottom")+
  labs(x = "Number of edges selected")+scale_color_manual(values=color_codes)+
  facet_grid(case~lambda2,labeller = labeller(label_both,lambda2 = lam2.labs))
dev.off()

#------Plotting PPV---------------
png(paste0(folder_path, "Plots/S1_S2_PPV_",
           "p_",p,"_n_",n,".png"),
    height = 1200, width = 1400,res = 220)
ggplot(PPV3, aes(x=Recall, y=Precision,color = Method)) + 
  scale_color_manual(values=color_codes)+
  geom_point(size=0, shape=0) + geom_line(aes(linetype=Method)) + 
  theme(legend.position="bottom", axis.text.x = element_text(size = 8, 
  angle = 90, vjust = 0, hjust=0))+
  facet_grid(case~lambda2,labeller = labeller(label_both,lambda2 = lam2.labs))
dev.off()
