#---------Code to run RFGL and RCFGL on simulated datasets under scenario S1-------------

import os
import numpy as np
import pandas as pd
import time
import sys
sys.path.insert(0, '')
from scipy.sparse import csr_matrix
from scipy.sparse import save_npz

from sklearn.preprocessing import StandardScaler

RCFGL_path = '/Users/seals/Documents/GitHub/RCFGL' #address of the RCFGL package folder
os.chdir(RCFGL_path)
sys.path.insert(0, 'Python_functions')
from RCFGL import RFGL, RCFGL
os.chdir(RCFGL_path)
from Dstream_functions import*

import pyreadr
from sklearn.preprocessing import StandardScaler


folder_path = "/Users/seals/Desktop/CSPH/CFGL/31Aug/Codes_used_in_paper/" #address of the paper data/codes folder


lambdas = np.array([0.01,0.05,0.01,0.1,0.01,0.2,0.03,0.05,
          0.03,0.1,0.03,0.2,0.04,0.05,0.04,0.1,
          0.04,0.2,0.05,0.05,0.05,0.1,0.05,0.2,0.1,0.05,
          0.1,0.1,0.1,0.2,0.15,0.05,
          0.15,0.1,0.15,0.2,0.2,0.05,0.2,0.1,
          0.2,0.2,0.25,0.05,0.25,0.1,
          0.25,0.2]).reshape([24,2]) #different combinations of lambda1 and lambda2 for which the methods were run

#-----------Running RFGL for different combinations of lambda1 and lambda2------------

p = 500; n = 100;
for sim in range(10):
  sim = sim + 1
  for j in range(24): 
    lambda1 = lambdas[j,0]; lambda2 = lambdas[j,1];
    B = pd.read_csv(folder_path + 'Simulated_datasets/3condition/twosame_third_partially_same_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(1)+'_sim_'+str(sim)+'.csv');

    H = pd.read_csv(folder_path + 'Simulated_datasets/3condition/twosame_third_partially_same_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(2)+'_sim_'+str(sim)+'.csv');

    M = pd.read_csv(folder_path + 'Simulated_datasets/3condition/twosame_third_partially_same_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(3)+'_sim_'+str(sim)+'.csv');

    scalerB = StandardScaler(with_std=False); scalerB.fit(B)
    B = scalerB.transform(B)

    scalerH = StandardScaler(with_std=False);  scalerH.fit(H)
    H = scalerH.transform(H)

    scalerM = StandardScaler(with_std=False); scalerM.fit(M)
    M = scalerM.transform(M)
    
    A = []
    A.append(B);
    A.append(H);
    A.append(M);
    print(j)

    RFGL_output = RFGL(lambda1 = lambda1, lambda2 = lambda2, A = A, ADMMmaxiter = 100, admmtol = 0.001)

    for k in range(3):
     save_npz(folder_path+'Simulated_datasets/Results/Python_res/3condition/RFGL_twosame_third_partially_same_new_morediff_theta'+
                  'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(k)+'_sim_'+
                  str(sim)+'_condition_'+str(k)+'_lambda1_'+str(lambda1)+
                  '_lambda2_'+str(lambda2)+'.npz',csr_matrix(RFGL_output[0][:,:,k]))
  print(sim)

#-----------Running RCFGL for different combinations of lambda1 and lambda2------------

p = 500; n = 100;
for sim in range(10):
  sim = sim + 1
  for j in range(24): 
    lambda1 = lambdas[j,0]; lambda2 = lambdas[j,1];
    B = pd.read_csv(folder_path + 'Simulated_datasets/3condition/twosame_third_partially_same_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(1)+'_sim_'+str(sim)+'.csv');

    H = pd.read_csv(folder_path + 'Simulated_datasets/3condition/twosame_third_partially_same_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(2)+'_sim_'+str(sim)+'.csv');

    M = pd.read_csv(folder_path + 'Simulated_datasets/3condition/twosame_third_partially_same_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(3)+'_sim_'+str(sim)+'.csv');

    scalerB = StandardScaler(with_std=False); scalerB.fit(B)
    B = scalerB.transform(B)

    scalerH = StandardScaler(with_std=False);  scalerH.fit(H)
    H = scalerH.transform(H)

    scalerM = StandardScaler(with_std=False); scalerM.fit(M)
    M = scalerM.transform(M)
    
    A = []
    A.append(B);
    A.append(H);
    A.append(M);
    print(j)

    RCFGL_output = RCFGL(lambda1 = lambda1, lambda2 = lambda2, A = A, ADMMmaxiter = 100, admmtol = 0.001)

    for k in range(3):
     save_npz(folder_path+'Simulated_datasets/Results/Python_res/3condition/RCFGL_twosame_third_partially_same_new_morediff_theta'+
                  'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(k)+'_sim_'+
                  str(sim)+'_condition_'+str(k)+'_lambda1_'+str(lambda1)+
                  '_lambda2_'+str(lambda2)+'.npz',csr_matrix(RCFGL_output[0][:,:,k]))
  print(sim)

