
#---------Code to run RFGL and RCFGL on simulated datasets under scenario S4-------------

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


lambdas = np.array([0.01,0.05,0.01,0.1,0.01,0.2,
                   0.02,0.05,0.02,0.1,0.02,0.2,
                   0.03,0.05,0.03,0.1,0.03,0.2,
                   0.04,0.05,0.04,0.1,0.04,0.2,
                   0.05,0.05,0.05,0.1,0.05,0.2,
                   0.06,0.05,0.06,0.1,0.06,0.2,
                   0.07,0.05,0.07,0.1,0.07,0.2,
                   0.08,0.05,0.08,0.1,0.08,0.2,
                   0.1,0.05, 0.1,0.1,0.1,0.2,
                   0.15,0.05, 0.15,0.1,0.15,0.2,
                   0.2,0.05, 0.2,0.1,0.2,0.2]).reshape([33,2]) #different combinations of lambda1 and lambda2 for which the methods were run

#-----------Running RFGL for different combinations of lambda1 and lambda2------------

p = 500; n = 100;
for sim in range(10):
  sim = sim + 1
  for j in range(33): 
    lambda1 = lambdas[j,0]; lambda2 = lambdas[j,1];
    B = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(1)+'_sim_'+str(sim)+'.csv');

    H = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(2)+'_sim_'+str(sim)+'.csv');

    M = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(3)+'_sim_'+str(sim)+'.csv');

    R = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(4)+'_sim_'+str(sim)+'.csv');


    scalerB = StandardScaler(with_std=False); scalerB.fit(B)
    B = scalerB.transform(B)

    scalerH = StandardScaler(with_std=False); scalerH.fit(H)
    H = scalerH.transform(H)

    scalerM = StandardScaler(with_std=False); scalerM.fit(M)
    M = scalerM.transform(M)

    scalerR = StandardScaler(with_std=False); scalerR.fit(R)
    R = scalerR.transform(R)
    
    A = []
    A.append(B);
    A.append(H);
    A.append(M);
    A.append(R);
    print(j)

    RFGL_output = RFGL(lambda1 = lambda1, lambda2 = lambda2, A = A, ADMMmaxiter = 200, admmtol = 0.001)

    for k in range(4):
     save_npz(folder_path+'Simulated_datasets/Results/Python_res/4condition/RFGL_twosame_morediff_theta'+
                  'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(k)+'_sim_'+
                  str(sim)+'_condition_'+str(k)+'_lambda1_'+str(lambda1)+
                  '_lambda2_'+str(lambda2)+'.npz',csr_matrix(RFGL_output[0][:,:,k]))
  print(sim)

#-----------Running RCFGL for different combinations of lambda1 and lambda2------------

p = 500; n = 100;
for sim in range(10):
  sim = sim + 1
  for j in range(33): 
    lambda1 = lambdas[j,0]; lambda2 = lambdas[j,1];
    B = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(1)+'_sim_'+str(sim)+'.csv');

    H = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(2)+'_sim_'+str(sim)+'.csv');

    M = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                    'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(3)+'_sim_'+str(sim)+'.csv');

    R = pd.read_csv(folder_path + 'Simulated_datasets/4condition/twosame_morediff'+
                'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(4)+'_sim_'+str(sim)+'.csv');


    scalerB = StandardScaler(with_std=False); scalerB.fit(B)
    B = scalerB.transform(B)

    scalerH = StandardScaler(with_std=False); scalerH.fit(H)
    H = scalerH.transform(H)

    scalerM = StandardScaler(with_std=False); scalerM.fit(M)
    M = scalerM.transform(M)

    scalerR = StandardScaler(with_std=False); scalerR.fit(R)
    R = scalerR.transform(R)
    
    A = []
    A.append(B);
    A.append(H);
    A.append(M);
    A.append(R);
    print(j)

    RCFGL_output = RCFGL(lambda1 = lambda1, lambda2 = lambda2, A = A, ADMMmaxiter = 200, admmtol = 0.001)

    for k in range(4):
     save_npz(folder_path+'Simulated_datasets/Results/Python_res/4condition/RCFGL_twosame_morediff_theta'+
                  'p_'+str(p)+'_n_'+str(n)+'_condition_'+str(k)+'_sim_'+
                  str(sim)+'_condition_'+str(k)+'_lambda1_'+str(lambda1)+
                  '_lambda2_'+str(lambda2)+'.npz',csr_matrix(RCFGL_output[0][:,:,k]))
  print(sim)
 
