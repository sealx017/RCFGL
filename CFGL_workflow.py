
import os
import numpy as np
import pandas as pd
os.chdir("/Users/seals/Documents/Github/RCFGL/Python_functions")
import ADMM_py_function_new as AP

os.chdir("/Users/seals/Documents/Github/RCFGL/Python_functions")
import CFGL_ADMM as CFGL_AP

os.chdir("/Users/seals/Documents/Github/RCFGL/Python_functions")
import get_screening as scr

os.chdir("/Users/seals/Desktop/CSPH/CFGL/PyModule")

which_dat = '500';
Acbc = pd.read_csv('/Users/seals/Desktop/CSPH/CFGL/MATLAB_data/'+which_dat+'_Acbc_better.csv');
Acbc = Acbc.iloc[:,1:Acbc.shape[1]]

IL = pd.read_csv('/Users/seals/Desktop/CSPH/CFGL/MATLAB_data/'+which_dat+'_IL_better.csv');
IL = IL.iloc[:,1:IL.shape[1]]

LHB = pd.read_csv('/Users/seals/Desktop/CSPH/CFGL/MATLAB_data/'+which_dat+'_LHB_better.csv');
LHB = LHB.iloc[:,1:LHB.shape[1]]

p = Acbc.shape[1];
lambda1 = 0.01
lambda2 = 0.05;
rho = 1;
ADMMmaxiter = 250;
#pen_diag = "True";
admmtol = 1e-4;
difftol = 1e-4;
params = []
params.extend((lambda1, lambda2, rho, p, ADMMmaxiter, admmtol, difftol))

A = []
A.append(Acbc);
A.append(IL);
A.append(LHB);

K = len(A);

S = np.zeros((p,p,K));
P = np.zeros((p,p,K));
n = np.zeros(K);
trueSparsity  = 0;

for k in np.array(range(K)):
 n[k] = A[k].shape[0];   
 S[:,:,k] = np.dot(np.cov((A[k]).T),(n[k]-1)/n[k]);
 P[:,:,k] = np.diag(1/np.diag(S[:,:,k]));




#Simple FGL ADMM

[P_ADMM, funVal] = AP.FGL_ADMM(params, S, P, diff_tol = False)


#CFGL ADMM

#Pairwise screening matrix computation

Weight = np.zeros((2,p,p))

for k in range(K-1):
    Weight[k] = scr.get_scr_mat(np.array(A[k]),np.array(A[k+1]))

Sum_Weight = np.sum(Weight,axis = 0)
lower_indices = np.tril_indices(p,k = -1)
Sum_Weight[lower_indices] = -1
which_K = np.where(Sum_Weight==K-1)
which_not_K = np.where((Sum_Weight!=K-1) & (Sum_Weight!=-1))

#Fitting CFGL with the computed screening matrices

[P_CFGL_ADMM, funVal] = CFGL_AP.CFGL_ADMM(params, S, P, Weight, which_K, which_not_K, diff_tol = False)



