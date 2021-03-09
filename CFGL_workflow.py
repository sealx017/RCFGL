
import os
import numpy as np
import pandas as pd
os.chdir("/Users/seals/Documents/Github/RCFGL/Python_functions")
import ADMM_py_function as AP
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
lambda1 = 0.05
lambda2 = 0.005;
rho = 1;
ADMMmaxiter = 250;
pen_diag = "True"
admmtol = 0.00001
params = []
params.extend((lambda1, lambda2, rho, p, ADMMmaxiter, pen_diag, admmtol))

A = []
A.append(Acbc);
A.append(IL);
A.append(LHB);

K = len(A);

S = np.zeros((p,p,K));
P = np.zeros((p,p,K));

trueSparsity  = 0;

for k in np.array(range(K)):
 S[:,:,k] = np.cov((A[k]).T);
 P[:,:,k] = np.diag(1/np.diag(S[:,:,k]));

#Simple FGL ADMM

[P_ADMM, funVal] = AP.ADMM(params,S,P)


#Pairwise screening matrix computation

expr1 = np.array(A[0])
expr2 = np.array(A[1])

W = scr.get_diff_W(expr1,expr2,s=2)