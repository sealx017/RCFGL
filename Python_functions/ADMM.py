#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 17:39:24 2021

@author: seals
"""

import cppyy
import os
os.chdir("/Users/seals/Documents/Github/RCFGL")
cppyy.include('fmgl_n.h')


from cppyy.gbl import fmgl_subfusedLasso_n
import numpy as np
import pandas as pd
import os
from scipy.linalg import blas

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
#sigma = 1e-3;
#maxlineiter = 50;
#maxiter = 50;
#paramSPGmaxiter = 250;
ADMMmaxiter = 250;
pen_diag = "True"
admmtol = 0.00001
params = []
params.extend((lambda1, lambda2, rho, p, ADMMmaxiter, pen_diag, admmtol))

A = []
A.append(Acbc);
A.append(IL);
#A.append(LHB);

K = len(A);

S = np.zeros((p,p,K));
P = np.zeros((p,p,K));

trueSparsity  = 0;

for k in np.array(range(K)):
 S[:,:,k] = np.cov((A[k]).T);
 P[:,:,k] = np.diag(1/np.diag(S[:,:,k]));


def computLogDet(P, S, K, lam1, lam2):

    td = 0;
    for i in range(K):
       td = td - np.log(np.linalg.det(P[:,:,i])) + np.sum(np.multiply(S[:,:,i],P[:,:,i]));
       
    td = td + np.dot(lam1,np.sum(abs(P)));
    P = P[:,:,range(K-2)] - P[:,:,range(1,K-1)];
    td = td + np.dot(lam2,np.sum(abs(P)));
    return td;


lambda1 = params[0];
lambda2 = params[1];
rho = params[2];
p = int(params[3]);
maxiter = params[4]
pen_diag = params[5];
admmtol = params[6];
K = S.shape[2]
U = np.zeros((p,p,K));
Z = np.zeros((p,p,K));
fx = np.zeros((int(np.multiply((p+1),p)/2),1));
fy = np.zeros((int(np.multiply((p+1),p)/2),1));

iter = 0;
   
for i in range(p):
 for j in range(i,p):
  fx[iter,:] = np.array([i]);
  fy[iter,:] = np.array([j]);
  iter = iter+1;
  
funVal = np.zeros([3, int(maxiter)])
#iter = 0
for iter in range(maxiter):
 W = -S + np.multiply(rho,(Z - U));
 for k in range(K):
  wd, V = np.linalg.eigh(W[:,:,k]);
  D = np.diag((wd + np.sqrt(wd**2 + np.dot(4,rho)))/np.dot(2,rho));
  P[:,:,k] = blas.sgemm(1,blas.sgemm(1,V,D),V,trans_b = 1);
   
 W = P + U;
 #W = W.reshape(W.shape,order='F')
 Wt = W.reshape(-1,order='C').astype("float32")  
  
 fmgl_subfusedLasso_n(Z,W,fx,fy,lambda1/rho, lambda2/rho, p, K, int(np.dot((p+1),p/2)));
 #Z = Z.reshape(Z.shape,order='F')

 U = U + (P - Z);
   
 funVal[0,iter] = computLogDet( P, S, K, lambda1, lambda2);
 funVal[1,iter] =  np.sum(abs(P-Z));      
 if (funVal[0,iter]) == np.float64('-inf') or abs(funVal[1,iter] - funVal[1,iter-1]) < admmtol:
  break;
 #Zo = Z;
 print(iter); 
 
funVal = funVal[:,range(iter+1)]