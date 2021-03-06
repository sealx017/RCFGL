"""
Created on Sun Mar  7 17:39:24 2021

@author: seals
"""

import cppyy
import os
os.chdir("/Users/seals/Documents/Github/RCFGL/C_functions")
cppyy.include('fmgl_n.h')
from cppyy.gbl import fmgl_subfusedLasso_n
import numpy as np
import os
from scipy.linalg import blas

def ADMM(params,S,P,diff_tol):

    lambda1 = params[0];
    lambda2 = params[1];
    rho = params[2];
    p = int(params[3]);
    maxiter = params[4]
    #pen_diag = params[5];
    admmtol = params[6];
    difftol = params[7];
    K = S.shape[2]
    U = np.zeros((p,p,K));
    Z = np.zeros((p,p,K));
    fx = np.zeros((int(np.multiply((p+1),p)/2),1));
    fy = np.zeros((int(np.multiply((p+1),p)/2),1));
    iter = 0
    for i in range(p):
     for j in range(i,p):
      fx[iter,:] = np.array([i]);
      fy[iter,:] = np.array([j]);
      iter = iter+1;
         
    iter = 0;      
    funVal = np.zeros([3, int(maxiter)])
    for iter in range(maxiter):
     print(iter);
     P_prev = np.copy(P);
     W = -S + np.multiply(rho,(Z - U));
     for k in range(K):
      wd, V = np.linalg.eigh(W[:,:,k]);
      D = np.diag((wd + np.sqrt(wd**2 + np.dot(4,rho)))/np.dot(2,rho));
      P[:,:,k] = blas.sgemm(1,blas.sgemm(1,V,D),V,trans_b = 1);
       
     W = P + U;      
     fmgl_subfusedLasso_n(Z,W,fx,fy,lambda1/rho, lambda2/rho, p, K, int(np.dot((p+1),p/2)));    
     U = U + (P - Z);
       
     funVal[0,iter] = computLogDet( P, S, K, lambda1, lambda2);
     #funVal[1,iter] =  np.sum(abs(P-Z));  
     #funVal[1,iter] =  np.sum(abs(P-P_prev));      
     #if (funVal[0,iter]) == np.float64('-inf') or abs(funVal[1,iter] - funVal[1,iter-1]) < admmtol:
      #break;
     #if abs(funVal[1,iter] - funVal[1,iter-1]) < admmtol:
     diff_value = 0
     for k in range(K):
       diff_value += np.sum(abs(P[:,:,k] - P_prev[:,:,k]))/np.sum(abs(P_prev[:,:,k]))
     
     funVal[1,iter] = diff_value
     funVal[2,iter] =  np.sum(abs(P-Z))
     if diff_tol == True:
      if abs(funVal[1,iter] - funVal[1,iter-1]) < difftol:
       break;  
     else:
      if abs(funVal[1,iter]) < admmtol:
       break;
     
    funVal = funVal[:,range(iter)]
    result = [];
    result.append(P)
    result.append(funVal)
    return(result)
    

def computLogDet(P, S, K, lam1, lam2):

    td = 0;
    for i in range(K):
       td -= np.log(np.linalg.det(P[:,:,i])) + np.sum(np.multiply(S[:,:,i],P[:,:,i]));
       
    td = td + np.dot(lam1,np.sum(abs(P)));
    P = P[:,:,range(K-2)] - P[:,:,range(1,K-1)];
    td = td + np.dot(lam2,np.sum(abs(P)));
    return td;
   
