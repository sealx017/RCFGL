#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 20:39:38 2021

@author: seals
"""

    lambda1 = params[0];
    lambda2 = params[1];
    rho = params[2];
    p = int(params[3]);
    maxiter = params[4]
    admmtol = params[5];
    difftol = params[6];
    K = S.shape[2]
    U = np.zeros((p,p,K));
    Z = np.zeros((p,p,K));
    fx = np.zeros(len(which_K[0]))
    fy = np.zeros(len(which_K[0]))
    for i in range(len(which_K[0])):
        fx[i] = which_K[0][i]
        fy[i] = which_K[1][i]

    fx_not = which_not_K[0]
    fy_not = which_not_K[1]

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
     fmgl_subfusedLasso_nc(Z,W,fx,fy,lambda1/rho, lambda2/rho, p, K, len(which_K[0]));    

     
     for i in fx_not:
      for j in fy_not:
       Z[i,j,:] = ptv_z.ptv_z(W[i,j,:],Weight[:,i,j],lambda1/rho, lambda2/rho, K)
       Z[j,i,:] = Z[i,j,:] 


     U = U + (P - Z);
       
     funVal[0,iter] = computLogDet( P, S, K, lambda1, lambda2);
     diff_value = 0
     for k in range(K):
       diff_value += np.sum(abs(P[:,:,k] - P_prev[:,:,k]))/np.sum(abs(P_prev[:,:,k]))
     
     funVal[1,iter] = diff_value
     funVal[2,iter] =  np.sum(abs(P-Z))