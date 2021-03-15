#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 20:39:08 2021

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
    for iter in range(240):
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
     diff_value = 0
     for k in range(K):
       diff_value += np.sum(abs(P[:,:,k] - P_prev[:,:,k]))/np.sum(abs(P_prev[:,:,k]))
     
     funVal[1,iter] = diff_value
     funVal[2,iter] =  np.sum(abs(P-Z))