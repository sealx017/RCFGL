#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 23:09:25 2021

@author: seals
"""


import cppyy
import os
os.chdir("/home/seals/CFGL/Codes/RCFGL/C_functions")
cppyy.include('fmgl_n.h')
from cppyy.gbl import fmgl_subfusedLasso_n
from cppyy.gbl import fmgl_subfusedLasso_diagonal

import numpy as np
import os
from scipy.linalg import blas


def FGL_ADMM(params,S,P,diff_tol):

    lambda1 = params[0];
    lambda2 = params[1];
    rho = params[2];
    p = S.shape[1]
    maxiter = params[4]
    admmtol = params[5];
    difftol = params[6];
    n = params[7];
    n1 = params[7];
    n = n/np.sum(n)
    K = S.shape[2]
    U = np.zeros((p,p,K));
    Z = np.zeros((p,p,K));
    W = np.zeros((p,p,K));
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
     #print(iter);
     P_prev = np.copy(P);

     for k in range(K):
      W[:,:,k] = -S[:,:,k] + np.multiply(rho/n[k],(Z[:,:,k] - U[:,:,k]));
      wd, V = np.linalg.eigh(W[:,:,k]);
      D = np.diag((wd + np.sqrt(wd**2 + np.dot(4,rho/n[k])))/np.dot(2,rho/n[k]));
      P[:,:,k] = blas.sgemm(1,blas.sgemm(1,V,D),V,trans_b = 1);

     W = P + U;
     fmgl_subfusedLasso_n(Z,W,fx,fy,lambda1/rho, lambda2/rho, p, K, int(np.dot((p+1),p/2)));
     U = U + (P - Z);

     funVal[0,iter] = computLogDet( P, S, K, n1, lambda1, lambda2);
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




def FGL_ADMM_unconnected(params,S,P,diff_tol):

    lambda1 = params[0];
    lambda2 = params[1];
    rho = params[2];
    p = S.shape[1]
    maxiter = params[4]
    admmtol = params[5];
    difftol = params[6];
    n = params[7];
    n1 = params[7];
    n = n/np.sum(n)
    K = S.shape[2]
    U = np.zeros((p,p,K));
    Z = np.zeros((p,p,K));
    W = np.zeros((p,p,K));
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
     #print(iter);
     P_prev = np.copy(P);

     for k in range(K):
      W[:,:,k] = -S[:,:,k] + np.multiply(rho/n[k],(Z[:,:,k] - U[:,:,k]));
      wd, V = np.linalg.eigh(W[:,:,k]);
      D = np.diag((wd + np.sqrt(wd**2 + np.dot(4,rho/n[k])))/np.dot(2,rho/n[k]));
      P[:,:,k] = blas.sgemm(1,blas.sgemm(1,V,D),V,trans_b = 1);

     W = P + U;
     fmgl_subfusedLasso_n(Z,W,fx,fy,lambda1/rho, lambda2/rho, p, K, int(np.dot((p+1),p/2)));
     U = U + (P - Z);

     funVal[0,iter] = computLogDet( P, S, K, n1, lambda1, lambda2);
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


def FGL_ADMM_unconnected_n(params,S,P,diff_tol):

    lambda1 = params[0];
    lambda2 = params[1];
    rho = params[2];
    p = S.shape[1]
    maxiter = params[4]
    admmtol = params[5];
    difftol = params[6];
    n = params[7];
    n1 = params[7];
    n = n/np.sum(n)
    K = S.shape[2]
    U = np.zeros((p,p,K));
    Z = np.zeros((p,p,K));
    W = np.zeros((p,p,K));
    fx1 = np.zeros((p,1));
    fy1 = np.zeros((p,1));
    iter = 0
    for i in range(p):
     fx1[iter,:] = np.array([i]);
     fy1[iter,:] = np.array([i]);
     iter = iter+1;


    iter = 0;
    funVal = np.zeros([3, int(maxiter)])
    for iter in range(maxiter):
     #print(iter);
     P_prev = np.copy(P);

     for k in range(K):
      W[:,:,k] = -np.diag(S[:,:,k]) + np.multiply(rho/n[k],(Z[:,:,k] - U[:,:,k]));
      wd = np.copy(np.diag(W[:,:,k]))
      P[:,:,k] = np.diag((wd + np.sqrt(wd**2 + np.dot(4,rho/n[k])))/np.dot(2,rho/n[k]));

     W = P + U;
     fmgl_subfusedLasso_diagonal(Z,W,fx1,fy1,lambda1/rho, lambda2/rho, p, K, p);
     U = U + (P - Z);

     #funVal[0,iter] = computLogDet_diag( P, S, K, n1, lambda1, lambda2);
     funVal[0,iter] = 0;

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


def computLogDet(P, S, K, n, lam1, lam2):

    td = 0;
    for i in range(K):
       td = td + n[i]*(np.log(np.linalg.det(P[:,:,i])) - np.sum(np.multiply(S[:,:,i],P[:,:,i])));

    td = td - np.dot(lam1,np.sum(abs(P)));
    P = P[:,:,range(K-2)] - P[:,:,range(1,K-1)];
    td = td - np.dot(lam2,np.sum(abs(P)));
    td = -td
    return td;


def computLogDet_diag(P, S, K, n, lam1, lam2):

    td = 0;
    for i in range(K):
       td = td + n[i]*(np.log(np.dot(np.diag(P[:,:,i]))) - np.sum(np.multiply(np.diag(S[:,:,i]),np.diag(P[:,:,i]))));

    td = td - np.dot(lam1,np.sum(abs(P)));
    P = P[:,:,range(K-2)] - P[:,:,range(1,K-1)];
    td = td - np.dot(lam2,np.sum(abs(P)));
    td = -td
    return td;
