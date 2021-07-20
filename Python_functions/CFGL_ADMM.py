import cppyy
import os
os.chdir("/Users/seals/Documents/Github/RCFGL/C_functions")
cppyy.include('fmgl_n_copy.h')
from cppyy.gbl import fmgl_subfusedLasso_nc
os.chdir("/Users/seals/Documents/Github/RCFGL/Python_functions")
import numpy as np
from scipy.linalg import blas
import prox_tv_Z_solver as ptv_z

def CFGL_ADMM(params,S,P,Weight,which_K,which_not_K,diff_tol):

    lambda1 = params[0];
    lambda2 = params[1];
    rho = params[2];
    p = S.shape[1];
    maxiter = params[4];
    admmtol = params[5];
    difftol = params[6];
    n = params[7];
    n1 = params[7];
    n = n/np.sum(n)
    K = S.shape[2]
    U = np.zeros((p,p,K));
    Z = np.zeros((p,p,K));
    W = np.zeros((p,p,K));
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
     #W = -S + np.multiply(rho,(Z - U));
     for k in range(K):
      W[:,:,k] = -S[:,:,k] + np.multiply(rho/n[k],(Z[:,:,k] - U[:,:,k]));
      wd, V = np.linalg.eigh(W[:,:,k]);
      D = np.diag((wd + np.sqrt(wd**2 + np.dot(4,rho/n[k])))/np.dot(2,rho/n[k]));
      P[:,:,k] = blas.sgemm(1,blas.sgemm(1,V,D),V,trans_b = 1);
       
     W = P + U;    
     fmgl_subfusedLasso_nc(Z,W,fx,fy,lambda1/rho, lambda2/rho, p, K, len(which_K[0]));    

     
     for m in range(len(fx_not)):
       i = fx_not[m]; j = fy_not[m];  
       Z[i,j,:] = ptv_z.ptv_z(W[i,j,:],Weight[:,i,j],lambda1/rho, lambda2/rho, K)
       Z[j,i,:] = Z[i,j,:] 


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
    


def computLogDet(P, S, K, n, lam1, lam2):

    td = 0;
    for i in range(K):
       td = td + n[i]*(np.log(np.linalg.det(P[:,:,i])) - np.sum(np.multiply(S[:,:,i],P[:,:,i])));
       
    td = td - np.dot(lam1,np.sum(abs(P)));
    P = P[:,:,range(K-2)] - P[:,:,range(1,K-1)];
    td = td - np.dot(lam2,np.sum(abs(P)));
    td = -td
    return td;
   
   
