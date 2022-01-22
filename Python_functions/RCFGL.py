import os
import numpy as np
import time
import sys
#sys.path.insert(0, '')

#os.chdir("/Users/seals/Documents/GitHub/RCFGL/Python_functions")
import screening_by_thm as scthm

#os.chdir("/Users/seals/Documents/GitHub/RCFGL/Python_functions")
from RCFGL_RFGL_base_functions import CFGL_ADMM, FGL_ADMM, FGL_ADMM_unconnected_n

#os.chdir("/Users/seals/Documents/GitHub/RCFGL/Python_functions")
import get_screening_2 as scr


def RFGL(lambda1, lambda2, A, ADMMmaxiter = 100,trueSparsity  = 0, truncate=1e-05, rho = 1, admmtol = 1e-4, difftol = 1e-4):
    K = len(A);
    p = A[0].shape[1];
    S = np.zeros((p,p,K));
    P = np.zeros((p,p,K));
    n = np.zeros(K);
    for k in np.array(range(K)):
         n[k] = A[k].shape[0];
         S[:,:,k] = np.dot(np.cov((A[k].T)),(n[k]-1)/n[k]);
         P[:,:,k] = np.diag(1/np.diag(S[:,:,k]));

    weights = n/np.sum(n)
    start_time = time.time()
    screening = scthm.thm_screening(lambda1,S,weights)
    unconnected = np.zeros((1,1))
    blocklist = []
    for i in range(screening[0]):
      if np.sum(screening[1]==i)==1:
              unconnected = np.hstack((unconnected,np.where(screening[1]==i)))
      if np.sum(screening[1]==i)>1:
          blocklist.append(np.where(screening[1]==i))
    
    unconnected = unconnected[0,range(1,unconnected.shape[1])].astype(int)
    t1 = (time.time() - start_time)
    start_time = time.time()
    params = []
    params.extend((lambda1, lambda2, rho, p, ADMMmaxiter, admmtol, difftol, n))
    
    theta_conn = []
    for i in range(len(blocklist)):
     blockS = np.zeros((len(blocklist[i][0]),len(blocklist[i][0]),K))
     blockP = np.zeros((len(blocklist[i][0]),len(blocklist[i][0]),K))
     for k in np.array(range(K)):
      blockS[:,:,k] = S[:,:,k][np.ix_(blocklist[i][0],blocklist[i][0])]
      blockP[:,:,k] = P[:,:,k][np.ix_(blocklist[i][0],blocklist[i][0])]
     [conn, fun_conn] = FGL_ADMM(params,blockS,blockP,diff_tol = False)
     theta_conn.append(conn)
    
    
    del conn;
    if len(unconnected)>0:
     S_unconn = np.zeros((len(unconnected),len(unconnected),K))
     P_unconn = np.zeros((len(unconnected),len(unconnected),K))
     for k in np.array(range(K)):
       S_unconn[:,:,k] = S[:,:,k][np.ix_(unconnected,unconnected)]
       P_unconn[:,:,k] = P[:,:,k][np.ix_(unconnected,unconnected)]
     [unconn, fun_unconn] = FGL_ADMM_unconnected_n(params,S_unconn,P_unconn,diff_tol = False)
    
    final_theta = np.zeros((p,p,K));
    #print(len(unconnected))
    for k in np.array(range(K)):
     if len(unconnected)>0:
      final_theta[:,:,k][np.ix_(unconnected,unconnected)] = unconn[:,:,k]
     for i in range(len(blocklist)):
      final_theta[:,:,k][np.ix_(blocklist[i][0],blocklist[i][0])] = theta_conn[i][:,:,k]
    
    roundoff = (np.abs(final_theta)<truncate)
    final_theta[roundoff] = 0
    logdet = np.zeros((K,1))

    #for k in np.array(range(K)):
    # nonz[k] = len(np.where(np.triu(final_theta[:,:,k],0)!=0)[0])

    AIC = 0
    for k in range(K):
     logdet[k] = np.linalg.slogdet(final_theta[:,:,k])[1]
     nonz = (np.nonzero(final_theta[:,:,k])[0].shape[0]-p)/2 + p
     AIC = AIC + weights[k]*(np.sum(np.multiply(S[:,:,k], final_theta[:,:,k]))-logdet[k])+2*nonz
     
    t2 = (time.time() - start_time)
    list = []
    list.append(final_theta)
    list.append(AIC)
    list.append(t1+t2)
    #list.append(logdet)
    return list




def RCFGL(lambda1, lambda2, A, ADMMmaxiter = 100,trueSparsity  = 0, truncate=1e-05, rho = 1, admmtol = 1e-4, difftol = 1e-4):
    K = len(A);
    p = A[0].shape[1];
    S = np.zeros((p,p,K));
    P = np.zeros((p,p,K));
    n = np.zeros(K);
    for k in np.array(range(K)):
         n[k] = A[k].shape[0];
         S[:,:,k] = np.dot(np.cov((A[k].T)),(n[k]-1)/n[k]);
         P[:,:,k] = np.diag(1/np.diag(S[:,:,k]));

    weights = n/np.sum(n)
    start_time = time.time()
    screening = scthm.thm_screening(lambda1,S,weights)
    unconnected = np.zeros((1,1))
    blocklist = []
    for i in range(screening[0]):
      if np.sum(screening[1]==i)==1:
              unconnected = np.hstack((unconnected,np.where(screening[1]==i)))
      if np.sum(screening[1]==i)>1:
          blocklist.append(np.where(screening[1]==i))

    newblocklist = []
    for i in range(len(blocklist)):
      if len(blocklist[i][0])<5 and i == (len(blocklist)-1):
          newblocklist.insert(i-1, np.append(newblocklist[i-1],blocklist[i][0]))
          del newblocklist[i]
      else:
          newblocklist.insert(i, blocklist[i][0])

    unconnected = unconnected[0,range(1,unconnected.shape[1])].astype(int)
    t1 = (time.time() - start_time)
    start_time = time.time()
    params = []
    params.extend((lambda1, lambda2, rho, p, ADMMmaxiter, admmtol, difftol, n))

    theta_conn = []
    for i in range(len(newblocklist)):
     len_block = len(newblocklist[i])
     blockS = np.zeros((len_block,len_block,K))
     blockP = np.zeros((len_block,len_block,K))
     blockWeight = np.zeros((K-1,len_block,len_block))
     blockA = []
     for k in np.array(range(K)):
         blockA.append(A[k][:,newblocklist[i]]);
         
     for k in range(K-1):
      blockWeight[k] = scr.get_scr_mat(np.array(blockA[k]),np.array(blockA[k+1]),alpha = 0.6)


     for k in np.array(range(K)):
      blockS[:,:,k] = S[:,:,k][np.ix_(newblocklist[i],newblocklist[i])]
      blockP[:,:,k] = P[:,:,k][np.ix_(newblocklist[i],newblocklist[i])]


     Sum_Weight = np.sum(blockWeight,axis = 0)
     lower_indices = np.tril_indices(len_block,k = -1)
     Sum_Weight[lower_indices] = -1
     which_K_block = np.where(Sum_Weight==K-1)
     which_not_K_block = np.where((Sum_Weight!=K-1) & (Sum_Weight!=-1))
     [conn, fun_conn] = CFGL_ADMM(params,blockS,blockP,blockWeight,which_K_block,which_not_K_block,diff_tol = False)
     theta_conn.append(conn)
     del conn;

    if len(unconnected)>0:
     S_unconn = np.zeros((len(unconnected),len(unconnected),K))
     P_unconn = np.zeros((len(unconnected),len(unconnected),K))

     for k in np.array(range(K)):
       S_unconn[:,:,k] = S[:,:,k][np.ix_(unconnected,unconnected)]
       P_unconn[:,:,k] = P[:,:,k][np.ix_(unconnected,unconnected)]

     [unconn, fun_unconn] =  FGL_ADMM_unconnected_n(params,S_unconn,P_unconn,diff_tol = False)

    final_theta = np.zeros((p,p,K));
    for k in np.array(range(K)):
     if len(unconnected)>0:
      final_theta[:,:,k][np.ix_(unconnected,unconnected)] = unconn[:,:,k]
     for i in range(len(newblocklist)):
      final_theta[:,:,k][np.ix_(newblocklist[i],newblocklist[i])] = theta_conn[i][:,:,k]

    roundoff = (np.abs(final_theta)<truncate)
    final_theta[roundoff] = 0
    logdet = np.zeros((K,1))

    #for k in np.array(range(K)):
    # nonz[k] = len(np.where(np.triu(final_theta[:,:,k],0)!=0)[0])+p

    AIC = 0
    for k in range(K):
     logdet[k] = np.linalg.slogdet(final_theta[:,:,k])[1]
     nonz = (np.nonzero(final_theta[:,:,k])[0].shape[0]-p)/2 + p
     #print(nonz)
     AIC = AIC + weights[k]*(np.sum(np.multiply(S[:,:,k], final_theta[:,:,k]))-logdet[k])+2*nonz
     
    t2 = (time.time() - start_time)
    list = []
    list.append(final_theta)
    list.append(AIC)
    list.append(t1+t2)
    #list.append(logdet)
    return list

