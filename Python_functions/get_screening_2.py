

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 18:37:35 2021

@author: seals
"""

import numpy as np
import cppyy
import os
from scipy.stats import norm
cppyy.include('C_functions/for_loop.h')
from cppyy.gbl import screening_loop
from sklearn.linear_model import ElasticNet
 

def s_selection(expr1,expr2,ss=np.arange(0.1,2.1,0.1)):
  p = expr1.shape[1]
  temp1 = norm.cdf(np.sqrt(np.log(p)))
  d = np.zeros((ss.shape[0],10))
  
  for i in range(ss.shape[0]):
    W = get_diff_W(expr1,expr2,s = ss[i])
    for j in range(10):
      temp2 = (1-temp1)*(j+1)/10
      check = abs(W)>norm.ppf(1-temp2)
      check = check.reshape(-1)
      nomi = np.sum(check)
      deno = temp2*p*(p-1)
      d[i,j] = (nomi/deno-1)**2
  s_selected = ss[np.argmin(np.sum(d,axis = 1))]
  return(s_selected)



def get_diff_W(expr1,expr2,s=2):
     
  n1 = expr1.shape[0]
  n2 = expr2.shape[0]
  p = expr1.shape[1]

  expr1t = expr1.T
  expr2t = expr2.T
  covm1 = np.cov(expr1t)
  covm2 = np.cov(expr2t)
  sigma1 = np.diag(covm1)
  sigma2 = np.diag(covm2)
  
  
  b1 = np.zeros([(p-1),p])
  b2 = np.zeros([(p-1),p])
  c1 = np.zeros([n1,p])
  c2 = np.zeros([n2,p])
  r1 = np.zeros([p,p])
  r2 = np.zeros([p,p])
  W = np.zeros([p,p])
  #Num = int(np.multiply((p+1),p)/2)
  
  for i in range(p):
    y1 = expr1[:,i]
    x1 = np.delete(expr1,i,axis = 1)
    lam1 = ((s*np.sqrt(sigma1[i]*np.log(p)/n1)))    
    x1_s = (x1-np.mean(x1,axis = 0)).T*(1/(np.std(x1, axis = 0)))[:, np.newaxis]
    y1_s = (y1-np.mean(y1,axis = 0))
    regr = ElasticNet(alpha = lam1,l1_ratio = 1,fit_intercept = False,random_state=0, normalize = False)
    regr.fit(x1_s.T,y1_s)
    b1[:,i] = np.multiply(1/(np.std(x1, axis = 0)),regr.coef_)
    

    y2 = expr2[:,i]
    x2 = np.delete(expr2,i,axis = 1)
    lam2 = (s*np.sqrt(sigma2[i]*np.log(p)/n2))
    x2_s = (x2-np.mean(x2,axis = 0)).T*(1/(np.std(x2, axis = 0)))[:, np.newaxis]
    y2_s = (y2-np.mean(y2,axis = 0))
    regr = ElasticNet(alpha = lam2,l1_ratio = 1,fit_intercept = False,random_state=0, normalize = False)
    regr.fit(x2_s.T,y2_s)
    b2[:,i] = np.multiply(1/(np.std(x2, axis = 0)),regr.coef_)

    c1[:,i] = y1-np.mean(y1)-np.dot((x1-np.mean(x1,axis = 0)),b1[:,i])    
    c2[:,i] = y2-np.mean(y2)-np.dot((x2-np.mean(x2,axis = 0)),b2[:,i])

  r1 = np.dot(c1.T,c1)/n1
  r2 = np.dot(c2.T,c2)/n2
  s1 = np.mean(np.power(c1,2), axis = 0).reshape((p,1))
  s2 = np.mean(np.power(c2,2), axis = 0).reshape((p,1))
  screening_loop(W, r1, s1, b1, r2, s2, b2, p, n1, n2)
  W = W+W.T
  return(W)


def get_W_theshold(W,alpha):
    
  p = W.shape[0]
  t_upper = 2*np.sqrt(np.log(p))
  t0 = abs(W[np.triu_indices_from(W, k=1)])
  t1 = t0[np.where(t0<t_upper)]
  t2 = np.sort(t1)[::-1]
  temp = (p**2-p)/2
  
  thes = None
  use_t_upper = False
  x = np.zeros(t2.shape[0])
  for i in range(t2.shape[0]):
   x[i] = 2*(1-norm.cdf(t2[i]))*temp/(i+1)
   if x[i]>=alpha:  
     if i>0:
        thes = t2[i-1]
     if i==0:
        thes = t_upper; use_t_upper = True
     break 

  W_thes = abs(W)<thes
  res = []
  res.append(W_thes)
  res.append(thes)
  res.append(use_t_upper)
  
  return(res)


def get_scr_mat(expr1,expr2,s = None,s_seq=np.arange(0.2,2.2,0.2),alpha=0.4):
  if s == None:
    s_selected = s_selection(expr1,expr2,ss = s_seq)

  W = get_diff_W(expr1 = expr1,expr2 = expr2,s = s_selected)
  temp = get_W_theshold(W,alpha = alpha)
  W_diff_m = temp[0]
  return(W_diff_m)


 