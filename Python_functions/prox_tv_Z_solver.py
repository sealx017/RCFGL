#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 17:57:32 2021

@author: seals
"""
import numpy as np
import prox_tv as ptv
import pywt

def ptv_z(x,Weightij,lambda1, lambda2, K):

 ws = np.zeros(K-1)
 for k in range(K-1):
  if Weightij[k] == 0:
      ws[k] = 0
  else:
      ws[k] = lambda2
 solution = ptv.tv1w_1d(x, ws);
 return pywt.threshold(solution, lambda1, 'soft')
 
 

