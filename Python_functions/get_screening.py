import numpy as np
import cppyy
import os
os.chdir("/Users/seals/Documents/Github/RCFGL/C_functions")
cppyy.include('for_loop.h')
from cppyy.gbl import screening_loop
os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
from rpy2.robjects.packages import importr
import rpy2 


base = importr('base');
glmnet = importr('glmnet');
rpy2.robjects.r('''
        # create glmnet function `f`
        f <- function(x1,y1,lam1) {

        temp1 <- glmnet(x1,y1,family = "gaussian", alpha = 1, lambda = lam1)
        return(as.matrix(coefficients(temp1)[-1]))
    
        }
        
        ''')
        
r_f = rpy2.robjects.r['f']   


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
    x1_r = rpy2.robjects.r.matrix(rpy2.robjects.FloatVector(np.ascontiguousarray(x1.reshape(-1))), nrow = n1, ncol = (p-1), byrow = True)
    y1_r = rpy2.robjects.r.matrix(rpy2.robjects.FloatVector(np.ascontiguousarray(y1)), nrow = n1, ncol = 1, byrow = True)
    lam1 = rpy2.robjects.FloatVector(np.ascontiguousarray(s*np.sqrt(sigma1[i]*np.log(p)/n1)))
    reg = r_f(x1_r,y1_r,lam1) 
    b1[:,i] = np.array(reg).reshape(-1)


    y2 = expr2[:,i]
    x2 = np.delete(expr2,i,axis = 1)
    x2_r = rpy2.robjects.r.matrix(rpy2.robjects.FloatVector(np.ascontiguousarray(x2.reshape(-1))), nrow = n2, ncol = (p-1), byrow = True)
    y2_r = rpy2.robjects.r.matrix(rpy2.robjects.FloatVector(np.ascontiguousarray(y2)), nrow = n2, ncol = 1, byrow = True)
    lam1 = rpy2.robjects.FloatVector(np.ascontiguousarray(s*np.sqrt(sigma2[i]*np.log(p)/n2)))
    reg = r_f(x2_r,y2_r,lam1) 
    b2[:,i] = np.array(reg).reshape(-1)
   
    c1[:,i] = y1-np.mean(y1)-np.dot((x1-np.mean(x1,axis = 0)),b1[:,i])
    c2[:,i] = y2-np.mean(y2)-np.dot((x2-np.mean(x2,axis = 0)),b2[:,i])

  r1 = np.dot(c1.T,c1)/n1
  r2 = np.dot(c2.T,c2)/n2
  s1 = np.mean(np.power(c1,2), axis = 0).reshape((p,1))
  s2 = np.mean(np.power(c2,2), axis = 0).reshape((p,1))
  screening_loop(W, r1, s1, b1, r2, s2, b2, p, n1, n2)
  return(W)
