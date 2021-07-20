
#import igraph
#os.chdir("/Users/seals/Documents/Github/RCFGL/C_functions")
#cppyy.include('thm_screening.h')
#from cppyy.gbl import globalScreening
import igl
import numpy as np
from scipy.sparse import csr_matrix
def thm_screening(lambda1,S,weights):
    p = S.shape[1];
    K = S.shape[2]
    #adj = np.zeros((p,p));
    #globalScreening(S,adj,p,K,lambda1,lambda2)
    #g = igraph.Graph.Adjacency(adj.tolist())
    #g = igl.connected_components(csr_matrix(adj,dtype=np.int32)

    crit1 = np.zeros((p,p,K))
    #crit2 = np.zeros((p,p));
    for k in range(K):
        crit1[:,:,k] =  abs(S[:,:,k]) > lambda1/weights[k] 

    critall = np.sum(crit1, axis=2)
    critall = critall != 0
    g = igl.connected_components(csr_matrix(critall.astype(int),dtype=np.int32))
    return(g)
