import pandas as pd
import sys
import os
from sklearn.preprocessing import StandardScaler

RCFGL_path = '/Users/seals/Documents/GitHub/RCFGL'
os.chdir(RCFGL_path)
sys.path.insert(0, 'Python_functions')
from RCFGL import RFGL, RCFGL
os.chdir(RCFGL_path)
from Dstream_functions import*


B = pd.read_csv('Data/Data_1.csv');
H = pd.read_csv('Data/Data_2.csv');
M = pd.read_csv('Data/Data_3.csv');

scalerB = StandardScaler(with_std=False); scalerB.fit(B)
B = scalerB.transform(B)

scalerH = StandardScaler(with_std=False);  scalerH.fit(H)
H = scalerH.transform(H)

scalerM = StandardScaler(with_std=False); scalerM.fit(M)
M = scalerM.transform(M)

A = []; A.append(B); A.append(H); A.append(M);


RFGL_result = RFGL(lambda1 = 0.1, lambda2 = 0.1, A = A, 
        ADMMmaxiter = 200)


RCFGL_result = RCFGL(lambda1 = 0.1, lambda2 = 0.1, A = A, 
        ADMMmaxiter = 200)


truncation_value = 0.05; top_N = 100;

    

theta1 = RCFGL_result[0][:,:,0]
names = np.array(range(0, theta1.shape[0]))
Adjacency1 = MakeAdjMatrix(theta1, truncation_value = truncation_value, top_N = top_N, names = names)


theta2 = RCFGL_result[0][:,:,1]
Adjacency2 = MakeAdjMatrix(theta2, truncation_value = truncation_value, top_N = top_N, names = names)


theta3 = RCFGL_result[0][:,:,2]
Adjacency3 = MakeAdjMatrix(theta3, truncation_value = truncation_value, top_N = top_N, names = names)





Adjacency_all = MakeAdjMatrix_all(RCFGL_result, truncation_value = 0.05,
top_N = 75)

PairNetworkPlotter(Adj_all, pair = [2, 3])
AllNetworkPlotter(Adjacency_all)
NetworkPlotter(Adj_all, which = 1)


Adj_all = MakeAdjMatrix_all(RCFGL_result, truncation_value = 0.1, top_N = 100)

theta = RCFGL_result[0][:,:,0]
SummaryTable(Adj_all)

