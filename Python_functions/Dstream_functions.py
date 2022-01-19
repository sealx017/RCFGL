import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
from matplotlib.lines import Line2D
#from tabulate import tabulate
import networkx as nx

def CovtoCor(covariance):
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    return correlation

def MakeAdjMatrix(theta, truncation_value, top_N, names):
    theta = abs(CovtoCor(theta)); theta[theta < truncation_value] = 0
    if names == 'default':
       names = np.array(range(0, theta.shape[0]))
    arr = np.triu(theta, 1)
    if top_N == 'all':
       top_N = len(np.where(arr != 0)[0])
    
    idx = np.argpartition(arr, arr.size - top_N, axis=None)[-top_N:]
    locations = np.column_stack(np.unravel_index(idx, arr.shape))
    adj = np.zeros(arr.shape)
    for k in locations:
     adj[k[0], k[1]] = adj[k[1], k[0]] = 1
    
    nonz = np.where(np.sum(adj, axis = 0) != 0)
    nonz_names = names[nonz[0]]
    adj = adj[np.ix_(nonz[0],nonz[0])]
    Adjacency = []
    Adjacency.append(adj)
    Adjacency.append(nonz_names)
    return Adjacency

def MakeAdjMatrix_all(result, truncation_value = 0.05, top_N = 'all', names = 'default'):
    Adjacency_all = []
    K = (result[0]).shape[2]
    for k in range(K):
     theta = result[0][:,:,k]
     Adjacency1 = MakeAdjMatrix(theta = theta, 
     truncation_value = truncation_value, top_N = top_N, names = names)
     A21 = pd.DataFrame(Adjacency1[0], index=Adjacency1[1], 
     columns=Adjacency1[1])
     Adjacency_all.append(A21)
    return Adjacency_all



def NetworkPlotter(Adjacency_all, which = 1):
 A2 = Adjacency_all[which-1]
 G = nx.from_pandas_adjacency(A2)
 plt.rcParams['figure.dpi'] = 500
 pos = nx.circular_layout(G)
 options = {
    "node_color": "#A0CBE2",
    "edge_color": "lightblue",
    "width": 1,
    "with_labels": True,
    "font_family":'sans-serif',
    "font_size" : 5,
    "node_size" : 0.5,
 }
 nx.draw_networkx(G, pos, **options)
 plt.show()

def PairNetworkPlotter(Adjacency_all, pair = [1, 2]):
    A21= Adjacency_all[pair[0]-1]; A22 = Adjacency_all[pair[1]-1]
    all_ind = np.union1d(A21.columns, A22.columns)
    A21_Full = pd.DataFrame(np.zeros([len(all_ind), len(all_ind)]), 
                            index=all_ind, columns=all_ind)
    A22_Full = pd.DataFrame(np.zeros([len(all_ind), len(all_ind)]), 
                            index=all_ind, columns=all_ind)
    A21_locs = [A21_Full.columns.get_loc(col) for col in A21.columns]
    A22_locs = [A22_Full.columns.get_loc(col) for col in A22.columns]
    
    A21_Full.iloc[A21_locs, A21_locs]  = A21
    A22_Full.iloc[A22_locs, A22_locs]  = A22
    A_fin = A21_Full + A22_Full
    A_diff = A21_Full - A22_Full; 
    A_diff = A_diff.replace(1, 0); A_diff = A_diff.replace(-1, -2)
    A2 = A_fin + A_diff; 
    plt.rcParams['figure.dpi'] = 500
    G = nx.from_pandas_adjacency(A2)
    #pos = nx.drawing.nx_agraph.graphviz_layout(G, prog='neato')
    #pos = nx.spring_layout(G, k = 0.1, iterations = 10)
    pos = nx.circular_layout(G)
    
    #cmap = cm.get_cmap('Set3') 
    unique_clrs = ['lightblue', 'lightgreen', 'orange']
    edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())
    
    clrs = []
    weights_low = []
    for i in range(len(weights)):
     if weights[i] == 2:
         clrs.append(unique_clrs[2])
         weights_low.append(2)
     elif weights[i] == 1:
         clrs.append(unique_clrs[1])
         weights_low.append(0.75)
     else:
         clrs.append(unique_clrs[0])
         weights_low.append(0.75)
    
    plt.ion()
    
    plt.figure(figsize = (10, 7), num=1); plt.clf()
    # draw the graph in several steps
    h1 = nx.draw_networkx_nodes(G, pos=pos, node_color = "#A0CBE2",
    alpha = 0.5, node_size = 1, linewidths=1)
    # we need the LineCollection of the edges to produce the legend (h2)
    h2 = nx.draw_networkx_edges(G, pos=pos, width=weights_low, edge_color=clrs)
    
    # and just show the node labels to check the labels are right!
    h3 = nx.draw_networkx_labels(G, pos=pos, font_size=5, font_color='black')
    
        
    def make_proxy(clr, mappable, **kwargs):
        return Line2D([0, 1], [0, 1], color=clr, **kwargs)
    
    # generate proxies with the above function
    labels = ["Only in cond."+ str(pair[0]), "Only in cond."+ str(pair[1]), "Common in both"]
    proxies = [make_proxy(clr, h2, lw=2) for clr in unique_clrs]
    if len(set(clrs)) == 1:
        x = unique_clrs.index(list(set(clrs))[0])
        labels = [labels[x]]
        proxies = [proxies[x]]
    elif len(set(clrs)) == 2:
        x = unique_clrs.index(list(set(clrs))[0])
        y = unique_clrs.index(list(set(clrs))[1])
        labels = [labels[x], labels[y]]
        proxies = [proxies[x], proxies[y]]
    
    legend = plt.legend(proxies, labels, title="Edge type",
    loc=1, fontsize='small', fancybox=True)
    frame = legend.get_frame() #sets up for color, edge, and transparency
    frame.set_facecolor('gainsboro') #color of legend
    frame.set_edgecolor('black') #edge color of legend
    frame.set_alpha(1) #deals with transparency
    
    plt.show()



def AllNetworkPlotter(Adjacency_all): 
    all_ind = np.zeros(1)
    num_condn = len(Adjacency_all)
    for k in range(num_condn):
     A21 = Adjacency_all[k]
     all_ind = np.union1d(all_ind, A21.columns)
     if k == 0: 
         int_ind = A21.columns
     else: 
         int_id = np.intersect1d(int_ind, A21.columns)
    
    Adj_all = pd.DataFrame(np.zeros((len(int_id), len(int_id))), 
                           index=int_id, columns=int_id)
    for k in range(len(Adjacency_all)):
     A2 = Adjacency_all[k]
     A2_locs = np.where(np.isin(A2.columns, int_id))[0].tolist()
     Adj_all = Adj_all + A2.iloc[A2_locs, A2_locs]
    
    Adj_all = Adj_all.replace(num_condn-1, 0); Adj_all = Adj_all.replace(num_condn-2, 0)
    Adj_all = Adj_all.loc[(Adj_all.sum(axis=1) != 0), (Adj_all.sum(axis=0) != 0)]
    plt.rcParams['figure.dpi'] = 500
    G = nx.from_pandas_adjacency(Adj_all)
    #pos = nx.spring_layout(G, k=0.3, iterations=50, seed = 63)
    pos = nx.circular_layout(G)
    
    options = {
       "node_color": "#A0CBE2",
       "edge_color": "orange",
       "width": 2,
       "with_labels": True,
       "font_size" : 5,
       "node_size" : 0.5,
    }
    nx.draw_networkx(G, pos, **options)    
    plt.show()


def SummaryTable(Adjacency_all): 
    K = len(Adjacency_all)
    S = list()
    for k1 in range(K):
     A21 = Adjacency_all[k1]; 
     for k2 in range(K):
      A22 = Adjacency_all[k2]
      all_ind = np.union1d(A21.columns, A22.columns)
      A21_Full = pd.DataFrame(np.zeros([len(all_ind), len(all_ind)]), 
                            index=all_ind, columns=all_ind)
      A22_Full = pd.DataFrame(np.zeros([len(all_ind), len(all_ind)]), 
                            index=all_ind, columns=all_ind)
      A21_locs = [A21_Full.columns.get_loc(col) for col in A21.columns]
      A22_locs = [A22_Full.columns.get_loc(col) for col in A22.columns]
      A21_Full.iloc[A21_locs, A21_locs]  = A21
      A22_Full.iloc[A22_locs, A22_locs]  = A22
      A_fin = A21_Full + A22_Full
      A_diff = A22_Full - A21_Full
      A_fin = np.triu(A_fin.values, 1)
      A_diff = np.triu(A_diff.values, 1)
      both_loc = np.where(A_fin == 2)
      all_loc = np.where(A_diff == -1)
      if k1 != k2:
       S.append({'Pair': (k1+1, k2+1), 
       'Prop.': round(len(both_loc[0])/(len(all_loc[0])+len(both_loc[0])),3),
       '# edges': (len(all_loc[0])+len(both_loc[0]))})
    return(S)
    