import torch
from torch_geometric.data import Data
from torch_geometric.data import InMemoryDataset
from torch_geometric.loader import DataLoader
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, GCNConv, SAGEConv # Graph Neural Network 
from torch_geometric.nn import global_mean_pool 
import rpy2.robjects as robjects # load R object 
from rpy2.robjects import pandas2ri # load R object 
from tqdm import tqdm # print progress bar 
import pickle # save object 
import matplotlib.pyplot as plt
import numpy as np
import random as rd 



pandas2ri.activate()

fname_graph = "trees-dataset/ntrees-10009-ntaxa-100-1000-lambda-0.1-1-q-0.01-0.1-sscheck-TRUE-df.rds"
fname_param = "trees-dataset/ntrees-10009-ntaxa-100-1000-lambda-0.1-1-q-0.01-0.1-sscheck-TRUE-param.rds"

readRDS = robjects.r['readRDS']
df_graph = readRDS(fname_graph)
df_graph = pandas2ri.rpy2py(df_graph)
df_param = readRDS(fname_param)
df_param = pandas2ri.rpy2py(df_param)


n_param = len(df_param)
n_trees = len(df_graph)
print(n_trees, n_param)

def convert_df_to_tensor(df_node, df_edge, params):

    """
    Convert the data frames containing node and edge information 
    to a torch tensor that can be used to feed neural 
    """

    n_node, n_edge = df_node.shape[0], df_edge.shape[0]

    l1, l2 = [], []
    
    for i in range(n_edge):
        edge = df_edge.iloc[i]
        u, v = edge[0]-1, edge[1]-1
        l1 = l1 + [u,v]
        l2 = l2 + [v,u]

    edge_index = torch.tensor([l1,l2], dtype=torch.long)

    x = []

    for i in range(n_node):
        node_attr = list(df_node.iloc[i])
        x.append(node_attr)

    x = torch.tensor(x, dtype = torch.float)

    y = torch.tensor(params, dtype = torch.float)

    data = Data(x = x, edge_index = edge_index, y = y)

    return(data)


batch_size_max = 64
data_list  = []

for n in tqdm(range(n_trees)):
    df_node, df_edge = df_graph[n][0], df_graph[n][1]
    params = [df_param[i][n] for i in range(n_param)]
    data = convert_df_to_tensor(df_node, df_edge, params)
    data_list.append(data)


fname = fname_graph[:-6] + "geomtensor" + ".obj" # file name 
print("Save")
file = open(fname, "wb") # file handler 
pickle.dump(data_list, file) # save data_list
