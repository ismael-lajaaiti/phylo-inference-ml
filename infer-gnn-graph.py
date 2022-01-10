# -----------------------------------------------------------
# INFERING MACROEVOLUTIONARY RATES WITH GRAPH NEURAL NETWORKS
# -----------------------------------------------------------


# Importing libraries 

import torch
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, ChebConv, SAGEConv, global_mean_pool 
import rpy2.robjects as robjects # load R object 
from rpy2.robjects import pandas2ri # load R object 
from tqdm import tqdm # print progress bar 
import pickle # save object 
import matplotlib.pyplot as plt # plot
import numpy as np
import random as rd 


# Global parameters
load_data = True # if data is already saved, don't compute just load it
device = "cuda:2" # which GPU to use 
batch_size_max = 64 # max. number of trees per batch 
n_train = 9000 # size of training set 
n_valid = 500  # size of validation set 
n_test  = 500  # size of test set 


# Loading trees and their corresponding parameters

pandas2ri.activate()
fname_graph = "trees-dataset/ntrees-10000-ntaxa-100-1000-lambda-0.1-1-epsilon-0-0.9-sscheck-TRUE-df.rds"
fname_param = "trees-dataset/ntrees-10000-ntaxa-100-1000-lambda-0.1-1-epsilon-0-0.9-sscheck-TRUE-param.rds"
readRDS = robjects.r['readRDS']
df_graph = readRDS(fname_graph)
df_graph = pandas2ri.rpy2py(df_graph) # data.frame containing tree information
df_param = readRDS(fname_param)
df_param = pandas2ri.rpy2py(df_param) # data.frame containing target parameters

n_param = len(df_param) # number of parameters to guess for each tree 
n_trees = len(df_graph) # total number of trees of the dataset 


# Format data 

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


fname = fname_graph[:-6] + "geomtensor" + ".obj" # file name 

if (not load_data):

    data_list  = []
    print("Formating data...")
    for n in tqdm(range(n_trees)):
        df_node, df_edge = df_graph[n][0], df_graph[n][1]
        params = [df_param[i][n] for i in range(n_param)]
        data = convert_df_to_tensor(df_node, df_edge, params)
        data_list.append(data)
    print("Formating data... Done.")
    
    file = open(fname, "wb") # file handler 
    pickle.dump(data_list, file) # save data_list
    print("Formated data saved.")

else:

    file = open(fname, "rb")
    data_list = pickle.load(file)
    print("Formated data loaded.")


# Creating train, valid and test set 

# Choosing the tree indices for training, validation and test randomly 
ind = np.arange(0, n_trees) 
rd.shuffle(ind) 
train_ind = ind[0:n_train]  
valid_ind = ind[n_train:n_train + n_valid]  
test_ind  = ind[n_train + n_valid:] 

# Splitting the dataset between training, validation and test. 
train_data = [data_list[i].to(device=device) for i in train_ind]
valid_data = [data_list[i].to(device=device) for i in valid_ind]
test_data  = [data_list[i].to(device=device) for i in test_ind]

# Converting the list to DataLoader
train_dl = DataLoader(train_data, batch_size = batch_size_max, shuffle = True)
valid_dl = DataLoader(valid_data, batch_size = batch_size_max, shuffle = True)
test_dl  = DataLoader(test_data , batch_size = 1)


# Creating the GNN architecture

class GCN(torch.nn.Module):
    def __init__(self, n_in, n_hidden, n_out):
        super().__init__()
        self.conv1 = SAGEConv(n_in, n_hidden)
        self.conv2 = SAGEConv(n_hidden, n_hidden)
        self.conv3 = SAGEConv(n_hidden, 2*n_hidden)
        self.lin1  = torch.nn.Linear(2*n_hidden, n_hidden)
        self.lin2  = torch.nn.Linear(n_hidden, n_out)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p = 0.001, training=self.training)
        x = self.conv2(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p = 0.001, training=self.training)
        x = self.conv3(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p = 0.001, training=self.training)
        x = global_mean_pool(x, batch)
        x = self.lin1(x)
        x = self.lin2(x)
        return x


# Defining training and validation loop 

def train(model, batch):
    optimizer.zero_grad()
    out = model(batch)
    batch_size = int(max(data.batch) + 1) # number of trees in the batch 
    loss = F.mse_loss(out, data.y.reshape([batch_size, n_out])) # compute loss 
    loss.backward() # backward propagation 
    optimizer.step()
    return(loss)

def valid(model, batch):
    out = model(batch)
    batch_size = int(max(data.batch) + 1) # number of trees in the batch 
    loss = F.mse_loss(out, data.y.reshape([batch_size, n_out])) # compute loss
    return(loss)


# Setting up the training 
n_in = data_list[0].num_node_features
n_out = len(data_list[0].y)
n_hidden = 100
n_epochs = 10
model = GCN(n_in, n_hidden, n_out).to(device=device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)

# Training loop 
for epoch in range(n_epochs):

    # Training 
    model.train()
    train_loss = []
    for data in tqdm(train_dl):
        loss = train(model, data) # train model and get loss
        loss = float(loss.to(device = "cpu"))
        train_loss.append(loss)
    mean_loss = np.mean(train_loss)
    print("Epoch %d - Train Loss %.4f" % (epoch, float(mean_loss))) # print progression 

    # Validation 
    model.eval()
    valid_loss = []
    for data in tqdm(valid_dl):
        loss = valid(model, data) # train model and get loss
        loss = float(loss.to(device = "cpu"))
        valid_loss.append(loss)
    mean_loss = np.mean(valid_loss)
    print("Epoch %d - Valid Loss %.4f" % (epoch, float(mean_loss))) # print progression 