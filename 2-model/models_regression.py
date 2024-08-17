import random
import torch 
import torch.nn as nn
import torch.nn.functional as F
from sklearn import metrics
from torch.utils.data import Dataset, DataLoader, BatchSampler
import torch.utils.data as torch_data
from matplotlib import pyplot as plt
import seaborn as sn
import pandas as pd 
import numpy as np
import time as time
import argparse as ap
import data_aug as aug
import math

# Sacrifices spatial relationships, but probably a decent baseline to use  
class FeedForwardNN(nn.Module):
    def __init__(self): 
        super(FeedForwardNN, self).__init__()
        self.layer1 = nn.Linear(in_features=2400, out_features=1000, bias=True)
        self.layer2 = nn.Linear(in_features=1000, out_features=10, bias=True)
        self.layer3 = nn.Linear(in_features=10, out_features=1)
    def forward(self, x): 
        x = x.reshape((1,2400))
        x = self.layer1(x)
        x = F.relu(x)
        x = self.layer2(x)
        x = F.relu(x)
        x = self.layer3(x)
        return x

# testing cnn, based on very simple architecture (like deepchrome)
# conv -> maxpool -> dropout -> linear -> linear -> linear layers
# note i'll treat each nucleotide as a feature for the 1d CNN 
class SimpleCNN(nn.Module):
    def __init__(self, simple_config): 
        super(SimpleCNN, self).__init__()
        self.conv_one = nn.Conv1d(in_channels=simple_config['conv1']['Cin'], out_channels=simple_config['conv1']['Cout'], kernel_size=simple_config['conv1']['kernel_size']) # output dim: (1, 300, 582)
        self.max_pool = nn.MaxPool1d(kernel_size=simple_config['mp1']['kernel_size']) 
        self.dropout = nn.Dropout1d(simple_config['dropout']['probability'])
        self.linear_1 = nn.Linear(in_features=15000, out_features=simple_config['fc1']['out'])
        self.linear_2 = nn.Linear(in_features=simple_config['fc1']['out'], out_features=simple_config['fc2']['out'])
        self.linear_3 = nn.Linear(in_features=simple_config['fc2']['out'], out_features=simple_config['fc3']['out'])
    def forward(self, x): 
        x = self.conv_one(x)
        x = F.relu(x)
        x = self.max_pool(x)
        x = self.dropout(x)
        x = x.view(x.shape[0], x.shape[1] * x.shape[2])
        x = self.linear_1(x)
        x = F.relu(x)
        x = self.linear_2(x)
        x = F.relu(x)
        x = self.linear_3(x)
        return x
    
# this CNN is derived from the publicized architecture of the Basset paper. Purely a re-implementation 
# just to see how it performs on this type of data. 
class BassetCNN(nn.Module): 
    def __init__(self): 
        super(BassetCNN, self).__init__()
        self.conv_one = nn.Conv1d(in_channels=4, out_channels=300, kernel_size=19)
        self.batchnorm_one = nn.BatchNorm1d(num_features=300) # input = num channels 
        self.pool_one = nn.MaxPool1d(kernel_size=3) 
        self.conv_two = nn.Conv1d(in_channels=300, out_channels=200, kernel_size=11)
        self.batchnorm_two = nn.BatchNorm1d(num_features=200) 
        self.pool_two = nn.MaxPool1d(kernel_size=4)
        self.conv_three = nn.Conv1d(in_channels=200, out_channels=200, kernel_size=7)
        self.batchnorm_three = nn.BatchNorm1d(num_features=200) 
        self.pool_three = nn.MaxPool1d(kernel_size=4)
        self.fc1 = nn.Linear(2000, 1000) 
        self.dropout_one = nn.Dropout1d(0.3) 
        self.fc2 = nn.Linear(1000, 1000) # unsure of specific rationale why they kept the layer the same size, guess it was just optimal? 
        self.dropout_two = nn.Dropout1d(0.3)
        self.fc3 = nn.Linear(1000, 1) # output dim should be [1x164] over all cells BUT since i'm only working with 1, should just be 1
    def forward(self, x): 
        x = self.conv_one(x)
        x = self.batchnorm_one(x)
        x = F.relu(x)
        x = self.pool_one(x)
        x = self.conv_two(x)
        x = self.batchnorm_two(x)
        x = F.relu(x)
        x = self.pool_two(x)
        x = self.conv_three(x)
        x = self.batchnorm_three(x)
        x = F.relu(x)
        x = self.pool_three(x)
        x = x.reshape(x.shape[0], -1)
        x = self.fc1(x)
        x = F.relu(x)
        x = self.dropout_one(x)
        x = self.fc2(x)
        x = F.relu(x)
        x = self.dropout_two(x)
        x = self.fc3(x)
        x = torch.squeeze(x) # flatten since it is [64,1] by the end but just want [64]
        return x
