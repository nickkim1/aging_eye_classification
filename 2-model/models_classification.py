import random
import torch 
import torch.nn as nn
import torch.nn.functional as F
import collections
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

# TODOs in this section: 
    # 1. Read through some of the papers that Ananya sent, experiment w their architectures 
        # - atp have given up on original results
    # 2. Finish up using the Basset architecture 
    # 3. Try a regression architecture instead of binary classification <- priority
    # 4. Check up on ATAC-seq preprocessing and do a merged model 
    # 5. Hyperparameter tuning (just do random and/or grid search for now)


# sacrifices spatial relationships, but probably a decent baseline to use  
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
        x = F.sigmoid(torch.tensor(torch.flatten(x), dtype=torch.float32, requires_grad=True))
        return x

# testing cnn, based on very simple architecture (like deepchrome)
# conv -> maxpool -> dropout -> linear -> linear -> linear layers
# note i'll treat each nucleotide as a feature for the 1d CNN 

class SimpleCNN(nn.Module):
    def __init__(self, simple_config): 
        super(SimpleCNN, self).__init__()
        model_architecture = simple_config['model_architecture']
        self.conv_one = nn.Conv1d(in_channels=simple_config['conv1']['Cin'], out_channels=simple_config['conv1']['Cout'], kernel_size=simple_config['conv1']['kernel_size']) # output dim: (1, 300, 582)
        self.max_pool = nn.MaxPool1d(kernel_size=simple_config['mp1']['kernel_size']) 
        self.dropout = nn.Dropout1d(simple_config['dropout']['probability'])
        self.linear_1 = nn.Linear(in_features=self.calculate_flattened_size(simple_config, model_architecture), out_features=simple_config['fc1']['out'])
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
        x = torch.sigmoid(x)
        x = torch.squeeze(x) 
        return x
    
    def calculate_flattened_size(self, simple_config:dict, model_architecture): 
        list_of_layers = [a.strip() for a in model_architecture.split(",")]
        L_in = 0
        # this only works since the essential formula for mp and also conv are like the same 
        layer_inout_size_queue = collections.deque([])
        for layer_number, layer in enumerate(list_of_layers): 
            if simple_config[layer]: 
                # if you hit this condition you know that you're at the first convolutional layer 
                if layer_number == 0 and 'Lin' in simple_config[layer]: 
                    L_in = simple_config[layer]['Lin']
                else:
                    L_in = layer_inout_size_queue.popleft()
                # you should've hit a convolutional layer before this so pop left
                L_out_new = self.calculate_cnn_mp_output_size(L_in, simple_config, layer)
                layer_inout_size_queue.append(L_out_new)
        # for the purposes of this prj ig this -2 assumption is safe. Cout should be prev convolutional layer's output size 
        final_input = simple_config[list_of_layers[len(list_of_layers)-2]]['Cout']
        final_output = layer_inout_size_queue.popleft()
        return final_input * final_output
    
    def calculate_cnn_mp_output_size(self, L_in, simple_config, layer):
        padding, dilation, kernel_size, stride = simple_config[layer]['padding'], simple_config[layer]['dilation'], simple_config[layer]['kernel_size'], simple_config[layer]['stride']
        L_out_numerator = L_in + 2 * padding - dilation * (kernel_size - 1) - 1
        L_out_denominator = stride
        return math.floor(L_out_numerator / L_out_denominator) + 1 
    
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
        x = F.sigmoid(x)
        x = torch.squeeze(x) # flatten since it is [64,1] by the end but just want [64]
        return x
