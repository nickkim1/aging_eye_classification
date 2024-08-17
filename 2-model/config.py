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

# TODO (1): make this more extensible - don't want to have to mannually change what script I'm importing from to specify type of task
# TODO (2): convert this to some sort of JSON format or something idk 

from models_regression import SimpleCNN, FeedForwardNN, BassetCNN

parser = ap.ArgumentParser(prog='This program runs the models')
parser.add_argument('-primarymodel', type=str, help='Input model type')
parser.add_argument('-trainfile', type=str, help='Input filepath to the training data', nargs=1)
parser.add_argument('-testfile', type=str, help='Input filepath to the testing data', nargs=1)
parser.add_argument('-seqlen', type=int, help='Input length of the sequences in your data')
parser.add_argument('-optim', type=str, help='Input optimizer type')
parser.add_argument('-loss', type=str, help='Input loss type')
parser.add_argument('-batchsize',type=int, help='Input batch size')
args = vars(parser.parse_args())

# This is the female config dictionary
simple_female_config = {
    'model_architecture': 'conv1, mp1', 
    'conv1': {
        'Cin': 4,
        'Lin': 100, 
        'Cout': 500, 
        'kernel_size': 20, 
        'stride': 1, 
        'dilation': 1, 
        'padding': 0
    }, 
    'mp1': {
        'kernel_size': 19, 
        'stride': 19, 
        'dilation': 1, 
        'padding': 0
    },
    'dropout': {
        'probability': 0.4
    }, 
    'fc1': {
        'out': 1500
    }, 
    'fc2': {
        'out': 300
    }, 
    'fc3': {
        'out': 8
    }
}

# This is the male config dictionary
simple_male_config = {
    'model_architecture': 'conv1, mp1', 
    'conv1': {
        'Cin': 4,
        'Lin': 100, 
        'Cout': 500, 
        'kernel_size': 20, 
        'stride': 1, 
        'dilation': 1, 
        'padding': 0
    }, 
    'mp1': {
        'kernel_size': 19, 
        'stride': 19, 
        'dilation': 1, 
        'padding': 0
    },
    'dropout': {
        'probability': 0.4
    }, 
    'fc1': {
        'out': 1500
    }, 
    'fc2': {
        'out': 300
    }, 
    'fc3': {
        'out': 8
    }
}

# This is the overall male configuration 
model_config = {
    'DeepChrome_F': SimpleCNN(simple_female_config),  
    'DeepChrome_M': SimpleCNN(simple_male_config),  
    'MLP': FeedForwardNN(),
    'Basset': BassetCNN()
}

# This is a function that populates your settings for you 
def populate_settings(train_filepath, test_filepath, seqlen, loss, optim, model_type, batch_size): 
    data_config = {
        'train': train_filepath, 
        'test': test_filepath, 
        'seqlen': seqlen
    }
    hyperparameter_config = {
        'learning_rate': 0.001, 
        'num_epochs': 50, 
        'batch_size': batch_size
    }
    optim_config = {
        'Adam': torch.optim.Adam(model_config[model_type].parameters(), lr=hyperparameter_config['learning_rate']), 
        'RMSprop': torch.optim.RMSprop(model_config[model_type].parameters(), lr=hyperparameter_config['learning_rate']), 
        'SGD': torch.optim.SGD(model_config[model_type].parameters(), lr=hyperparameter_config['learning_rate']) 
    }
    loss_config = {
        'BCE': nn.BCELoss(), 
        'MSE': nn.MSELoss()
    }
    settings = {
        'data': data_config,
        'optim_name': optim,
        'model_name': model_type,
        'model': model_config[model_type], 
        'hyperparameters': hyperparameter_config, 
        'optimizer': optim_config[optim], 
        'loss': loss_config[loss], 
    }
    return settings

