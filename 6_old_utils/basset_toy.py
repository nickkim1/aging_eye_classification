import csv
import torch 
import torch.nn as nn
import torch.nn.functional as F
import random
from sklearn import metrics
from torch.utils.data import Dataset, DataLoader, BatchSampler
import torch.utils.data as torch_data
from matplotlib import pyplot as plt
import seaborn as sn
import pandas as pd 
import numpy as np
import math
import time as time
import os
import argparse as ap
import data_preprocessing as p

### NOTE ### 
# 1. this is not the full implementation of the Basset architecture
# 2. no need to preprocess because this is just bulk data, not single cell (can just use raw counts)

### COMMANDS FOR IMPLEMENTATION ### 

# first set the device used
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print("===================================")
print(f"0. Using {device} device!")

### toy dataset for use ### 
class ToyDataset(Dataset):
    def __init__(self, n_samples): 
        self.n_samples = n_samples
        # create a list of 100 sequences (samples) to play with 
        self.features = np.zeros([self.n_samples * 600, 4])
        # print(np.shape(self.features))
        # print('original: ', self.features[:,0:3])
        for i in range(self.n_samples):
            s = self.features[600*i:(600*i)+600,:]
            for j in range(600): # iterate over each row 
                s[j][np.random.randint(4)] = 1 # randomly set a base -> 1 (indicate it exists @ that position in the sequence)
            self.features[600*i:(600*i)+600,:] = s # set the sequence @ the range of a given gene -> the specified seequence 
        # print('exploration of what feats look like, ', self.features[0:50,:])
        self.features = torch.from_numpy(self.features).float()
        # print('shape of features vector ', self.features.shape)
        self.labels = np.random.negative_binomial(20, 0.75, (self.n_samples * 600, 1, 164)) 
        self.labels = torch.from_numpy(self.labels).float()
        # print('exploration of what labels look like, ', self.labels[0:50, :])
        # print('shape of labels vector, ', self.labels.shape)
        self.gene_ids = torch.from_numpy(np.arange(0, self.n_samples * 600, 1))
        # print('shape of gene_ids, ', self.gene_ids.shape)
        self.indices = {}
        self.classes = []
        for idx, gene_id in enumerate(self.gene_ids): # since we sample by idx
            if (idx+600) % 600 == 0:
                self.classes.append(int(gene_id.item()))
                self.indices[int(gene_id.item())] = np.arange(idx, idx+600, 1) # 4 bp
        # print('length of indices vector', len(self.indices))
        # print('length of classes ', len(self.classes))
    def __len__(self):
        return self.n_samples
    def __getitem__(self, idx):
        return self.features[idx], self.labels[idx]
    def gene_ids_and_indices(self):
        return self.classes, self.indices

### custom batch sampler class (for generating the 4 x 600 inputs) etc ###
class BSampler(BatchSampler): 
    def __init__(self, num_rows, gene_ids, indices, batch_size):
        super(BSampler, self).__init__(train_dataset, batch_size, drop_last=False) # i forget why dataset is needed here 
        self.gene_ids = gene_ids
        self.num_rows = num_rows
        self.indices = indices 
        self.n_batches = int(num_rows / batch_size) # batch size is 4; num_rows is 400 -> 100 batches (samples)
        self.batch_size = batch_size
        print('num batches: ', self.n_batches)
    def __iter__(self):
        batches = []
        for i in range(self.n_batches):
            batch = []
            batch_class = random.choice(self.gene_ids) # randomly choose an idx (sample)
            # print(i, batch_class)
            vals = torch.from_numpy(self.indices[batch_class]) # get all 600 rows for that idx 
            # print(i, vals.shape)
            for ignore, idx in enumerate(vals):
                batch.append(idx) # should be a list of indices (single element tensors)
            batches.append(batch) # append each batch to the batches list
        # print('exploration of what batch entries look like, ', batches[1][len(batches[1])-3:len(batches[1])-1])
        # print('length of batches list, ', len(batches))
        return iter(batches)
    def __len__(self):
        # this doesn't return anything informative unless i change the num_rows into constructor param
        return self.num_rows

#### create the datasets from toy dataset class ####
#### just worked with a 60:20:20 split under the assumption all the data would be loaded in separately #### 
#### ideally i would've been able to randomly split on the fake training data i generated, but i couldn't get it to work ### 

# train dataset 
train_dataset = ToyDataset(n_samples=60)
train_dataloader = DataLoader(dataset=train_dataset, batch_sampler=BSampler(num_rows=(len(train_dataset) * 600), gene_ids=train_dataset.gene_ids_and_indices()[0], indices=train_dataset.gene_ids_and_indices()[1], batch_size=600))
train_feat, train_label = next(iter(train_dataloader))
n_total_steps = len(train_dataloader) / 600

# validation dataset 
valid_dataset = ToyDataset(n_samples=20)
valid_dataloader = DataLoader(dataset=valid_dataset, batch_sampler=BSampler(num_rows=(len(valid_dataset) * 600), gene_ids=valid_dataset.gene_ids_and_indices()[0], indices=valid_dataset.gene_ids_and_indices()[1], batch_size=600))
valid_feat, valid_label = next(iter(valid_dataloader))

# test dataset 
test_dataset = ToyDataset(n_samples=20)
test_dataloader = DataLoader(dataset=test_dataset, batch_sampler=BSampler(num_rows=(len(test_dataset) * 600), gene_ids=test_dataset.gene_ids_and_indices()[0], indices=test_dataset.gene_ids_and_indices()[1], batch_size=600))
test_feat, test_label = next(iter(valid_dataloader))

### define the CNN architecture for basset ###
class BassetCNN(nn.Module):
    def __init__(self): 
        super(BassetCNN, self).__init__()
        self.conv_one = nn.Conv1d(in_channels=4, out_channels=300, kernel_size=19)
        self.batchnorm_one = nn.BatchNorm1d(582) # input = num channels 
        self.pool_one = nn.MaxPool1d(3) 
        self.conv_two = nn.Conv1d(in_channels=300, out_channels=200, kernel_size=11)
        self.batchnorm_two = nn.BatchNorm1d(184) # input = num channels 
        self.pool_two = nn.MaxPool1d(4)
        self.conv_three = nn.Conv1d(in_channels=200, out_channels=200, kernel_size=7)
        self.batchnorm_three = nn.BatchNorm1d(40) # input = num channels 
        self.pool_three = nn.MaxPool1d(4)
        self.fc1 = nn.Linear(2000, 1000) # 1000 unit linear layer (same as paper)
        self.dropout_one = nn.Dropout1d(0.3) # for some reason the 1D dropout doesn't work  
        self.fc2 = nn.Linear(1000, 1000) # unsure of specific rationale why they kept the layer the same size, guess it was just optimal? 
        self.dropout_two = nn.Dropout1d(0.3)
        self.fc3 = nn.Linear(1000, 164) # output dim should be [1x164] since i unsqueezed @ flattening above
    def forward(self, x): 
        x = self.conv_one(x)
        # print('first conv layer: ', x.shape)
        x = self.batchnorm_one(x)
        # print('first batchnorm layer: ', x.shape)
        x = F.relu(x)
        # print('first relu: ', x.shape)
        x = self.pool_one(x)
        # print('first pooling: ', x.shape)
        x = self.conv_two(x)
        # print('second conv layer: ', x.shape)
        x = self.batchnorm_two(x)
        # print('second batchnorm layer: ', x.shape)
        x = F.relu(x)
        # print('second relu: ', x.shape)
        x = self.pool_two(x)
        # print('second pooling: ', x.shape)
        x = self.conv_three(x)
        # print('third conv layer: ', x.shape)
        x = self.batchnorm_three(x)
        # print('third batchnorm layer: ', x.shape)
        x = F.relu(x)
        # print('third relu ', x.shape)
        x = self.pool_three(x)
        # print('third pooling: ', x.shape)
        x = x.view(2000) # view for the output layer -- should just be the product of the first and escond dimension
        # print('flattened layer: ', x.shape)
        x = self.fc1(x.unsqueeze(0)) # unsqueeze after flattening
        # print('first fc layer: ', x.shape)
        x = F.relu(x)
        # print('fourth relu ', x.shape)
        x = self.dropout_one(x)
        # print('first dropout ', x.shape)
        x = self.fc2(x)
        # print('second fc layer: ', x.shape)
        x = self.dropout_two(x)
        # print('second dropout layer: ', x.shape)
        x = self.fc3(x)
        # print('final fc layer: ', x.shape)
        x = F.sigmoid(x)
        # print('output sigmoid layer: ', x.shape)
        return x

### initialize and send the model to the device above ### 
model = BassetCNN().to(device)

### define hyperparameters as described in the paper ### 
### of note: they did bayesian search w/ Spearmint for theirs - tbd here ### 

learning_rate = 0.001 # my own arbitrary choice for the lr 
criterion = nn.BCELoss() # this is taken from the paper 
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)

### define some other parameters up here for use downstream in training/testing ### 
num_epochs = 50
n_batches = 6 # re-defined here purely for the purpose of training

### train the model ### 
### i didn't include any measures of accuracy of anything here -- just loss -- because i didn't know what the ground truth data should really look like? ### 
def train_model():
    ### initialize the weights. unsure if they did random initialization however ### 
    def init_weights(m):
        if isinstance(m, nn.Linear) or isinstance(m, nn.Conv1d):
            torch.nn.init.xavier_uniform_(m.weight)
            m.bias.data.fill_(0.01)
        model.apply(init_weights)

    ### initialize log dictionary for keeping track of training loss, accuracy, etc.. ### 
    train_log = {'training_loss_per_epoch':[]}

    ### for x number of epochs, run through all mini-batches
    for epoch in range(num_epochs):
        # switch the model to training mode - for incorporating stuff like dropout layers 
        model.train() 
        # init log of training loss 
        t_loss_per_batch = []
        ### optimize for each mini-batch ### 
        for i, (samples, labels) in enumerate(train_dataloader): 
            
            samples = samples.permute(1, 0) # tranpose the sample so that it is (4x600)
            # print('minibatch sample dims', samples.shape)
            samples = samples.to(device) # send the samples to the device
            # print(labels.type()) <- in case some inconsistency in typing from data in tensor format
            labels = labels[0] # assume that the first row is just the labels (i'm not sure how the real input data will be formatted, but this definitely isn't a great assumption)
            labels = labels.to(device) # send the labels to the device
            # print('minibatch label dims', labels.shape)
            
            model.to(device) # send the model to the device
            predicted = model(samples) # get the predicted output 
            # calculate the loss 
            loss = criterion(predicted, labels) 
            t_loss_per_batch.append(loss.item())# append the loss to the tracker of the loss per batch 

            optimizer.zero_grad() # zero accumulated gradients 
            loss.backward() # backprop
            optimizer.step() # step on params 

            # print message for every batch 
            if (i+1) % n_batches == 0: # -- this is where the i term above is used in for loop
                print (f'Epoch [{epoch+1}/{num_epochs}], Step [{i+1}/{int(n_total_steps)}], Batch Loss: {loss.item():.4f}')
        
        # insert the AVERAGE batch training loss -> log 
        t_loss = sum(t_loss_per_batch) / len(t_loss_per_batch)
        train_log['training_loss_per_epoch'].append(t_loss)

    # save the model's params after fully training it 
    PATH = './basset_params.pth'
    torch.save(model.state_dict(), PATH) # -- save the model's params 

    return train_log # return the log of loss values 

### test the model; load in saved params from specified PATH ### 
def test_model(PATH):
    # load in the state dict from the passed in path 
    model.load_state_dict(torch.load(PATH)) 
    model.eval() # switch to eval mode to switch off layers like dropout 

    ### initialize log dictionary for keeping track of training loss, accuracy, etc.. ### 
    test_log = {'testing_loss_per_epoch':[]}
    
    with torch.no_grad():
        for epoch in range(num_epochs):
            testing_loss_per_batch = []
            for i, (samples, labels) in enumerate(test_dataloader):
                # set testing samples 
                samples = samples.permute(1, 0)
                samples = samples.to(device)
                # set testing labels 
                labels = labels[0]
                labels = labels.to(device)
                # send the model to the device
                model.to(device)
                predicted = model(samples)
                # calculate the loss 
                loss = criterion(predicted, labels)
                testing_loss_per_batch.append(loss.item())
            # calculate the average test loss 
            test_loss = sum(testing_loss_per_batch) / len(testing_loss_per_batch)
            test_log['testing_loss_per_epoch'].append(test_loss)
    
    return test_log

# call the methods to train and test the model 
train_log = train_model()
test_log = test_model('./basset_params.pth')

### quick function to get a visual of the loss ###
### of note: strange results for testing loss -- constant -- but i still have to debug and see what the specific issue is ### 
def plot_loss(num_epochs, train_loss, test_loss):
    f, a = plt.subplots(figsize=(10,7.5), layout='constrained') # don't need to specify 2x2 or anything here, bc i'm just going to plot the loss 
    f.suptitle('Calculated Loss')
    a.plot(num_epochs, train_loss, label='Training Loss')
    a.plot(num_epochs, test_loss, label=f'Testing Loss')
    a.set_xlabel('Number of Epochs')
    a.set_ylabel('Average Loss')
    a.set_title(f'Training and Testing Loss')
    a.legend()
    plt.show() # just show it - i guess i'll save it later once i have a more functional model running 
 
# call the loss function plotter 
plot_loss(np.linspace(1, num_epochs, num=num_epochs).astype(int), train_log['training_loss_per_epoch'], test_log['testing_loss_per_epoch'])













##### OLD CODE ~ IGNORE ~ #####

# then do training, valid, and testing splits (80:20) on the dataset 
# train_size = int(0.8 * len(dataset))
# call on the radnom split function 
# train_dataset, valid_dataset = torch.utils.data.random_split(dataset, [train_size, len(dataset) - train_size])


# feat, label = next(iter(dataloader))

# i'm not too sure if this is right bc if randomly splitting dataset -> dk correspondence w idx vals? 
# print(len(train_dataset))
# print(len(valid_dataset))
# print('shape ', len(train_dataset))
# print('length', len(dataset.gene_ids_and_indices()[0]))
# print('length keys', len(dataset.gene_ids_and_indices()[1].keys()))
# train_dataloader = DataLoader(dataset=train_dataset, batch_sampler=BSampler(num_rows=(len(train_dataset) * 600), gene_ids=train_dataset.dataset.gene_ids_and_indices()[0], indices=train_dataset.dataset.gene_ids_and_indices()[1], batch_size=600))
# print('dataloader length', len(train_dataloader.dataset))
# valid_dataloader = DataLoader(dataset=valid_dataset, batch_sampler=BSampler(num_rows=(len(valid_dataset) * 600), indices=valid_dataset.dataset.gene_ids_and_indices()[1][0:20], batch_size=600))

# print(feat.shape, label.shape)
# print(label)

#---for one liner---# 
# sampler = {}
# sampler['train'] = train_data
# sampler['valid'] = valid_data
# sampler['test'] = test_data

# set a seed for reproducible results (splits) in this context 
# generator = torch.Generator().manual_seed(42)
# randomly split into training/val/test splits (only 100 samples so should be 60/20/20 split)
# train_data, valid_data, test_data = torch_data.random_split(dataset, [0.6, 0.2, 0.2], generator=generator)
# print(len(train_data), len(valid_data), len(test_data)) 

# num_seq = 100 
# create a list of 100 sequences to play with 
# features = np.zeros([num_seq * 4, 600])
# print(np.shape(features))
# print('original: ', features[0:3,:])
# for i in range(num_seq):
#     s = features[4*i:(4*i)+4,:]
#     for j in range(600):
#         s[np.random.randint(4)][j] = 1
#     features[4*i:(4*i)+4,:] = s
# print('new: ', features[0:3,:])
# labels = np.random.negative_binomial(20, 0.75, (num_seq, 600)) 
# print(np.shape(labels))

# classification problem -- create peak calls base by base 
# label = np.random.negative_binomial(20, 0.75, (1, 600)) 
# data[seq] = (feat, label)

# each input should be a 4x600 matrix and labels should be 1x600
# labels should be 100,600 cuz no need for basepairs