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
import time as time
import os
import sys
import argparse as ap
import load_data as p
import data_aug as aug

# TODO 1: 
# * this is not the full implementation of the Basset architecture
# * no need to preprocess because this is just bulk data, not single cell (can just use raw counts)
# - incorporate in the control data and predict with that? should have two outputs in the model: one for control, one for experimental X 
# - see where the loss curves are messing up X 
# - make each figure take in a description field automatically filled in w/ model's relevant params at time of running <-- 

# TODO 2: 
# -- augment the data by reverse complementing the data, inverting it, shifting it 

# specify model type 
model_type = sys.argv[1]

# first set the device used
device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print(f"0. Using {device} device!")
print('===============================================================================')

# for the dataset 
class CustomDataset(Dataset): 
    def __init__(self, dataset, n_samples): 
        # xy = np.loadtxt(file_path, delimiter=",", dtype=np.float32) # data loaded in will be comma delimited
        self.n_samples = n_samples
        self.features = np.empty([self.n_samples, 600, 4]) # (13854x1) 
        self.labels = np.empty([self.n_samples, 1])
        for i in range(self.n_samples): # n_samples is the number of rows, 13854
            s = np.zeros([600,4]) # <- (13854x600)x4
            matching_row_sequence = dataset[i][1]
            matching_row_label = dataset[i][2]
            for j in range(600): 
                # one hot encode sequences as: A T C G (top -> bottom)
                idx_to_mark = 0
                if matching_row_sequence[j] == 'A': 
                    idx_to_mark = 0
                elif matching_row_sequence[j] == 'T':
                    idx_to_mark = 1
                elif matching_row_sequence[j] == 'C': 
                    idx_to_mark = 2
                elif matching_row_sequence[j] == 'G':
                    idx_to_mark = 3
                s[j][idx_to_mark] = 1
            self.features[i] = s # set the sequence @ the range of a given gene -> the specified seequence 
            self.labels[i] = matching_row_label
        self.features = torch.from_numpy(self.features).float() 
        self.gene_ids = np.arange(0, self.n_samples, 1)
        # self.gene_ids = []
        self.classes = list(self.gene_ids)
    def __len__(self):
        return self.n_samples
    def __getitem__(self, index):
        return self.features[index], self.labels[index]
    def gene_ids_and_indices(self):
        return self.classes, self.features

#  custom batch sampler class (for generating the 4 x 600 inputs) et
class BSampler(BatchSampler): 
    def __init__(self, num_rows, gene_ids, indices, batch_size):
        super(BSampler, self).__init__(train_dataset, batch_size, drop_last=False) # i forget why dataset is needed here 
        self.gene_ids = gene_ids
        self.num_rows = num_rows
        self.indices = indices 
    def __iter__(self):
        # np.random.seed(0)
        np.random.shuffle(self.gene_ids)
        batches = []
        for ignore in range(self.num_rows):
            # randomly choose an idx (sample), from 0 to 13584 (full dataset)
            batch = [random.choice(self.gene_ids)]
            batches.append(batch)
        return iter(batches)
    def __len__(self):
        # this doesn't return anything informative unless i change the num_rows into constructor param
        return self.num_rows

# train dataset 
training, testing, validation = p.training, p.testing, p.validation
train_dataset = CustomDataset(dataset=training, n_samples=np.shape(training)[0])
train_dataloader = DataLoader(dataset=train_dataset, batch_sampler=BSampler(num_rows=len(train_dataset), gene_ids=train_dataset.gene_ids_and_indices()[0], indices=train_dataset.gene_ids_and_indices()[1], batch_size=1))
train_feat, train_label = next(iter(train_dataloader))
n_total_steps = len(train_dataloader) 

# test dataset 
test_dataset = CustomDataset(dataset=testing, n_samples=np.shape(testing)[0])
test_dataloader = DataLoader(dataset=test_dataset, batch_sampler=BSampler(num_rows=len(test_dataset), gene_ids=test_dataset.gene_ids_and_indices()[0], indices=test_dataset.gene_ids_and_indices()[1], batch_size=1))
test_feat, test_label = next(iter(test_dataloader))

# define the CNN architecture for basset
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
        self.dropout_one = nn.Dropout1d(0.3)  
        self.fc2 = nn.Linear(1000, 164) # 
        self.dropout_two = nn.Dropout1d(0.3)
        self.fc3 = nn.Linear(164, 1) # output dim should be [1x164] since i unsqueezed @ flattening above
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
        x = x.view(2000) 
        x = self.fc1(x.unsqueeze(0)) # unsqueeze after flattening
        x = F.relu(x)
        x = self.dropout_one(x)
        x = self.fc2(x)
        x = self.dropout_two(x)
        x = self.fc3(x)
        x = F.sigmoid(torch.tensor(torch.flatten(x), dtype=torch.float32, requires_grad=True))
        return x

# sacrifices spatial relationships 
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

# initialize and send the model to the device above
if model_type == "cnn": 
    model = BassetCNN().to(device)
elif model_type == "nn":
    model = FeedForwardNN().to(device)

# init data augmentation class
data_aug = aug.DataAugmentation

# define hyperparameters as described in the paper
learning_rate = 0.001 
criterion = nn.BCELoss() 
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
# optimizer = torch.optim.RMSprop(model.parameters(), lr=learning_rate)
# define some other parameters up here for use downstream in training/testing
num_epochs = 20
n_batches = len(train_dataloader)

# train the model: 
def train_model():
    
    # initialize the weights. unsure if they did random initialization however 
    def init_weights(m):
        if isinstance(m, nn.Linear) or isinstance(m, nn.Conv1d):
            torch.nn.init.xavier_uniform_(m.weight)
            m.bias.data.fill_(0.01)
        model.apply(init_weights)
    
    # keep track of the parameters with this dictionary 
    train_log = {'training_loss_per_epoch':[], 'training_accuracy_per_epoch':[]} 

    for epoch in range(num_epochs):  

        # switch the model to training mode 
        model.train() 
        t_loss_per_batch = []
        t_accuracy_per_batch = []

        for i, (samples, labels) in enumerate(train_dataloader): 

            # tranpose the sample so that it is (4x600)
            samples = torch.reshape(samples,(samples.shape[0]*samples.shape[1],samples.shape[2])).permute(1, 0)
            # every other epoch, reverse complement the sample 
            # if epoch % 2 == 1: 
            #     # samples = data_aug.reverse_complement_sequence(data_aug, samples)
            #     samples = data_aug.reverse_complement_sequence(data_aug, samples)
            samples = samples.to(device) # send the samples to the device

            # assume the first row is the labels? 
            labels = labels[0].type(torch.float32)
            labels = labels.to(device) 
            model.to(device) 
            predicted = model(samples) 
            loss = criterion(predicted, labels) 
            t_loss_per_batch.append(loss.item()) 

            if (predicted.item() > 0.5 and labels.item() == 1) or (predicted.item() <= 0.5 and labels.item() == 0): 
                t_accuracy_per_batch.append(1)
            else: 
                t_accuracy_per_batch.append(0)

            # zero accumulated gradients, backprop, and step on params 
            optimizer.zero_grad()
            loss.backward() 
            optimizer.step() 

            # print message for every batch 
            if (i+1) % n_batches == 0: # -- this is where the i term above is used in for loop
                print (f'Epoch [{epoch+1}/{num_epochs}], Step [{i+1}/{int(n_total_steps)}], Loss: {loss.item():.4f}')
        
        # insert the AVERAGE batch training loss -> log. technically could use total_batches here but opted not to
        t_loss = sum(t_loss_per_batch) / len(t_loss_per_batch)
        train_log['training_loss_per_epoch'].append(t_loss)
        t_accuracy = sum(t_accuracy_per_batch) / len(t_accuracy_per_batch)
        print('Training Accuracy: ', t_accuracy)
        train_log['training_accuracy_per_epoch'].append(t_accuracy)

    # save the model's params after fully training it 
    PATH = '../5_saved_models/basset_params.pth'
    torch.save(model.state_dict(), PATH)

    # return the log of loss values 
    return train_log 

# test the model; load in saved params from specified PATH
def test_model(PATH):
    model.load_state_dict(torch.load(PATH)) 
    model.eval() # switch to eval mode to switch off layers like dropout

    test_log = {'testing_loss_per_epoch':[], 'testing_accuracy_per_epoch': []}
    
    with torch.no_grad():
        for epoch in range(num_epochs):
            testing_loss_per_batch = []
            testing_accuracy_per_batch = []
            for i, (samples, labels) in enumerate(test_dataloader):
                
                samples = torch.reshape(samples, (samples.shape[0]*samples.shape[1], samples.shape[2])).permute(1, 0)
                samples = samples.to(device)
                labels = labels[0].type(torch.float32)
                labels = labels.to(device)
                model.to(device)
                predicted = model(samples)
                loss = criterion(predicted, labels)
                testing_loss_per_batch.append(loss.item())

                if (predicted.item() > 0.5 and labels.item() == 1) or (predicted.item() <= 0.5 and labels.item() == 0): 
                    testing_accuracy_per_batch.append(1)
                else: 
                    testing_accuracy_per_batch.append(0)

            test_loss = sum(testing_loss_per_batch) / len(testing_loss_per_batch)
            test_log['testing_loss_per_epoch'].append(test_loss)
            test_accuracy = sum(testing_accuracy_per_batch) / len(testing_accuracy_per_batch)
            print('Test Accuracy', test_accuracy)
            test_log['testing_accuracy_per_epoch'].append(test_accuracy)

    return test_log

# call the methods to train and test the model 
train_log = train_model()
test_log = test_model('../5_saved_models/basset_params.pth')

# quick function to get a visual of the loss
def plot_loss(id, num_epochs, train_loss, test_loss):
    f, a = plt.subplots(figsize=(10,7.5), layout='constrained') # don't need to specify 2x2 or anything here, bc i'm just going to plot the loss 
    f.suptitle('Calculated Loss')
    a.plot(num_epochs, train_loss, label='Training Loss')
    a.plot(num_epochs, test_loss, label=f'Testing Loss')
    a.set_xlabel('Number of Epochs')
    a.set_ylabel('Average Loss')
    a.set_title(f'Training and Testing Loss')
    a.legend()
    # plt.show() 
    plt.savefig(f'../4_viz/loss_curves_{id}')
 
def plot_accuracy(id, num_epochs, train_accuracy, test_accuracy): 
    f, a = plt.subplots(figsize=(10,7.5), layout='constrained') # don't need to specify 2x2 or anything here, bc i'm just going to plot the loss 
    f.suptitle('Calculated Loss')
    a.plot(num_epochs, train_accuracy, label='Training Accuracy')
    a.plot(num_epochs, test_accuracy, label=f'Testing Accuracy')
    a.set_xlabel('Number of Epochs')
    a.set_ylabel('Accuracy')
    a.set_title(f'Training and Testing Accuracy')
    a.legend()
    # plt.show() 
    plt.savefig(f'../4_viz/accuracy_curves_{id}')


# call the loss function plotter 
plot_loss("03", np.linspace(1, num_epochs, num=num_epochs).astype(int), train_log['training_loss_per_epoch'], test_log['testing_loss_per_epoch'])
plot_accuracy("03", np.linspace(1, num_epochs, num=num_epochs).astype(int), train_log['training_accuracy_per_epoch'], test_log['testing_accuracy_per_epoch'])