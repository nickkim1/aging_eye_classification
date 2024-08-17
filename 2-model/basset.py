import random
import torch 
import torch.nn as nn
import torch.nn.functional as F
import math
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

# TODO 1: (x meeting)
# * this is not the full implementation of the Basset architecture
# * no need to preprocess because this is just bulk data, not single cell (can just use raw counts)
# * idk what more to do besides just using a different architecture 
# - control data vs. exp data - two separate models 
# - make each figure take in a description field automatically filled in w/ model's relevant params at time of running <-- 

# TODO 2: (x meeting)
# augment the data by reverse complementing the data, inverting it, shifting it 

# TODO 3: set up a config dict for downstream use 

# TODO 4 (Jan 29 meeting): 
# 1. First check preprocessing 
    # - check if getFasta and the subsequent up/down regions on the chromosome are fetching the right things 
    # - check if one hot encoding is correct
    # - check if i'm fetching the correct TSS sites (use the experimentally validated TSS sites for testing)
# 2. Next check if the model is convolving correctly 
    # - try with smaller kernel size and smaller input size (ex: 4 x 10 input with kernel size of 5, and manually calculate to see if output shape is correct)
# 3. Try different model architecture, i.e. maybe just one dropout and one pooling layer 
# 4. If the RNAseq data fails, then try to make the ATAC-seq training work 
    # - implement above suggestions for the ATAC-seq data trained model, but in terms of input data, input the sequences corresponding to whole genes (not binned intervals)
# > If all else fails, try the CAGE dataset instead of RNAseq data 
# 5. IFF 1-4/> have been tried and done then scale the output of the model so that it oeprates for >1 dataset input 

# TODO: 
# plot the distribution of the counts 


# parse in the arguments 
# parser = ap.ArgumentParser(prog='This program runs the models')
# parser.add_argument('-primarymodel', type=str, help='Input model type')
# parser.add_argument('-trainfile', type=str, help='Input filepath to the training data', nargs=1)
# parser.add_argument('-testfile', type=str, help='Input filepath to the testing data', nargs=1)
# args = vars(parser.parse_args())

# for running with the debugger 
model_type = 'testing'
# train_filepath = '../0-data/2-final/1-finalcsv/agingeye/cts_bed_train.csv'
train_filepath = '/Users/nicolaskim/Desktop/research/singh_lab/basset_basenji/aging_eye_classification/0-data/2-final/1-finalcsv/agingeye/cts_bed_train.csv'
# test_filepath = '../0-data/2-final/1-finalcsv/agingeye/cts_bed_test.csv'
test_filepath = '/Users/nicolaskim/Desktop/research/singh_lab/basset_basenji/aging_eye_classification/0-data/2-final/1-finalcsv/agingeye/cts_bed_test.csv'

# load in model type parameter, training filepath, testing filepath 
# model_type = args['primarymodel']
# train_filepath = args['trainfile'][0]
# test_filepath = args['testfile'][0]

train_data = pd.read_csv(train_filepath, header=None)
test_data = pd.read_csv(test_filepath, header=None )
train_data.to_numpy()
sequence_length = 1000

# print dimensions for reference
print('===============================================================================')
print("This is the shape of the training data: ", train_data.shape)
print("This is a brief sample of the training data: ", train_data[:2])
test_data.to_numpy()

# set the device 
device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print(f"0. Using {device} device!")
print('===============================================================================')

# define the custom dataset for use 
class CustomDataset(Dataset): 
    def __init__(self, dataset, n_samples): 
        self.n_samples = n_samples
        self.features = []
        self.labels = []
        # iterate over all the training samples
        for i in range(self.n_samples): 
            s = np.zeros([sequence_length, 4]) 
            sequence = dataset[0][i]
            label = dataset[1][i]
            # iterate over each base within the sequence for that training example 
            for j in range(sequence_length): 
                base_at_jth_position = self.one_hot_encode(0, sequence[j]) # one hot encode sequences as: A T C G (top -> bottom)
                s[j][base_at_jth_position] = 1
            self.features.append(s)
            self.labels.append(label)
    def __len__(self):
        return self.n_samples
    def __getitem__(self, index):
        return self.features[index], self.labels[index]
    def one_hot_encode(self, marker, base):
        if base.upper() == 'A': marker = 0
        elif base.upper() == 'T': marker = 1
        elif base.upper() == 'C': marker = 2
        elif base.upper() == 'G': marker = 3
        return marker

class CustomDataset2(Dataset): 
    def __init__(self, dataset, n_samples): 
        self.n_samples = n_samples
        self.features = []
        self.labels = []
        # iterate over all the rows of the training dataset (training examples)
        for i in range(self.n_samples): 
            s = np.zeros([sequence_length, 4]) 
            sequence = dataset[0][i]
            label = dataset[1][i]
            # iterate over each base within the sequence for that training example 
            for j in range(sequence_length): 
                base_at_jth_position = self.one_hot_encode(0, sequence[j]) # one hot encode sequences as: A T C G (top -> bottom)
                s[j][base_at_jth_position] = 1
            self.features.append(s)
            self.labels.append(label)
    def __len__(self):
        return self.n_samples
    def __getitem__(self, index):
        return self.features[index], self.labels[index]
    def one_hot_encode(self, marker, base):
        if base.upper() == 'A': marker = 0
        elif base.upper() == 'T': marker = 1
        elif base.upper() == 'C': marker = 2
        elif base.upper() == 'G': marker = 3
        return marker


class BSampler(BatchSampler): 
    def __init__(self, num_rows_in_train, gene_ids, indices, batch_size):
        super(BSampler, self).__init__(train_dataset, batch_size, drop_last=False) # i forget why dataset is needed here 
        self.gene_ids = gene_ids
        self.indices = indices 
        self.num_batches = int(num_rows_in_train / batch_size)
        self.batch_size = batch_size
    def __iter__(self):
        batches = []
        for _ in range(self.num_batches):
            batch = []
            batch_gene = random.choice(self.gene_ids) 
            batches.append(batch)
        return iter(batches)
    def __len__(self):
        # this doesn't return anything informative unless i change the num_rows into constructor param
        return self.num_rows


# train dataset 
train_dataset = CustomDataset(dataset=train_data, n_samples=np.shape(train_data)[0])
train_dataloader = DataLoader(dataset=train_dataset, batch_size=256, shuffle=True)
train_feat, train_label = next(iter(train_dataloader))
n_total_steps = len(train_dataloader) 

# test dataset 
test_dataset = CustomDataset(dataset=test_data, n_samples=np.shape(test_data)[0])
test_dataloader = DataLoader(dataset=test_dataset, batch_size=256, shuffle=True)
test_feat, test_label = next(iter(test_dataloader))

basset_cnn_config = {
}

ff_nn_config = {
}
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
    
simple_config = {
    'conv_1_in': 4, 
    'conv_1_out': 100, 
    'conv_1_kernel_size': 50, 
    'pool_kernel_size': 10, 
    'dropout_probability': 0.5, 
    'num_samples_per_batch': 600, 
    'batch_size': 1,
    'output_size':1, 
    'num_features': 4, 
    'fc1_out': 2500, 
    'fc2_out': 600, 
}

class SimpleCNN(nn.Module):
    def __init__(self): 
        super(SimpleCNN, self).__init__()
        self.conv_one = nn.Conv1d(in_channels=simple_config['conv_1_in'], out_channels=simple_config['conv_1_out'], kernel_size=simple_config['conv_1_kernel_size']) # output dim: (1, 300, 582)
        self.max_pool = nn.MaxPool1d(kernel_size=simple_config['pool_kernel_size']) # output dim: 
        self.dropout = nn.Dropout1d(simple_config['dropout_probability'])
        self.linear_1 = nn.Linear(in_features=int((((simple_config['num_samples_per_batch']-simple_config['conv_1_kernel_size'])*simple_config['batch_size'])/simple_config['pool_kernel_size'])*simple_config['conv_1_out']), out_features=simple_config['fc1_out'])
        self.linear_2 = nn.Linear(in_features=simple_config['fc1_out'], out_features=simple_config['fc2_out'])
        self.linear_3 = nn.Linear(in_features=simple_config['fc2_out'], out_features=simple_config['output_size'])
    def forward(self, x): 
        x = self.conv_one(x)
        # print('first conv layer: ', x.shape)
        x = F.relu(x)
        x = self.max_pool(x)
        # print('first pooling layer: ', x.shape)
        x = self.dropout(x)
        # print('dropout: ', x.shape)
        x = x.view(int((((simple_config['num_samples_per_batch']-simple_config['conv_1_kernel_size'])*simple_config['batch_size'])/simple_config['pool_kernel_size'])*simple_config['conv_1_out']))
        # print('post flattening: ', x.shape)
        x = self.linear_1(x)
        # print('first linear: ' ,x.shape)
        x = F.relu(x)
        x = self.linear_2(x)
        # print('second linear: ', x.shape)
        x = F.relu(x)
        x = self.linear_3(x)
        # print('third linear: ', x.shape)
        x = torch.sigmoid(x)
        # print('simgoid', x.shape)
        return x


# initialize and send the model to the device above
# if model_type == "cnn": 
    # model = BassetCNN().to(device)
if model_type == "nn":
    model = FeedForwardNN().to(device)
if model_type == "testing": 
    model = SimpleCNN().to(device)

# other params 
data_aug = aug.DataAugmentation
learning_rate = 0.001 
criterion = nn.BCELoss() 
# optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
# optimizer = torch.optim.RMSprop(model.parameters(), lr=learning_rate)
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# try Adam as well
num_epochs = 50
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

            # tranpose the sample so its (1, 4, 600) for input into the model
            samples = samples.permute(0, 2, 1)

            # every other epoch, reverse complement the sample 
            # if epoch % 2 == 1: 
            #     # samples = data_aug.reverse_complement_sequence(data_aug, samples)
            #     samples = data_aug.reverse_complement_sequence(data_aug, samples)

            samples = samples.type(torch.float32) # send the samples to the device
            samples = samples.to(device)

            # assume the first row is the labels? 
            labels = labels.type(torch.float32)
            labels = labels.to(device) 
            model.to(device) 
            predicted = model(samples) 
            loss = criterion(predicted, labels) 
            t_loss_per_batch.append(loss.item()) 

            if (predicted.item() > 0.5 and labels.item() == 1) or (predicted.item() <= 0.5 and labels.item() == 0): 
                t_accuracy_per_batch.append(1)
            elif (predicted.item() > 0.5 and labels.item() == 0) or (predicted.item() <= 0.5 and labels.item() == 1): 
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
    PATH = '../4-saved/basset_params.pth'
    torch.save(model.state_dict(), PATH)

    # return the log of loss values 
    return train_log 

# test the model - load in saved params from specified PATH
def test_model(PATH):
    model.load_state_dict(torch.load(PATH)) 

    # switch to eval mode to switch off layers like dropout
    model.eval() 

    test_log = {'testing_loss_per_epoch':[], 'testing_accuracy_per_epoch': []}
    
    with torch.no_grad():
        for epoch in range(num_epochs):
            testing_loss_per_batch = []
            testing_accuracy_per_batch = []
            for i, (samples, labels) in enumerate(test_dataloader):
                
                samples = torch.reshape(samples, (samples.shape[0]*samples.shape[1], samples.shape[2])).permute(1, 0)
                samples = samples.type(torch.float32)
                samples = samples.to(device)
                labels = labels.type(torch.float32)
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
# test_log = test_model('../4-saved/basset_params.pth')

# plots the loss 
def plot_loss(id, num_epochs, train_loss, test_loss):
    f, a = plt.subplots(figsize=(10,7.5), layout='constrained') # don't need to specify 2x2 or anything here, bc i'm just going to plot the loss 
    f.suptitle('Calculated Loss')
    a.plot(num_epochs, train_loss, label='Training Loss')
    a.plot(num_epochs, test_loss, label=f'Testing Loss')
    a.set_xlabel('Number of Epochs')
    a.set_ylabel('Average Loss')
    a.set_title(f'Training and Testing Loss')
    a.legend()
    plt.savefig(f'../3-viz/loss_curves_{id}')
 
# plots the accuracy
def plot_accuracy(id, num_epochs, train_accuracy, test_accuracy): 
    f, a = plt.subplots(figsize=(10,7.5), layout='constrained') # don't need to specify 2x2 or anything here, bc i'm just going to plot the loss 
    f.suptitle('Calculated Loss')
    a.plot(num_epochs, train_accuracy, label='Training Accuracy')
    a.plot(num_epochs, test_accuracy, label=f'Testing Accuracy')
    a.set_xlabel('Number of Epochs')
    a.set_ylabel('Accuracy')
    a.set_title(f'Training and Testing Accuracy')
    a.legend()
    plt.savefig(f'../3-viz/accuracy_curves_{id}')

# call the loss function plotter 
# plot_loss("04", np.linspace('ff_nn_50_epoch', num_epochs, num=num_epochs).astype(int), train_log['training_loss_per_epoch'], test_log['testing_loss_per_epoch'])
# plot_accuracy("04", np.linspace('ff_nn_50_epoch', num_epochs, num=num_epochs).astype(int), train_log['training_accuracy_per_epoch'], test_log['testing_accuracy_per_epoch'])

