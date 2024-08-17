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
from config import populate_settings

# parse in the arguments 
parser = ap.ArgumentParser(prog='This program runs the models')
parser.add_argument('-primarymodel', type=str, help='Input model type')
parser.add_argument('-trainfile', type=str, help='Input filepath to the training data', nargs=1)
parser.add_argument('-testfile', type=str, help='Input filepath to the testing data', nargs=1)
parser.add_argument('-seqlen', type=int, help='Input length of the sequences in your data')
parser.add_argument('-optim', type=str, help='Input optimizer type')
parser.add_argument('-loss', type=str, help='Input loss type')
parser.add_argument('-batchsize',type=int, help='Input batch size')
args = vars(parser.parse_args())

# load in various parameters from the command line 
model_type = args['primarymodel']
train_filepath = args['trainfile'][0]
test_filepath = args['testfile'][0]
optim = args['optim']
loss = args['loss']
seqlen = args['seqlen']
batch_size = args['batchsize']
settings = populate_settings(train_filepath, test_filepath, seqlen, loss, optim, model_type, batch_size)

# read in the data 
train_data = pd.read_csv(settings['data']['train'], header=None)
test_data = pd.read_csv(settings['data']['test'], header=None)
sequence_length = settings['data']['seqlen']
train_data.to_numpy()
test_data.to_numpy()

# set the device 
device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print(f"0. Using {device} device!")

# define custom dataset for regression 
class RegressionDataset(): 
    def __init__(self, dataset, n_samples): 
        self.n_samples = n_samples
        self.features = []
        self.labels = []
        self.num_age_pts = 8
        for i in range(self.n_samples): # iterate over all the training samples
            s = np.zeros([sequence_length, 4]) 
            sequence = dataset[0][i]
            label = []
            for l in range(1, self.num_age_pts+1): # Skip the sequence column
                label.append(dataset[l][i]) # Get the age labels @ our current sample 
            label = np.array(label)
            for j in range(sequence_length): # iterate over each base within the sequence for that training example 
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

# train dataset 
train_dataset = RegressionDataset(dataset=train_data, n_samples=train_data.shape[0])
train_dataloader = DataLoader(dataset=train_dataset, batch_size=settings['hyperparameters']['batch_size'], shuffle=True)
train_feat, train_label = next(iter(train_dataloader))
n_total_steps = len(train_dataloader) 

# test dataset 
test_dataset = RegressionDataset(dataset=test_data, n_samples=test_data.shape[0])
test_dataloader = DataLoader(dataset=test_dataset, batch_size=settings['hyperparameters']['batch_size'], shuffle=True)
test_feat, test_label = next(iter(test_dataloader))

# initialize and send the model to the device above
model = settings['model'].to(device)

# other params that you load in from the config file itself 
data_aug = aug.DataAugmentation
num_epochs = settings['hyperparameters']['num_epochs']
learning_rate = settings['hyperparameters']['learning_rate']
criterion = settings['loss']
optimizer = settings['optimizer']
n_batches = len(train_dataloader)

# train the model: 
def train_model():  
    
    # initialize the weights. unsure if they did random initialization however 
    def init_weights(m):
        if isinstance(m, nn.Linear) or isinstance(m, nn.Conv1d):
            torch.nn.init.xavier_uniform_(m.weight) 
            m.bias.data.fill_(0.01)
    
    # apply weight + bias initialization for the model
    model.apply(init_weights)
    
    # keep track of the parameters with this dictionary 
    train_log = {'training_loss_per_epoch':[], 'training_accuracy_per_epoch':[]} 
    test_log = {'testing_loss_per_epoch':[], 'testing_accuracy_per_epoch': []}

    for epoch in range(num_epochs):  

        # switch the model to training mode 
        model.train() 
        t_loss_per_batch = []
        r2_per_batch = []

        for i, (samples, labels) in enumerate(train_dataloader): 
            samples = samples.permute(0, 2, 1)  # tranpose the sample so its (batch_size, 4, 600) for input into the model
            samples = samples.type(torch.float32) # send the samples to the device
            samples = samples.to(device)

            # assume the first row is the labels? 
            labels = labels.type(torch.float32)
            labels = labels.to(device) 
            predicted = model(samples) 

            r2_score = metrics.r2_score(labels.detach().cpu().numpy(), predicted.detach().cpu().numpy())
            loss = criterion(predicted, labels) 
            t_loss_per_batch.append(loss.item()) 
            r2_per_batch.append(r2_score)

            # zero accumulated gradients, backprop, and step on params 
            optimizer.zero_grad()
            loss.backward() 
            optimizer.step() 

            # print message for every batch 
            if (i+1) % n_batches == 0: # -- this is where the i term above is used in for loop
                print (f'Epoch [{epoch+1}/{num_epochs}], Step [{i+1}/{int(n_total_steps)}]')
        
        # insert the AVERAGE batch training loss -> log. technically could use total_batches here but opted not to
        t_loss = sum(t_loss_per_batch) / len(t_loss_per_batch)
        train_log['training_loss_per_epoch'].append(t_loss)
        print('Train Loss: ', t_loss)

        # avg r2 score
        r2_score = sum(r2_per_batch) / len(r2_per_batch)
        print('Train R2: ', r2_score)

        # save the model's params after fully training it 
        PATH = '../../4-saved/basset_params.pth'
        torch.save(model.state_dict(), PATH)

        # run the model on the testing data after each training epoch  
        test_model(PATH, test_log)

    # return the log of loss values 
    return train_log, test_log

# test the model - load in saved params from specified PATH
def test_model(PATH, test_log):
    model.load_state_dict(torch.load(PATH)) 

    # switch to eval mode to switch off layers like dropout
    model.eval() 
    
    with torch.no_grad():
        testing_loss_per_batch = []
        r2_per_batch = []

        for i, (samples, labels) in enumerate(test_dataloader):
            
            samples = samples.permute(0, 2, 1)
            samples = samples.type(torch.float32)
            samples = samples.to(device)

            labels = labels.type(torch.float32)
            labels = labels.to(device)

            predicted = model(samples)

            r2_score = metrics.r2_score(labels.detach().cpu().numpy(), predicted.detach().cpu().numpy())
            r2_per_batch.append(r2_score)

            loss = criterion(predicted, labels)
            testing_loss_per_batch.append(loss.item())

        test_loss = sum(testing_loss_per_batch) / len(testing_loss_per_batch)
        print('Test Loss', test_loss)

        r2_score = sum(r2_per_batch) / len(r2_per_batch)
        print('Test R2: ', r2_score)

        test_log['testing_loss_per_epoch'].append(test_loss)

# call the methods to train and test the model 
train_log, test_log = train_model()

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
    info_text = '\n'.join([f'{key}: {value}' for key, value in settings.items()])
    plt.text(1, 4.5, info_text, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    plt.savefig(f'../../3-viz/debug/4-22-debug/loss_curves_{id}')
 
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
    info_text = '\n'.join([f'{key}: {value}' for key, value in settings.items()])
    plt.text(1, 4.5, info_text, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    plt.savefig(f'../../3-viz/debug/4-22-debug/accuracy_curves_{id}')

# call the loss function plotter 
model_name, optim_name, batch_size = settings['model_name'], settings['optim_name'], settings['hyperparameters']['batch_size']
loss_message = f'{model_name}-loss-|{optim_name}_optimizer_{sequence_length}_seqlen_{batch_size}_batchsize_{num_epochs}_epochs|'
accuracy_message = f'{model_name}-accuracy-|{optim_name}_optimizer_{sequence_length}_seqlen_{batch_size}_batchsize_{num_epochs}_epochs|'
plot_loss(loss_message, np.linspace(1, num_epochs, num=num_epochs).astype(int), train_log['training_loss_per_epoch'], test_log['testing_loss_per_epoch'])
plot_accuracy(accuracy_message, np.linspace(1, num_epochs, num=num_epochs).astype(int), train_log['training_accuracy_per_epoch'], test_log['testing_accuracy_per_epoch'])