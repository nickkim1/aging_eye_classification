

# ====================================== OLD CODE  ======================================
# train_dataloader = DataLoader(dataset=train_dataset, batch_sampler=BSampler(num_rows=len(train_dataset), gene_ids=train_dataset.gene_ids_and_indices()[0], indices=train_dataset.gene_ids_and_indices()[1], batch_size=1))
# test_dataloader = DataLoader(dataset=test_dataset, batch_sampler=BSampler(num_rows=len(test_dataset), gene_ids=test_dataset.gene_ids_and_indices()[0], indices=test_dataset.gene_ids_and_indices()[1], batch_size=1))

#  custom batch sampler class (for generating the 4 x 600 inputs) et
# class BSampler(BatchSampler): 
#     def __init__(self, num_rows, gene_ids, indices, batch_size):
#         super(BSampler, self).__init__(train_dataset, batch_size, drop_last=False) # i forget why dataset is needed here 
#         self.gene_ids = gene_ids
#         self.num_rows = num_rows
#         self.indices = indices 
#     def __iter__(self):
#         np.random.seed(0)
#         # np.random.shuffle(self.gene_ids)
#         batches = []
#         for ignore in range(self.num_rows):
#             # randomly choose an idx (sample) from 0 to the full dataset 
#             batch = [random.choice(self.gene_ids)]
#             batches.append(batch)
#         return iter(batches)
#     def __len__(self):
#         # this doesn't return anything informative unless i change the num_rows into constructor param
#         return self.num_rows

# define the CNN architecture for basset
# class BassetCNN(nn.Module):
#     def __init__(self): 
#         super(BassetCNN, self).__init__()
#         self.conv_one = nn.Conv1d(in_channels=4, out_channels=300, kernel_size=19)
#         self.batchnorm_one = nn.BatchNorm1d(300) # input = num channels 
#         self.pool_one = nn.MaxPool1d(kernel_size=3) 
#         self.conv_two = nn.Conv1d(in_channels=300, out_channels=200, kernel_size=11)
#         self.batchnorm_two = nn.BatchNorm1d(200) # input = num channels 
#         self.pool_two = nn.MaxPool1d(4)
#         self.conv_three = nn.Conv1d(in_channels=200, out_channels=200, kernel_size=7)
#         self.batchnorm_three = nn.BatchNorm1d(200) # input = num channels 
#         self.pool_three = nn.MaxPool1d(kernel_size=4)
#         self.fc1 = nn.Linear(2000, 1000) # 1000 unit linear layer (same as paper)
#         self.dropout_one = nn.Dropout1d(0.3)  
#         self.fc2 = nn.Linear(1000, 1000) # 
#         self.dropout_two = nn.Dropout1d(0.3)
#         self.fc3 = nn.Linear(1000, 1) # output dim should be [1x164] since i unsqueezed @ flattening above
#     def forward(self, x): 
#         x = self.conv_one(x)
#         x = self.batchnorm_one(x)
#         x = F.relu(x)
#         x = self.pool_one(x)
#         x = self.conv_two(x)
#         x = self.batchnorm_two(x)
#         x = F.relu(x)
#         x = self.pool_two(x)
#         x = self.conv_three(x)
#         x = self.batchnorm_three(x)
#         x = F.relu(x)
#         x = self.pool_three(x)
#         x = self.fc1(x.unsqueeze(0)) # unsqueeze after flattening
#         x = F.relu(x)
#         x = self.dropout_one(x)
#         x = self.fc2(x)
#         x = self.dropout_two(x)
#         x = self.fc3(x)
#         x = F.sigmoid(torch.tensor(torch.flatten(x), dtype=torch.float32, requires_grad=True))
#         return x
