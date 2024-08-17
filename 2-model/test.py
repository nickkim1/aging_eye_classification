import torch 
import torch.nn as nn

# initialize input
seq = torch.randn(1, 4, 600)
conv = nn.Conv1d(in_channels=4, out_channels=300, kernel_size=19)
output = conv(seq)
print(output.shape)
