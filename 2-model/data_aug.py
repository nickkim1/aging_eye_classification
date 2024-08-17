import torch
import numpy as np

# class for augmenting the initial dataset, because the model isn't learning anything 
class DataAugmentation():
    # don't want to do anything more  
    def __init__(self):
        pass 
    def invert_sequence(self, sequence_tensor:torch.tensor):
        return torch.flip(sequence_tensor, [1])
    
    def reverse_complement_sequence(self, sequence_tensor:torch.tensor): 
        # complement, then reverse 
        for i in range(sequence_tensor.shape[1]): 
            # assume A - T - C - G
            base_vector = sequence_tensor[:, i]
            # print('before', base_vector)
            if base_vector[0] == 1: # A
                base_vector[0] = 0
                base_vector[1] = 1 # T
            elif base_vector[1] == 1: 
                base_vector[1] = 0 
                base_vector[0] = 1 # A 
            elif base_vector[2] == 1: # C
                base_vector[2] = 0
                base_vector[3] = 1
            elif base_vector[3] == 1: # G
                base_vector[3] = 0 
                base_vector[2] = 1
        return torch.flip(sequence_tensor, [1])

# brief testing 
test_class = DataAugmentation()
seq = torch.tensor([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]]).reshape(4, 3)
reverse_complemented_seq = test_class.reverse_complement_sequence(seq) # <- [A, T, C] => [G, A, T]
# print('after', reverse_complemented_seq)