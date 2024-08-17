from matplotlib import pyplot as plt

'''
This script is meant for getting working visualizations/metrics for the model
'''

# Loss plotter 
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
 
# Accuracy plotter 
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
