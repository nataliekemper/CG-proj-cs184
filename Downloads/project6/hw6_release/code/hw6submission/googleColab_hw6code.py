# -*- coding: utf-8 -*-
"""Copy of CS189_HW6_NN_final.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1lK4w1gaarDk_rTvYyXNYPe6ClUY4ZMr9

# CS 189 HW 6: Neural Networks
**Note:** before starting this notebook, please make a copy of it, otherwise your changes will not persist.

This part of the assignment is designed to get you familiar with how engineerings in the real world train neural network systems. It isn't designed to be difficult. In fact, everything you need to complete the assignment is available directly on the pytorch website [here](https://pytorch.org/tutorials/beginner/blitz/neural_networks_tutorial.html). This note book will have the following components:

1. Understanding the basics of Pytorch (no deliverables)
2. Training a simple neural network on MNIST (Deliverable = training graphs)
3. Train a model on CIFAR-10 for Kaggle (Deliverable = kaggle submission and explanation of methods)

The last part of this notebook is left open for you to explore as many techniques as you want to do as well as possible on the dataset.

You will also get practice being an ML engineer by reading documentation and using it to implement models. The first section of this notebook will cover an outline of what you need to know -- we are confident that you can find the rest on your own.

Note that like all other assignments, you are free to use this notebook or not. You just need to complete the deliverables and turn in your code. If you want to run everything outside of the notebook, make sure to appropriately install pytorch to download the datasets and copy out the code for kaggle submission. If you don't want to use pytorch and instead want to use Tensorflow, feel free, but you may still need to install pytorch to download the datasets.
"""

# Imports for pytorch
import numpy as np
import torch
import torchvision
from torch import nn
import matplotlib
from matplotlib import pyplot as plt
import tqdm

"""# 1. Understanding Pytorch

Pytorch is based on the "autograd" paradigm. Essentially, you perform operations on multi-dimensional arrays like in numpy, except pytorch will automatically handle gradient tracking. In this section you will understand how to use pytorch.

This section should help you understand the full pipeline of creating and training a model in pytorch. Feel free to re-use code from this section in the assigned tasks.

Content in this section closely follows this pytorch tutorial: https://pytorch.org/tutorials/beginner/basics/intro.html

## Tensors

Tensors can be created from numpy data or by using pytorch directly.
"""

data = [[1, 2],[3, 4]]
x_data = torch.tensor(data)

np_array = np.array(data)
x_np = torch.from_numpy(np_array)

shape = (2,3,)
rand_tensor = torch.rand(shape)
np_rand_array = rand_tensor.numpy()

print(f"Tensor from np: \n {x_np} \n")
print(f"Rand Tensor: \n {rand_tensor} \n")
print(f"Rand Numpy Array: \n {np_rand_array} \n")

"""They also support slicing and math operations very similar to numpy. See the examples below:"""

# Slicing
tensor = torch.ones(4, 4)
print('First row: ',tensor[0])
print('First column: ', tensor[:, 0])

# Matrix Operations
y1 = tensor @ tensor.T
y2 = tensor.matmul(tensor.T)

# Getting a single item
scalar = torch.sum(y1) # sums all elements
item = scalar.item()
print("Sum as a tensor:", scalar, ", Sum as an item:", item)

"""## Autograd
This small section shows you how pytorch computes gradients. When we create tenors, we can set `requires_grad` to be true to indicate that we are using gradients. For most of the work that you actually do, you will use the `nn` package, which automatically sets all parameter tensors to have `requires_grad=True`.
"""

# Below is an example of computing the gradient for a single data point in logistic regression using pytorch's autograd.

x = torch.ones(5)  # input tensor
y = torch.zeros(1) # label
w = torch.randn(5, 1, requires_grad=True)
b = torch.randn(1, requires_grad=True) 
pred = torch.sigmoid(torch.matmul(x, w) + b)
loss = torch.nn.functional.binary_cross_entropy(pred, y)
loss.backward() # Computers gradients
print("W gradient:", w.grad)
print("b gradient:", b.grad)

# when we want to actually take an update step, we can use optimizers:
optimizer = torch.optim.SGD([w, b], lr=0.1)
print("Weight before", w)
optimizer.step() # use the computed gradients to update 
# Print updated weights
print("Updated weight", w)

# Performing operations with gradients enabled is slow...
# You can disable gradient computation using the following enclosure:
with torch.no_grad():
    # Perform operations without gradients
    ...

"""## Devices
Pytorch supports accelerating computation using GPUs which are available on google colab. To use a GPU on google colab, go to runtime -> change runtime type -> select GPU.

Note that there is some level of strategy for knowing when to use which runtime type. Colab will kick users off of GPU for a certain period of time if you use it too much. Thus, its best to run simple models and prototype to get everything working on CPU, then switch the instance type over to GPU for training runs and parameter tuning.

Its best practice to make sure your code works on any device (GPU or CPU) for pytorch, but note that numpy operations can only run on the CPU. Here is a standard flow for using GPU acceleration:
"""

# Determine the device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device", device)
# Next create your tensors
tensor = torch.zeros(4, 4, requires_grad=True)
# Move the tensor to the device you want to use
tensor = tensor.to(device)

# Perform whatever operations you want.... (often this will involve gradients)
# These operations will be accelerated by GPU.
tensor = 10*(tensor + 1)

# bring the tensor back to CPU, first detaching it from any gradient computations
tensor = tensor.detach().cpu()

tensor_np = tensor.numpy() # Convert to numpy if you want to perform numpy operations.

"""## The NN Package
Pytorch implements composable blocks in `Module` classes. All layers and modules in pytorch inherit from `nn.Module`. When you make a module you need to implement two functions: `__init__(self, *args, **kwargs)` and `foward(self, *args, **kwargs)`. Modules also have some nice helper functions, namely `parameters` which will recursively return all of the parameters. Here is an example of a logistic regression model:
"""

class Perceptron(nn.Module):
  def __init__(self, in_dim):
    super().__init__()
    self.layer = nn.Linear(in_dim, 1) # This is a linear layer, it computes Xw + b

  def forward(self, x):
    return torch.sigmoid(self.layer(x)).squeeze(-1)

perceptron = Perceptron(10)
perceptron = perceptron.to(device) # Move all the perceptron's tensors to the device
print("Parameters", list(perceptron.parameters()))

"""## Datasets

Pytorch has nice interfaces for using datasets. Suppose we create a logistic regression dataset as follows:
"""

c1_x1, c1_x2 = np.random.multivariate_normal([-2.5,3], [[1, 0.3],[0.3, 1]], 100).T
c2_x1, c2_x2 = np.random.multivariate_normal([1,1], [[2, 1],[1, 2]], 100).T
c1_X = np.vstack((c1_x1, c1_x2)).T
c2_X = np.vstack((c2_x1, c2_x2)).T
train_X = np.concatenate((c1_X, c2_X))
train_y = np.concatenate((np.zeros(100), np.ones(100)))
# Shuffle the data
permutation = np.random.permutation(train_X.shape[0])
train_X = train_X[permutation, :]
train_y = train_y[permutation]
# Plot the data
plt.plot(c1_x1, c1_x2, 'x')
plt.plot(c2_x1, c2_x2, 'o')
plt.axis('equal')
plt.show()

"""We can then create a pytorch dataset object as follows. Often times, the default pytorch datasets will create these objects for you. Then, we can apply dataloaders to iterate over the dataset in batches."""

dataset = torch.utils.data.TensorDataset(torch.from_numpy(train_X), torch.from_numpy(train_y))
# We can create a dataloader that iterates over the dataset in batches.
dataloader = torch.utils.data.DataLoader(dataset, batch_size=10, shuffle=True)
for x, y in dataloader:
    print("Batch x:", x)
    print("Batch y:", y)
    break

# Clean up the dataloader as we make a new one later
del dataloader

"""## Training Loop Example
Here is an example of training a full logistic regression model in pytorch. Note the extensive use of modules -- modules can be used for storing networks, computation steps etc.
"""

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device", device)

epochs = 10
batch_size = 10
learning_rate = 0.01

num_features = dataset[0][0].shape[0]
model = Perceptron(num_features).to(device)
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
criterion = torch.nn.BCELoss()
dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)

model.train() # Put model in training mode
for epoch in range(epochs):
    training_losses = []
    for x, y in tqdm.tqdm_notebook(dataloader, unit="batch"):
        x, y = x.float().to(device), y.float().to(device)
        optimizer.zero_grad() # Remove the gradients from the previous step
        pred = model(x)
        loss = criterion(pred, y)
        loss.backward()
        optimizer.step()
        training_losses.append(loss.item())
    print("Finished Epoch", epoch + 1, ", training loss:", np.mean(training_losses))
        
# We can run predictions on the data to determine the final accuracy.
with torch.no_grad():
    model.eval() # Put model in eval mode
    num_correct = 0
    for x, y in dataloader:
        x, y = x.float().to(device), y.float().to(device)
        pred = model(x)
        num_correct += torch.sum(torch.round(pred) == y).item()
    print("Final Accuracy:", num_correct / len(dataset))
    model.train() # Put model back in train mode

"""# Task 1: MLP For FashionMNIST
Earlier in this course you trained SVMs and GDA models on MNIST. Now you will train a multi-layer perceptron model on an MNIST-like dataset. Your deliverables are as follows:

1. Code for training an MLP on MNIST (can be in code appendix, tagged in your submission).
2. A plot of the training loss and validation loss for each epoch of training after trainnig for at least 8 epochs.
3. A plot of the training and validation accuracy, showing that it is at least 82% for validation by the end of training. 

Below we will create the training and validation datasets for you, and provide a very basic skeleton of the code. Please leverage the example training loop from above.

Some pytorch components you should definetily use:
1. `nn.Linear`
2. Some activation: `nn.ReLU`, `nn.Tanh`, `nn.Sigmoid`, etc.
3. `nn.CrossEntropyLoss`

Here are challenges you will need to overcome:
1. The data is default configured in image form ie (28 x 28), versus one feature vector. You will need to reshape it somewhere to feed it in as vector to the MLP. There are many ways of doing this.
2. You need to write code for plotting.
3. You need to find appropriate hyper-parameters to achieve good accuracy.

Your underlying model must be fully connected or dense, and may not have convolutions etc., but you can use anything in torch.optim or any layers in torch.nn besides nn.Linear that do not have weights. 
"""

# Creating the datasets
transform = torchvision.transforms.ToTensor() # feel free to modify this as you see fit.

training_data = torchvision.datasets.FashionMNIST(
    root="data",
    train=True,
    download=True,
    transform=transform,
)

validation_data = torchvision.datasets.FashionMNIST(
    root="data",
    train=False,
    download=True,
    transform=transform,
)

class MLP(nn.Module):
    def __init__(self, in_size = 784, out_size = 10, layers = [120, 84]):
        super().__init__()
        self.fc1 = nn.Linear(in_size, layers[0])
        self.fc2 = nn.Linear(layers[0], layers[1])
        self.fc3 = nn.Linear(layers[1], out_size)
        self.relu = nn.ReLU()
        self.logsoftmax = nn.LogSoftmax()
    
    def forward(self, X):
        X = self.relu(self.fc1(X))
        X = self.relu(self.fc2(X))
        X = self.fc3(X)
        X = self.logsoftmax(X)
        return X

from re import X
### YOUR CODE HERE ###
epochs = 10
batch_size = 10
learning_rate = 0.01
torch.manual_seed(12)

model = MLP()
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
criterion = torch.nn.CrossEntropyLoss()
 
train_loader = torch.utils.data.DataLoader(training_data, batch_size=batch_size, shuffle=True)
test_loader = torch.utils.data.DataLoader(validation_data, batch_size=batch_size, shuffle=False)

model.train() # Put model in training mode

training_losses = []
validation_losses = []

training_accuracy = []
validation_accuracy = []

for epoch in range(epochs):

    epoch_losses = []
    train_epoch_accuracy = []

    for x, y in tqdm.notebook.tqdm(train_loader, unit = "batch"):
        # flattening data
        x_train = x.view(batch_size, -1)
        optimizer.zero_grad() # Remove the gradients from the previous step
        pred = model(x_train)
        loss = criterion(pred, y)

        loss.backward()
        optimizer.step()
        epoch_losses.append(loss.item())

    print("Finished Epoch", epoch + 1, ", training loss:", np.mean(epoch_losses))
    training_losses.append(np.mean(epoch_losses))
    
    with torch.no_grad():
      model.eval() # Put model in eval mode
      num_correct = 0
      for x, y in train_loader:
          x_train = x.view(batch_size, -1)
          pred = model(x_train)

          num_correct += torch.sum(torch.round(pred) == y).item()
          train_epoch_accuracy.append(num_correct)
      print("Accuracy for Epoch", epoch + 1, ":", num_correct / len(train_loader))
      training_accuracy.append(num_correct / len(train_loader))

      model.train() # Put model back in train mode

    # run through validation data

    val_epoch_losses = []
    val_epoch_accuracy = []

    for x, y in tqdm.notebook.tqdm(test_loader, unit = "batch"):
      optimizer.zero_grad()
      x_val = x.view(batch_size, -1)
      pred = model(x_val)
      loss = criterion(pred, y)

      loss.backward()
      optimizer.step()
      val_epoch_losses.append(loss.item())

    print("Finished Epoch", epoch + 1, ", validation loss:", np.mean(val_epoch_losses))
    validation_losses.append(np.mean(val_epoch_losses))

    with torch.no_grad():
      model.eval() # Put model in eval mode
      num_correct = 0
      for x, y in test_loader:
          x_train = x.view(batch_size, -1)
          pred = model(x_train)

          num_correct += torch.sum(torch.round(pred) == y).item()
          val_epoch_accuracy.append(num_correct)
      print("Accuracy for Epoch", epoch + 1, ":", num_correct / len(test_loader))
      validation_accuracy.append(num_correct / len(test_loader))

      model.train() # Put model back in train mode

training_accuracy = []
validation_accuracy = []

for epoch in range(epochs):
  train_epoch_accuracy = []
  val_epoch_accuracy = []
# We can run predictions on the data to determine the final accuracy.
  with torch.no_grad():
    model.eval() # Put model in eval mode
    num_correct = 0
    for x, y in train_loader:
        x_train = x.view(batch_size, -1)
        pred = model(x_train)

        num_correct += torch.sum(torch.round(pred) == y).item()
        train_epoch_accuracy.append(num_correct)
    print("Accuracy for Epoch", epoch + 1, ":", num_correct / len(train_loader))
    training_accuracy.append(num_correct / len(train_loader))

model.train() # Put model back in train mode

# We can run predictions on the data to determine the final accuracy.
with torch.no_grad():
    model.eval() # Put model in eval mode
    num_correct = 0
    for x, y in test_loader:
        x_train = x.view(batch_size, -1)
        pred = model(x_train)

        num_correct += torch.sum(torch.round(pred) == y).item()
    print("Final Accuracy:", num_correct / len(test_loader))
    model.train() # Put model back in train mode

plt.plot(training_losses, label = 'training loss')
plt.plot(validation_losses, label = 'validation loss')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Loss of Training and Validation Per Training Epcoh")
plt.legend();

plt.plot(training_accuracy, label = 'training accuracy')
plt.plot(validation_accuracy, label = 'validation accuracy')
plt.xlabel("Epoch")
plt.ylabel("Accuracy")
plt.title("Accuracy of Training and Validation Per Training Epcoh")
plt.legend();

"""# Task 2: CNNs for CIFAR-10

In this section, you will create a CNN for the CIFAR dataset, and submit your predictions to Kaggle. It is recommended that you use GPU acceleration for this part.

Here are some of the components you should consider using:
1. `nn.Conv2d`
2. `nn.ReLU`
3. `nn.Linear`
4. `nn.CrossEntropyLoss`
5. `nn.MaxPooling2d` (though many implementations without it exist)

We encourage you to explore different ways of improving your model to get higher accuracy. Here are some suggestions for things to look into:
1. CNN architectures: AlexNet, VGG, ResNets, etc.
2. Different optimizers and their parameters (see torch.optim)
3. Image preprocessing / data augmentation (see torchvision.transforms)
4. Regularization or dropout (see torch.optim and torch.nn respectively)
5. Learning rate scheduling: https://pytorch.org/docs/stable/optim.html#how-to-adjust-learning-rate
6. Weight initialization: https://pytorch.org/docs/stable/nn.init.html

Though we encourage you to explore, there are some rules:
1. You are not allowed to install or use packages not included by default in the Colab Environment.
2. You are not allowed to use any pre-defined architectures or feature extractors in your network.
3. You are not allowed to use **any** pretrained weights, ie no transfer learning.
4. You cannot train on the test data.

Otherwise everything is fair game.

Your deliverables are as follows:
1. Submit to Kaggle and include your test accuracy in your report.
2. Provide at least (1) training curve for your model, depicting loss per epoch or step after training for at least 8 epochs.
3. Explain the components of your final model, and how you think your design choices contributed to it's performance.

After you write your code, we have included skeleton code that should be used to submit predictions to Kaggle. **You must follow the instructions below under the submission header**. Note that if you apply any processing or transformations to the data, you will need to do the same to the test data otherwise you will likely achieve very low accuracy. 

It is expected that this task will take a while to train. Our simple solution achieves a training accuracy of 90.2% and a test accuracy of 74.8% after 10 epochs (be careful of overfitting!). This easily beats the best SVM based CIFAR10 model submitted to the HW 1 Kaggle! It is possible to achieve 95% or higher test accuracy on CIFAR 10 with good model design and tuning.
"""

# Creating the datasets, feel free to change this as long as you do the same to the test data.
# You can also modify this to split the data into training and validation.
# See https://pytorch.org/docs/stable/data.html#torch.utils.data.random_split

transform = torchvision.transforms.ToTensor()

training_data = torchvision.datasets.CIFAR10(
    root="data",
    train=True,
    download=True,
    transform=transform,
)
# If you make a train-test partition it is up to you.

batch_size = 10

cifar_train, cifar_validation = torch.utils.data.random_split(training_data, [0.75, 0.25])

train_loader = torch.utils.data.DataLoader(cifar_train, batch_size = batch_size, shuffle=True)

test_loader = torch.utils.data.DataLoader(cifar_validation, batch_size = batch_size, shuffle=True)

class CIFAR10CNN(nn.Module):
    def __init__(self):
        super(CIFAR10CNN, self).__init__()
        
        # First convolutional layer
        self.conv1 = nn.Conv2d(3, 32, 3, padding=1)
        self.relu1 = nn.ReLU()
        self.pool1 = nn.MaxPool2d(2, 2)
        
        # Second convolutional layer
        self.conv2 = nn.Conv2d(32, 64, 3, padding=1)
        self.relu2 = nn.ReLU()
        self.pool2 = nn.MaxPool2d(2, 2)
        
        # Third convolutional layer
        self.conv3 = nn.Conv2d(64, 128, 3, padding=1)
        self.relu3 = nn.ReLU()
        self.pool3 = nn.MaxPool2d(2, 2)
        
        # Fully connected layers
        self.fc1 = nn.Linear(4 * 4 * 128, 512)
        self.relu4 = nn.ReLU()
        self.fc2 = nn.Linear(512, 10)
        
    def forward(self, x):
        x = self.conv1(x)
        x = self.relu1(x)
        x = self.pool1(x)
        
        x = self.conv2(x)
        x = self.relu2(x)
        x = self.pool2(x)
        
        x = self.conv3(x)
        x = self.relu3(x)
        x = self.pool3(x)
        
        x = x.view(x.size(0), -1)
        
        x = self.fc1(x)
        x = self.relu4(x)
        x = self.fc2(x)
        
        return x

from re import X
### YOUR CODE HERE ###
epochs = 10
batch_size = 10
learning_rate = 0.01
torch.manual_seed(12)

model = CIFAR10CNN()
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
criterion = torch.nn.CrossEntropyLoss()
 
train_loader = torch.utils.data.DataLoader(cifar_train, batch_size = batch_size, shuffle=True)
test_loader = torch.utils.data.DataLoader(cifar_validation, batch_size = batch_size, shuffle=False)

model.train() # Put model in training mode

training_losses = []
validation_losses = []

training_accuracy = []
validation_accuracy = []

for epoch in range(epochs):

    epoch_losses = []
    train_epoch_accuracy = []

    for x, y in tqdm.tqdm_notebook(train_loader, unit = "batch"):
        # flattening data
        #x_train = x.view(batch_size, -1)
        optimizer.zero_grad() # Remove the gradients from the previous step
        pred = model(x)
        loss = criterion(pred, y)

        loss.backward()
        optimizer.step()
        epoch_losses.append(loss.item())

    print("Finished Epoch", epoch + 1, ", training loss:", np.mean(epoch_losses))
    training_losses.append(np.mean(epoch_losses))
    
    with torch.no_grad():
      model.eval() # Put model in eval mode
      num_correct = 0
      for x, y in train_loader:
          #x_train = x.view(batch_size, -1)
          pred = model(x)

          num_correct += torch.sum(np.round(pred) == y).item()
          #num_correct = (torch.argmax(pred) == y).sum().item()
          #num_correct += torch.sum(torch.round(pred) == y).item()
          train_epoch_accuracy.append(num_correct)
      print("Accuracy for Training Epoch", epoch + 1, ":", num_correct / len(cifar_train))
      training_accuracy.append(num_correct / len(cifar_train))

      model.train() # Put model back in train mode

    # run through validation data

    val_epoch_losses = []
    val_epoch_accuracy = []

    for x, y in tqdm.notebook.tqdm(test_loader, unit = "batch"):
      optimizer.zero_grad()
      #x_val = x.view(batch_size, -1)
      pred = model(x)
      loss = criterion(pred, y)

      loss.backward()
      optimizer.step()
      val_epoch_losses.append(loss.item())

    print("Finished Epoch", epoch + 1, ", validation loss:", np.mean(val_epoch_losses))
    validation_losses.append(np.mean(val_epoch_losses))

    with torch.no_grad():
      model.eval() # Put model in eval mode
      num_correct = 0
      for x, y in test_loader:
          #x_train = x.view(batch_size, -1)
          pred = model(x)

          num_correct += torch.sum(torch.round(pred) == y).item()
          #num_correct = (torch.argmax(pred) == y).sum().item()
          #num_correct += torch.sum(torch.round(pred) == y).item()
          val_epoch_accuracy.append(num_correct)
      print("Accuracy for Validation Epoch", epoch + 1, ":", num_correct / len(cifar_validation))
      validation_accuracy.append(num_correct / len(cifar_validation))

      model.train() # Put model back in train mode

plt.plot(training_losses, label = 'training loss')
plt.plot(validation_losses, label = 'validation loss')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Loss of Training and Validation Per Training Epcoh")
plt.legend();

"""### Kaggle Submission
The following code is for you to make your submission to kaggle. Here are the steps you must follow:

1. Upload `cifar_test_data.npy` to the colab notebook by going to files on the right hand pane, then hitting "upload". 
2. Run the following cell to generate the dataset object for the test data. Feel free to modify the code to use the same transforms that you use for the training data. By default, this will re-use the `transform` variable.
3. In the second cell, write code to run predictions on the testing dataset and store them into an array called `predictions`.
4. Run the final cell which will convert your predictions array into a CSV for kaggle.
5. Go to the files pane again, and download the file called `submission.csv` by clicking the three dots and then download.

"""

from PIL import Image
import os

class CIFAR10Test(torchvision.datasets.VisionDataset):
    
    def __init__(self, transform=None, target_transform=None):
        super(CIFAR10Test, self).__init__(None, transform=transform,
                                      target_transform=target_transform)
        assert os.path.exists("cifar_test_data.npy"), "You must upload the test data to the file system."
        self.data = [np.load("cifar_test_data.npy", allow_pickle=False)]

        self.data = np.vstack(self.data).reshape(-1, 3, 32, 32)
        self.data = self.data.transpose((0, 2, 3, 1))  # convert to HWC

    def __getitem__(self, index: int):
        img = self.data[index]
        img = Image.fromarray(img)
        if self.transform is not None:
            img = self.transform(img)
        return img

    def __len__(self) -> int:
        return len(self.data)

# Create the test dataset
testing_data = CIFAR10Test(
    transform=transform, # NOTE: Make sure transform is the same as used in the training dataset.
)

### YOUR CODE HERE ###

# Recommendation: create a `test_dataloader` from torch.utils.data.DataLoader with `shuffle=False` to iterate over the test data in batches.

# Store a numpy vector of the predictions for the test set in the variable `predictions`.
predictions = None

# This code below will generate kaggle_predictions.csv file. Please download it and submit to kaggle.
import pandas as pd

if isinstance(predictions, np.ndarray):
    predictions = predictions.astype(int)
else:
    predictions = np.array(predictions, dtype=int)
assert predictions.shape == (len(testing_data),), "Predictions were not the correct shape"
df = pd.DataFrame({'Category': predictions})
df.index += 1  # Ensures that the index starts at 1. 
df.to_csv('submission.csv', index_label='Id')

# Now download the submission.csv file to submit.

"""Congrats! You made it to the end."""