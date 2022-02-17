from random import gauss
from csv import writer
import sys
from tqdm import tqdm
import os
import csv
import numpy as np



file = open("dataset.csv", "w")
dim = 10000000
seque = 1000
writer = csv.writer(file)
x = np.random.randint(1,seque,dim)
y = np.random.randint(1,seque,dim)
for w in range(dim):
    writer.writerow([x[w], y[w]])

file.close()
