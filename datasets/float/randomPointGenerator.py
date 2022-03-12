from random import gauss
from csv import writer
import sys
from tqdm import tqdm
import os
import csv
import numpy as np		
		
file = open("dataset.csv", "w")
dim = int(input("dataset size: "))
writer = csv.writer(file)
for w in range(dim):
    writer.writerow([round(np.random.uniform(low = 0, high = 1), 3), round(np.random.uniform(low = 0, high = 1), 3)])

file.close()
