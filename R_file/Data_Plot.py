import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

with open('Tomato_data_det.csv', newline='') as File:
    reader = csv.reader(File)
    for row in reader:
       print(row)
print(reader)

df=pd.read_csv('Tomato_data_det.csv', sep=',',header=None)
