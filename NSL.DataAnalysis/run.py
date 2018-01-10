import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
if __name__ == '__main__':
	# xl = pd.ExcelFile('./Dataset/a.csv')
	df = {}
	dataset = pd.read_csv('./Dataset/a.csv')
	df['RXE']  = dataset[dataset.LogRecType=='RXE']
	df['STE'] = dataset[dataset.LogRecType == 'STE']
	df['TXE']= dataset[dataset.LogRecType == 'TXE']
	for key in df.keys() :
		df[key] = df[key].dropna(axis  = 1 ,how='all')
		print(df[key].keys())
	print('end')
	df[key]