
import matplotlib.pyplot as plt
import numpy as np

from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

import pickle

test_chrom_filename = 'test_chrom_submission.txt'
test_data_filename = 'test_data_submission.txt'

test_chrom = np.loadtxt(test_chrom_filename, dtype=int, delimiter='\t', skiprows=1)
test_data = np.loadtxt(test_data_filename, dtype=float, delimiter='\t', skiprows=1)

print('Shape of test_data:{}'.format(test_data.shape))


model_filename = 'random_forest_model.sav'
rf = pickle.load(open(model_filename, 'rb'))

test_predicted = rf.predict(test_data)

print('Shape of test_predicted:{}'.format(test_predicted.shape))

ind = np.lexsort((test_chrom[:,1],test_chrom[:,0]))


test_predicted_filename = 'test_predicted_submission.bed'

with open(test_predicted_filename,'w') as f_test_predicted_filename:
	for i in range(len(test_predicted)):
		index = ind[i]
		if test_predicted[index] == 1:

			if test_chrom[index,0] == 23:
				chrom_name = 'X'
			elif test_chrom[index,0] == 24:
				chrom_name = 'Y'
			else:
				chrom_name = str(test_chrom[index,0])

			pos = str(test_chrom[index,1])

			f_test_predicted_filename.write('{}\t{}\t{}\n'.format(chrom_name,pos,pos))


print('Number of positive predictions: {} out of {} candidate positions'.format(np.sum(test_predicted),len(test_predicted)))














