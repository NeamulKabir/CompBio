
import matplotlib.pyplot as plt
import numpy as np

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

def precision_recall_fscore(truth_label, predicted_label):
	TP = np.sum(np.where((truth_label+predicted_label==2),1,0))
	FP = np.sum(np.where((truth_label-predicted_label==-1),1,0))
	FN = np.sum(np.where((truth_label-predicted_label==1),1,0))

	precision = TP/(TP+FP)
	recall = TP/(TP+FN)
	fscore = (2*precision*recall)/(precision+recall)

	return precision, recall, fscore


# tumor_names = ['real1','real2','syn1','syn2','syn3','syn4','syn5']
tumor_names = ['syn1']
# tumor_names = ['syn5']

M = 99999.0

# read datasets of selected tumor samples

firt_tumor_flag = True

for tumor_name in tumor_names:
	print('Data:{}'.format(tumor_name))

	train_data_file = tumor_name + '/' + tumor_name + '_train_data.txt'
	train_data_mask_file = tumor_name + '/' + tumor_name + '_train_data_mask.txt'
	train_truth_file = tumor_name + '/' + tumor_name + '_train_truth.txt'

	val_data_file = tumor_name + '/' + tumor_name + '_val_data.txt'
	val_data_mask_file = tumor_name + '/' + tumor_name + '_val_data_mask.txt'
	val_truth_file = tumor_name + '/' + tumor_name + '_val_truth.txt'

	temp_train_data = np.loadtxt(train_data_file, dtype=float, delimiter='\t', skiprows=1)
	temp_train_data_mask = np.loadtxt(train_data_mask_file, dtype=bool, delimiter='\t', skiprows=1)
	temp_train_truth = np.loadtxt(train_truth_file, dtype=int, delimiter='\t', skiprows=1)

	temp_val_data = np.loadtxt(val_data_file, dtype=float, delimiter='\t', skiprows=1)
	temp_val_data_mask = np.loadtxt(val_data_mask_file, dtype=bool, delimiter='\t', skiprows=1)
	temp_val_truth = np.loadtxt(val_truth_file, dtype=int, delimiter='\t', skiprows=1)

	# print('Shape of temp_train_data:{}'.format(temp_train_data.shape))
	# print('Shape of temp_train_data_mask:{}'.format(temp_train_data_mask.shape))
	# print('Shape of temp_train_truth:{}'.format(temp_train_truth.shape))

	# print('Shape of temp_val_data:{}'.format(temp_val_data.shape))
	# print('Shape of temp_val_data_mask:{}'.format(temp_val_data_mask.shape))
	# print('Shape of temp_val_truth:{}'.format(temp_val_truth.shape))

	temp_train_data[temp_train_data_mask] = M
	temp_val_data[temp_val_data_mask] = M

	# print('First row of train data after imputation:{}'.format(temp_train_data[0,:]))
	# print('First row of validation data after imputation:{}'.format(temp_val_data[0,:]))

	if firt_tumor_flag:
		train_data = temp_train_data
		train_truth = temp_train_truth

		val_data = temp_val_data
		val_truth = temp_val_truth

		firt_tumor_flag = False
	else:
		train_data = np.vstack((train_data,temp_train_data))
		train_truth = np.hstack((train_truth,temp_train_truth))

		val_data = np.vstack((val_data,temp_val_data))
		val_truth = np.hstack((val_truth,temp_val_truth))


print('Shape of train_data:{}'.format(train_data.shape))
print('Shape of train_truth:{}'.format(train_truth.shape))

print('Shape of val_data:{}'.format(val_data.shape))
print('Shape of val_truth:{}'.format(val_truth.shape))

del temp_train_data
del temp_train_data_mask
del temp_train_truth

del temp_val_data
del temp_val_data_mask
del temp_val_truth


print('############## DECISION TREE ##############')
clf = DecisionTreeClassifier()
clf = clf.fit(train_data, train_truth)

train_predicted = clf.predict(train_data)
precision, recall, f1_score = precision_recall_fscore(train_truth, train_predicted)
print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision,recall,f1_score))

val_predicted = clf.predict(val_data)
precision, recall, f1_score = precision_recall_fscore(val_truth, val_predicted)
print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision,recall,f1_score))


print('############## RANDOM FOREST ##############')
clf = RandomForestClassifier(n_estimators=10,min_samples_leaf=10)
clf.fit(train_data, train_truth)

train_predicted = clf.predict(train_data)
precision, recall, f1_score = precision_recall_fscore(train_truth, train_predicted)
print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision,recall,f1_score))

val_predicted = clf.predict(val_data)
precision, recall, f1_score = precision_recall_fscore(val_truth, val_predicted)
print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision,recall,f1_score))
