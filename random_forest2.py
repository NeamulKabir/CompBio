
import matplotlib.pyplot as plt
import numpy as np

from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

from subprocess import call

def precision_recall_fscore(truth_label, predicted_label):
	TP = np.sum(np.where((truth_label+predicted_label==2),1,0))
	FP = np.sum(np.where((truth_label-predicted_label==-1),1,0))
	FN = np.sum(np.where((truth_label-predicted_label==1),1,0))

	precision = TP/(TP+FP)
	recall = TP/(TP+FN)
	fscore = (2*precision*recall)/(precision+recall)

	return precision, recall, fscore


# tumor_names = ['real1','real2','syn1','syn2','syn3','syn4','syn5']

tumor_names_train = ['real1','syn1','syn2','syn3','syn4','syn5']
tumor_names_val = ['real2']

selected_features = ['freebayes_pred','mutect_pred','vardict_pred','varscan_pred','F_MQ','F_DPRA','F_QD','F_GT_GQ','F_ODDS','F_QA','F_GT_DP','F_HaplotypeScore','F_RPP','F_PAIRED','F_QR','F_SAF','F_AB','F_PAIREDR','F_GT_PL1','F_GT_PL3','F_RPR','M_GT_FREQ','M_GT_AD_ALT','M_GT_BQ','M_GT_AD_REF','VD_GT_MQ','VD_GT_AD_REF','VD_SSF','VD_GT_HIAF','VD_GT_AF','VD_GT_NM','VD_GT_RD2','VD_GT_PMEAN','VD_GT_DP','VD_GT_VD','VD_MSI','VS_SPV','VS_MQ','VS_MQ0','VS_GT_DP','VS_GT_FREQ','VS_GT_DP4_1','VS_HaplotypeScore','VS_BaseQRankSum','VS_MQRankSum']

M = 99999.0

# read datasets of selected tumor samples

firt_tumor_flag = True
for tumor_name in tumor_names_train:
	print('Training Data:{}'.format(tumor_name))

	train_data_file = tumor_name + '/' + tumor_name + '_train_data.txt'
	train_data_mask_file = tumor_name + '/' + tumor_name + '_train_data_mask.txt'
	train_truth_file = tumor_name + '/' + tumor_name + '_train_truth.txt'
	
	temp_train_data = np.loadtxt(train_data_file, dtype=float, delimiter='\t', skiprows=1)
	temp_train_data_mask = np.loadtxt(train_data_mask_file, dtype=bool, delimiter='\t', skiprows=1)
	temp_train_truth = np.loadtxt(train_truth_file, dtype=int, delimiter='\t', skiprows=1)

	# print('Shape of temp_train_data:{}'.format(temp_train_data.shape))
	# print('Shape of temp_train_data_mask:{}'.format(temp_train_data_mask.shape))
	# print('Shape of temp_train_truth:{}'.format(temp_train_truth.shape))
	
	temp_train_data[temp_train_data_mask] = M

	# print('First row of train data after imputation:{}'.format(temp_train_data[0,:]))
	
	if firt_tumor_flag:
		train_data = temp_train_data
		train_truth = temp_train_truth

		firt_tumor_flag = False
	else:
		train_data = np.vstack((train_data,temp_train_data))
		train_truth = np.hstack((train_truth,temp_train_truth))

firt_tumor_flag = True
for tumor_name in tumor_names_val:
	print('Validation Data:{}'.format(tumor_name))

	val_data_file = tumor_name + '/' + tumor_name + '_val_data.txt'
	val_data_mask_file = tumor_name + '/' + tumor_name + '_val_data_mask.txt'
	val_truth_file = tumor_name + '/' + tumor_name + '_val_truth.txt'

	temp_val_data = np.loadtxt(val_data_file, dtype=float, delimiter='\t', skiprows=1)
	temp_val_data_mask = np.loadtxt(val_data_mask_file, dtype=bool, delimiter='\t', skiprows=1)
	temp_val_truth = np.loadtxt(val_truth_file, dtype=int, delimiter='\t', skiprows=1)

	# print('Shape of temp_val_data:{}'.format(temp_val_data.shape))
	# print('Shape of temp_val_data_mask:{}'.format(temp_val_data_mask.shape))
	# print('Shape of temp_val_truth:{}'.format(temp_val_truth.shape))

	temp_val_data[temp_val_data_mask] = M

	# print('First row of validation data after imputation:{}'.format(temp_val_data[0,:]))

	if firt_tumor_flag:
		val_data = temp_val_data
		val_truth = temp_val_truth

		firt_tumor_flag = False
	else:
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


# print('############## DECISION TREE ##############')
clf = DecisionTreeClassifier(min_samples_leaf=15)
clf = clf.fit(train_data, train_truth)

train_predicted = clf.predict(train_data)
precision_train, recall_train, f1_score_train = precision_recall_fscore(train_truth, train_predicted)

val_predicted = clf.predict(val_data)
precision_val, recall_val, f1_score_val = precision_recall_fscore(val_truth, val_predicted)

# print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_train, recall_train, f1_score_train))
# print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_val, recall_val, f1_score_val))

print('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(precision_train, recall_train, f1_score_train, precision_val, recall_val, f1_score_val))

importances = np.array(clf.feature_importances_)
# print(importances)
sorted_ind = (importances.argsort())[::-1]
# print(sorted_ind.shape)
sorted_importances = importances[sorted_ind]
sorted_temp_features = np.array(selected_features)[sorted_ind]
# print('Features:{}'.format(sorted_temp_features))
# print('Feature Importances:{}'.format(sorted_importances)) 

tree.export_graphviz(clf, out_file='train_{}_val_{}_dt.dot'.format(tumor_names_train[0],tumor_names_val[0]), feature_names=selected_features) 

call(["dot", "-Tpng", "train_{}_val_{}_dt.dot".format(tumor_names_train[0],tumor_names_val[0]), "-o", "train_{}_val_{}_dt.png".format(tumor_names_train[0],tumor_names_val[0]) ])

# print('############## RANDOM FOREST ##############')
clf = RandomForestClassifier(n_estimators=20, min_samples_leaf=10)
clf.fit(train_data, train_truth)

train_predicted = clf.predict(train_data)
precision_train, recall_train, f1_score_train = precision_recall_fscore(train_truth, train_predicted)

val_predicted = clf.predict(val_data)
precision_val, recall_val, f1_score_val = precision_recall_fscore(val_truth, val_predicted)

# print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_train, recall_train, f1_score_train))
# print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_val, recall_val, f1_score_val))

print('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(precision_train, recall_train, f1_score_train, precision_val, recall_val, f1_score_val))
