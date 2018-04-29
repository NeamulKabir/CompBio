
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


tumor_names_train = ['real1','real2','syn1','syn2','syn3','syn4','syn5']
tumor_names_val = ['real2']

selected_features = ['freebayes_pred','mutect_pred','vardict_pred','varscan_pred','F_MQ','F_DPRA','F_QD','F_GT_GQ','F_ODDS','F_QA','F_GT_DP','F_HaplotypeScore','F_RPP','F_PAIRED','F_QR','F_SAF','F_AB','F_PAIREDR','F_GT_PL1','F_GT_PL3','F_RPR','M_GT_FREQ','M_GT_AD_ALT','M_GT_BQ','M_GT_AD_REF','VD_GT_MQ','VD_GT_AD_REF','VD_SSF','VD_GT_HIAF','VD_GT_AF','VD_GT_NM','VD_GT_RD2','VD_GT_PMEAN','VD_GT_DP','VD_GT_VD','VD_MSI','VS_SPV','VS_MQ','VS_MQ0','VS_GT_DP','VS_GT_FREQ','VS_GT_DP4_1','VS_HaplotypeScore','VS_BaseQRankSum','VS_MQRankSum']

selected_features_str = 'freebayes_pred\tmutect_pred\tvardict_pred\tvarscan_pred\tF_MQ\tF_DPRA\tF_QD\tF_GT_GQ\tF_ODDS\tF_QA\tF_GT_DP\tF_HaplotypeScore\tF_RPP\tF_PAIRED\tF_QR\tF_SAF\tF_AB\tF_PAIREDR\tF_GT_PL1\tF_GT_PL3\tF_RPR\tM_GT_FREQ\tM_GT_AD_ALT\tM_GT_BQ\tM_GT_AD_REF\tVD_GT_MQ\tVD_GT_AD_REF\tVD_SSF\tVD_GT_HIAF\tVD_GT_AF\tVD_GT_NM\tVD_GT_RD2\tVD_GT_PMEAN\tVD_GT_DP\tVD_GT_VD\tVD_MSI\tVS_SPV\tVS_MQ\tVS_MQ0\tVS_GT_DP\tVS_GT_FREQ\tVS_GT_DP4_1\tVS_HaplotypeScore\tVS_BaseQRankSum\tVS_MQRankSum'

# Big M value to impute missing data
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


print('############## DECISION TREE ##############')
clf = DecisionTreeClassifier(min_samples_leaf=15)
clf = clf.fit(train_data, train_truth)

train_predicted = clf.predict(train_data)
precision_train, recall_train, f1_score_train = precision_recall_fscore(train_truth, train_predicted)

val_predicted = clf.predict(val_data)
precision_val, recall_val, f1_score_val = precision_recall_fscore(val_truth, val_predicted)

# print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_train, recall_train, f1_score_train))
# print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_val, recall_val, f1_score_val))

print('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(precision_train, recall_train, f1_score_train, precision_val, recall_val, f1_score_val))


print('############## RANDOM FOREST ##############')
clf = RandomForestClassifier(n_estimators=20, min_samples_leaf=10)
clf.fit(train_data, train_truth)

train_predicted = clf.predict(train_data)
precision_train, recall_train, f1_score_train = precision_recall_fscore(train_truth, train_predicted)

val_predicted = clf.predict(val_data)
precision_val, recall_val, f1_score_val = precision_recall_fscore(val_truth, val_predicted)

# print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_train, recall_train, f1_score_train))
# print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_val, recall_val, f1_score_val))

print('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(precision_train, recall_train, f1_score_train, precision_val, recall_val, f1_score_val))





