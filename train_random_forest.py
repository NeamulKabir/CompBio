
import matplotlib.pyplot as plt
import numpy as np

from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

from subprocess import call

import pickle



def precision_recall_fscore(truth_label, predicted_label):
	TP = np.sum(np.where((truth_label+predicted_label==2),1,0))
	FP = np.sum(np.where((truth_label-predicted_label==-1),1,0))
	FN = np.sum(np.where((truth_label-predicted_label==1),1,0))

	precision = TP/(TP+FP)
	recall = TP/(TP+FN)
	fscore = (2*precision*recall)/(precision+recall)

	return precision, recall, fscore

selected_features = ['freebayes_pred','mutect_pred','vardict_pred','varscan_pred','F_MQ','F_DPRA','F_QD','F_GT_GQ','F_ODDS','F_QA','F_GT_DP','F_HaplotypeScore','F_RPP','F_PAIRED','F_QR','F_SAF','F_AB','F_PAIREDR','F_GT_PL1','F_GT_PL3','F_RPR','M_GT_FREQ','M_GT_AD_ALT','M_GT_BQ','M_GT_AD_REF','VD_GT_MQ','VD_GT_AD_REF','VD_SSF','VD_GT_HIAF','VD_GT_AF','VD_GT_NM','VD_GT_RD2','VD_GT_PMEAN','VD_GT_DP','VD_GT_VD','VD_MSI','VS_SPV','VS_MQ','VS_MQ0','VS_GT_DP','VS_GT_FREQ','VS_GT_DP4_1','VS_HaplotypeScore','VS_BaseQRankSum','VS_MQRankSum']

selected_features_str = 'freebayes_pred\tmutect_pred\tvardict_pred\tvarscan_pred\tF_MQ\tF_DPRA\tF_QD\tF_GT_GQ\tF_ODDS\tF_QA\tF_GT_DP\tF_HaplotypeScore\tF_RPP\tF_PAIRED\tF_QR\tF_SAF\tF_AB\tF_PAIREDR\tF_GT_PL1\tF_GT_PL3\tF_RPR\tM_GT_FREQ\tM_GT_AD_ALT\tM_GT_BQ\tM_GT_AD_REF\tVD_GT_MQ\tVD_GT_AD_REF\tVD_SSF\tVD_GT_HIAF\tVD_GT_AF\tVD_GT_NM\tVD_GT_RD2\tVD_GT_PMEAN\tVD_GT_DP\tVD_GT_VD\tVD_MSI\tVS_SPV\tVS_MQ\tVS_MQ0\tVS_GT_DP\tVS_GT_FREQ\tVS_GT_DP4_1\tVS_HaplotypeScore\tVS_BaseQRankSum\tVS_MQRankSum'

train_data_filename = 'train_data_submission.txt'
train_truth_filename = 'train_truth_submission.txt'
train_data = np.loadtxt(train_data_filename, dtype=float, delimiter='\t', skiprows=1)
train_truth = np.loadtxt(train_truth_filename, dtype=int, delimiter='\t', skiprows=1)

val_data_filename = 'val_data_submission.txt'
val_truth_filename = 'val_truth_submission.txt'
val_data = np.loadtxt(val_data_filename, dtype=float, delimiter='\t', skiprows=1)
val_truth = np.loadtxt(val_truth_filename, dtype=int, delimiter='\t', skiprows=1)

print('Shape of train_data:{}'.format(train_data.shape))
print('Shape of train_truth:{}'.format(train_truth.shape))

print('Shape of val_data:{}'.format(val_data.shape))
print('Shape of val_truth:{}'.format(val_truth.shape))


print('############## DECISION TREE ##############')
dt = DecisionTreeClassifier(min_samples_leaf=15)
dt = dt.fit(train_data, train_truth)

train_predicted = dt.predict(train_data)
precision_train, recall_train, f1_score_train = precision_recall_fscore(train_truth, train_predicted)

val_predicted = dt.predict(val_data)
precision_val, recall_val, f1_score_val = precision_recall_fscore(val_truth, val_predicted)

# print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_train, recall_train, f1_score_train))
# print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_val, recall_val, f1_score_val))
print('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(precision_train, recall_train, f1_score_train, precision_val, recall_val, f1_score_val))

importances = np.array(dt.feature_importances_)
sorted_ind = (importances.argsort())[::-1]
sorted_importances = importances[sorted_ind]
sorted_temp_features = np.array(selected_features)[sorted_ind]

print('Feature\tImportance')
for i in range(train_data.shape[1]):
	print('{}\t{:.3f}'.format(sorted_temp_features[i],sorted_importances[i]))

# tree.export_graphviz(dt, out_file='trained_dt.dot', feature_names=selected_features)
# call(["dot", "-Tpng", "trained_dt.dot", "-o", "trained_dt.png"])

# Plot the feature importances of the forest
fig1 = plt.figure(1,figsize=(12, 7.5), dpi=100)
plt.title("Decision Tree - Feature importances")
# plt.bar(range(train_data.shape[1]), sorted_importances, color="r", yerr=std[sorted_ind], align="center")
plt.bar(range(train_data.shape[1]), sorted_importances, color="r", align="center")
plt.xticks(np.arange(len(sorted_importances)), sorted_temp_features, rotation='vertical')
plt.xlim([-1, train_data.shape[1]])
plt.grid()
plt.tight_layout()

fig_path = 'decision_tree_feature_importance_scores.png'
fig1.savefig(fig_path)


print('############## RANDOM FOREST ##############')
rf = RandomForestClassifier(n_estimators=20, min_samples_leaf=10)
rf.fit(train_data, train_truth)

train_predicted = rf.predict(train_data)
precision_train, recall_train, f1_score_train = precision_recall_fscore(train_truth, train_predicted)

val_predicted = rf.predict(val_data)
precision_val, recall_val, f1_score_val = precision_recall_fscore(val_truth, val_predicted)

# print('Training set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_train, recall_train, f1_score_train))
# print('Validation set - precision:{:.3f}, recall:{:.3f}, f1_score:{:.3f}'.format(precision_val, recall_val, f1_score_val))
print('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(precision_train, recall_train, f1_score_train, precision_val, recall_val, f1_score_val))

# save random forest model
model_filename = 'random_forest_model.sav'
pickle.dump(rf, open(model_filename, 'wb'))

importances = np.array(rf.feature_importances_)
sorted_ind = (importances.argsort())[::-1]
sorted_importances = importances[sorted_ind]
sorted_temp_features = np.array(selected_features)[sorted_ind]

print('Feature\tImportance')
for i in range(train_data.shape[1]):
	print('{}\t{:.3f}'.format(sorted_temp_features[i],sorted_importances[i]))

# std = np.std([tr.feature_importances_ for tr in rf.estimators_], axis=0)

# count=0
# for tr in rf.estimators_:
# 	tree.export_graphviz(tr, out_file='random_forest_tree_{}.dot'.format(count), feature_names=selected_features)
# 	call(["dot", "-Tpng", "random_forest_tree_{}.dot".format(count), "-o", "random_forest_tree_{}.png".format(count)])

# 	count += 1


# Plot the feature importances of the forest
fig2 = plt.figure(2,figsize=(12, 7.5), dpi=100)
plt.title("Random Forest - Feature importances")
# plt.bar(range(train_data.shape[1]), sorted_importances, color="r", yerr=std[sorted_ind], align="center")
plt.bar(range(train_data.shape[1]), sorted_importances, color="r", align="center")
plt.xticks(np.arange(len(sorted_importances)), sorted_temp_features, rotation='vertical')
plt.xlim([-1, train_data.shape[1]])
plt.grid()
plt.tight_layout()

fig_path = 'random_forest_feature_importance_scores.png'
fig2.savefig(fig_path)

plt.show()



