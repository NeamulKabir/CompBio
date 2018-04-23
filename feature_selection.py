
import matplotlib.pyplot as plt
import numpy as np

from sklearn.feature_selection import SelectFromModel
from sklearn.tree import DecisionTreeClassifier


tumor_names = ['real1','real2','syn1','syn2','syn3','syn4','syn5']
# tumor_names = ['syn5']

for tumor_name in tumor_names:
	print('Data:{}'.format(tumor_name))

	feature_names = ['chrom','loc','truth','freebayes_pred','mutect_pred','vardict_pred','varscan_pred','F_GT_DP','F_GT_GQ','F_GT_PL1','F_GT_PL2','F_GT_PL3','F_GT_RO','F_AB','F_ABP','F_AC','F_AF','F_AN','F_AO','F_BaseQRankSum','F_DP','F_DPB','F_DPRA','F_EPP','F_FS','F_GC','F_GTI','F_Hrun','F_HaplotypeScore','F_LEN','F_MEANALT','F_MQ','F_MQ0','F_MQM','F_MQMR','F_MQRankSum','F_NS','F_NUMALT','F_ODDS','F_PAIRED','F_PAIREDR','F_PAO','F_PQA','F_PQR','F_PRO','F_QA','F_QD','F_QR','F_RO','F_RPL','F_RPP','F_RPPR','F_RPR','F_RUN','F_ReadPosRankSum','F_SAF','F_SAP','F_SAR','F_SRF','F_SRP','F_SRR','M_GT_AD_REF','M_GT_AD_ALT','M_GT_BQ','M_GT_DP','M_GT_FREQ','M_GT_SS','VD_GT_AD_REF','VD_GT_AD_ALT','VD_GT_ADJAF','VD_GT_AF','VD_GT_ALD_FOR','VD_GT_ALD_REV','VD_GT_BIAS_REF','VD_GT_BIAS_ALT','VD_GT_DP','VD_GT_HIAF','VD_GT_MQ','VD_GT_NM','VD_GT_ODDRATIO','VD_GT_PMEAN','VD_GT_PSTD','VD_GT_QSTD','VD_GT_QUAL','VD_GT_RD1','VD_GT_RD2','VD_GT_SBF','VD_GT_SN','VD_GT_VD','VD_MSI','VD_MSILEN','VD_SOR','VD_SSF','VS_GT_DP','VS_GT_DP4_1','VS_GT_DP4_2','VS_GT_DP4_3','VS_GT_DP4_4','VS_GT_FREQ','VS_AC','VS_AF','VS_AN','VS_BaseQRankSum','VS_DP','VS_FS','VS_GC','VS_GPV','VS_Hrun','VS_HaplotypeScore','VS_MQ','VS_MQ0','VS_MQRankSum','VS_ReadPosRankSum','VS_SPV','VS_SS','VS_SSC']

	methods = ['freebayes', 'mutect', 'vardict', 'varscan']
	# methods = ['freebayes']

	ind_offset = 7
	ind_dict = {'freebayes':np.array([0,54]),
				'mutect':np.array([54,60]),
				'vardict':np.array([60,86]),
				'varscan':np.array([86,109])
	}

	pred_dict = {'freebayes':3,
				'mutect':4,
				'vardict':5,
				'varscan':6
	}


	dataset_file = tumor_name + '/' + tumor_name + '_dataset.txt'
	mask_file = tumor_name + '/' + tumor_name + '_dataset_mask.txt'

	dataset = np.loadtxt(dataset_file, dtype=float, delimiter='\t', skiprows=1)
	mask = np.loadtxt(mask_file, dtype=np.int, delimiter='\t', skiprows=1)

	# print('Shape of dataset:{}'.format(dataset.shape))
	# print('Shape of mask:{}'.format(mask.shape))

	for method in methods:
		temp_ind = np.where(dataset[:,pred_dict[method]] == 1)[0]

		temp_features = feature_names[ind_offset + ind_dict[method][0]:ind_offset + ind_dict[method][1]]
		# print('Features of "{}" method:{}'.format(method,temp_features))

		temp_data = dataset[temp_ind,ind_offset + ind_dict[method][0]:ind_offset + ind_dict[method][1]]
		temp_mask = mask[temp_ind,ind_dict[method][0]:ind_dict[method][1]]
		temp_truth = dataset[temp_ind,2]

		# print(temp_data.shape)
		# print(temp_mask.shape)
		# print(temp_truth.shape)

		mask_row_sum = np.sum(temp_mask,axis=1)
		keep_data_rows = np.where(mask_row_sum == 0)[0]

		filtered_data = temp_data[keep_data_rows,:]
		filtered_truth = temp_truth[keep_data_rows]

		# print(filtered_data.shape)
		# print(filtered_truth.shape)

		clf = DecisionTreeClassifier(max_depth=5)
		clf = clf.fit(filtered_data, filtered_truth)
		importances = np.array(clf.feature_importances_)
		# print(importances)
		sorted_ind = (importances.argsort())[::-1]
		# print(sorted_ind.shape)
		sorted_importances = importances[sorted_ind]
		sorted_temp_features = np.array(temp_features)[sorted_ind]
		print('Features of "{}" method:{}'.format(method,sorted_temp_features))
		print('Feature Importances of "{}" method:{}'.format(method, sorted_importances)) 

		# model = SelectFromModel(clf, prefit=True)
		# selected_data = model.transform(filtered_data)
		# print(selected_data.shape)


		num_selected_features = len(np.where(importances >= importances.mean())[0])
		# print(selected_ind)
		print('Selected Features of "{}" method:{}'.format(method,sorted_temp_features[:num_selected_features]))
		print('Selected Features Importances of "{}" method:{}\n'.format(method, sorted_importances[:num_selected_features]))

		fig = plt.figure(1,figsize=(12, 7.5), dpi=100)
		plt.bar(np.arange(len(sorted_importances)),sorted_importances)
		plt.xticks(np.arange(len(sorted_importances)), sorted_temp_features, rotation='vertical')
		plt.ylabel('Importance Score')
		plt.title('Feature Importance Scores for "{}/{}"'.format(tumor_name,method))
		plt.grid()
		plt.tight_layout()
		# plt.show()
		fig_path = tumor_name + '/' + tumor_name + '_feature_importance_scores_' + method + '.png'
		fig.savefig(fig_path)

		fig.clf()



