
import numpy as np

from sklearn.feature_selection import SelectFromModel
from sklearn.tree import DecisionTreeClassifier

feature_names = ['chrom','loc','truth','freebayes_pred','mutect_pred','vardict_pred','varscan_pred','F_GT_DP','F_GT_GQ','F_GT_PL1','F_GT_PL2','F_GT_PL3','F_GT_RO','F_AB','F_ABP','F_AC','F_AF','F_AN','F_AO','F_BaseQRankSum','F_DP','F_DPB','F_DPRA','F_EPP','F_FS','F_GC','F_GTI','F_Hrun','F_HaplotypeScore','F_LEN','F_MEANALT','F_MQ','F_MQ0','F_MQM','F_MQMR','F_MQRankSum','F_NS','F_NUMALT','F_ODDS','F_PAIRED','F_PAIREDR','F_PAO','F_PQA','F_PQR','F_PRO','F_QA','F_QD','F_QR','F_RO','F_RPL','F_RPP','F_RPPR','F_RPR','F_RUN','F_ReadPosRankSum','F_SAF','F_SAP','F_SAR','F_SRF','F_SRP','F_SRR','M_GT_AD_REF','M_GT_AD_ALT','M_GT_BQ','M_GT_DP','M_GT_FREQ','M_GT_SS','VD_GT_AD_REF','VD_GT_AD_ALT','VD_GT_ADJAF','VD_GT_AF','VD_GT_ALD_FOR','VD_GT_ALD_REV','VD_GT_BIAS_REF','VD_GT_BIAS_ALT','VD_GT_DP','VD_GT_HIAF','VD_GT_MQ','VD_GT_NM','VD_GT_ODDRATIO','VD_GT_PMEAN','VD_GT_PSTD','VD_GT_QSTD','VD_GT_QUAL','VD_GT_RD1','VD_GT_RD2','VD_GT_SBF','VD_GT_SN','VD_GT_VD','VD_MSI','VD_MSILEN','VD_SOR','VD_SSF','VS_GT_DP','VS_GT_DP4_1','VS_GT_DP4_2','VS_GT_DP4_3','VS_GT_DP4_4','VS_GT_FREQ','VS_AC','VS_AF','VS_AN','VS_BaseQRankSum','VS_DP','VS_FS','VS_GC','VS_GPV','VS_Hrun','VS_HaplotypeScore','VS_MQ','VS_MQ0','VS_MQRankSum','VS_ReadPosRankSum','VS_SPV','VS_SS','VS_SSC']

selected_features = 'freebayes_pred\tmutect_pred\tvardict_pred\tvarscan_pred\tF_MQ\tF_DPRA\tF_QD\tF_GT_GQ\tF_ODDS\tF_QA\tF_GT_DP\tF_HaplotypeScore\tF_RPP\tF_PAIRED\tF_QR\tF_SAF\tF_AB\tF_PAIREDR\tF_GT_PL1\tF_GT_PL3\tF_RPR\tM_GT_FREQ\tM_GT_AD_ALT\tM_GT_BQ\tM_GT_AD_REF\tVD_GT_MQ\tVD_GT_AD_REF\tVD_SSF\tVD_GT_HIAF\tVD_GT_AF\tVD_GT_NM\tVD_GT_RD2\tVD_GT_PMEAN\tVD_GT_DP\tVD_GT_VD\tVD_MSI\tVS_SPV\tVS_MQ\tVS_MQ0\tVS_GT_DP\tVS_GT_FREQ\tVS_GT_DP4_1\tVS_HaplotypeScore\tVS_BaseQRankSum\tVS_MQRankSum'

feature_col_dict = {'freebayes':3,
			'mutect':4,
			'vardict':5,
			'varscan':6
}

mask_col_dict = {'F_BaseQRankSum':5,
				'F_ODDS':9,
				'VS_BaseQRankSum':36
}

col_features = [0,1,2,3,4,5,6,7,8,9,11,13,22,28,31,38,39,40,45,46,47,50,52,55,61,62,63,65,67,70,75,76,77,78,80,85,88,89,92,93,94,98,102,108,109,110,111,113]
col_mask = [0,1,2,4,6,12,15,21,24,31,32,33,38,39,40,43,45,48,54,55,56,58,60,63,68,69,70,71,73,78,81,82,85,86,87,91,95,101,102,103,104,106]

tumor_names = ['real1','real2','syn1','syn2','syn3','syn4','syn5']
# tumor_names = ['real1','real2']

for tumor_name in tumor_names:
	print('Data:{}'.format(tumor_name))

	dataset_file = tumor_name + '/' + tumor_name + '_dataset.txt'
	mask_file = tumor_name + '/' + tumor_name + '_dataset_mask.txt'

	dataset = np.loadtxt(dataset_file, dtype=float, delimiter='\t', skiprows=1, usecols=col_features)
	print('Shape of dataset:{}'.format(dataset.shape))

	mask = np.loadtxt(mask_file, dtype=np.int, delimiter='\t', skiprows=1, usecols=col_mask)
	print('Shape of mask:{}'.format(mask.shape))
	
	# get F_BaseQRankSum mask column
	F_BaseQRankSum_mask = mask[:,mask_col_dict['F_BaseQRankSum']]
	# get F_ODDS mask column
	F_ODDS_mask = mask[:,mask_col_dict['F_ODDS']]
	# get VS_BaseQRankSum mask column
	VS_BaseQRankSum_mask = mask[:,mask_col_dict['VS_BaseQRankSum']]

	mask = np.delete(mask, mask_col_dict['F_BaseQRankSum'], axis=1)
	print('Shape of mask:{}'.format(mask.shape))

	F_ind = np.where(dataset[:,feature_col_dict['freebayes']] == 1)[0]
	print('Shape of F_ind:{}'.format(F_ind.shape))
	F_BaseQRankSum_missing_ind = np.where(F_BaseQRankSum_mask[F_ind] == 1)[0]
	print('Shape of F_BaseQRankSum_missing_ind:{}'.format(F_BaseQRankSum_missing_ind.shape))
	F_ODDS_missing_ind = np.where(F_ODDS_mask[F_ind] == 1)[0]
	print('Shape of F_ODDS_missing_ind:{}'.format(F_ODDS_missing_ind.shape))

	VS_ind = np.where(dataset[:,feature_col_dict['varscan']] == 1)[0]
	print('Shape of VS_ind:{}'.format(VS_ind.shape))
	VS_BaseQRankSum_missing_ind = np.where(VS_BaseQRankSum_mask[VS_ind] == 1)[0]
	print('Shape of VS_BaseQRankSum_missing_ind:{}'.format(VS_BaseQRankSum_missing_ind.shape))

	missing_ind = np.array(list(set(F_BaseQRankSum_missing_ind) | set(F_ODDS_missing_ind) | set(VS_BaseQRankSum_missing_ind)))
	print('Shape of missing_ind:{}'.format(missing_ind.shape))


	FPs_filename = tumor_name + '/' + tumor_name + '_FPs_detected_by_using_missing_data.txt'
	np.savetxt(FPs_filename, dataset[missing_ind,0:3], fmt='%d', delimiter='\t', header='chrom\tloc\ttruth', comments='#')
	

	# Remove missing data rows
	clean_dataset = np.delete(dataset, missing_ind, axis=0)
	clean_mask = np.hstack((np.zeros((len(clean_dataset),4),dtype=np.int),np.delete(mask, missing_ind, axis=0)))
	print('Shape of clean_dataset:{}'.format(clean_dataset.shape))
	print('Shape of clean_mask:{}'.format(clean_mask.shape))

	# Do not remove missing data rows
	# clean_dataset = dataset
	# clean_mask = np.hstack((np.zeros((len(clean_dataset),4),dtype=np.int),mask))
	# print('Shape of clean_dataset:{}'.format(clean_dataset.shape))
	# print('Shape of clean_mask:{}'.format(clean_mask.shape))

	positive_ind = np.where(clean_dataset[:,2] == 1)[0]
	negative_ind = np.where(clean_dataset[:,2] == 0)[0]

	np.random.shuffle(positive_ind)
	np.random.shuffle(negative_ind)

	num_train_pos = int(len(positive_ind)*0.8)
	num_train_neg = int(len(negative_ind)*0.8)
	
	train_pos_ind = positive_ind[:num_train_pos]
	train_neg_ind = negative_ind[:num_train_neg]
	print('Number of positive samples in training dataset:{}'.format(len(train_pos_ind)))
	print('Number of negative samples in training dataset:{}'.format(len(train_neg_ind)))

	val_pos_ind = positive_ind[num_train_pos:]
	val_neg_ind = negative_ind[num_train_neg:]
	print('Number of positive samples in validation dataset:{}'.format(len(val_pos_ind)))
	print('Number of negative samples in validation dataset:{}'.format(len(val_neg_ind)))

	train_ind = np.concatenate((train_pos_ind,train_neg_ind),axis=0)
	val_ind = np.concatenate((val_pos_ind,val_neg_ind),axis=0)
	print('Number of samples in training dataset:{}'.format(len(train_ind)))
	print('Number of samples in validation dataset:{}'.format(len(val_ind)))
	
	np.random.shuffle(train_ind)
	np.random.shuffle(val_ind)


	train_chrom_filename = tumor_name + '/' + tumor_name + '_train_chrom_cleaned.txt'
	train_truth_filename = tumor_name + '/' + tumor_name + '_train_truth_cleaned.txt'
	train_data_filename = tumor_name + '/' + tumor_name + '_train_data_cleaned.txt'
	train_data_mask_filename = tumor_name + '/' + tumor_name + '_train_data_mask_cleaned.txt'

	val_chrom_filename = tumor_name + '/' + tumor_name + '_val_chrom_cleaned.txt'
	val_truth_filename = tumor_name + '/' + tumor_name + '_val_truth_cleaned.txt'
	val_data_filename = tumor_name + '/' + tumor_name + '_val_data_cleaned.txt'
	val_data_mask_filename = tumor_name + '/' + tumor_name + '_val_data_mask_cleaned.txt'

	np.savetxt(train_chrom_filename, clean_dataset[train_ind,0:2], fmt='%d', delimiter='\t', header='chrom\tloc', comments='#')
	np.savetxt(train_truth_filename, clean_dataset[train_ind,2], fmt='%d', delimiter='\t', header='truth', comments='#')
	np.savetxt(train_data_filename, clean_dataset[train_ind,3:], fmt='%.3f', delimiter='\t', header=selected_features, comments='#')
	np.savetxt(train_data_mask_filename, clean_mask[train_ind,:], fmt='%d', delimiter='\t', header=selected_features, comments='#')

	np.savetxt(val_chrom_filename, clean_dataset[val_ind,0:2], fmt='%d', delimiter='\t', header='chrom\tloc', comments='#')
	np.savetxt(val_truth_filename, clean_dataset[val_ind,2], fmt='%d', delimiter='\t', header='truth', comments='#')
	np.savetxt(val_data_filename, clean_dataset[val_ind,3:], fmt='%.3f', delimiter='\t', header=selected_features, comments='#')
	np.savetxt(val_data_mask_filename, clean_mask[val_ind,:], fmt='%d', delimiter='\t', header=selected_features, comments='#')
	



