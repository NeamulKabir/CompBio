import numpy as np


# DATASET_TYPE = 'TRAIN'
DATASET_TYPE = 'TEST'

print('DATASET_TYPE:{}'.format(DATASET_TYPE))

tumor_names_train = ['real1','real2','syn1','syn2','syn3','syn4','syn5']
tumor_names_val = ['real2']

selected_features = ['freebayes_pred','mutect_pred','vardict_pred','varscan_pred','F_MQ','F_DPRA','F_QD','F_GT_GQ','F_ODDS','F_QA','F_GT_DP','F_HaplotypeScore','F_RPP','F_PAIRED','F_QR','F_SAF','F_AB','F_PAIREDR','F_GT_PL1','F_GT_PL3','F_RPR','M_GT_FREQ','M_GT_AD_ALT','M_GT_BQ','M_GT_AD_REF','VD_GT_MQ','VD_GT_AD_REF','VD_SSF','VD_GT_HIAF','VD_GT_AF','VD_GT_NM','VD_GT_RD2','VD_GT_PMEAN','VD_GT_DP','VD_GT_VD','VD_MSI','VS_SPV','VS_MQ','VS_MQ0','VS_GT_DP','VS_GT_FREQ','VS_GT_DP4_1','VS_HaplotypeScore','VS_BaseQRankSum','VS_MQRankSum']

selected_features_str = 'freebayes_pred\tmutect_pred\tvardict_pred\tvarscan_pred\tF_MQ\tF_DPRA\tF_QD\tF_GT_GQ\tF_ODDS\tF_QA\tF_GT_DP\tF_HaplotypeScore\tF_RPP\tF_PAIRED\tF_QR\tF_SAF\tF_AB\tF_PAIREDR\tF_GT_PL1\tF_GT_PL3\tF_RPR\tM_GT_FREQ\tM_GT_AD_ALT\tM_GT_BQ\tM_GT_AD_REF\tVD_GT_MQ\tVD_GT_AD_REF\tVD_SSF\tVD_GT_HIAF\tVD_GT_AF\tVD_GT_NM\tVD_GT_RD2\tVD_GT_PMEAN\tVD_GT_DP\tVD_GT_VD\tVD_MSI\tVS_SPV\tVS_MQ\tVS_MQ0\tVS_GT_DP\tVS_GT_FREQ\tVS_GT_DP4_1\tVS_HaplotypeScore\tVS_BaseQRankSum\tVS_MQRankSum'

# Big M value to impute missing data
M = 99999.0

if DATASET_TYPE == 'TRAIN':

	# read datasets of selected tumor samples
	first_tumor_flag = True
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
		
		if first_tumor_flag:
			train_data = temp_train_data
			train_truth = temp_train_truth

			first_tumor_flag = False
		else:
			train_data = np.vstack((train_data,temp_train_data))
			train_truth = np.hstack((train_truth,temp_train_truth))

	first_tumor_flag = True
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

		if first_tumor_flag:
			val_data = temp_val_data
			val_truth = temp_val_truth

			first_tumor_flag = False
		else:
			val_data = np.vstack((val_data,temp_val_data))
			val_truth = np.hstack((val_truth,temp_val_truth))

	train_data_filename = 'train_data_submission.txt'
	train_truth_filename = 'train_truth_submission.txt'
	np.savetxt(train_data_filename, train_data, fmt='%.3f', delimiter='\t', header=selected_features_str, comments='#')
	np.savetxt(train_truth_filename, train_truth, fmt='%d', delimiter='\t', header='truth', comments='#')

	val_data_filename = 'val_data_submission.txt'
	val_truth_filename = 'val_truth_submission.txt'
	np.savetxt(val_data_filename, val_data, fmt='%.3f', delimiter='\t', header=selected_features_str, comments='#')
	np.savetxt(val_truth_filename, val_truth, fmt='%d', delimiter='\t', header='truth', comments='#')

elif DATASET_TYPE == 'TEST':

	test_chrom_file = 'test/test_chrom.txt'
	test_data_file = 'test/test_data.txt'
	test_data_mask_file = 'test/test_data_mask.txt'
	
	test_chrom = np.loadtxt(test_chrom_file, dtype=int, delimiter='\t', skiprows=1)
	test_data = np.loadtxt(test_data_file, dtype=float, delimiter='\t', skiprows=1)
	test_data_mask = np.loadtxt(test_data_mask_file, dtype=bool, delimiter='\t', skiprows=1)

	test_data[test_data_mask] = M

	test_chrom_filename = 'test_chrom_submission.txt'
	test_data_filename = 'test_data_submission.txt'

	np.savetxt(test_chrom_filename, test_chrom, fmt='%d', delimiter='\t', header='chrom\tloc', comments='#')
	np.savetxt(test_data_filename, test_data, fmt='%.3f', delimiter='\t', header=selected_features_str, comments='#')








