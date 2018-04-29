
# import matplotlib.pyplot as plt
import numpy as np

feature_names = ['GT_AD_REF,','GT_AD_ALT','GT_ADJAF','GT_AF','GT_ALD_FOR','GT_ALD_REV','GT_BIAS_REF','GT_BIAS_ALT','GT_BQ','GT_DP','GT_DP4_1','GT_DP4_2','GT_DP4_3','GT_DP4_4','GT_FREQ','GT_GQ','GT_HIAF','GT_MQ','GT_NM','GT_ODDRATIO','GT_PL1','GT_PL2','GT_PL3','GT_PMEAN','GT_PSTD','GT_QSTD','GT_QUAL','GT_RD1','GT_RD2','GT_RO','GT_SBF','GT_SN','GT_SS','GT_VD','QUAL','AB','ABP','AC','AF','AN','AO','BaseQRankSum','DP','DPB','DPRA','EPP','EPPR','FS','GC','GPV','GTI','Hrun','HaplotypeScore','LEN','MEANALT','MQ','MQ0','MQM','MQMR','MQRankSum','MSI','MSILEN','NS','NUMALT','ODDS','PAIRED','PAIREDR','PAO','PQA','PQR','PRO','QA','QD','QR','RO','RPL','RPP','RPPR','RPR','RUN','ReadPosRankSum','SAF','SAP','SAR','SOR','SPV','SRF','SRP','SRR','SS','SSC','SSF','CHROM','POS']

tumor_names = ['test']
# tumor_names = ['syn5']

for tumor_name in tumor_names:

	print('Data:{}'.format(tumor_name))

	# truth_file = tumor_name + '/' + tumor_name + '_truth.bed'
	freebayes_file = tumor_name + '/' + tumor_name + '_freebayes_features.txt'
	freebayes_mask_file = tumor_name + '/' + tumor_name + '_freebayes_features_mask.txt'
	mutect_file = tumor_name + '/' + tumor_name + '_mutect_features.txt'
	mutect_mask_file = tumor_name + '/' + tumor_name + '_mutect_features_mask.txt'
	vardict_file = tumor_name + '/' + tumor_name + '_vardict_features.txt'
	vardict_mask_file = tumor_name + '/' + tumor_name + '_vardict_features_mask.txt'
	varscan_file = tumor_name + '/' + tumor_name + '_varscan_features.txt'
	varscan_mask_file = tumor_name + '/' + tumor_name + '_varscan_features_mask.txt'


	############## TRUTH ##############
	truth_list = list()
	# with open(truth_file,'r') as f_truth_file:
	# 	for line in f_truth_file:
	# 		if line[0] == '#':
	# 			continue

	# 		data = line.split('\t')

	# 		chrom_name = data[0]

	# 		if chrom_name == 'X':
	# 			chrom_name = '23'
	# 		elif chrom_name =='Y':
	# 			chrom_name = '24'

	# 		# if chrom_name == 'X' or chrom_name =='Y':
	# 		# 	continue
	# 		pos = data[1]
	# 		key = chrom_name + '_' + pos

	# 		truth_list.append(key)

	truth_set = set(truth_list)

	############## FREEBAYES ##############
	mask = np.loadtxt(freebayes_mask_file, dtype=np.int, delimiter='\t')
	mask_col_sum = np.sum(mask,axis=0)
	valid_cols_freebayes = np.where(mask_col_sum<len(mask))[0]
	num_valid_cols_freebayes = len(valid_cols_freebayes) 

	freebayes_list = list()
	freebayes_data_dict = dict()
	freebayes_mask_dict = dict()
	count = 0
	with open(freebayes_file,'r') as f_freebayes_file:
		next(f_freebayes_file)
		for line in f_freebayes_file:
			data = line.split('\t')

			chrom_name = data[-3]
			pos = data[-2]
			key = chrom_name + '_' + pos
			# print('key:{}'.format(key))

			freebayes_list.append(key)

			freebayes_data_dict[key] = np.array(data)[valid_cols_freebayes]
			# print(freebayes_data_dict[key])
			freebayes_mask_dict[key] = mask[count,valid_cols_freebayes]
			# print(freebayes_mask_dict[key])

			count += 1

	freebayes_set = set(freebayes_list)


	# ############## MUTECT ##############
	mask = np.loadtxt(mutect_mask_file, dtype=np.int, delimiter='\t')
	mask_col_sum = np.sum(mask,axis=0)
	valid_cols_mutect = np.where(mask_col_sum<len(mask))[0]
	num_valid_cols_mutect = len(valid_cols_mutect) 

	mutect_list = list()
	mutect_data_dict = dict()
	mutect_mask_dict = dict()
	count = 0
	with open(mutect_file,'r') as f_mutect_file:
		next(f_mutect_file)
		for line in f_mutect_file:
			data = line.split('\t')

			chrom_name = data[-3]
			pos = data[-2]
			key = chrom_name + '_' + pos
			# print('key:{}'.format(key))

			mutect_list.append(key)

			mutect_data_dict[key] = np.array(data)[valid_cols_mutect]
			# print(mutect_data_dict[key])
			mutect_mask_dict[key] = mask[count,valid_cols_mutect]
			# print(mutect_mask_dict[key])

			count += 1

	mutect_set = set(mutect_list)

	# ############## VARDICT ##############
	mask = np.loadtxt(vardict_mask_file, dtype=np.int, delimiter='\t')
	mask_col_sum = np.sum(mask,axis=0)
	valid_cols_vardict = np.where(mask_col_sum<len(mask))[0]
	num_valid_cols_vardict = len(valid_cols_vardict) 

	vardict_list = list()
	vardict_data_dict = dict()
	vardict_mask_dict = dict()
	count = 0
	with open(vardict_file,'r') as f_vardict_file:
		next(f_vardict_file)
		for line in f_vardict_file:
			data = line.split('\t')

			chrom_name = data[-3]
			pos = data[-2]
			key = chrom_name + '_' + pos
			# print('key:{}'.format(key))

			vardict_list.append(key)

			vardict_data_dict[key] = np.array(data)[valid_cols_vardict]
			# print(vardict_data_dict[key])
			vardict_mask_dict[key] = mask[count,valid_cols_vardict]
			# print(vardict_mask_dict[key])

			count += 1

	vardict_set = set(vardict_list)

	# ############## VARSCAN ##############
	mask = np.loadtxt(varscan_mask_file, dtype=np.int, delimiter='\t')
	mask_col_sum = np.sum(mask,axis=0)
	valid_cols_varscan = np.where(mask_col_sum<len(mask))[0]
	num_valid_cols_varscan = len(valid_cols_varscan) 

	varscan_list = list()
	varscan_data_dict = dict()
	varscan_mask_dict = dict()
	count = 0
	with open(varscan_file,'r') as f_varscan_file:
		next(f_varscan_file)
		for line in f_varscan_file:
			data = line.split('\t')

			chrom_name = data[-3]
			pos = data[-2]
			key = chrom_name + '_' + pos
			# print('key:{}'.format(key))

			varscan_list.append(key)

			varscan_data_dict[key] = np.array(data)[valid_cols_varscan]
			# print(varscan_data_dict[key])
			varscan_mask_dict[key] = mask[count,valid_cols_varscan]
			# print(varscan_mask_dict[key])

			count += 1

	varscan_set = set(varscan_list)


	set_list=[truth_set,freebayes_set,mutect_set,vardict_set,varscan_set]
	m = len(set_list)

	union_set = truth_set | freebayes_set | mutect_set | vardict_set | varscan_set
	sorted_union_list = sorted(list(union_set))
	num_preds = len(sorted_union_list)

	# Column names: truth,freebayes,mutect,vardict,varscan 
	labels = np.zeros((num_preds,m)) 
	for i in range(num_preds):
		key = sorted_union_list[i]
		for j in range(m):
			if key in set_list[j]:
				labels[i,j] = 1

	data_file = tumor_name + '/' + tumor_name +'_dataset.txt'
	with open(data_file,'w') as f_data_file:
		f_data_file.write('#chrom\tloc\ttruth\tfreebayes_pred\tmutect_pred\tvardict_pred\tvarscan_pred\t')
		for i in range(num_valid_cols_freebayes):
			ftr_name = 'F_' + str(feature_names[valid_cols_freebayes[i]])
			f_data_file.write(ftr_name + '\t')
		for i in range(num_valid_cols_mutect):
			ftr_name = 'M_' + str(feature_names[valid_cols_mutect[i]])
			f_data_file.write(ftr_name + '\t')
		for i in range(num_valid_cols_vardict):
			ftr_name = 'VD_' + str(feature_names[valid_cols_vardict[i]])
			f_data_file.write(ftr_name + '\t')
		for i in range(num_valid_cols_varscan-1):
			ftr_name = 'VS_' + str(feature_names[valid_cols_varscan[i]])
			f_data_file.write(ftr_name + '\t')
		ftr_name = 'VS_' + str(feature_names[valid_cols_varscan[-1]])
		f_data_file.write(ftr_name + '\n')

		for i in range(num_preds):
			key = sorted_union_list[i]
			chrom = key.split('_')[0]
			loc = key.split('_')[1]
			f_data_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chrom,loc,labels[i,0],labels[i,1],labels[i,2],labels[i,3],labels[i,4]))
			
			if key in freebayes_data_dict:
				temp_data = freebayes_data_dict[key]
			else:
				temp_data = np.zeros(num_valid_cols_freebayes,dtype=np.int)
			for j in range(len(temp_data)):
				f_data_file.write('{}\t'.format(temp_data[j]))

			if key in mutect_data_dict:
				temp_data = mutect_data_dict[key]
			else:
				temp_data = np.zeros(num_valid_cols_mutect,dtype=np.int)
			for j in range(len(temp_data)):
				f_data_file.write('{}\t'.format(temp_data[j]))

			if key in vardict_data_dict:
				temp_data = vardict_data_dict[key]
			else:
				temp_data = np.zeros(num_valid_cols_vardict,dtype=np.int)
			for j in range(len(temp_data)):
				f_data_file.write('{}\t'.format(temp_data[j]))

			if key in varscan_data_dict:
				temp_data = varscan_data_dict[key]
			else:
				temp_data = np.zeros(num_valid_cols_varscan,dtype=np.int)
			for j in range(len(temp_data)-1):
				f_data_file.write('{}\t'.format(temp_data[j]))

			f_data_file.write('{}\n'.format(temp_data[-1]))

	mask_file = tumor_name + '/' + tumor_name +'_dataset_mask.txt'
	with open(mask_file,'w') as f_mask_file:
		f_mask_file.write('#')
		for i in range(num_valid_cols_freebayes):
			ftr_name = 'F_' + str(feature_names[valid_cols_freebayes[i]])
			f_mask_file.write(ftr_name + '\t')
		for i in range(num_valid_cols_mutect):
			ftr_name = 'M_' + str(feature_names[valid_cols_mutect[i]])
			f_mask_file.write(ftr_name + '\t')
		for i in range(num_valid_cols_vardict):
			ftr_name = 'VD_' + str(feature_names[valid_cols_vardict[i]])
			f_mask_file.write(ftr_name + '\t')
		for i in range(num_valid_cols_varscan-1):
			ftr_name = 'VS_' + str(feature_names[valid_cols_varscan[i]])
			f_mask_file.write(ftr_name + '\t')
		ftr_name = 'VS_' + str(feature_names[valid_cols_varscan[-1]])
		f_mask_file.write(ftr_name + '\n')



		for i in range(num_preds):
			key = sorted_union_list[i]

			if key in freebayes_mask_dict:
				temp_mask = freebayes_mask_dict[key]
			else:
				temp_mask = np.ones(num_valid_cols_freebayes,dtype=np.int)
			for j in range(len(temp_mask)):
				f_mask_file.write('{}\t'.format(temp_mask[j]))

			if key in mutect_mask_dict:
				temp_mask = mutect_mask_dict[key]
			else:
				temp_mask = np.ones(num_valid_cols_mutect,dtype=np.int)
			for j in range(len(temp_mask)):
				f_mask_file.write('{}\t'.format(temp_mask[j]))

			if key in vardict_mask_dict:
				temp_mask = vardict_mask_dict[key]
			else:
				temp_mask = np.ones(num_valid_cols_vardict,dtype=np.int)
			for j in range(len(temp_mask)):
				f_mask_file.write('{}\t'.format(temp_mask[j]))

			if key in varscan_mask_dict:
				temp_mask = varscan_mask_dict[key]
			else:
				temp_mask = np.ones(num_valid_cols_varscan,dtype=np.int)
			for j in range(len(temp_mask)-1):
				f_mask_file.write('{}\t'.format(temp_mask[j]))

			f_mask_file.write('{}\n'.format(temp_mask[-1]))






