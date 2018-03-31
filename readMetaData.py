import vcf
import numpy as np

tumor_name = 'syn3'
file_name = 'freebayes'
vcf_reader = vcf.Reader(open(tumor_name+'/'+tumor_name+'_'+file_name+'.vcf','r'))

counter = 0
missingCount = dict()
feature = list()
maskedList = list()
uniqueGT = set()

for record in vcf_reader:
	templist = list()
	tempMaskedList = list()

	GTflag = False
	for sample in record.samples:
		gtValue = sample['GT']
		uniqueGT.add(gtValue)
		if gtValue == '.':
			GTflag = True
		# elif gtValue == ('0/1' | '1/1' | '1' | '1|1' | '0|1'):
		elif gtValue in ['0/1','1/1','1','1|1','0|1','1|0','1/0']:
			#print('{}\tneamul'.format(gtValue))
			GTflag = False
			try:
				data = sample['AF']
				templist.append(data)
				tempMaskedList.append(0)
			except:
				templist.append(0)
				tempMaskedList.append(1)

			try:
				data = sample['BIAS']
				templist.append(data)
				tempMaskedList.append(0)
			except:
				templist.append(0)
				tempMaskedList.append(1)

			try:
				data = sample['DP']
				templist.append(data)
				tempMaskedList.append(0)
			except:
				templist.append(0)
				tempMaskedList.append(1)

			try:
				data = sample['GQ']
				templist.append(data)
				tempMaskedList.append(0)

			except:
				templist.append(0)
				tempMaskedList.append(1)

			break

		# else:
			#print('{}\tmustafa'.format(gtValue))
	if GTflag :
		#print('{}\tHerty'.format(gtValue))
		continue

	
	try:
		data = record.INFO['AB']
		templist.append(data[0])
		tempMaskedList.append(0)

	except KeyError:
		templist.append(0)
		tempMaskedList.append(1)
		if 'AB' not in missingCount:
			missingCount['AB'] = 1
		else:
			missingCount['AB'] += 1
	
	try:
		data = record.INFO['AC']
		templist.append(data[0])
		tempMaskedList.append(0)
	except KeyError:
		templist.append(0)
		tempMaskedList.append(1)
		if 'AC' not in missingCount:
			missingCount['AC'] = 1
		else:
			missingCount['AC'] += 1
	

	try:
		data = record.INFO['AF']
		templist.append(data[0])
		tempMaskedList.append(0)
	except KeyError:
		templist.append(0)
		tempMaskedList.append(1)
		if 'AF' not in missingCount:
			missingCount['AF'] = 1
		else:
			missingCount['AF'] += 1
	

	try:
		data = record.INFO['BaseQRankSum']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'BaseQRankSum' not in missingCount:
			missingCount['BaseQRankSum'] = 1
		else:
			missingCount['BaseQRankSum'] += 1
	templist.append(data)


	try:
		data = record.INFO['DP']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'DP' not in missingCount:
			missingCount['DP'] = 1
		else:
			missingCount['DP'] += 1
	templist.append(data)


	try:
		data = record.INFO['FS']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'FS' not in missingCount:
			missingCount['FS'] = 1
		else:
			missingCount['FS'] += 1
	templist.append(data)

	try:
		data = record.INFO['GC']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'GC' not in missingCount:
			missingCount['GC'] = 1
		else:
			missingCount['GC'] += 1
	templist.append(data)

	try:
		data = record.INFO['HaplotypeScore']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'HaplotypeScore' not in missingCount:
			missingCount['HaplotypeScore'] = 1
		else:
			missingCount['HaplotypeScore'] += 1
	templist.append(data)

	try:
		data = record.INFO['MQ']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'MQ' not in missingCount:
			missingCount['MQ'] = 1
		else:
			missingCount['MQ'] += 1
	templist.append(data)

	try:
		data = record.INFO['MQRankSum']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'MQRankSum' not in missingCount:
			missingCount['MQRankSum'] = 1
		else:
			missingCount['MQRankSum'] += 1
	templist.append(data)

	try:
		data = record.INFO['ReadPosRankSum']
		tempMaskedList.append(0)
	except KeyError:
		tempMaskedList.append(1)
		data = 0
		if 'ReadPosRankSum' not in missingCount:
			missingCount['ReadPosRankSum'] = 1
		else:
			missingCount['ReadPosRankSum'] += 1
	templist.append(data)

	templist.append(record.CHROM)
	templist.append(record.POS)

	feature.append(templist)
	maskedList.append(tempMaskedList)
	#print(len(tempMaskedList))

print uniqueGT

labels_file = tumor_name + '/' + tumor_name +'_'+file_name+'_features1.txt'
with open(labels_file,'w') as f_labels_file:
	f_labels_file.write('GT_AF\tGT_BIAS\tGT_DP\tGT_GQ\tAB\tAC\tAF\tBaseQRankSum\tFS\tGC\tHaplotypeScore\tMQ\tMQRanksum\tReadPosRakSum\tchrom\tloc\n')
	for i in range(len(feature)):
		for j in range(len(feature[0])):
			f_labels_file.write('{}\t'.format(feature[i][j]))
		f_labels_file.write('\n')

featureArr2 = np.array(feature)[:,:-2]
featureArr = np.array(featureArr2, dtype=np.float)

maskedArr = np.array(maskedList, dtype=np.bool)

mask_file = tumor_name + '/' + tumor_name +'_'+file_name+'_mask.txt'
np.savetxt(mask_file,maskedArr,delimiter = '\t', fmt='%d')

print('type of feature array:{}'.format(type(featureArr)))
print('shape of feature array:{}'.format(featureArr.shape))
print('shape of masked array:{}'.format(maskedArr.shape))

maskedFeature = np.ma.masked_array(featureArr, mask=maskedArr)

featureMin = np.amin(maskedFeature, axis=0)
featureMax = np.amax(maskedFeature, axis=0)
print('feature min:{}'.format(featureMin))
print('feature max:{}'.format(featureMax))

print('Missing counts:{}\n Total records: {}'.format(missingCount,len(feature)))














