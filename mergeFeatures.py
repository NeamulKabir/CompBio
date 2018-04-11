import vcf
import numpy as np 


tumor_name = 'syn2'
file_name = ['freebayes', 'mutect', 'vardict', 'varscan']

featureList = list()
for iterator in range(len(file_name)):
	freebayesList = list()
	with open(tumor_name+"/Data/"+tumor_name+'_'+file_name[iterator]+'_features1.txt', 'r') as f:
		next(f)
		for line in f:
			tempList = list()
			data = line.split("\t")
			for i in range(len(data) - 1):
				if data[i] == 'X':
					tempList.append(25)
				elif data[i] == 'Y':
					tempList.append(26)
				else:
					tempList.append(data[i])
			freebayesList.append(tempList)
	featureList.append(freebayesList)

labelsList = list()
with open(tumor_name+"/Data/"+tumor_name+'_modified_labels.txt', 'r') as f:
	next(f)
	for line in f:
		tempList = list()
		data = line.split("\t")
		for i in range(len(data)):
			tempList.append(data[i])
		labelsList.append(tempList)

m = len(labelsList)
n = len(featureList[0][0])

print m,n
#print featureList[0][0]

finalFeatures = np.zeros((m,n))
iterator1 = 0
iterator2 = 0
iterator3 = 0
iterator4 = 0
flag = [False, False, False, False ]

'''for i in range(10):
	print labelsList[i][0], labelsList[i][1]
	if labelsList[i][3] == '1':
		if flag[0] == False:
			flag[0] = True 
		else:
			print "freeBayes"
			print featureList[0][iterator1][15],featureList[0][iterator1][16]
			iterator1 += 1

	if labelsList[i][4] == '1':
		if flag[1] == False:
			flag[1] = True 
		else:
			print "mutect" 
			print featureList[1][iterator2][15], featureList[1][iterator2][16]
			iterator2 += 1

	if labelsList[i][5] == '1':
		if flag[2] == False:
			flag[2] = True 
		else:
			print "vardict"
			print featureList[2][iterator3][15] , featureList[2][iterator3][16] 
			iterator3 += 1

	if labelsList[i][6] == '1':
		if flag[2] == False:
			flag[2] = True 
		else:
			print "varscan" 
			print featureList[3][iterator4][15], featureList[3][iterator4][16]
			iterator4 += 1

'''
print len(featureList[0]),len(featureList[1]),len(featureList[2]),len(featureList[3])
length = [len(featureList[0]),len(featureList[1]),len(featureList[2]),len(featureList[3])]



for i in range(m):
	#print labelsList[i][0], labelsList[i][1]
	if labelsList[i][3] == '1'and iterator1 < length[0]:
		if flag[0] == False:
			flag[0] = True 
		else:
			for j in range(n):
				if finalFeatures[i][j] < featureList[0][iterator1][j]:
					finalFeatures[i][j] = featureList[0][iterator1][j]
			iterator1 += 1

	if labelsList[i][4] == '1'and iterator2 < length[1]:
		if flag[1] == False:
			flag[1] = True 
		else:
			for j in range(n):
				if finalFeatures[i][j] < featureList[1][iterator2][j]:
					finalFeatures[i][j] = featureList[1][iterator2][j]
			iterator2 += 1

	if labelsList[i][5] == '1' and iterator3 < length[2]:
		if flag[2] == False:
			flag[2] = True 
		else:
			for j in range(n):
				#print iterator3,labelsList[i][0], labelsList[i][1], labelsList[i][0], labelsList[i][5]
				if finalFeatures[i][j] < featureList[2][iterator3][j]:
					finalFeatures[i][j] = featureList[2][iterator3][j] 
			iterator3 += 1

	if labelsList[i][6] == '1' and iterator4 < length[3]:
		if flag[2] == False:
			flag[2] = True 
		else:
			for j in range(n):
				if finalFeatures[i][j] < featureList[3][iterator4][j]:
					finalFeatures[i][j] = featureList[3][iterator4][j]
			iterator4 += 1

labels_file = tumor_name + '/Data/' + tumor_name +'_merged_features.txt'
with open(labels_file,'w') as f_labels_file:
	f_labels_file.write('GT_AF\tGT_BIAS\tGT_DP\tGT_GQ\tAB\tAC\tAF\tBaseQRankSum\tDP\tFS\tGC\tHaplotypeScore\tMQ\tMQRanksum\tReadPosRakSum\tchrom\tloc\n')
	for i in range(m):
		for j in range(n):
			f_labels_file.write('{}\t'.format(finalFeatures[i][j]))
		f_labels_file.write('\n')



target_file = tumor_name + '/Data/' + tumor_name +'_target.txt'
with open(target_file,'w') as f_labels_file:
	f_labels_file.write('truth\n')
	for i in range(m):
		f_labels_file.write('{}\n'.format(labelsList[i][2]))














