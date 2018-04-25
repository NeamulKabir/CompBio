from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

tumorName = "real2"
fType = "train"
ensmebleFile = tumorName+ "/Ensemble_"+ tumorName +"_train_retained.csv"
posFile = tumorName+ "/" + tumorName + "_"+fType+"_chrom.txt"
dataFile = tumorName+ "/" + tumorName + "_"+fType+"_data.txt"
maskFile = tumorName+ "/" + tumorName + "_"+fType+"_data_mask.txt"
print('Data File: {}\n'.format(tumorName))
####################################################################################################
posArr = [33,38,56,57,81,82,83,98]
sizePos = len(posArr)

ensembleDict = dict()
cropFeature = list()
ensPos = list()
counter = 0
with open(ensmebleFile) as f_file:
	featureNames = f_file.readline()
	for line in f_file:
		tempList = list()
		tempPos = list()

		data = line.split(',')
		if data[0] == 'X':
			tempPos.append('23')
		elif data[0] == 'Y':
			tempPos.append('24')						#read additional feature from ensemble
		else:
			tempPos.append(data[0])
		tempPos.append(data[1])
		temp = tempPos[0] +'_'+ tempPos[1]
		#if temp in ensembleDict.keys():
		#	print temp

		for i in xrange(sizePos):
			tempList.append(float(data[posArr[i]]))
		ensembleDict[temp] = tempList
		#print ensembleDict[temp]

		counter += 1

print len(ensembleDict)
print counter

selection = featureNames.split(',')
newFeatureName = list()
for i in xrange(sizePos):
	newFeatureName.append(selection[posArr[i]])
	print selection[posArr[i]]
##################################################################################################
missingData = list()
for i in xrange(sizePos):
	missingData.append('0.0')

dataPos = list()
with open(posFile, 'r') as p_file:
	next(p_file)
	for line in p_file:
		tempData = line.split('\n')
		data = tempData[0].split('\t')
		temp = data[0] + '_' + data[1]
		#print temp
		dataPos.append(temp)
print len(dataPos), len(data)
##################################################################################################
mergedData = list()
lineCount = 0
found = 0
with open(dataFile, 'r') as d_file:
	tempDataFeatures = d_file.readline()
	dataFeatures = tempDataFeatures.split('\n')
	for line in d_file:
		tempList = list()
		data = line.split('\t')
		size = len(data)
		for i in xrange(size):
			tempList.append(float(data[i]))				#merge ensemble feature with existing

		if dataPos[lineCount] in ensembleDict.keys():
			found += 1
			addData = ensembleDict[dataPos[lineCount]]
			for iterator in xrange(sizePos):
				tempList.append(float(addData[iterator]))
		else:
			for iterator in xrange(sizePos):
				tempList.append(float(missingData[iterator]))
		mergedData.append(tempList)
		lineCount += 1
mDataSize = len(mergedData)
mFeatureSize = len(mergedData[0])
print lineCount, mDataSize, mFeatureSize
print found
##################################################################################################
merged_data_file = tumorName+ "/" + tumorName +"_merged_ensemble_"+ fType+ "_data.txt"
with open(merged_data_file,'w') as md_file:
	md_file.write(dataFeatures[0])
	for i in xrange(sizePos):
		md_file.write('\t{}'.format(selection[posArr[i]]))
	md_file.write('\n')
	for i in xrange(mDataSize):
		for j in xrange(mFeatureSize-1):
			md_file.write('{}\t'.format(mergedData[i][j]))
		md_file.write('{}'.format(mergedData[i][mFeatureSize-1]))
		md_file.write('\n')
##################################################################################################
mergedMask = list()
lineCount = 0
found = 0
with open(maskFile, 'r') as m_file:
	next(m_file)
	for line in m_file:
		tempList = list()
		data = line.split('\t')
		size = len(data)
		for i in xrange(size):
			tempList.append(int(data[i]))
		if dataPos[lineCount] in ensembleDict.keys():		# update mask file for ensemble
			found += 1
			for iterator in xrange(sizePos):
				tempList.append(0)
		else:
			for iterator in xrange(sizePos):
				tempList.append(1)
		mergedMask.append(tempList)

		lineCount += 1
print lineCount, len(mergedMask), len(mergedMask[0])
##################################################################################################
merged_mask_file =  tumorName+ "/" +tumorName +"_merged_ensemble_"+fType+"_data_mask.txt"
with open(merged_mask_file,'w') as mask_merge_file:
	mask_merge_file.write(dataFeatures[0])
	for i in xrange(sizePos):
		mask_merge_file.write('\t{}'.format(newFeatureName[i]))
	mask_merge_file.write('\n')
	for i in xrange(mDataSize):
		for j in xrange(mFeatureSize-1):
			mask_merge_file.write('{}\t'.format(mergedMask[i][j]))
		mask_merge_file.write('{}'.format(mergedMask[i][mFeatureSize-1]))
		mask_merge_file.write('\n')


		




