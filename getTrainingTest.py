import vcf
import numpy as np
from random import randint
from sklearn.linear_model import LogisticRegression

tumor_name = 'syn3'
tumor_names = ['syn2', 'syn3', 'syn4', 'syn5', 'real1', 'real2_chr1to5']

trainTarget = list()
testTarget = list()
trainingList = list()
testList = list()

for it in range(6):
	featureList = list()
	targetList = list()

	with open(tumor_names[it]+"/Data/"+tumor_names[it]+'_merged_features.txt', 'r') as f:
		next(f)
		for line in f:
			tempList = list()
			data = line.split("\t")
			for i in range(len(data)-1):
				tempList.append(float(data[i]))
			featureList.append(tempList)

	with open(tumor_names[it]+"/Data/"+tumor_names[it]+'_target.txt', 'r') as f:
		next(f)
		for line in f:
			targetList.append(line)

	featureSize = len(featureList)
	trainingSize = int(0.8 * featureSize)
	testSize = featureSize - trainingSize
	trainCount = 0
	testCount = 0

	for iterator in range(featureSize):
		a = randint(0,9)
		if a>4 and testCount < testSize:
			testList.append(featureList[iterator])
			testTarget.append(targetList[iterator])
			testCount += 1
		else:
			trainingList.append(featureList[iterator])
			trainTarget.append(targetList[iterator])
			trainCount += 1


print len(trainingList), len(trainTarget)
print len(testList), len(testTarget)

logisticRegr = LogisticRegression()

logisticRegr.fit(trainingList,trainTarget)
testSample = np.array(testList)
predictions = logisticRegr.predict(testSample)

size = len(predictions)
match = 0
mismatch  = 0

tp,fp,tn,fn = 0,0,0,0

for iterate in range(size):
	if(testTarget[iterate] == predictions[iterate]):
		match += 1
		if(testTarget[iterate] == 0):
			tn += 1
		else:
			tp += 1
	else:
		mismatch += 1
		if(testTarget[iterate] == 1):
			fn +=1
		else:
			fp += 1

mPer = match * 100 / size
mmPer = mismatch * 100 /size
precision = float(tp)/float((tp+fp))
recall = float(tp)/float((tp+fn))
accuracy = float(float(tp+tn)/float(tp+tn+fp+fn))

print match, mismatch, mPer, mmPer
print precision, recall, accuracy


#print len(predictions)
'''m = len(trainingList)
n = len(trainingList[0])

labels_file = 'training.txt'
with open(labels_file,'w') as f_labels_file:
	for i in range(m):
		for j in range(n):
			f_labels_file.write('{}\t'.format(trainingList[i][j]))
		f_labels_file.write('\n')

m = len(testList)
n = len(testList[0])
labels_file = 'test.txt'
with open(labels_file,'w') as f_labels_file:
	for i in range(m):
		for j in range(n):
			f_labels_file.write('{}\t'.format(testList[i][j]))
		f_labels_file.write('\n')

m = len(trainTarget)
labels_file = 'trainingTarget.txt'
with open(labels_file,'w') as f_labels_file:
	for i in range(m):
		f_labels_file.write('{}\n'.format(trainTarget[i]))

m = len(testTarget)
labels_file = 'testTarget.txt'
with open(labels_file,'w') as f_labels_file:
	for i in range(m):
		f_labels_file.write('{}\n'.format(testTarget[i]))


'''

