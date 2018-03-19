
import matplotlib.pyplot as plt
import numpy as np

tumor_name = 'syn1'

truth_file = tumor_name + '/' + tumor_name + '_truth.bed'
freebayes_file = tumor_name + '/' + tumor_name + '_freebayes.vcf'
mutect_file = tumor_name + '/' + tumor_name + '_mutect.vcf'
vardict_file = tumor_name + '/' + tumor_name + '_vardict.vcf'
varscan_file = tumor_name + '/' + tumor_name + '_varscan.vcf'


############## TRUTH ##############
truth_list = list()
with open(truth_file,'r') as f_truth_file:
	for line in f_truth_file:
		if line[0] == '#':
			continue

		data = line.split('\t')

		chrom_name = data[0]
		pos = data[1]
		key = chrom_name + '_' + pos

		truth_list.append(key)

truth_set = set(truth_list)

############## FREEBAYES ##############
freebayes_list = list()
with open(freebayes_file,'r') as f_freebayes_file:
	for line in f_freebayes_file:
		if line[0] == '#':
			continue

		data = line.split('\t')

		chrom_name = data[0]
		pos = data[1]
		key = chrom_name + '_' + pos

		freebayes_list.append(key)

freebayes_set = set(freebayes_list)

############## MUTECT ##############
mutect_list = list()
with open(mutect_file,'r') as f_mutect_file:
	for line in f_mutect_file:
		if line[0] == '#':
			continue

		data = line.split('\t')

		chrom_name = data[0]
		pos = data[1]
		key = chrom_name + '_' + pos

		mutect_list.append(key)

mutect_set = set(mutect_list)

############## VARDICT ##############
vardict_list = list()
with open(vardict_file,'r') as f_vardict_file:
	for line in f_vardict_file:
		if line[0] == '#':
			continue

		data = line.split('\t')

		chrom_name = data[0]
		pos = data[1]
		key = chrom_name + '_' + pos

		vardict_list.append(key)

vardict_set = set(vardict_list)

############## VARSCAN ##############
varscan_list = list()
with open(varscan_file,'r') as f_varscan_file:
	for line in f_varscan_file:
		if line[0] == '#':
			continue

		data = line.split('\t')

		chrom_name = data[0]
		pos = data[1]
		key = chrom_name + '_' + pos

		varscan_list.append(key)

varscan_set = set(varscan_list)


set_list=[truth_set,freebayes_set,mutect_set,vardict_set,varscan_set]
m = len(set_list)
# intersection_count = np.zeros((m,m))
# for i in range(m):
# 	first_set = set_list[i]
# 	intersection_count[i,i] = len(first_set)
# 	for j in range(i+1,m):
# 		second_set = set_list[j]
# 		intersection_count[i,j] = len(first_set.intersection(second_set))

# print(intersection_count)

union_set = freebayes_set | mutect_set | vardict_set | varscan_set
sorted_union_list = sorted(list(union_set))
num_preds = len(sorted_union_list)

# Column names: truth,freebayes,mutect,vardict,varscan 
labels = np.zeros((num_preds,m)) 
for i in range(num_preds):
	key = sorted_union_list[i]
	for j in range(m):
		if key in set_list[j]:
			labels[i,j] = 1

labels_file = tumor_name + '/' + tumor_name +'_labels.txt'
with open(labels_file,'a') as f_labels_file:
	f_labels_file.write('chrom\tloc\ttruth\tfreebayes\tmutect\tvardict\tvarscan\n')
	for i in range(num_preds):
		key = sorted_union_list[i]
		chrom = key.split('_')[0]
		loc = key.split('_')[1]
		f_labels_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom,loc,labels[i,0],labels[i,1],labels[i,2],labels[i,3],labels[i,4]))











