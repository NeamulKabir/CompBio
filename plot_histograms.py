import numpy as np
import matplotlib.pyplot as plt

feature_names = ['GT_AF','GT_BIAS_REF','GT_BIAS_ALT','GT_DP','GT_GQ','AB','AC','AF','BaseQRankSum','DP','FS','GC','HaplotypeScore','MQ','MQRanksum','ReadPosRakSum','chrom','loc']

# tumor_names = ['real1','syn1','syn2','syn3','syn4','syn5']
# file_names = ['freebayes','mutect','vardict','varscan']

tumor_names = ['real1']
file_names = ['varscan']

for tumor_name in tumor_names:
	for file_name in file_names:

		print('tumor_name:{} - file_name:{}'.format(tumor_name,file_name))

		features_file = tumor_name + '/' + tumor_name +'_' + file_name + '_features.txt'
		features_mask_file = tumor_name + '/' + tumor_name +'_' + file_name + '_features_mask.txt'
		
		features = np.loadtxt(features_file, delimiter='\t', skiprows=1, usecols = np.arange(15))
		features_mask = np.loadtxt(features_mask_file, delimiter='\t', skiprows=0)

		for i in range(len(feature_names)):
			ft_name = feature_names[i]

			ft_mask = features_mask[:,i]
			if np.sum(ft_mask) == len(ft_mask):
				continue

			ft_data = features[:,i]

			masked_ft_data = ft_data[ft_mask==0]
			print('Missing data for {}:{}'.format(ft_name,len(ft_data)-len(masked_ft_data)))

			n, bins, patches = plt.hist(masked_ft_data, 10, density=True, facecolor='g', alpha=0.75)


			plt.xlabel('Value')
			plt.ylabel('Probability')
			plt.title('Histogram of {} - {}/{}'.format(ft_name,tumor_name,file_name))
			plt.xlim((np.amin(masked_ft_data), np.amax(masked_ft_data)))
			# plt.yscale('log')
			plt.grid(True)
			plt.show()


