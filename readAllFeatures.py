import vcf
import numpy as np

def getGenData(sample, index):
	try:
		genData = sample[index]
		flag = 0
	except:
		genData = 0
		flag = 1
	return genData,flag


tumor_name = 'real1'
file_name = 'vardict'
vcf_reader = vcf.Reader(open(tumor_name+'/'+tumor_name+'_'+file_name+'.vcf','r'))

counter = 0
missingCount = dict()
feature = list()
maskedList = list()
uniqueGT = set()

for record in vcf_reader:
	templist = list()
	tempMaskedList = list()

	if record.CHROM != 'X' or record.CHROM !='Y':

		GTflag = False
		for sample in record.samples:
			data = getGenData(sample,'AD')
			templist.append(data[0])
			tempMaskedList.append(data[1])
			
			data = getGenData(sample,'ADJAF')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'AF')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'ALD')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'AO')
			templist.append(data[0])
			tempMaskedList.append(data[1])
				
			data = getGenData(sample,'BIAS')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'BQ')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'DP')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'DP4')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'FREQ')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			tdata = getGenData(sample,'GQ')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'HIAF')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'MQ')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'NM')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'ODDRATIO')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'PL')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'PMEAN')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'PSTD')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'QA')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'QR')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'QSTD')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'QUAL')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'RD')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'RO')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'SBF')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'SN')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'SS')
			templist.append(data[0])
			tempMaskedList.append(data[1])

			data = getGenData(sample,'VD')
			templist.append(data[0])
			tempMaskedList.append(data[1])

		try:
			data = record.INFO['QUAL']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'QUAL' not in missingCount:
				missingCount['QUAL'] = 1
			else:
				missingCount['QUAL'] += 1
		
		try:
			data = record.INFO['AB']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'AB' not in missingCount:
				missingCount['AB'] = 1
			else:
				missingCount['AB'] += 1

		try:
			data = record.INFO['ABP']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'ABP' not in missingCount:
				missingCount['ABP'] = 1
			else:
				missingCount['ABP'] += 1
		
		try:
			data = record.INFO['AC']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
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
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'AF' not in missingCount:
				missingCount['AF'] = 1
			else:
				missingCount['AF'] += 1
		
		try:
			data = record.INFO['AN']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'AN' not in missingCount:
				missingCount['AN'] = 1
			else:
				missingCount['AN'] += 1
		templist.append(data)

		try:
			data = record.INFO['AO']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'AO' not in missingCount:
				missingCount['AO'] = 1
			else:
				missingCount['AO'] += 1

		try:
			data = record.INFO['BaseQRankSum']
			tempMaskedList.append(0)
		except:
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
		except:
			tempMaskedList.append(1)
			data = 0
			if 'DP' not in missingCount:
				missingCount['DP'] = 1
			else:
				missingCount['DP'] += 1
		templist.append(data)

		try:
			data = record.INFO['DPB']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'DPB' not in missingCount:
				missingCount['DPB'] = 1
			else:
				missingCount['DPB'] += 1
		templist.append(data)

		try:
			data = record.INFO['DPRA']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'DPRA' not in missingCount:
				missingCount['DPRA'] = 1
			else:
				missingCount['DPRA'] += 1


		try:
			data = record.INFO['EPP']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'EPP' not in missingCount:
				missingCount['EPP'] = 1
			else:
				missingCount['EPP'] += 1

		try:
			data = record.INFO['EPPR']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'EPPR' not in missingCount:
				missingCount['EPPR'] = 1
			else:
				missingCount['EPPR'] += 1

		try:
			data = record.INFO['FS']
			tempMaskedList.append(0)
		except:
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
		except:
			tempMaskedList.append(1)
			data = 0
			if 'GC' not in missingCount:
				missingCount['GC'] = 1
			else:
				missingCount['GC'] += 1
		templist.append(data)

		try:
			data = record.INFO['GPV']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'GPV' not in missingCount:
				missingCount['GPV'] = 1
			else:
				missingCount['GPV'] += 1
		templist.append(data)

		try:
			data = record.INFO['GTI']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'GTI' not in missingCount:
				missingCount['GTI'] = 1
			else:
				missingCount['GTI'] += 1
		templist.append(data)

		try:
			data = record.INFO['HRun']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'HRun' not in missingCount:
				missingCount['HRun'] = 1
			else:
				missingCount['HRun'] += 1
		templist.append(data)

		try:
			data = record.INFO['HaplotypeScore']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'HaplotypeScore' not in missingCount:
				missingCount['HaplotypeScore'] = 1
			else:
				missingCount['HaplotypeScore'] += 1
		templist.append(data)

		try:
			data = record.INFO['LEN']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'LEN' not in missingCount:
				missingCount['LEN'] = 1
			else:
				missingCount['LEN'] += 1

		try:
			data = record.INFO['MEANALT']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'MEANALT' not in missingCount:
				missingCount['MEANALT'] = 1
			else:
				missingCount['MEANALT'] += 1

		try:
			data = record.INFO['MQ']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'MQ' not in missingCount:
				missingCount['MQ'] = 1
			else:
				missingCount['MQ'] += 1
		templist.append(data)

		try:
			data = record.INFO['MQ0']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'MQ0' not in missingCount:
				missingCount['MQ0'] = 1
			else:
				missingCount['MQ0'] += 1
		templist.append(data)

		try:
			data = record.INFO['MQM']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'MQM' not in missingCount:
				missingCount['MQM'] = 1
			else:
				missingCount['MQM'] += 1

		try:
			data = record.INFO['MQMR']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'MQMR' not in missingCount:
				missingCount['MQMR'] = 1
			else:
				missingCount['MQMR'] += 1
		templist.append(data)

		try:
			data = record.INFO['MQRankSum']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'MQRankSum' not in missingCount:
				missingCount['MQRankSum'] = 1
			else:
				missingCount['MQRankSum'] += 1
		templist.append(data)

		try:
			data = record.INFO['MSI']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'MSI' not in missingCount:
				missingCount['MSI'] = 1
			else:
				missingCount['MSI'] += 1
		templist.append(data)

		try:
			data = record.INFO['MSILEN']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'MSILEN' not in missingCount:
				missingCount['MSILEN'] = 1
			else:
				missingCount['MSILEN'] += 1
		templist.append(data)

		try:
			data = record.INFO['NS']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'NS' not in missingCount:
				missingCount['NS'] = 1
			else:
				missingCount['NS'] += 1
		templist.append(data)

		try:
			data = record.INFO['NUMALT']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'NUMALT' not in missingCount:
				missingCount['NUMALT'] = 1
			else:
				missingCount['NUMALT'] += 1
		templist.append(data)

		try:
			data = record.INFO['ODDS']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'ODDS' not in missingCount:
				missingCount['ODDS'] = 1
			else:
				missingCount['ODDS'] += 1
		templist.append(data)

		try:
			data = record.INFO['PAIRED']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'PAIRED' not in missingCount:
				missingCount['PAIRED'] = 1
			else:
				missingCount['PAIRED'] += 1

		try:
			data = record.INFO['PAIREDR']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'PAIREDR' not in missingCount:
				missingCount['PAIREDR'] = 1
			else:
				missingCount['PAIREDR'] += 1
		templist.append(data)

		try:
			data = record.INFO['PAO']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'PAO' not in missingCount:
				missingCount['PAO'] = 1
			else:
				missingCount['PAO'] += 1
		templist.append(data)

		try:
			data = record.INFO['PQA']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'PQA' not in missingCount:
				missingCount['PQA'] = 1
			else:
				missingCount['PQA'] += 1
		templist.append(data)

		try:
			data = record.INFO['PQR']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'PQR' not in missingCount:
				missingCount['PQR'] = 1
			else:
				missingCount['PQR'] += 1
		templist.append(data)

		try:
			data = record.INFO['PRO']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'PRO' not in missingCount:
				missingCount['PRO'] = 1
			else:
				missingCount['PRO'] += 1
		templist.append(data)

		try:
			data = record.INFO['QA']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'QA' not in missingCount:
				missingCount['QA'] = 1
			else:
				missingCount['QA'] += 1

		try:
			data = record.INFO['QD']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'QD' not in missingCount:
				missingCount['QD'] = 1
			else:
				missingCount['QD'] += 1
		templist.append(data)

		try:
			data = record.INFO['QR']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'QR' not in missingCount:
				missingCount['QR'] = 1
			else:
				missingCount['QR'] += 1
		templist.append(data)

		try:
			data = record.INFO['RO']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'RO' not in missingCount:
				missingCount['RO'] = 1
			else:
				missingCount['RO'] += 1
		templist.append(data)

		try:
			data = record.INFO['RPL']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'RPL' not in missingCount:
				missingCount['RPL'] = 1
			else:
				missingCount['RPL'] += 1

		try:
			data = record.INFO['RPP']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'RPP' not in missingCount:
				missingCount['RPP'] = 1
			else:
				missingCount['RPP'] += 1

		try:
			data = record.INFO['RPPR']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'RPPR' not in missingCount:
				missingCount['RPPR'] = 1
			else:
				missingCount['RPPR'] += 1
		templist.append(data)

		try:
			data = record.INFO['RPR']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'RPR' not in missingCount:
				missingCount['RPR'] = 1
			else:
				missingCount['RPR'] += 1

		try:
			data = record.INFO['RUN']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'RUN' not in missingCount:
				missingCount['RUN'] = 1
			else:
				missingCount['RUN'] += 1

		try:
			data = record.INFO['ReadPosRankSum']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'ReadPosRankSum' not in missingCount:
				missingCount['ReadPosRankSum'] = 1
			else:
				missingCount['ReadPosRankSum'] += 1
		templist.append(data)

		try:
			data = record.INFO['SAF']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'SAF' not in missingCount:
				missingCount['SAF'] = 1
			else:
				missingCount['SAF'] += 1

		try:
			data = record.INFO['SAP']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'SAP' not in missingCount:
				missingCount['SAP'] = 1
			else:
				missingCount['SAP'] += 1

		try:
			data = record.INFO['SAR']
			templist.append(data[0])
			tempMaskedList.append(0)
		except:
			templist.append(0)
			tempMaskedList.append(1)
			if 'SAR' not in missingCount:
				missingCount['SAR'] = 1
			else:
				missingCount['SAR'] += 1
		try:
			data = record.INFO['SOR']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SOR' not in missingCount:
				missingCount['SOR'] = 1
			else:
				missingCount['SOR'] += 1
		templist.append(data)

		try:
			data = record.INFO['SPV']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SPV' not in missingCount:
				missingCount['SPV'] = 1
			else:
				missingCount['SPV'] += 1
		templist.append(data)

		try:
			data = record.INFO['SRF']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SRF' not in missingCount:
				missingCount['SRF'] = 1
			else:
				missingCount['SRF'] += 1
		templist.append(data)

		try:
			data = record.INFO['SRP']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SRP' not in missingCount:
				missingCount['SRP'] = 1
			else:
				missingCount['SRP'] += 1
		templist.append(data)

		try:
			data = record.INFO['SRR']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SRR' not in missingCount:
				missingCount['SRR'] = 1
			else:
				missingCount['SRR'] += 1
		templist.append(data)

		try:
			data = record.INFO['SS']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SS' not in missingCount:
				missingCount['SS'] = 1
			else:
				missingCount['SS'] += 1
		templist.append(data)

		try:
			data = record.INFO['SSC']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SSC' not in missingCount:
				missingCount['SSC'] = 1
			else:
				missingCount['SSC'] += 1
		templist.append(data)

		try:
			data = record.INFO['SSF']
			tempMaskedList.append(0)
		except:
			tempMaskedList.append(1)
			data = 0
			if 'SSF' not in missingCount:
				missingCount['SSF'] = 1
			else:
				missingCount['SSF'] += 1
		templist.append(data)


		templist.append(record.CHROM)
		templist.append(record.POS)

		feature.append(templist)
		maskedList.append(tempMaskedList)
	#print(len(tempMaskedList))

print uniqueGT
feature.sort(key = lambda x: x[15])

labels_file = tumor_name + '/' + tumor_name +'_'+file_name+'_total_features_list.txt'
with open(labels_file,'w') as f_labels_file:
	f_labels_file.write('GT_AD1\tGT_ADJAF1\tGT_AF1\tGT_ALD1\tGT_AO1\tGT_BIAS1\tGT_BQ1\tGT_DP1\tGT_DP41\tGT_FREQ1\tGT_GQ1\tGT_HIAF1\tGT_MQ1\tGT_NM1\tGT_ODDRATIO1\tGT_PL1\tGT_PMEAN1\tGT_PSTD1\tGT_QA1\tGT_QR1\tGT_QSTD1\tGT_QUAL1\tGT_RD1\tGT_RO1\tGT_SBF1\tGT_SN1\tGT_SS1\tGT_VD1\tGT_AD2\tGT_ADJAF2\tGT_AF2\tGT_ALD2\tGT_AO2\tGT_BIAS2\tGT_BQ2\tGT_DP2\tGT_DP42\tGT_FREQ2\tGT_GQ2\tGT_HIAF2\tGT_MQ2\tGT_NM2\tGT_ODDRATIO2\tGT_PL2\tGT_PMEAN2\tGT_PSTD2\tGT_QA2\tGT_QR2\tGT_QSTD2\tGT_QUAL2\tGT_RD2\tGT_RO2\tGT_SBF2\tGT_SN2\tGT_SS2\tGT_VD2\tQUAL\tAB\tABP\tAC\tAF\tAN\tAO\tBaseQRankSum\tDP\tDPB\tDPRA\tEPP\tEPPR\tFS\tGC\tGPV\tGTI\tHRun\tHaplotypeScore\tLEN\tMEANALT\tMQ\tMQ0\tMQM\tMQMR\tMQRanksum\tMSI\tMSILEN\tNS\tNUMALT\tODDS\tPAIRED\tPAIREDR\tPAO\tPQA\tPQR\tPRO\tQA\tQD\tQR\tRO\tRPL\tRPP\tRPPR\tRPR\tRUN\tReadPosRakSum\tSAF\tSAP\tSAR\tSOR\tSPV\tSRF\tSS\tSSC\tSSF\tchrom\tloc\n')
	for i in range(len(feature)):
		for j in range(len(feature[0])):
			f_labels_file.write('{}\t'.format(feature[i][j]))
		f_labels_file.write('\n')

#featureArr2 = np.array(feature)[:,:-2]
#featureArr = np.array(featureArr2, dtype=np.float)

maskedArr = np.array(maskedList, dtype=np.bool)

mask_file = tumor_name + '/' + tumor_name +'_'+file_name+'_total_features_mask.txt'
np.savetxt(mask_file,maskedArr,delimiter = '\t', fmt='%d')

#print('type of feature array:{}'.format(type(featureArr)))
#print('shape of feature array:{}'.format(featureArr.shape))
#print('shape of masked array:{}'.format(maskedArr.shape))

#maskedFeature = np.ma.masked_array(featureArr, mask=maskedArr)

#featureMin = np.amin(maskedFeature, axis=0)
#featureMax = np.amax(maskedFeature, axis=0)
#print('feature min:{}'.format(featureMin))
#print('feature max:{}'.format(featureMax))

print('Missing counts:{}\n Total records: {}'.format(missingCount,len(feature)))














