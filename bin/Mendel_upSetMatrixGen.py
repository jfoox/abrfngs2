#! /usr/bin/env python3
## July 13, 2020
## Create UpSet matrices of each Mendelian Violation stratum

import sys
import os

strata = ['SNP', 'INS_5', 'INS_6to15', 'INS_15', 'DEL_5', 'DEL_6to15', 'DEL_15']

samples = ['BGISEQ500_LAB01',
           'HiSeq2500_LAB01', 'HiSeq2500_LAB02', 'HiSeq2500_LAB03',  
           'HiSeq4000_LAB01', 'HiSeq4000_LAB02', 'HiSeq4000_LAB03',  
           'HiSeqX10_LAB01',  'HiSeqX10_LAB02',  'HiSeqX10_LAB03',
           'MGISEQ2000_LAB01','NovaSeq2x150_LAB01', 'NovaSeq2x150_LAB02', 'NovaSeq2x250_LAB02' 
          ]
          
for stratum in strata:
	global_matrix = {}
	with open('forUpset__'+stratum+'.csv', 'w') as outfile:
		outfile.write('chr_pos')
		for sample in samples:
			with open('ViolationRegions__'+sample+'_'+stratum+'.bed', 'r') as infile:
				for line in infile:
					chrpos = line.strip().split('\t')[0] + '_' + line.strip().split('\t')[1]
					if chrpos not in global_matrix:
						global_matrix[chrpos] = {}
					global_matrix[chrpos][sample] = "1"
			outfile.write(','+sample)
		outfile.write('\n')

		for locus in global_matrix:
			outfile.write(locus)
			for sample in samples:
				if sample in global_matrix[locus]:
					outfile.write(',1')
				else:
					outfile.write(',0')
			outfile.write('\n')

