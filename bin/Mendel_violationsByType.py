#! /usr/bin/env python3
## June 24, 2020

import sys
import os

inputfile = sys.argv[1]
samplename = inputfile.replace('_trio.vcf', '')

with open('ViolationRegions__' + samplename + '_SNP.bed', 'w') as outfile_SNP, \
     open('ViolationRegions__' + samplename + '_INS_5.bed', 'w') as outfile_INS_5, \
     open('ViolationRegions__' + samplename + '_INS_6to15.bed', 'w') as outfile_INS_6to15, \
     open('ViolationRegions__' + samplename + '_INS_15.bed', 'w') as outfile_INS_15, \
     open('ViolationRegions__' + samplename + '_DEL_5.bed', 'w') as outfile_DEL_5, \
     open('ViolationRegions__' + samplename + '_DEL_6to15.bed', 'w') as outfile_DEL_6to15, \
     open('ViolationRegions__' + samplename + '_DEL_15.bed', 'w') as outfile_DEL_15:
     with open(inputfile, 'r') as infile:
     	for line in infile:
     		#print(line)
     		if not line.startswith('#') and 'MD=2' in line:
     			parts = line.strip().split('\t')
     			chr = parts[0]
     			pos = parts[1]
     			site = chr+'\t'+pos+'\t'+pos+'\n'
     			
     			ref = parts[3]
     			refLen = len(ref)
     			
     			alt = parts[4]
     			if ',' in alt:
     				# this takes first genotype from son as "alt";
     				# probably needs to be refined based on how VBT calculates
     				son = int(parts[11].split('/')[0])
     				if son == 0:
     					altLen = refLen
     				else:
	     				altLen = len(alt.split(',')[son-1])
	     		else:
     			    altLen = len(alt)
     			size = altLen - refLen
     			if size > 0 and size < 6:
     				outfile_INS_5.write(site)
     			elif size >= 6 and size < 15:
     				outfile_INS_6to15.write(site)
     			elif size >= 15:
     				outfile_INS_15.write(site)
     			elif size < 0 and size > -6:
     				outfile_DEL_5.write(site)
     			elif size <= -6 and size > -15:
     				outfile_DEL_6to15.write(site)
     			elif size <= -15:
     				outfile_DEL_15.write(site)
     			elif size == 0:
     				outfile_SNP.write(site)
     			else:
     				print('Something funky happened at:\n' + site)
     				sys.exit()
