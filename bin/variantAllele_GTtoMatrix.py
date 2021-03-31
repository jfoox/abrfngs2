#! /usr/bin/env python

# convert GT matrix into numeric matrix for heatmap viz

# get files
import sys
import gzip
rawfile  = sys.argv[1]
max_numtaxa_absent = sys.argv[2]
loci = {}   # account for the same locus reported more than once

# execute
with open(rawfile+'.toMatrix', 'w') as outfile:
	with open(rawfile, 'r') as infile:
		for line in infile:
			if line.startswith('Chr'):
				outfile.write(line.replace('Chr\tPos\tRef\tAlt','chr_pos_ref_alt'))
			else:
				if line.count('./.') < int(max_numtaxa_absent):
					splitline = line.strip().split('\t')
					if not ',' in splitline[3]:
						# deal with chr_locus appearing more than once
						locus = splitline[0]+'_'+splitline[1]+'_'+splitline[2]+'_'+splitline[3]
						if locus not in loci:
							loci[locus] = 1
						else:
							loci[locus] += 1
						if loci[locus] == 1:
							outfile.write(locus)
						else:
							outfile.write(locus+'_'+str(loci[locus]))
					
						# export per-sample GT
						for i in splitline[4:]:
							i = i.replace('|', '/')
							if i == '0/0' or i == '1/0':
								outfile.write('\t0.0')
							elif i == '0/1':
								outfile.write('\t0.5')
							elif i == '1/1':
								outfile.write('\t1.0')
							else:
								outfile.write('\tNA')
						outfile.write('\n')


