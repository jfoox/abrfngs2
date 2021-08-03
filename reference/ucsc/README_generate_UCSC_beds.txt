# get UCSC RepeatMasker tracks from Table Browser:
http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1028261817_PokSqZrEpX6gHaLLd38TEoh0tZIt&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chrX%3A15%2C560%2C138-15%2C602%2C945&hgta_outputType=bed&hgta_outFileName=

# generate tracks
zcat UCSC_RepeatMasker_GRCh38.bed.gz | tail -n+2 | awk '{ if ($13 == "Alu") print $6,$7,$8 }' | tr ' ' '\t' | sort -k1,1V -k2,2n > ucsc.Alu.bed
zcat UCSC_RepeatMasker_GRCh38.bed.gz | tail -n+2 | awk '{ if ($13 == "L1")  print $6,$7,$8 }' | tr ' ' '\t' | sort -k1,1V -k2,2n > ucsc.L1.bed
zcat UCSC_RepeatMasker_GRCh38.bed.gz | tail -n+2 | awk '{ if ($13 == "Low_complexity") print $6,$7,$8 }' | tr ' ' '\t' | sort -k1,1V -k2,2n > ucsc.lowcomplexity.bed
zcat UCSC_RepeatMasker_GRCh38.bed.gz | tail -n+2 | awk '{ if ($12 == "LTR") print $6,$7,$8 }' | tr ' ' '\t' | sort -k1,1V -k2,2n > ucsc.LTR.bed
zcat UCSC_RepeatMasker_GRCh38.bed.gz | tail -n+2 | awk '{ if ($13 == "Satellite") print $6,$7,$8 }' | tr ' ' '\t' | sort -k1,1V -k2,2n > ucsc.satellite.bed
zcat UCSC_RepeatMasker_GRCh38.bed.gz | tail -n+2 | awk '{ if ($12 == "Simple_repeat") print $6,$7,$8 }' | tr ' ' '\t' | sort -k1,1V -k2,2n > ucsc.simplerepeat.bed

