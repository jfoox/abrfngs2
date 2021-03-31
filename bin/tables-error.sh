### get bbmap mhist outputs
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/error/outs/mhists
for i in $(ls *.txt | grep -v PromethION | grep -v PacBio | grep -v Genapsys); do 
  sample=$(echo $i | cut -d'.' -f2)
  gc=$(echo $i | cut -d'.' -f5 | sed -e "s/to/,/")
  tail -n +2 $i | cut -f1-7 | sed -r "s/\t+/,/g" | sed -e "s/^/${sample},${gc},/g"
done | sed -e "1s/^/sample,GCmin,GCmax,base,Match1,Sub1,Del1,Ins1,N1,Other1\n/" > ../../../tables/mhist_v2.csv

for i in $(ls *Genapsys*); do
  sample=$(echo $i | cut -d'.' -f2)
  gc=$(echo $i | cut -d'.' -f3 | sed -e "s/to/,/")
  tail -n +2 $i | cut -f1-7 | sed -r "s/\t+/,/g" | sed -e "s/^/${sample},${gc},/g"
done >> ../../../tables/mhist_v2.csv

for i in $(ls *Pro* | grep -v sorted | grep -v calmd | grep -v forbbmap); do 
  sample=$(echo $i | cut -d'.' -f2)
  gc=$(echo $i | cut -d'.' -f3 | sed -e "s/to/,/")
  tail -n +2 $i | cut -f1-7 | sed -r "s/\t+/,/g" | sed -e "s/^/${sample},${gc},/g"
done >> ../../../tables/mhist_v2.csv

for i in $(ls *Pac*ds25*); do 
  sample=$(echo $i | cut -d'.' -f2)
  gc=$(echo $i | cut -d'.' -f4 | sed -e "s/to/,/")
  tail -n +2 $i | cut -f1-7 | sed -r "s/\t+/,/g" | sed -e "s/^/${sample},${gc},/g"
done >> ../../../tables/mhist_v2.csv


### get samstats error in UCSC contexts
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/samstats/outs/ds_ucsc
for i in *.samstats; do
  sample=$(echo $i | cut -d'.' -f1)
  context=$(echo $i | rev | cut -d'.' -f2 | rev)
  rate=$(grep -m1 'error rate' $i | cut -f3)
  echo ${sample},${context},${rate}
done | sed -e "1s/^/sample,context,error\n/" > /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/mismatch_UCSC.csv
