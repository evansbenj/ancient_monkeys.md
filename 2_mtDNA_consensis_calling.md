# Making mtDNA consensus

I first ran the sediment pipeline to get reads that map preferentially to Cercopithicines over other mammals.  These reads are in this directory:
```
/mnt/scratch/ben_evans/ancient_macaques/analyzed_runs/180502_G5925_Macaca/out/blast/Cercopithecidae
```
The Cercopithicines mtDNA genome in the sediment db is from a baboon. Now I'll map these reads to a M. fascicularis reference to generate whole genome assemblies.

# map 
bwa bam2bam -n 0.01 -o 2 -l 16500 -g /mnt/scratch/ben_evans/ancient_frogz/xenTr_MT.fasta -f R7931_mapped_to_XT.bam R7931.bam

# sort
samtools sort R7931_mapped_to_XT.bam -o R7931_mapped_to_XT_sorted.bam

# index
samtools index R7931_mapped_to_XT_sorted.bam

# Remove unmapped, non-merged, filter-flagged sequences, remove duplicates, create summary statistic
/home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl -qual 25 -paired R7931_mapped_to_XT_sorted.bam 

# consensus caller
/home/mmeyer/perlscripts/solexa/analysis/consensus_from_bam.pl -ref /mnt/scratch/ben_evans/ancient_frogz/xenTr_MT.fasta R7931_mapped_to_XT_sorted.uniq.L35MQ25.bam


This was sufficient using Matthias' scripts:
* bwa bam2bam -n 0.01 -o 2 -l 16500 -g /mnt/scratch/ben_evans/ancient_frogz/xenXL_MT.fasta -f R7935_mapped_to_XL.bam R7935.bam
* samtools sort R7935_mapped_to_XL.bam -o R7935_mapped_to_XL_sorted.bam
* /home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl -qual 25 -paired R7935_mapped_to_XL_sorted.bam 
* /home/mmeyer/perlscripts/solexa/analysis/consensus_from_bam.pl -ref /mnt/scratch/ben_evans/ancient_frogz/xenXL_MT.fasta R7935_mapped_to_XL_sorted.uniq.L35MQ25.bam
