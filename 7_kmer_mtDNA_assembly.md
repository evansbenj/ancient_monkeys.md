# De novo assembly using kmers

One concern is that the gene conversion tracks may have low coverage.  I am going to use RepBase to assemble mtDNA using high count kmers.

For now I am working in this directory:
```
/scratch/ben/SEAsian_macaques_original_rawdata/maura_PF615/reads_from_bam
```

and using reads in bam files from here:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam
```

and I used this sbatch script (in this directory: `/scratch/ben/SEAsian_macaques_original_rawdata/maura_PF615`)to pull out the reads:
```
#!/bin/sh
#SBATCH --job-name=gatk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=16gb
#SBATCH --output=gatk.%J.out
#SBATCH --error=gatk.%J.err
#SBATCH --account=def-ben

/home/ben/project/ben/bin/samtools-1.9/samtools bam2fq -1 ./reads_from_bam/maura_PF615sorted_ddedup_rg_realigned.bamBS
QR_paired1.fq -2 ./reads_from_bam/maura_PF615sorted_ddedup_rg_realigned.bamBSQR_paired2.fq -0 /dev/null -s /dev/null -
n -F 0x900 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/maura_PF615sorted_ddedup_
rg_realigned.bamBSQR.bam
```