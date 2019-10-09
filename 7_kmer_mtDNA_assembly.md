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

# Count kmers
In this directory: 
```
/scratch/ben/SEAsian_macaques_original_rawdata/maura_PF615/
```
After copying over the RepPark script and modifying the directories to be "", I am running RepPark using this sbatch command, which loads jellyfish and trinity modules and hopefully works.  The -t sets the threshold count to 80, so we should be retaining only high coverage kmers:

```
#!/bin/sh
#SBATCH --job-name=reppark
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=3gb
#SBATCH --output=reppark.%J.out
#SBATCH --error=reppark.%J.err
#SBATCH --account=def-ben

module load rsem/1.3.0
module load samtools
#module load nixpkgs/16.09
#module load intel/2016.4
#module load bowtie2/2.3.4.1
module load nixpkgs/16.09 gcc/5.4.0 openmpi/2.1.1
module load salmon/0.9.1
module load transdecoder/3.0.1
module load jellyfish/2.2.6
module load trinity/2.6.5
module load velvet/1.2.10
./RepARK.pl -l $1 -l $2 -k 31 -t 80 -o $3
```
Here is an example of the command:
```
sbatch RepPark_sbatch.sh /scratch/ben/SEAsian_macaques_original_rawdata/maura_PF615/reads_from_bam/maura_PF615sorted_ddedup_rg_realigned.bamBSQR_paired1.fq /scratch/ben/SEAsian_macaques_original_rawdata/maura_PF615/reads_from_bam/maura_PF615sorted_ddedup_rg_realigned.bamBSQR_paired2.fq maura_PF615_kmer_31
```

# Making a blast db out of the kmer contigs
```
module load nixpkgs/16.09
module load gcc/7.3.0
module load blast+/2.7.1
```
make a blast db in this directory:
```
/scratch/ben/SEAsian_macaques_original_rawdata/maura_PF615/maura_PF615_kmer_31/velvet_repeat_lib
```
like this:
```
makeblastdb -in contigs.fa -dbtype nucl -out contigs.fa_blastable
```

now blast a mtDNA genome against this blast database like this:
```
blastn -query maura_PF615_mapped_to_rhesus.fa -db contigs.fa_blastable -outfmt 6 -out maura_PF615_mapped_to_rhesus_to_highabundancecontigs.out
```
I can check the coverage of the high abundance contigs like this:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa | sort -rn | head -n 10
```

This produced nothing.  Maybe there is a problem with using the reads from the bam file?  I will try doing this now with the raw data, but it will take a while to upload it probably.


