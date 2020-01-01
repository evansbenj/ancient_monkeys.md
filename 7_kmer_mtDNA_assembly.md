# Coverage per site in WGS mtDNA

After identifying evidence of gene conversion in mtDNA genomes mapped to the rhesus reference, I checked the coverage per site using samtools:
```
samtools depth -a -d 0 -r chrM:1-16564 maura_PF615sorted_ddedup_rg_realigned.bamBSQR.bam > maura_PF615sorted_ddedup_rg_realigned.bamBSQR.depth_per_base
```
In this individual, 318 sites had zero coverage and probably got the reference sequence as a result. This is a huge problem and probably explains why there was evidence of gene conversion. To remedy this, I have two ideas.  First, I can try to assemble the mtDNA from high coverage kmers extracted from the bam file or the original trimmed reads. The former may be bad if the unmapped reads are not in the bam file.  I need to look into that. And second, I could try to use a more closely related reference sequence.  There are some Sulawesi macaque sequences on Genbank but I think they also have the same issue with reference sequence mixed in with read sequences.  So maybe the kmer approach will work?  See below for more.

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

This produced nothing.  Maybe there is a problem with using the reads from the bam file?  I will try doing this now with the raw data. When I did this, there were many mtDNA hits. But the problem is how do I distinguish them from numts?  Also the real mtDNA hits seem to be fragmented (i.e. there isn't one contig that is ~16,000 bp long).

Amazingly, it appears that I have 23 of 29 individually trimmed data on graham (safetly) here:
```
/home/ben/projects/rrg-ben/ben/SEAsian_macaques_rawdata_MPIexpressions/
```

So I should be able to run this again like this:
```
sbatch RepPark_sbatch.sh /home/ben/projects/rrg-ben/ben/SEAsian_macaques_rawdata_MPIexpressions/maura_PF615/PF615_all_R1scythe_and_trimm_paired.cor.fastq.gz /home/ben/projects/rrg-ben/ben/SEAsian_macaques_rawdata_MPIexpressions/maura_PF615/PF615_all_R2scythe_and_trimm_paired.cor.fastq.gz maura_PF615_kmer_31
```

# Extracting fastq files from bam

I did not use this one because it does not extract unmapped reads (bedtools.sh):
```
#!/bin/sh
#SBATCH --job-name=bedtools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=16gb
#SBATCH --output=bedtools.%J.out
#SBATCH --error=bedtools.%J.err
#SBATCH --account=def-ben
module load samtools/1.9
samtools sort -n $1 -o aln.qsort
# the aln.qsort.bam file will be overwritten every time I run this
module load bedtools/2.27.1
bedtools bamtofastq -i aln.qsort.bam -fq bru_PF707.end1.fq -fq2 bru_PF707.end2.fq
```

But this one does (fastqs_from_bam.sh):
```
#!/bin/sh
#SBATCH --job-name=fastqs_from_bam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mem=16gb
#SBATCH --output=fastqs_from_bam.%J.out
#SBATCH --error=fastqs_from_bam.%J.err
#SBATCH --account=def-ben

# this script should extract mapped and unmapped reads from a bam file
# and generate perfect paired fastq files
# load samtools
module load samtools/1.9
# directions are here
# https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd 
# R1 mapped, R2 mapped
 samtools view -u -f 1 -F 12 $1 > lib_002_map_map.bam
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 $1 > lib_002_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 $1 > lib_002_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 $1 > lib_002_unmap_unmap.bam

# merge
samtools merge -u lib_002_unmapped.bam lib_002_unmap_map.bam lib_002_map_unmap.bam lib_002_unmap_unmap.bam

# sort
samtools sort -n lib_002_map_map.bam -o lib_002_mapped.sort.bam
samtools sort -n lib_002_unmapped.bam -o lib_002_unmapped.sort.bam

# get some summary stats
samtools flagstat lib_002.sorted.md.bam
samtools view -c lib_002_mapped.sort.bam
samtools view -c lib_002_unmapped.sort.bam


module load bedtools/2.27.1
#extract the FASTQ reads into two paired read files
bamToFastq -i lib_002_mapped.sort.bam -fq lib_002_mapped.1.fastq -fq2 lib_002_mapped.2.fastq 
bamToFastq -i lib_002_unmapped.sort.bam -fq lib_002_unmapped.1.fastq -fq2 lib_002_unmapped.2.fastq 

# combine both the first and paired reads together from the mapped and unmapped files
cat lib_002_mapped.1.fastq lib_002_unmapped.1.fastq > lib_002.1.fastq
cat lib_002_mapped.2.fastq lib_002_unmapped.2.fastq > lib_002.2.fastq

# You might need to run [BBTools](https://sourceforge.net/projects/bbmap/) [Repair.sh](https://jgi.doe.gov/data-and-tools/b
btools/bb-tools-user-guide/repair-guide/) to repair the read1 and read2 pairing

module load bbmap/37.36
repair.sh in=lib_002.1.fastq in2=lib_002.2.fastq out=lib_002.1.fastq.gz out2=lib_002.2.fastq.gz repair
```
I executed the above on graham like this:
```
sbatch fastqs_from_bam.sh /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/nigrescens_PM654sorted_ddedup_rg_realigned.bamBSQR.bam
```


# Starting to look better

When I used the raw trimmed data to feed into RepArk, the resulting kmer contigs were very fragmented but matched pretty much the entire mtDNA.  I think some of these are probably numts.

Next step is to try just mapping the raw data to the consensus seq from the map to rhemac mtDNA.  I think this might work.

Here are some commands I used to get the highest coverage and longest kmer contigs:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa | sort -rn | head -n 10
grep -o -P '(?<=length_).*(?=)' contigs.fa | sort -rn | head -n 10
grep '1666.413818' contigs.fa
awk -v seq="NODE_693971_length_29_cov_1666.413818" -v RS='>' '$1 == seq {print RS $0}' contigs.fa
```

# De novo assembly using abyss

Another approach is to do a de novo assembly and hope that the entire mtDNA is assembled.  I useed the trimmeeed raw data and abyss like this (for maura_PF615 in this directory: /scratch/ben/SEAsian_macaques_original_rawdata/maura_PF615):
```
#!/bin/sh
#SBATCH --job-name=abyss
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=96:00:00
#SBATCH --mem=256gb
#SBATCH --output=abyss.%J.out
#SBATCH --error=abyss.%J.err
#SBATCH --account=def-ben

module load abyss/2.0.2

abyss-pe np=8 name=$1 k=64 in="/home/ben/projects/rrg-ben/ben/SEAsian_macaques_rawdata_MPIexpressions/maura_PF615/PF615_
all_R1scythe_and_trimm_paired.cor.fastq.gz /home/ben/projects/rrg-ben/ben/SEAsian_macaques_rawdata_MPIexpressions/maura_
PF615/PF615_all_R2scythe_and_trimm_paired.cor.fastq.gz"
```
This seems to be working ok. Hopefully a 4 day timee wall will suffice.

# Norgal

I am also trying another approach that tries to assemble mtDNA genomes using high coverage kmers.  

```
#!/bin/sh                                                                                                                                          
#SBATCH --job-name=norgal                                                                                                               #SBATCH --nodes=1                                                                                                                       #SBATCH --ntasks-per-node=8                                                                                                             
#SBATCH --time=12:00:00                                                                                                                 
#SBATCH --mem=128gb                                                                                                                     
#SBATCH --output=norgal.%J.out                                                                                                           
#SBATCH --error=norgal.%J.err                                                                                                           
#SBATCH --account=def-ben                                                                                                               

module load python/2.7.14
module load scipy-stack/2019b

python ../norgal/norgal.py -i /home/ben/projects/rrg-ben/ben/SEAsian_macaques_rawdata_MPIexpressions/hecki_PF505/PF505_all_R1scythe_and_trimm_paired.cor.fastq.gz /home/ben/projects/rrg-ben/ben/SEAsian_macaques_rawdata_MPIexpressions/hecki_PF505/PF505_all_R2scythe_and_trimm_paired.cor.fastq.gz -t 8 -m 1000 -r hecki_PF505-3.fa -o norgal_output --blast
```
# Making a subset for testing

```
module load seqtk/1.2
seqtk sample -s100 PF505_trim1_repaired.fq 10000000 > PF505_trim1_repaired1000.fq
seqtk sample -s100 PF505_trim2_repaired.fq 10000000 > PF505_trim2_repaired1000.fq
```
