# Senkenburg redo
I realize now that the reference I used for the Senkenburg nemestrina was not great.  So instead I am redoing them using the denovo mtDNA genomes I made from the WGS data.

I am working here on graham:
```
/home/ben/scratch/SEAsian_macaques_original_rawdata/Senckenburg/G10637
```
I ran this sbatch  script (`maps_leipzig_from_bam.sh`):
```
#!/bin/sh
#SBATCH --job-name=leipbam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=16gb
#SBATCH --output=leipbam.%J.out
#SBATCH --error=leipbam.%J.err
#SBATCH --account=def-ben
module load samtools/1.9
samtools sort -n $1 -o aln.qsort
# the aln.qsort.bam file will be overwritten every time I run this
module load bedtools/2.27.1
bedtools bamtofastq -i aln.qsort.bam -fq $1_r1.corfixed.fastq -fq2 $1_r2.corfixed.fastq
# map the reads to the mt genome
module load bwa/0.7.17
bwa mem $2 $1_r1.corfixed.fastq $1_r2.corfixed.fastq | samtools view -Shu - | samtools sort - -o $1_to_$2_sorted.bam
```
using this command:
```
sbatch maps_leipzig_from_bam.sh /home/ben/projects/rrg-ben/ben/ancient_macaques/analyzed_runs/Senckenburg_monkeyz/out/blast/Cercopithecidae/G10637_extractedReads-Cercopithecidae.bam ../Circularized_assembly_1_nsang_NOVO.fasta
```
Hopefully the result will be a well mapped genome.  I may also try to do NOVOplasty on these reads to denovo assemble them.
