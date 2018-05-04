# Using MPI sediment scripts to assess species composition of libraries

Because the samples are from old museum specimens or ancient specimens, it is valuable to check what species the libraries tend to match.  To accomplish this, I'm using some scripts that are at the MPI that assess what mammalian family reads that map to a set of concatenated mtDNA genomes from mammals actually are. 

There are two shotgun runs but I am not working with those now.
Plate 79 was shotgun sequenced on lane
180322_D00829_0128_lane5

Plate 80 was shotgun sequenced on lane
180322_D00829_0128_lane4


There are two capture runs - one is completed and the other is not.

The first step, as always here, is to merge the reads.  This has been done.I completed the processing for lane 4; a summary and details can be found here:
/home/mmeyer/solexa_data/180322_D00829_0128_lane4_SimaPhosphate_Macaques_gelex_others_WGS

I completed the processing for lane 4; a summary and details can be found here:
/home/mmeyer/solexa_data/180322_D00829_0128_lane5_Macaques_LukasUSel_others_WGS

The macaque libraries were then capture with mtDNA probes (capture pools are G5924 and G5925).

G5925 was sequenced on MiSeq 180502


This command takes the mapped reads and blasts them against individual mtDNA genomes and identifies which one they match best.
```
snakemake --jobs 48 -s ~frederic_romagne/MetaGen/metagen.p1.SnakeFile all --config bamfile=/mnt/ngs_data//180502_M02279_0260_000000000-BPCW5_BN_G5925/Bustard/BWA/proc1/s_1_sequence_ancient_TomiMtGENOMES.bam byfile=/mnt/scratch/ben_evans/ancient_macaques/analyzed_runs/180502_G5925_Macaca/180502_G5925_Macaca.txt && snakemake --jobs 48 -s ~frederic_romagne/MetaGen/metagen.p2.SnakeFile && snakemake --jobs 48 -s ~frederic_romagne/MetaGen/metagen.p3.SnakeFile
```
then this command queries these data bases; probably have to update the -db flag:
```
/home/mmeyer/perlscripts/solexa/analysis/sediment_summary.pl -db /mnt/sediments/sediment_database.txt -mammal
 -query "180424" -reprocess
 ```
