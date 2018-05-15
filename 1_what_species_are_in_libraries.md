# Using MPI sediment scripts to assess species composition of libraries

Because the samples are from old museum specimens or ancient specimens, it is valuable to check what species the libraries tend to match.  To accomplish this, I'm using some scripts that are at the MPI that assess what mammalian family reads that map to a set of concatenated mtDNA genomes from mammals actually are. 

There are two shotgun runs but I am not working with those now.
Plate 79 was shotgun sequenced on lane
```
180322_D00829_0128_lane5
```
Plate 80 was shotgun sequenced on lane
```
180322_D00829_0128_lane4
```

There are two capture runs - one is completed and the other is not.

The first step, as always here, is to merge the reads.  This has been done.I completed the processing for lane 4; a summary and details can be found here:
```
/home/mmeyer/solexa_data/180322_D00829_0128_lane4_SimaPhosphate_Macaques_gelex_others_WGS
```
Matthias completed the processing for lane 4; a summary and details can be found here:
```
/home/mmeyer/solexa_data/180322_D00829_0128_lane5_Macaques_LukasUSel_others_WGS
```
The macaque libraries were then capture with mtDNA probes (capture pools are G5924 and G5925).
```
G5925 was sequenced on MiSeq 180502
```

A next step for the sediment pipeline is to map the reads to a concatenated mammalian mtDNA db.  This can be done by requesting a job on this link:
```
http://bioweb/~sbsuser/webform/form.php
Launch the lauch
Choose ssDNA library (5nt deletion)
Choose second overlap merging option
Tick this box: Do not fail reads due to poor sample assignment
Copy paste the first three columns of the index file sheet, change header as suggested on the web page
Request mapping to TomiMtGENOMES with ancient parameters
```
This has been done already for `G5925` and still needs to be done for `G5924` when it is ready.  It generates a bam file that has mapped and unmapped reads to the `TomiMtGENOMES` database.

This bam file is here:
```
/mnt/ngs_data//180502_M02279_0260_000000000-BPCW5_BN_G5925/Bustard/BWA/proc1/s_1_sequence_ancient_TomiMtGENOMES.bam 
```


Then this command takes the mapped reads and blasts them against individual mtDNA genomes and identifies which one they match best.
```
snakemake --jobs 48 -s ~frederic_romagne/MetaGen/metagen.p1.SnakeFile all --config bamfile=/mnt/ngs_data//180502_M02279_0260_000000000-BPCW5_BN_G5925/Bustard/BWA/proc1/s_1_sequence_ancient_TomiMtGENOMES.bam byfile=/mnt/scratch/ben_evans/ancient_macaques/analyzed_runs/180502_G5925_Macaca/180502_G5925_Macaca.txt && snakemake --jobs 48 -s ~frederic_romagne/MetaGen/metagen.p2.SnakeFile && snakemake --jobs 48 -s ~frederic_romagne/MetaGen/metagen.p3.SnakeFile
```
This was executed in this folder:
```
/mnt/scratch/ben_evans/ancient_macaques/analyzed_runs/180502_G5925_Macaca
```

and the commandline references a file called `180502_G5925_Macaca.txt` which has index information.  This was sent to me by Matthias.


then this command queries these data bases; probably have to update the -db flag:
```
/home/mmeyer/perlscripts/solexa/analysis/sediment_summary.pl -db /mnt/sediments/sediment_database.txt -mammal -query "180424" -reprocess
 ```
