# Add HiSeq macaque data to mtDNA project

To make the mtDNA paper more interesting, it would be useful to add the mtDNA genomes from the 29 HiSeqX bam files I have.  I should be able to get the fastq files that map to the rhesus mtDNA like this: 
```
samtools view -b input.bam "chrM" > output.bam
```

The BQSR bam files are on goblin here: `/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam`
