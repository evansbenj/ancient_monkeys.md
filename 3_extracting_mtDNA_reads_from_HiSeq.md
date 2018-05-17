# Add HiSeq macaque data to mtDNA project

To make the mtDNA paper more interesting, it would be useful to add the mtDNA genomes from the 29 HiSeqX bam files I have.  I should be able to get the fastq files that map to the rhesus mtDNA like this: 
```
samtools view -b input.bam "chrM" > output.bam
```

The BQSR bam files are on goblin here: `/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam`

I just chatted with Janet about strategies to avoid numts.  Several good points were brought up.  First, numts inserted before divergence from rhesus and SEAsian macaques should be pulled out of the mtDNA mapping.  So mainly numts that were inserted afterwards in SEAsian macaques might be problematic.  Also, I can make a de novo assembly, remap the reads, genotype them, and look for heterozygous sites.  This could be due to numts.  But overall this should be fine anyhow because the mtDNA will have MUCH higher coverage than numts.  Woo hoo!
