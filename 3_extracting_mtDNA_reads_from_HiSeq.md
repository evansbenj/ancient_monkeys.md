# Add HiSeq macaque data to mtDNA project

To make the mtDNA paper more interesting, it would be useful to add the mtDNA genomes from the 29 HiSeqX bam files I have. 
The BQSR bam files are on goblin here: `/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam`

I just chatted with Janet about strategies to avoid numts.  Several good points were brought up.  First, numts inserted before divergence from rhesus and SEAsian macaques should be pulled out of the mtDNA mapping.  So mainly numts that were inserted afterwards in SEAsian macaques might be problematic.  Also, I can make a de novo assembly, remap the reads, genotype them, and look for heterozygous sites.  This could be due to numts.  But overall this should be fine anyhow because the mtDNA will have MUCH higher coverage than numts.  Woo hoo!

I should be able to get the fastq files that map to the rhesus mtDNA like this (30_samtools_get_mtDNA_reads.pl): 
```
#!/usr/bin/perl
# This script will make commandlines for sharcnet to make bam files with 
# only those reads that map to mtDNA

my $samtoolspath = "/work/ben/2017_SEAsian_macaques/bin/samtools-1.6/";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/*males/";

@files = glob($majorpath."*sorted_ddedup_rg_realigned.bamBSQR.bam");

foreach my $file (@files){
    my $commandline = "sqsub -r 2h --mpp 6G -o ".$file."\.log ";
    $commandline = $commandline.$samtoolspath."samtools view -b ".$file." chrM ";
    $commandline = $commandline."| samtools bam2fq - > ".substr($file,0,-4).".fastq\"";
    print $commandline,"\n";
#    $status = system($commandline);
}
```

