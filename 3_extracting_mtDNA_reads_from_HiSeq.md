# Add HiSeq macaque data to mtDNA project

To make the mtDNA paper more interesting, it would be useful to add the mtDNA genomes from the 29 HiSeqX bam files I have. 
The BQSR bam files are on goblin here: `/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam`

I just chatted with Janet about strategies to avoid numts.  Several good points were brought up.  First, numts inserted before divergence from rhesus and SEAsian macaques should be pulled out of the mtDNA mapping.  So mainly numts that were inserted afterwards in SEAsian macaques might be problematic.  Also, I can make a de novo assembly, remap the reads, genotype them, and look for heterozygous sites.  This could be due to numts.  But overall this should be fine anyhow because the mtDNA will have MUCH higher coverage than numts.  Woo hoo!

I can get a bam files with only mtDNA like this (30_samtools_get_mtDNA_reads.pl): 
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
    $commandline = $commandline." -o ".substr($file,0,-4)."_mtDNA_only.bam";
#    $commandline = $commandline."| samtools bam2fq - > ".substr($file,0,-4)."_mtDNA_only.fastq";
    print $commandline,"\n";
#    $status = system($commandline);
}
```
Then I then can use picard to clean up the paired sequences for downsteam assembly and remove mapping information (31.5_picard_get_mtDNA_reads.pl).
```
#!/usr/bin/perl
# This script will make commandlines for sharcnet to make bam files with 
# only those reads that map to mtDNA

my $picardpath = "/work/ben/2017_SEAsian_macaques/bin/picard/";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/*males/";

@files = glob($majorpath."*sorted_ddedup_rg_realigned.bamBSQR_mtDNA_only.bam");

foreach my $file (@files){
    my $commandline = "sqsub -r 2h --mpp 6G -o ".$file."\.log ";
    $commandline = $commandline."/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -jar ".$picardpath."picard.jar RevertSam ";
    $commandline = $commandline."I= ".$file." ";
    $commandline = $commandline."O= ".substr($file,0,-4)."_picard.bam ";
    $commandline = $commandline."SANITIZE=true ";
    $commandline = $commandline."SORT_ORDER=queryname ";
    $commandline = $commandline."RESTORE_ORIGINAL_QUALITIES=true ";
    $commandline = $commandline."REMOVE_ALIGNMENT_INFORMATION=true";
    print $commandline,"\n";
#    $status = system($commandline);
}
```
Then split up the reads into separate forward and reverse files
```
grep -A3 -P "@*/1" --no-group-separator /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR_mtDNA_only_picard.fastq >/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR_mtDNA_only_picard_1.fastq
```
```
grep -A3 -P "@*/2" --no-group-separator /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR_mtDNA_only_picard.fastq > /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR_mtDNA_only_picard_2.fastq
```


I installed `trinity` on goblin, which was a bit of a pain because I needed to have java version 1.8 in the $PATH and I also needed jellyfish to be installed.  I used trinity version 2.25 instead of the most recent because the latter required salmon, which did not install properly.  This script (32_trinity_assmple_mtDNA.pl) will execute commands to do assemblies:
```
#!/usr/bin/perl
# This script will make commandlines for sharcnet to make bam files with 
# only those reads that map to mtDNA

my $trinitypath = "/work/ben/2017_SEAsian_macaques/bin/trinityrnaseq-Trinity-v2.5.0/";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/*males/";

@files = glob($majorpath."*fastq");

foreach my $file (@files){
    my $commandline = "sqsub -r 2h --mpp 6G -o ".$file."\.log ";
    $commandline = $commandline.$trinitypath."Trinity --seqType fq --single ".$file." --no_normalize_reads --max_memory 10G --KMER_SIZE 29 --no_bowtie";
    print $commandline,"\n";
#    $status = system($commandline);
}
```
