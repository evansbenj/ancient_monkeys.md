# Making mtDNA consensus

I first ran the sediment pipeline to get reads that map preferentially to Cercopithicines over other mammals.  These reads are in this directory:
```
/mnt/scratch/ben_evans/ancient_macaques/analyzed_runs/180502_G5925_Macaca/out/blast/Cercopithecidae
```
I backed this up on graham (safely) here:
```
/home/ben/projects/rrg-ben/ben/ancient_macaques/analyzed_runs/180502_G5925_Macaca/out/blast/Cercopithecidae
```
and here
```
/home/ben/projects/rrg-ben/ben/ancient_macaques/analyzed_runs/Senckenburg_monkeyz/out/blast/Cercopithecidae
```

The Cercopithicines mtDNA genome in the sediment db is from a baboon. Now I'll map these reads to a M. fascicularis reference to generate whole genome assemblies.  I found out from Matthias that the reads also get mapped to 800 different mtDNA genomes, including M. fascicularis.  Nonetheless I redid the mapping to a mtDNA genome from Mauritus, which is more closely related to the Java+ Lesser Sunda Samples. For the nemestrina samples i am mapping to PM664 from Borneo, even though the names of the files incorrectly say PM655. This will be true of all nemestrina samples.  The mtDNA references are here:
```
/mnt/scratch/ben_evans/ancient_macaques/Mfasc_Mauritus_mtDNA_genome.fasta
```
```
/mnt/scratch/ben_evans/ancient_macaques/Mnem_PM664_mtDNAg_genome.fasta
```

Lineage assignment: Probably not worth it.  The only samples that have some low coverage but not ridiculously low are one M. fascicularis from East Java (A11226), one from west Flores (A11244) and we have two other complete genomes from this locality, and two M. nemestrina from Borneo.  The M. nemestrina may be helpful, especially if we add the mtDNA genomes from the HiSeqX, which I am inclined to do.  I need to get the sample ID information from Marie still to assess geographic overlap for those samples... 

# map 
bwa bam2bam -n 0.01 -o 2 -l 16500 -g /mnt/scratch/ben_evans/ancient_macaques/Mfasc_Mauritus_mtDNA_genome.fasta -f A11272_extractedReads-Cercopithecidae_mapped_to_MfascMAURITUS.bam A11272_extractedReads-Cercopithecidae.bam
# sort
samtools sort A11272_extractedReads-Cercopithecidae_mapped_to_MfascMAURITUS.bam -o A11272_extractedReads-Cercopithecidae_mapped_to_MfascMAURITUS_sorted.bam
# Remove unmapped, non-merged, filter-flagged sequences, remove duplicates, create summary statistic
/home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl -qual 25 -paired A11272_extractedReads-Cercopithecidae_mapped_to_MfascMAURITUS_sorted.bam
# consensus caller
/home/mmeyer/perlscripts/solexa/analysis/consensus_from_bam.pl -ref /mnt/scratch/ben_evans/ancient_macaques/Mfasc_Mauritus_mtDNA_genome.fasta A11272_extractedReads-Cercopithecidae_mapped_to_MfascMAURITUS_sorted.uniq.L35MQ25.bam


# Update for Senkenburg

I had to merge some bam files that were run with different library IDs.  In order to get the consensus calling to work properly, I had to replace the headers and then to the same pipeline as above:

```
samtools view -H G10633_G10639_extractedReads-Cercopithecidae.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' |  samtools reheader - G10632_G10638_extractedReads-Cercopithecidae.bam > G10633_G10639_extractedReads-Cercopithecidae_rg.bam

bwa bam2bam -n 0.01 -o 2 -l 16500 -g /mnt/scratch/ben_evans/ancient_macaques/Mfasc_Mauritus_mtDNA_genome.fasta -f G10633_G10639_extractedReads-Cercopithecidae_rg_mapped_to_MfascMAURITUS.bam G10633_G10639_extractedReads-Cercopithecidae_rg.bam

samtools sort G10633_G10639_extractedReads-Cercopithecidae_rg_mapped_to_MfascMAURITUS.bam -o G10633_G10639_extractedReads-Cercopithecidae_rg_mapped_to_MfascMAURITUS_sorted.bam

/home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl -qual 25 -paired G10633_G10639_extractedReads-Cercopithecidae_rg_mapped_to_MfascMAURITUS_sorted.bam

/home/mmeyer/perlscripts/solexa/analysis/consensus_from_bam.pl -ref /mnt/scratch/ben_evans/ancient_macaques/Mfasc_Mauritus_mtDNA_genome.fasta G10633_G10639_extractedReads-Cercopithecidae_rg_mapped_to_MfascMAURITUS_sorted.uniq.L35MQ25.bam

```
