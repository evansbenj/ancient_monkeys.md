# Ambiguous sections

I removed two regions of ambiguous homology: a portion of the non-genic region, positions 7697-7761 in KJ567054.1, and a portion of the d-loop, positions 15686-15741 in KJ567054.1

# Iqtree
I am using iqtree version 1.4.o to make an ML phylogeny and also to select a model of evolution. This model will be used for the time calibrated MrBayes analysis. I am in this directory on info:
```
/2/scratch/ben/2019_ancient_macaque_mtDNA/iqtree
```

```
iqtree -s All_data_align_for_GenBank_ambig_removed.nex -m TEST -nt 1 -pre All_data_align_for_GenBank_ambig_removed.nex_
```
```
iqtree -s All_data_align_for_GenBank_ambig_removed.nex -m TIM3+I+G4 -bb 1000
```


# BEAST

I'm using BEAST version 2 with the clade Age package as described here: https://taming-the-beast.org/tutorials/CladeAge-Tutorial/

I set the model that was selected by iqtree (TN+F+I+G4), which uses emprirical base frequencies. I generated a nj starting tree and scaled it by a factor of 25 in order for it to have the right dimensions (based on trials with Figtree). I did not use the carbon dates for the three LiangBua samples because they were all so recent (D1647 and D1541 were ~1600 years old and D1638 was ~200 years old). 

For the age of the macaques, I used 5.3-5.9 my based on 
* Alba DM, Delson E, Carnevale G, Colobero S, Delfino M, Giuntelli P, et al.
First joint record of Mesopithecus and cf. Macaca in the Miocene of Europe.
J Hum Evol. 2014;67:1–18.

For the age of the split between Theropithecus and Papio, I used 6.5-3.5 my based on:
* Leakey MG. Evolution of Theropithecus in the Turkana Basin. In: Jablonski NG, editor. Theropithecus, the Rise and Fall of a Primate Genus. Cambridge: Cambridge University Press; 1993. p. 85–124.
* Delson E. Cercopithecinae. In: Delson E, Tattersall I, Van Couvering JA,
Brooks AS, editors. Encyclopedia of Human Evolution and Prehistory. New
York: Garland Publishing Inc; 2000. p. 166–71.

I used a relaxed lon normal clock and estimated the clock rate. I used a Birth-Death model for the tree. I ran 100 million generations, sampling every 10,000 generations.

# summarizing trees and logs

I used tracer to check the posterior distributions.  ESS values still too low so I need to do some more runs

in info I am in this directory:
```
/2/scratch/ben/2019_ancient_macaque_mtDNA/beast/bin
```
I summarize the logs like this (for TN93):
```
./bin/logcombiner -burnin 45 -log TN93/run1/All_data_align_for_GenBank_ambig_removed.log -log TN93/run3/All_data_align_for_GenBank_ambig_removed.log -log TN93/run4/All_data_align_for_GenBank_ambig_removed.log -log TN93/run5/All_data_align_for_GenBank_ambig_removed.log -log TN93/run6/All_data_align_for_GenBank_ambig_removed.log -log TN93/run7/All_data_align_for_GenBank_ambig_removed.log -o combined.log

```
and (for GTR):
```
./bin/logcombiner -burnin 55 -log GTR/run1/All_data_align_for_GenBank_ambig_removed.log -log GTR/run2/All_data_align_for_GenBank_ambig_removed.log -log GTR/run3/All_data_align_for_GenBank_ambig_removed.log -log GTR/run4/All_data_align_for_GenBank_ambig_removed.log -log GTR/run5/All_data_align_for_GenBank_ambig_removed.log -log GTR/run6/All_data_align_for_GenBank_ambig_removed.log -log GTR/run7/All_data_align_for_GenBank_ambig_removed.log -log GTR/run9/All_data_align_for_GenBank_ambig_removed.log -log GTR/run10/All_data_align_for_GenBank_ambig_removed.log -o combined.log
```

and trees like this  (for TN93):

```

./bin/logcombiner -burnin 45 -log TN93/run1/All_data_align_for_GenBank_ambig_removed.trees -log  TN93/run3/All_data_align_for_GenBank_ambig_removed.trees -log TN93/run4/All_data_align_for_GenBank_ambig_removed.trees -log TN93/run5/All_data_align_for_GenBank_ambig_removed.trees -log TN93/run6/All_data_align_for_GenBank_ambig_removed.trees -log TN93/run7/All_data_align_for_GenBank_ambig_removed.trees -o combined.trees
```

and (for GTR):
```
./bin/logcombiner -burnin 55 -log GTR/run1/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run2/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run3/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run4/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run5/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run6/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run7/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run9/All_data_align_for_GenBank_ambig_removed.trees -log GTR/run10/All_data_align_for_GenBank_ambig_removed.trees -o combined.trees

```
I summarized the trees like this (for TN, after fixing cyclopsis and radiata):
```
./bin/treeannotator combined.trees TN_contree.tre
```
or this (for GTR):
```
./bin/treeannotator combined.trees GTR_contree.tre
```
