# Gblocks

I could have used gblocks here (or locally) but inspection of thew alignment indicates that there are no regions of ambiguous homology, even between Papio/Theropithecus and macaques.
```
http://phylogeny.lirmm.fr/phylo_cgi/one_task.cgi?task_type=gblocks
```

# Iqtree
I am using iqtree to make an ML phylogeny and also to select a model of evolution. This model will be used for the time calibrated MrBayes analysis. I am in this directory on sharcnet:
```
/work/ben/iqtree-1.6.8-Linux/bin/Ancient_macaque_mtDNA
```

```
../iqtree -s Liedigk_plus_new_genomez_align8_for_analysis_nogaps.nexus -m TEST -nt 1 -pre Liedigk_plus_new_genomez_align8_for_analysis_nogaps.nexus_
```
```
../iqtree -s Liedigk_plus_new_genomez_align8_for_analysis_nogaps.nexus -m TN+F+I+G4 -bb 1000```
```

I also did the analysis without the three weird Sulawesi mtDNA genomes from NCBI that had long branch lengths (despite the correct phylogenetic placement

```
../iqtree -s Liedigk_plus_new_genomez_align8_for_analysis_nogaps_noweirdsula.nexus -m TEST -nt 1 -pre Liedigk_plus_new_genomez_align8_for_analysis_nogaps_noweirdsula.nexus_
```
```
../iqtree -s Liedigk_plus_new_genomez_align8_for_analysis_nogaps_noweirdsula.nexus -m TN+F+I+G4 -bb 1000
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
I summarize the trees like this:
```
./treeannotator -b 40 Liedigk_plus_new_genomez_align8_for_analysis_nogaps_noweirdsula.trees Liedigk_plus_new_genomez_align8_for_analysis_nogaps_noweirdsula_con
```
