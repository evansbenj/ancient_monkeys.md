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
