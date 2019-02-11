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
../iqtree -s input.nxs -m TEST -nt 1 -pre input.nxs_
```
```
../iqtree -s input.nxs -m MODEL -bb 1000```
```