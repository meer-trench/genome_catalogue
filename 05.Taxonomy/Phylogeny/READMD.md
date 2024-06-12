# Phylogenetic Analysis for MMC  

## Identification and Multiple Sequence Alignment  
  
Firstly, markers from Species-level Genome Bins (SGBs) were identified and conducted multiple sequence alignments through `hmmalign`.

```bash
# All Species-level Genome Bins (SGBs)
# GTDB-Tk version 2.4.0
# GTDB release 220

gtdbtk identify --batchfile MEER_SGB.list --out_dir bin_ca_tetra_r220 --cpus 64 -x fa
gtdbtk align --identify_dir bin_ca_tetra_r220/ --out_dir bin_ca_tetra_r220/ --cpus 64
```
Generate trimmed msa for Archaea(ar53) and Bacteria(bac120).


## Build Archaea(ar53) and Bacteria(bac120) tree separately  

```bash
# Inferring Archaea Tree
gtdbtk infer --msa_file MEERv21_ar53_gtdbr220_msa.fasta.gz --out_dir bin_ca_tetra_r220_archaea --cpus 64

# Inferring Bacteria Tree
gtdbtk infer --msa_file MEERv21_bac120_gtdbr220_msa.fasta.gz --out_dir bin_ca_tetra_r220_bacteria --cpus 64
```

With `gtdbtk infer` module, `FastTree` was introduced to generate phylogenetic trees based on msa for Archaea and Bacteria, respectively. Trees were all need to be further rooted and decorated by taxonomic classification.

```bash
# Rooting Archaea Tree with 'p__Altiarchaeota'

gtdbtk root --input_tree bin_ca_tetra_r220_archaea/gtdbtk.unrooted.tree --outgroup_taxon p__Altiarchaeota --output_tree bin_ca_tetra_r220_archaea/gtdbtk.rooted.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_ar53_summary.tsv

gtdbtk decorate --input_tree bin_ca_tetra_r220_archaea/gtdbtk.rooted.tree --output_tree bin_ca_tetra_r220_archaea/gtdbtk.decorated.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_ar53_summary.tsv

# Rooting Bacteria Tree 'p__Patescibacteria'

gtdbtk root --input_tree bin_ca_tetra_r220_bacteria/gtdbtk.unrooted.tree --outgroup_taxon p__Patescibacteria --output_tree bin_ca_tetra_r220_bacteria/gtdbtk.rooted.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_bac120_summary.tsv

gtdbtk decorate --input_tree bin_ca_tetra_r220_bacteria/gtdbtk.rooted.tree --output_tree bin_ca_tetra_r220_bacteria/gtdbtk.decorated.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_bac120_summary.tsv
```

Taxa recommanded by GTDB-Tk in help message `p__Altiarchaeota` and `p__Patescibacteria` were used for rooting Tree of Archaea and Bacteria.

## Determine Ranks for SGBs and Nodes

Ranks and Relative Evolutionary Distances (REDs) of nodes were assigned and calculated by `gtdbtk infer_ranks` module.

```bash
# Inferring ranks and RED for Archaea Tree
gtdbtk infer_ranks --input_tree bin_ca_tetra_r220_archaea/gtdbtk.decorated.tree --ingroup_taxon d__Archaea --output_tree bin_ca_tetra_r220_archaea/gtdbtk.infer_RED_ranks.tree
# Inferring ranks and RED for Bacteria Tree
gtdbtk infer_ranks --input_tree bin_ca_tetra_r220_bacteria/gtdbtk.decorated.tree --ingroup_taxon d__Bacteria --output_tree bin_ca_tetra_r220_bacteria/gtdbtk.infer_RED_rankstree
```

## Pull out Monophyletic Clades from Trees

The rank and RED decorated trees were processed by `ete` for generating clades only composed by MMC SGBs.  
Detailed script can be found in `ete_parsetree.ipynb`.