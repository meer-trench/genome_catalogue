# Classification for MMC SGBs

## Generate Taxonomy Classification for SGBs with GTDB-Tk

```bash
# All Species-level Genome Bins (SGBs)
# GTDB-Tk version 2.4.0
# GTDB release 220

gtdbtk classify_wf --mash_db ./gtdbtk240/share/gtdbtk-2.4.0/db/release220/gtdb-tk-r220.msh \
--batchfile release_220/bin_ca_tetra/MEER_SGB.list -x fa \
--out_dir release_220/bin_ca_tetra/meerv2.1_ca_sgb_mash --cpus 64 --pplacer_cpus 16
```

Classification TSVs were renamed to `MEERv21_d95_gtdbr220_ar53_summary.tsv` and `MEERv21_d95_gtdbr220_bac120_summary.tsv`.

## Building Trees for MMC SGBs

Since a MASH and `ani-screen` mode were used in genome taxonomic classification, marker identification and alignment for ANI screened SGBs should be conducted.  

```bash
# All Species-level Genome Bins (SGBs)
# GTDB-Tk version 2.4.0
# GTDB release 220

gtdbtk identify --batchfile MEER_SGB.list --out_dir bin_ca_tetra_r220 --cpus 64 -x fa
gtdbtk align --identify_dir bin_ca_tetra_r220/ --out_dir bin_ca_tetra_r220/ --cpus 64
```

Then, we inferred trees from MSA with `FastTree`.

```bash
# MSA were renamed
# Inferring MMC Archaea Tree
gtdbtk infer --msa_file MEERv21_ar53_gtdbr220.user_msa.fasta.gz --out_dir bin_ca_tetra_user_arctree_r220 --cpus 64

# Inferring MMC Bacteria Tree
gtdbtk infer --msa_file MEERv21_bac120_gtdbr220.user_msa.fasta.gz --out_dir bin_ca_tetra_user_bactree_r220 --cpus 64
```

```bash
# Rooting Archaea Tree with 'p__Nanoarchaeota'

gtdbtk root --input_tree bin_ca_tetra_user_arctree_r220/gtdbtk.unrooted.tree --outgroup_taxon p__Nanoarchaeota --output_tree bin_ca_tetra_user_arctree_r220/gtdbtk.rooted.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_ar53_summary.tsv

gtdbtk decorate --input_tree bin_ca_tetra_user_arctree_r220/gtdbtk.rooted.tree --output_tree bin_ca_tetra_user_arctree_r220/gtdbtk.decorated.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_ar53_summary.tsv

# Rooting Bacteria Tree 'p__Patescibacteria'

gtdbtk root --input_tree bin_ca_tetra_user_bactree_r220/gtdbtk.unrooted.tree --outgroup_taxon p__Patescibacteria --output_tree bin_ca_tetra_user_bactree_r220/gtdbtk.rooted.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_bac120_summary.tsv

gtdbtk decorate --input_tree bin_ca_tetra_user_bactree_r220/gtdbtk.rooted.tree --output_tree bin_ca_tetra_user_bactree_r220/gtdbtk.decorated.tree --gtdbtk_classification_file MEERv21_d95_gtdbr220_bac120_summary.tsv
```

## Visualization

see Folder `Visualization`  
