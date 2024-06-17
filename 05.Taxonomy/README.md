# Taxonomic Classification of MMC SGBs and Phylogenetic Anaslysis  

This directory includes 4 directories of scripts and markdowns for different purposes.

- Classification: A markdown file indicates all commands used in taxonomic classification of MMC SGBs.  
- Phylogeny: A markdown file indicates all commands used in tree building of MMC SGBs with representatve SGBs of GTDB release 220. `ete_parsetree.ipynb` processes trees to locate monophyletic clades only constituted by MMC SGBs.
- Enrichment: A R markdown file processes the profile of MEER samples quantified by MMC SGBs to find out whether SGBs show enriched in certain geological location. Dunn-Test is applied.
- Visualization:  `tree-building.ipynb` visualized the MMC SGBs tree. The tree graph was annotated by Phylum-level taxonomic classification, enrichment results and SGB quality information. `Sankey-Taxa-gtdbr220.R` draws a sankey flow to describe the taxonomic source of underrepresented SGBs.  
