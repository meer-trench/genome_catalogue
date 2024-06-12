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

## Building Trees for MMC SGBs

