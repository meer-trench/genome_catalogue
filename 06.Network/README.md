# Code for network analysis
## Use fastspar to generate co-occurrence network
- https://github.com/scwatts/fastspar

### Step1 Calculate correlations
```bash
mkdir -p network_output
conda activate fastspar
fastspar --threshold 0.5 \
    --iterations 50 \
    --threads 50 \
    --otu_table network_input/meer.rpkm.tsv \
    --correlation network_output/meer.MAGs.correlation.tsv \
    --covariance network_output/meer.MAGs.covariance.tsv
```

### Step2 Calculate bootstrap counts
```bash
mkdir -p bootstrap_counts
fastspar_bootstrap --otu_table network_input/meer.rpkm.tsv \
    --number 1000 \
    --prefix bootstrap_counts/meer \
    --threads 64
```

### Step3 Calculate bootstrap correlations
```bash
mkdir -p bootstrap_correlation
parallel fastspar --otu_table {} \
    --correlation bootstrap_correlation/cor_{/} \
    --covariance bootstrap_correlation/cov_{/} \
    -i 5 ::: bootstrap_counts/*
```

### Step4 Calculate pvalues
```bash
mkdir -p network_output
fastspar_pvalues --otu_table network_input/meer.rpkm.tsv \
    --correlation network_output/meer.MAGs.correlation.tsv \
    --prefix bootstrap_correlation/cor_meer \
    --permutations 1000 \
    --outfile network_output/meer.MAGs.pvalues.tsv \
    --threads 64
```

## Use iDIRECT to scale co-occurrence network
- https://github.com/nxiao6gt/iDIRECT

To prevent an overload of edges that could impede the calculation by iDirect, the edge file should be scaled manually from the correlation matrix.

```bash
mkdir -p network_topology
python /data1/liliuyang/pylib/iDIRECT/idirect.py \
    network_topology/meer.to_idr_edge.txt \
    network_topology/meer.idr_edge.txt
```