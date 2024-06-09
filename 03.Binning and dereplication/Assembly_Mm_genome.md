# Predict rrna from assembled contigs
git clone https://gitlab.genomics.cn/meer/zetaSeq
cd zetaSeq
chmod +x biom_concat contister contister_gather keep_pair_aln lnster rrnaster solve_profile truncaster
export PATH="your/path/to/zetaSeq:$PATH"
export PYTHONPATH="your/path/to/zetaSeq:$PYTHONPATH"

mkdir -p rfam/cm/
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz;mv Rfam.tar.gz rfam/cm/;tar xzvf Rfam.tar.gz
rrnaster -g {sample}.contig.fa -d /path/to/rfam/cm/ -r all -o {sample}.rrna.tsv -c {sample}.rrna.fa -t 4 # rrnaster can search models for all Kingdoms at the same time, predicted region overlapped will be clustered later based on hit E-value.

## Combined predicted SSU or LSU with SLIVA sequences
cat expand({sample}.rrna.fa) > rrna_ssu.fa
cat expand({sample}.rrna.fa) > rrna_lsu.fa

cat SILVA_138.1_SSURef_NR99_tax_silva.fasta rrna_ssu.fa > ref_rrna_ssu.fna
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta rrna_lsu.fa > ref_rrna_lsu.fna

seqtk seq -L 1000 ref_rrna_ssu.fna > ref_rrna_ssu.1000.fna
seqtk seq -L 1000 ref_rrna_lsu.fna > ref_rrna_lsu.1000.fna

## Cluster the new ref at 0.97 similarity
vsearch --cluster_size ref_rrna_ssu.fna --id 0.97 --centroid ref_rrna_ssu.1000.97.fna --fast_width 0 --threads 12
vsearch --cluster_size ref_rrna_lsu.fna --id 0.97 --centroid ref_rrna_lsu.1000.97.fna --fast_width 0 --threads 12

## Classify using sintax
usearch -sintax ref_rrna_ssu.1000.97.fna -db 16s.udb -tabbedout ref_rrna_ssu.1000.97.sintax -strand both -sintax_cutoff 0.8
usearch -sintax ref_rrna_lsu.1000.97.fna -db 16s.udb -tabbedout ref_rrna_lsu.1000.97.sintax -strand both -sintax_cutoff 0.8

# Align reading onto new ref using BURST
## Build a burst index
Download BURST binary at https://github.com/knights-lab/BURST
burst12 -r ref_rrna_ssu.1000.97.fna -a databases/ref_rrna_ssu.1000.97.acx -o databases/ref_rrna_ssu.1000.97.edx -d name DNA qLen 600 -t 12
burst12 -r ref_rrna_lsu.1000.97.fna -a databases/ref_rrna_lsu.1000.97.acx -o databases/ref_rrna_lsu.1000.97.edx -d name DNA qLen 600 -t 12

## Align all sample reads onto the burst index (only SSU index was used)
snakemake -s snk_rrna_profile.smk -j {threads}

# Assmebly the fungal genome
From the combined 1194 sample profiles, the most abundant OTU turns out to be a mitochondria ssu of a fungus, by searching in NCBI.
We picked four samples with the highest abudnance with this OTU, combined them and assembled the combined dataset with megahit

for sn in {FDZ044-RBk24-26,FDZ039-Y14-16,FDZ081-YWG14-16,FDZ081-YWG16-18}
do
	cat data/samples/${sn}.clean.1.fq.gz >> data/cat/cat.1.fq.gz
	cat data/samples/${sn}.clean.2.fq.gz >> data/cat/cat.2.fq.gz
done

megahit -1 data/cat/cat.1.fq.gz -2 data/cat/cat.2.fq.gz -t 12 --presets meta-sensitive -o data/megahit/
seqtk seq -L 1000 data/megahit/final.contigs.fa > data/contigs_for_binning/cat.contigs.fa

snakemake -s snk_single_contig_depth.smk -j {threads}

metabat2 -i data/contigs_for_binning/cat.contigs.fa -o data/metabat2_bins/bin -a data/jgi_depth.depth

Metabat2 returned with 2 bins, one of them was a known Archaea, while another turns out to have a BUSCO score 98% as a Eukaryota and 90.8% as a Dothideomycota. A single cirular dna was also assmebled out and turns out to be a cirular mitochondria DNA. Both the fungal genome and the mtDNA point to a close relative of fungus Myriangium duriaei, an insect skin parasit.

The fungus discoved was given a code name Dasco (as for Deepsea asco), and annotate using funannotate

22 genomes of Dothideomycota in the () were downloaded from NCBI and annotate the same wat using funannotate.

The total of 23 genomes, including Dasco were compared using funannotate.

Based on the phylogentic distance of Dasco on the tree, it was given the name Ca. Myriangium mariana
