configfile: 'meer_snk_scripts/config.yaml'

import os
import sys
SAMPLES = {}
path = config['project']
if not path.endswith('/'): path += '/'
sample_path = path + 'data/samples/'

folder_checks = ['TMP/', path + 'data/contigs', 
                 path + 'data/binning/', 
                 path + 'data/bin_refinement/',
                 path + 'data/bin_passed_all/', 
                 path + 'data/gather_all_checkm/']
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

for file in os.listdir(sample_path):
    if file.endswith('.path'):
        sample_name = file.split('.')[0]
        SAMPLES[sample_name] = SAMPLES.get(sample_name, []) + [sample_path + file]

if len(SAMPLES) == 0:
    print('ERROR: no file under {0}.'.format(filepath))
    sys.exit('Pipeline terminated orz.')
else:
    for key, value in SAMPLES.items():
        if len(value) != 1:
            print('{0} does not have the right files.'.format(key))
            sys.exit('Pipeline terminated, please check.')
        else:
            SAMPLES[key] = tuple(sorted(value))


rule target:
    input:
        mk_gtdb = path + 'data/markers/gtdbtk_mags.mk',
        mk_drep = path + 'data/markers/drep_all_95.mk'

rule fastp_merged:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        ma = path + 'data/fastp_merged/{sample}.merged.fa.gz', 
        a1 = path + 'data/fastp_merged/{sample}.unmerged.1.fa.gz',
        a2 = path + 'data/fastp_merged/{sample}.unmerged.2.fa.gz'
    threads: 16
    params:
        sq = config['sequencer'],
        sn = '{sample}',
        raw = config['raw'],
        path = path + 'data/contigs/{sample}/'    
    log: path + 'logs/fastp_merged/{sample}.json'
    shell:
        """
        meer_snk_scripts/scripts/run_fastp.py -i {input[0]} -w {params.raw} -o {output.ma} -f {output.a1} -r {output.a2} -s {params.sq} -t {threads} -log {log}
        """
    
rule megahit:
    input:
        ma = path + 'data/fastp_merged/{sample}.merged.fa.gz',
        a1 = path + 'data/fastp_merged/{sample}.unmerged.1.fa.gz',
        a2 = path + 'data/fastp_merged/{sample}.unmerged.2.fa.gz'
    output:
        mk = path + 'data/markers/megahit_mk/{sample}.mk',
        fa = path + 'data/contigs_for_binning/{sample}.contigs.fa'
    threads: 48
    params:
        sn = '{sample}',
        path = path + 'data/contigs/{sample}/'
    shell:
        """
        if [ -d "{params.path}" ];then rm -r {params.path};else echo Folder not exist, MEGAHIT good to go.;fi
        megahit -1 {input.a1} -2 {input.a2} -r {input.ma} -t {threads} --presets meta-sensitive -o {params.path}
        seqtk rename {params.path}final.contigs.fa {params.sn}@ctg_ | seqtk seq -L 1000 - | seqtk seq -C - > {output.fa}
        touch {output.mk}
        """

rule fastp_not_merged:
    input:
        path = path + 'data/samples/{sample}.path'
    output:
        a1 = path + 'data/fastp_not_merged/{sample}.unmerged_1.fastq',
        a2 = path + 'data/fastp_not_merged/{sample}.unmerged_2.fastq'
    threads: 16
    params:
        sq = config['sequencer'],
        sn = '{sample}',
        raw = config['raw'],
        path = path + 'data/contigs/{sample}/'    
    log: path + 'logs/fastp_not_merged/{sample}.json'
    shell:
        """
        meer_snk_scripts/scripts/run_fastp.py -i {input.path} -w {params.raw} -f {output.a1} -r {output.a2} -s {params.sq} -t {threads} -log {log} -not_gz_out -fastq_out
        """


rule bwa_index:
    input:
        ctg = path + 'data/contigs_for_binning/{sample}.contigs.fa'
    output:
        mk = path + 'data/markers/bwa_index/{sample}.mk'
    threads: 48
    params:
        path = path + 'data/binning/{sample}_metabat2/'
    shell:
        """
        bwa index {input.ctg}
        touch {output.mk}
        """


rule random_draw_99:
	input:
        path = path + 'data/samples/{sample}.path'
    output:
        path = path + 'data/co_abundance/{sample}.list'
    threads: 1
    shell：
        """
        meer_snk_scripts/scripts/random_draw_99_samples.py {input.path} {output.path}
        """


rule bwa_mem:
    input:
        path = path + 'data/co_abundance/{sample}.list'
    output:
        mk = path + 'data/bwa_mem/{sample}/{sample}.mk'
    threads: 64
    shell:
        """
        meer_snk_scripts/scripts/bwa_mem_100_samples.py {input.path} {output.mk} {threads}
        """


rule jgi_depth:
    input:
        mk = path + 'data/bwa_mem/{sample}/{sample}.mk'
    output:
        jgi = path + 'data/jgi_depth/{sample}.depth'
    threads: 64
    params:
        path = path + 'data/bwa_mem/{sample}/*.bam'
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.jgi} {params.path}
        """


rule metabat2:
    input:
        jgi = path + 'data/jgi_depth/{sample}.depth',
        ctg = path + 'data/contigs_for_binning/{sample}.contigs.fa'
    output:
        mk = path + 'data/markers/metabat2/{sample}.mk'
    threads: 64
    params:
        path = path + 'data/metabat2_bins/{sample}/bin'
    shell:
        """
        metabat2 -i {input.ctg} -a {input.jgi} -o {params.path}
        touch {output.mk}
        """


rule checkm:   
    input:
        mk = path +'data/markers/metabat2/{sample}.mk'
    output:
        tsv = path +'data/checkm/{sample}.tsv
    threads: 64
    params:
        in_path = path + 'data/metabat2_bins/{sample}/',
        out_path = path +'data/checkm/{sample}/'
    shell:
        """
        checkm lineage_wf -t {threads} -x fa {params.in_path} {params.out_path}
        """
    

rule gather_all_checkm_csv:
    input:
        csv = expand(path + 'data/checkm/{sample}/{sample}.csv', sample=SAMPLES)
    output:
        csv = path + 'data/checkm_all_passed_bins.csv'
    threads: 1
    shell:
        """
        touch {output.csv}
        echo "genome,completeness,contamination" > {output.csv}
        cat {input.csv} >> {output.csv}
        """


rule gath_all_bins:
    input:
        mk = path + 'data/markers/metabat2/{sample}.mk'
    output:
        mk = path + 'data/markers/gath_all_bins/{sample}.mk'
    threads: 1
    params:
        in_path = path + 'data/metabat2_bins/{sample}/*.fa'
        out_path = path + 'data/all_passed_bins/'
    shell:
        """
        cp {params.in_path} {params.out_path}
        touch {output.mk}
        """


rule gtdbtk_mags:
    input:
        mk = expand(path + 'data/markers/metabat2/{sample}.mk', sample=SAMPLES)
    output:
        mk = path + 'data/markers/gtdbtk_mags.mk'
    params:
        in_path = path + 'data/bin_passed_all/',
        out_path = path + 'data/gtdb_passed_mags/'
    threads: 48
    conda:
        "gtdbtk207"
    shell:
        """
        gtdbtk classify_wf --genome_dir {params.in_path} --out_dir {params.out_path} -x fa --cpus {threads} --pplacer_cpus 8
        touch {output.mk}
        """


rule drep_all_99:
    input:
        csv = path + 'data/checkm_all_passed_bins.csv',
        mk = expand(path+ 'data/markers/gather_all_bins/{sample}.mk', sample=SAMPLES)
    output:
        mk = path + 'data/markers/drep_all_99.mk'
    threads: 48
    conda:
        "base"
    params:
        in_path = path + 'data/bin_passed_all/*.fa',
        out_path = path + 'data/drep_all_99/'
    shell:
        """
        dRep dereplicate {params.out_path} --genomeInfo {input.csv} -g {params.in_path} -p {threads} -comp 50 -con 10 -sa 0.99
        touch {output.mk}
        """

rule drep_all_95:
    input:
        csv = path + 'data/checkm_all_passed_bins.csv',
        mk = path + 'data/markers/drep_all_99.mk'
    output:
        mk = path + 'data/markers/drep_all_95.mk'
    threads: 48
    conda:
        "base"
    params:
        in_path = path + 'data/drep_all_99/dereplicated_genomes/*.fa',
        out_path = path + 'data/drep_all_95/'
    shell:
        """
        dRep dereplicate {params.out_path} --genomeInfo {input.csv} -g {params.in_path} -p {threads} -comp 50 -con 10 -sa 0.95
        touch {output.mk}
        """