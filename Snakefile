
import pandas as pd
from snakemake.utils import validate, min_version
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# Load configuration files
try:
    configfile_path = config['configfile_path']
except:
    configfile_path = "config.yaml"    
configfile: configfile_path

# Include the gtf biotypes yaml
configfile: config['META']['gtf_biotypes']

# experiment metadata
exp_mat = pd.read_table(config['INPUT']['exp_mat'])

# Define a few variables to make them easier to reference
ref_path = config['META']['reference-directory']
results_dir = config['LOCAL']['results']
logs_dir = config['LOCAL']['logs']
fastq_dir = config['LOCAL']['fastq']
temp_dir = config['LOCAL']['temp-directory']
flexbar_adapter_1 = config['FILTER']['FLEXBAR']['adapter_R1']
species=list(config['META']['species'])[0]
build=[config['META']['species'][species]['build']][0]
release=[config['META']['species'][species]['release']][0]
index_STAR = expand('{ref_path}/{species}_{build}_{release}/STAR_INDEX', ref_path=ref_path, species=species, build=build, release=release)[0]

# define source of fastq files
if config['INPUT']['fastq_method'] == 'ena':
    print('using ENA to download fastq files')
    srrs = exp_mat['run']
else:
    #srrs, = glob_wildcards(fastq_dir+"/{id}.fastq.gz")
    srrs = exp_mat['run']
    print('using fastq files in '+config['LOCAL']['fastq'])

# print some info
print("Genome info:")
print("species: "+species)
print("build: "+build)
print("release: "+str(release))

print("Using index: "+str(index_STAR))

print("These are the input run IDs:")
print(srrs)

###
### functions
###


### function to get ENA fasp file paths
# ascp -QT -l 300m -P33001 -i ~/.aspera/asperaweb_id_dsa.openssh  era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR644/005/SRR6441685/SRR6441685.fastq.gz
def get_ena_fasp(run):
  fasp=str('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/'+run[0:6]+'/00'+run[-1]+'/'+run+'/'+run+'.fastq.gz')
  return fasp
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR644/005/SRR6441685/SRR6441685.fastq.gz
def get_ena_ftp(run):
  fasp=str('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+run[0:6]+'/00'+run[-1]+'/'+run+'/'+run+'.fastq.gz')
  return fasp

### determine whether we should use only trimmed fastq files for STAR alignment
## use of .item() from: https://stackoverflow.com/a/36922103
def determine_input_filter(wildcards):
    if exp_mat[exp_mat.run == wildcards.sample].trimmed_only.item == bool(True):
        print('using only trimmed reads for '+wildcards.sample)
        return temp_dir+"/filter/{sample}.fastq"
    else:
        print('using all reads for '+wildcards.sample)
        return temp_dir+"/flexbar/{sample}.fastq"    


###
### rules
###

rule all:
    input:
        expand(
            # index
            ['{ref_path}/{species}_{build}_{release}/STAR_INDEX/SA',
            # fastq
            #fastq_dir+'/{sample}.fastq.gz',
            #qc
            '{logs_dir}/multiqc.html',
            #filter
            #mapping
            '{results_dir}/aln/{sample}.Aligned.sortedByCoord.out.bam',
            '{results_dir}/aln/{sample}.Aligned.sortedByCoord.out.bam.bai',
            #extract
            #merge
            ],
                sample=srrs,
                results_dir=results_dir,
                logs_dir=logs_dir,
                ref_path=ref_path,
                build=config['META']['reference-directory'],
                release=release,
                species=species)

 

rule download_ENA_fasp:
    input:
    output:
        fastq_dir+'/{sample}.fastq.gz'
    params:
        ascp="-q -QT -l 300m -P33001 -i ~/.aspera/asperaweb_id_dsa.openssh",
        fasp=lambda wildcards: get_ena_fasp(wildcards.sample),
        outDir=fastq_dir
    conda:
        config['LOCAL']['common_conda']
    resources: ascp_limit=1
    shell:
        '''
        ascp {params.ascp} {params.fasp} {params.outDir}
        '''    
                
#rule download_ENA_ftp:
#    input:
#        lambda wildcards: FTP.remote(get_ena_ftp('{sample}'.format(sample=wildcards.sample)), static=True, immediate_close=True)
#    output:
#        fastq_dir+'/{sample}.fastq.gz'
#    run:
#        shell("mv {input} {output}")
                
rule fastqc:
    input:
        fastq_dir+'/{sample}.fastq.gz'
    output:
        '{logs_dir}/fastqc/{sample}_fastqc.html'
    conda:
        config['LOCAL']['common_conda']
    params:
        outDir='{logs_dir}/fastqc'
    threads: 6
    shell:
        '''
        fastqc \
        --quiet \
        --outdir  {params.outDir} \
        --threads {threads} \
        --adapters ~/genomes/fastqc/adapter_list.txt \
        --contaminants ~/genomes/fastqc/contaminant_list.txt \
        {input}
        '''

rule head_fastq:
    input:
        fastq_dir+'/{sample}.fastq.gz'
    output:
        '{temp_dir}/fastqhead/{sample}.fastq.gz'
    shell:
        '''
        set +o pipefail ;
        zcat {input} \
        | head -n 400000 \
        | gzip \
        > {output}
        '''

    
rule flexbar_se:
    input:
        # temp_dir+'/fastqhead/{sample}.fastq.gz' ## head only first lines to test pipeline
        fastq_dir+'/{sample}.fastq.gz' ## full files
    output:
        temp_dir+'/flexbar/{sample}.fastq'
    log:
        logs_dir+'/flexbar/{sample}.log'
    params:
        adapter_1=flexbar_adapter_1,
        prefix=lambda wildcards: temp_dir+'/flexbar/'+wildcards.sample
    conda:
        config['LOCAL']['common_conda']
    threads: 4
    shell:
        '''
        flexbar \
        --reads {input} \
        --threads {threads} \
        --adapters {params.adapter_1} \
        --adapter-cycles 2 \
        --min-read-length 12 \
        --target {params.prefix} \
        --removal-tags &&
        mv {params.prefix}.log {log}
        '''

# rule flexbar_fix_logs:
#     input:
#         temp_dir+'/flexbar/{sample}.flexbar.log'
#     output:
#         logs_dir+'/flexbar/{sample}.flexbar.log'    
#     params:
#         run_id='{sample}'
#     shell:
#         '''
#         mv {input} {output} &&
#         sed -i -e 's/stdout/{params.run_id}/g'  {output}
#         '''
 
rule filter_nontrimmed:
    input:
        temp_dir+'/flexbar/{sample}.fastq'
    output:
        temp_dir+"/filter/{sample}.fastq"
    shell:
        '''
        awk '/Flexbar_removal/ {{{{print}} for(i=1; i<=3; i++) {{getline; print}}}}' {input} > {output}
        '''
        
rule align_STAR_SE:
    input:
        fastq=determine_input_filter,
        index=index_STAR+'/SA'
    output:
        results_dir+'/aln/{sample}.Aligned.sortedByCoord.out.bam'
    log:
        logs_dir+'/aln/{sample}.Log.final.out'
    params:
        run=lambda wildcards: wildcards.sample,
        prefix=lambda wildcards: results_dir+'/aln/'+wildcards.sample
    threads: 12
    shell:
        '''
        STAR \
        --runThreadN {threads} \
        --genomeDir {index_STAR} \
        --readFilesIn {input.fastq} \
        --outSAMattributes MD NH \
        --outFileNamePrefix {params.prefix}. \
        --outSAMtype BAM SortedByCoordinate &&
        mv {params.prefix}.Log.final.out {log}
        '''

rule index_bam:
    input:
        results_dir+'/aln/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        results_dir+'/aln/{sample}.Aligned.sortedByCoord.out.bam.bai'
    conda:
        config['LOCAL']['common_conda']
    shell:
        '''
        samtools index {input}
        '''

rule RnaSeqMetricsCollector:
    input:
        data=results_dir+'/aln/{sample}.Aligned.sortedByCoord.out.bam',
        refFlat=expand("{ref_path}/{species}_{build}_{release}/annotation.refFlat",
            ref_path=config['META']['reference-directory'],
            species=species,
            release=release,
            build=build),
        rRNA_intervals=expand("{ref_path}/{species}_{build}_{release}/annotation.rRNA.intervals",
            ref_path=config['META']['reference-directory'],
            species=species,
            release=release,
            build=build)
    params:     
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory']
    output:
        rna_metrics=logs_dir+'/aln/{sample}_rna_metrics.txt',
    conda: 'envs/preprocess.yaml'
    shell:
        """
        picard CollectRnaSeqMetrics -Xmx4g \
        INPUT={input.data}\
        OUTPUT={output}\
        STRAND=FIRST_READ_TRANSCRIPTION_STRAND\
        REF_FLAT={input.refFlat}\
        RIBOSOMAL_INTERVALS={input.rRNA_intervals}
        """
    
rule multiqc:
    input:
        expand(logs_dir+'/fastqc/{sample}_fastqc.html', sample=srrs),
        expand(logs_dir+'/flexbar/{sample}.log', sample=srrs),
        expand(results_dir+'/aln/{sample}.Aligned.sortedByCoord.out.bam',  sample=srrs),
        expand(logs_dir+'/aln/{sample}_rna_metrics.txt', sample=srrs)
    output:
        '{logs_dir}/multiqc.html'
    params:
        inDir=logs_dir,
        outDir=logs_dir,
        outFile="multiqc.html"
    conda:
        config['LOCAL']['common_conda']
    shell:
        '''
        multiqc \
        --outdir {params.outDir} \
        --filename {params.outFile} \
        --force \
        {params.inDir}
        '''
    
include: "rules/get_genomes.smk"
include: "rules/generate_meta.smk"
