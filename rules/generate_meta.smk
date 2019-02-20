import math
import platform
"""Generate all the meta data files"""

### from: https://github.com/Hoohm/dropSeqPipe/blob/master/rules/generate_meta.smk

#Which rules will be run on the host computer and not sent to nodes
localrules:
     create_dict,
     reduce_gtf,
     create_refFlat,
     create_intervals


rule create_dict:
    input:
        "{ref_path}/{species}_{build}_{release}/genome.fa"
    output:
        "{ref_path}/{species}_{build}_{release}/genome.dict"
    threads:1
    params:
        picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
        temp_directory=config['LOCAL']['temp-directory']
    conda: config['LOCAL']['common_conda']
    shell:
        """java -jar -Djava.io.tmpdir={params.temp_directory} {params.picard} CreateSequenceDictionary\
        REFERENCE={input}\
        OUTPUT={output}
        """

rule reduce_gtf:
    input:
        reference_dict="{ref_path}/{species}_{build}_{release}/genome.dict",
        annotation="{ref_path}/{species}_{build}_{release}/annotation.gtf"
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        "{ref_path}/{species}_{build}_{release}/reduced_annotation.gtf"
    conda: config['LOCAL']['common_conda']
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ReduceGtf -m {params.memory}\
        GTF={input.annotation}\
        OUTPUT={output}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        IGNORE_FUNC_TYPE='null'\
        ENHANCE_GTF='false'"""

rule create_refFlat:
    input:
        reference_dict="{ref_path}/{species}_{build}_{release}/genome.dict",
        annotation="{ref_path}/{species}_{build}_{release}/annotation.gtf"
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        "{ref_path}/{species}_{build}_{release}/annotation.refFlat"
    conda: config['LOCAL']['common_conda']
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ConvertToRefFlat -m {params.memory}\
        ANNOTATIONS_FILE={input.annotation}\
        OUTPUT={output}\
        SEQUENCE_DICTIONARY={input.reference_dict}
        """

rule create_intervals:
    input:
        annotation_reduced="{ref_path}/{species}_{build}_{release}/reduced_annotation.gtf",
        reference_dict="{ref_path}/{species}_{build}_{release}/genome.dict"
    params:
        memory=config['LOCAL']['memory'],
        reference_directory=config['META']['reference-directory'],
        temp_directory=config['LOCAL']['temp-directory'],
        prefix="{species}_{build}_{release}/annotation"
    output:
        intervals="{ref_path}/{species}_{build}_{release}/annotation.rRNA.intervals"
    conda: config['LOCAL']['common_conda']
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && CreateIntervalsFiles -m {params.memory}\
        REDUCED_GTF={input.annotation_reduced}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        O={params.reference_directory}\
        PREFIX={params.prefix}
        """

rule prep_star_index:
    input:
        reference_file="{ref_path}/{species}_genome.fa",
        config_file='config.yaml'
    output:
        '{reference_directory}/star_ref_config.txt'
    conda: config['LOCAL']['common_conda']
    script:
        '../scripts/prep_star.py'


rule create_star_index:
    input:
        reference_file="{ref_path}/{species}_{build}_{release}/genome.fa",
        annotation_file="{ref_path}/{species}_{build}_{release}/annotation.gtf"  
    params:
        genomeDir='{ref_path}/{species}_{build}_{release}/STAR_INDEX/'
    output:
        '{ref_path}/{species}_{build}_{release}/STAR_INDEX/'
    threads: 24
    conda: config['LOCAL']['common_conda']
    shell:
        """mkdir -p {params.genomeDir}; STAR\
        --runThreadN {threads}\
        --runMode genomeGenerate\
        --genomeDir {params.genomeDir}\
        --genomeFastaFiles {input.reference_file}\
        --sjdbGTFfile {input.annotation_file}
        """