from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()


localrules:
     download_genome,
     download_annotation

rule download_genome:
    input:
        FTP.remote(config['META']['genome_download_path'])
    output:
        "{ref_path}/{species}_{build}_{release}/genome.fa"
    shell:
        "gunzip -d -c {input} > {output}"

rule download_annotation:
    input:
        FTP.remote(config['META']['annotation_download_path'])
    output:
        "{ref_path}/{species}_{build}_{release}/annotation.gtf"
    shell:
        "gunzip -d -c {input} > {output}"

        
