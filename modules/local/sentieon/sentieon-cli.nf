process SENTIEON_CLI {
    tag "$meta.id"

    cpus = 32
    memory = 64.GB
    time = 6.h
    maxRetries = 3
    
    container "${params.sentieoncli_container_url}:${params.sentieoncli_container_tag}"

    input:
    tuple val(meta), val(read_group), path ( r1_fastq, stageAs: "?/*"), path ( r2_fastq, stageAs: "?/*")
    path index, stageAs: "index/*"
    tuple path(fasta), path(fai)
    path ml_model

    output:
    tuple val(meta), path("*_deduped.bam"), path("*_deduped.bam.bai"), emit: bam_bai
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*_metrics/*"), emit: metrics
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''

    def license = params.sentieon_license
    def prefix = meta.prefix ? "${meta.prefix}" : "${meta.id}"
    def read_group = meta.read_group ? "${meta.read_group}__${meta.sample}" : "RG"
    def platform = meta.platform ? "${meta.platform}" : "PL"
    def sample = meta.sample ? "${meta.sample}" : "SM"
    def model_option = ml_model ? "-x ${ml_model}/bwa.model" : ''
    def memory = task.memory.toString().replaceAll(' ', '').replaceAll('GB','G')

    //--readgroups "@RG\\tID:$read_group\\tSM:$sample\\tPL:$platform\\tLB:$sample" \\
    // Concatenate read_group into a space-separated string
    // --r1-fastq a_r1.fastq.gz b_r1.fastq.gz c_r1.fastq.gz --r2-fastq a_r2.fastq.gz b_r2.fastq.gz c_r2.fastq.gz --readgroups "rg1" "rg2" "rg3"
    
    def readgroups_string = read_group.collect { rg -> "@RG\\tID:${rg.read_group}__${rg.sample}\\tSM:${rg.sample}\\tPL:${rg.platform}\\tLB:${rg.sample}"}.join('" "') 

    """
    logfile=run.log
    exec > >(tee \$logfile)
    exec 2>&1

    LICENSE=$license
    export SENTIEON_LICENSE=$license

    echo "sentieon-cli dnascope -t $task.cpus -r $fasta --r1-fastq ${r1_fastq} --r2-fastq ${r2_fastq} --bam_format --readgroups "@RG\\tID:$read_group\\tSM:$sample\\tPL:$platform\\tLB:$sample" --model-bundle $ml_model --pcr-free --assay WGS ${prefix}.vcf.gz"

    sentieon-cli dnascope \\
        -t $task.cpus \\
        -r $fasta \\
        --r1-fastq ${r1_fastq} \\
        --r2-fastq ${r2_fastq} \\
        --bam_format \\
        --readgroups "$readgroups_string" \\
        --model-bundle $ml_model \\
        --pcr-free \\
        --assay WGS \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo /opt/sentieon/sentieon-genomics-* 2>&1 | sed 's/\\/opt\\/sentieon\\/sentieon-genomics-//')
    END_VERSIONS
    """
}
