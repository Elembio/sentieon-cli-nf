process SENTIEON_CLI_SE {
    tag "$meta.id"
    label 'process_high'
    
    container "${params.sentieoncli_container_url}:${params.sentieoncli_container_tag}"

    input:
    tuple val(meta), val(read_group), path(r1_fastq, stageAs: "?/*")
    path index, stageAs: "index/*"
    tuple path(fasta), path(fai)
    path ml_model
    val assay
    path target_region_bed

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
    def model_option = ml_model ? "-x ${ml_model}/bwa.model" : ''
    def pcr_option = params.pcr ? "" : "--pcr-free"
    def target_region_bed_option = target_region_bed ? "--bed ${target_region_bed}" : ''
    def memory = task.memory.toString().replaceAll(' ', '').replaceAll('GB','G')
    def readgroups_string = read_group.collect { rg -> "@RG\\tID:${rg.read_group}__${rg.sample}\\tSM:${rg.sample}\\tPL:${rg.platform}\\tLB:${rg.sample}"}.join('" "') 

    """
    logfile=run.log
    exec > >(tee \$logfile)
    exec 2>&1

    LICENSE=$license
    export SENTIEON_LICENSE=$license

    find -L $index/ -type f \\! -name "*.fa" -exec ln -s {} . \\;
    
    echo "sentieon-cli dnascope -t $task.cpus -r $fasta --r1-fastq ${r1_fastq} --bam_format --readgroups "$readgroups_string" --model-bundle $ml_model --assay $assay ${pcr_option} ${target_region_bed_option} ${prefix}.vcf.gz"

    sentieon-cli dnascope \\
        -t $task.cpus \\
        -r $fasta \\
        --r1-fastq ${r1_fastq} \\
        --bam_format \\
        --readgroups "$readgroups_string" \\
        --model-bundle $ml_model \\
        --assay $assay \\
        --skip-multiqc \\
        ${pcr_option} \\
        ${target_region_bed_option} \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo /opt/sentieon/sentieon-genomics-* 2>&1 | sed 's/\\/opt\\/sentieon\\/sentieon-genomics-//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def license = params.sentieon_license
    def prefix = meta.prefix ? "${meta.prefix}" : "${meta.id}"
    def model_option = ml_model ? "-x ${ml_model}/bwa.model" : ''
    def pcr_option = params.pcr ? "" : "--pcr-free"
    def bed_option = bed ? "--bed ${bed}" : ''
    def memory = task.memory.toString().replaceAll(' ', '').replaceAll('GB','G')
    def readgroups_string = read_group.collect { rg -> "@RG\\tID:${rg.read_group}__${rg.sample}\\tSM:${rg.sample}\\tPL:${rg.platform}\\tLB:${rg.sample}"}.join('" "') 
    """
    echo "stub for SENTIEON_CLI"
    touch ${prefix}_deduped.bam
    touch ${prefix}_deduped.bam.bai
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    mkdir ${prefix}_metrics/
    touch ${prefix}_metrics/metrics
    touch ${prefix}.log

    echo "sentieon-cli dnascope -t $task.cpus -r $fasta --r1-fastq ${r1_fastq} --bam_format --readgroups "$readgroups_string" --model-bundle $ml_model --assay $assay ${pcr_option} ${bed_option} ${prefix}.vcf.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: "stub version"
    END_VERSIONS
    """
}
