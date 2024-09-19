process SENTIEON_CLI {
    tag "$meta.id"
    label 'process_high'

    container "${params.sentieoncli_container_url}:${params.sentieoncli_container_tag}"

    input:
    tuple val(meta), path(r1_fastq), path(r2_fastq)
    path index, stageAs: "index/*"
    tuple path(fasta), path(fai)
    path ml_model

    output:
    tuple val(meta), path("**"), emit: output
    
    script:
    def args = task.ext.args ?: ''

    def license = params.sentieon_license
    def prefix = meta.prefix ? "${meta.prefix}" : "${meta.id}"
    def read_group = meta.read_group ? "${meta.read_group}__${meta.sample}" : "RG"
    def platform = meta.platform ? "${meta.platform}" : "PL"
    def sample = meta.sample ? "${meta.sample}" : "SM"
    def model_option = ml_model ? "-x ${ml_model}/bwa.model" : ''
    def memory = task.memory.toString().replaceAll(' ', '').replaceAll('GB','G')

    """
    logfile=run.log
    exec > >(tee \$logfile)
    exec 2>&1

    LICENSE=$license
    export SENTIEON_LICENSE=$license

    find -L $index/ -type f \\! -name "*.fa" -exec ln -s {} . \\;

    #INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    #FASTA=`find -L ./ -maxdepth 1 -name "*.fa"`

    echo "sentieon-cli dnascope -t $task.cpus -r $fasta --r1-fastq ${r1_fastq} --r2-fastq ${r2_fastq} --bam_format --readgroups "@RG\\tID:$read_group\\tSM:$sample\\tPL:$platform\\tLB:$sample" --model-bundle $ml_model --pcr-free --assay WGS ${prefix}.vcf.gz"

    sentieon-cli dnascope \\
        -t $task.cpus \\
        -r $fasta \\
        --r1-fastq ${r1_fastq} \\
        --r2-fastq ${r2_fastq} \\
        --bam_format \\
        --readgroups "@RG\\tID:$read_group\\tSM:$sample\\tPL:$platform\\tLB:$sample" \\
        --model-bundle $ml_model \\
        --pcr-free \\
        --assay WGS \\
        ${prefix}.vcf.gz

    touch ${prefix}.bam
    touch ${prefix}.bam.bai


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo /opt/sentieon/sentieon-genomics-* 2>&1 | sed 's/\\/opt\\/sentieon\\/sentieon-genomics-//')
    END_VERSIONS
    """
}
