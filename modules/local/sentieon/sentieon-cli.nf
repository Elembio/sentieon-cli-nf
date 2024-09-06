process SENTIEON_CLI {
    tag "$meta.id"
    label 'process_high'

    container "${params.sentieoncli_container_url}:${params.sentieoncli_container_tag}"

    input:
    tuple val(meta), path(fastq)
    path index
    path fasta 
    path ml_model

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_bai
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    def license    = params.sentieon_license
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

    sentieon bwa mem \\
        ${model_option} \\
        -M -R "@RG\\tID:$read_group\\tSM:$sample\\tPL:$platform" \\
        -t $task.cpus \\
        -K 10000000 \\
        \$INDEX \\
        $fastq \\
        | sentieon util sort \$bam_option -r $fasta -o ${prefix}.bam -t $task.cpus --sam2bam -i -

    # genome
    #sentieon-cli \
    #    -r GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta \
    #    --r1-fastq my_sample.r1.fastq.gz \
    #    --r2-fastq my_sample.r2.fastq.gz \
    #    --readgroups '@RG\tID:my_sample-1\tSM:my_sample\tPL:ELEMENT\tLB:my_sample-lib1' \
    #    --model-bundle DNAscopeElementBioWGS2.0.bundle \
    #    --bed canonical_contigs.bed \
    #    --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    #    --pcr-free \
    #    --assay WGS \
    #    my_sample.vcf.gz

    # exome
    #sentieon-cli \
    #    -r GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta \
    #    --r1-fastq my_sample.r1.fastq.gz \
    #    --r2-fastq my_sample.r2.fastq.gz \
    #    --readgroups '@RG\tID:my_sample-1\tSM:my_sample\tPL:ELEMENT\tLB:my_sample-lib1' \
    #    --model-bundle DNAscopeElementBioWGS2.0.bundle \
    #    --bed targets_hg38.bed \
    #      --interval_padding 200 \
    #    --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    #    --assay WES \
    #    my_sample.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo /opt/sentieon/sentieon-genomics-* 2>&1 | sed 's/\\/opt\\/sentieon\\/sentieon-genomics-//')
    END_VERSIONS
    """
}
