#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Pull in igenomes
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai = params.fasta+".fai"
params.bwa = WorkflowMain.getGenomeAttribute(params, 'bwa')

log.info """\
 ======================================================
            S E N T I E O N - C L I - N F 
 ======================================================
 input: ${params.input}
 outdir: ${params.outdir}
 ======================================================
 genome: ${params.genome}
 sentieon_ml_model: ${params.sentieon_ml_model}
 sentieon_license: ${params.sentieon_license}
 assay: ${params.assay}
 pcr: ${params.pcr}
 known_sites: ${params.known_sites}
 ======================================================
 fasta: ${params.fasta}
 fasta_fai: ${params.fasta_fai}
 bwa: ${params.bwa}
 ignore_samples: ${params.ignore_samples}
 igenomes_base: ${params.igenomes_base}
 igenomes_ignore: ${params.igenomes_ignore} 
 ======================================================
 custom_config_version: ${params.custom_config_version}
 custom_config_base: ${params.custom_config_base}
 ======================================================
 """

// Check mandatory parameters
def checkPathParamList = [ params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check input path parameters to see if they exist
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

include { SENTIEON_CLI } from './modules/local/sentieon/sentieon-cli'
include { MULTIQC } from './modules/nf-core/multiqc/main'

def model_file = params.sentieon_ml_model ? file(params.sentieon_ml_model, checkIfExists: true) : [] 
def target_region_bed = params.target_region_bed ? file(params.target_region_bed, checkIfExists: true) : [] 

workflow {  

    ch_versions = Channel.empty()

    ch_genome = [params.fasta, params.fasta_fai]
    
    Channel.value(ch_input)
        .splitCsv ( header:true, sep:',' )
        .set { sheet }

    ch_fastq = sheet
        .filter { row -> !params.ignore_samples.contains(row.sample) } // Skip matching samples
        .map { row -> [[row.sample], row] }
        .groupTuple()
        .map { meta, rows ->
            [rows, rows.size()]
        }
        .transpose()
        .map { row, numLanes ->
            create_fastq_channel(row + [num_lanes:numLanes])
        }
            
    ch_fastq
        .map{ meta, r1_fastq, r2_fastq ->
            grouped_id = meta.sample
            grouped_prefix = meta.id
            grouped_num_lanes = meta.num_lanes
            grouped_meta = [id:grouped_id, prefix: grouped_prefix, read_group: grouped_id, num_lanes: grouped_num_lanes]
            [grouped_meta, meta, r1_fastq, r2_fastq]
            }
        .groupTuple()
        .set { ch_grouped_fastq }
    
    // convert known sites to ch
    ch_known_sites = Channel.of(params.known_sites)

    // fastq -> bam (fq2bam)
    // turn off cli multiqc, pass outputs to below
    SENTIEON_CLI (
        ch_grouped_fastq,
        params.bwa,
        ch_genome,
        model_file,
        params.assay,
        target_region_bed
    )

    // MultiQC
    // get multiqc conf files
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.fromPath("$projectDir/assets/Element_Biosciences_Logo_Black_RGB.png", checkIfExists: true)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(SENTIEON_CLI.out.metrics.collect{it[1]}.ifEmpty([]))
    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

}

def create_fastq_channel(LinkedHashMap row) {

    def meta = [
        id: row.sample,
        sample: row.sample,
        prefix: row.sample + "__" + row.read_group,
        read_group: row.read_group,
        platform: row.platform,
        gender: row.gender,
        num_lanes: row.num_lanes,
        single_end: false  // Default to paired-end
    ]
    
    def fields = [
        'r1_fastq': ['meta': [:], 'read_num': 'R1'],
        'r2_fastq': ['meta': [:], 'read_num': 'R2'],
        'fastq_1': ['meta': [:], 'read_num': 'R1'],
        'fastq_2': ['meta': [:], 'read_num': 'R2']
    ]

    // Add paths of the fastq files to the meta map
    def fastq_files = []

    fields.each { key, value ->
        if (row[key]) {
            def file_path = file(row[key])
            if (!file_path.exists()) {
                error("ERROR: Please check input samplesheet -> ${value.read_num} FastQ file does not exist!\n${row[key]}")
            }
        }
    }

    // Set r1_fastq and r2_fastq explicitly
    def r1_fastq = null
    def r2_fastq = null
    
    // Validate R1 fastq file
    if (row.r1_fastq || row.fastq_1) {
        r1_fastq = file(row.r1_fastq ? row.r1_fastq : row.fastq_1)
        if (!r1_fastq.exists()) {
            error("ERROR: Please check input samplesheet -> R1 FastQ file does not exist!\n${r1_fastq}")
        }
    } else {
        error("ERROR: R1 FastQ file is required but not found in the samplesheet for sample ${row.sample}")
    }

    // Validate R1 fastq file
    if (row.r2_fastq || row.fastq_2) {
        r2_fastq = file(row.r2_fastq ? row.r2_fastq : row.fastq_2)
        if (!r2_fastq.exists()) {
            error("ERROR: Please check input samplesheet -> R2 FastQ file does not exist!\n${r2_fastq}")
        }
    } else {
        error("ERROR: R2 FastQ file is required but not found in the samplesheet for sample ${row.sample}")
    }

    // Return the meta and the explicit r1 and r2 fastq files
    return [meta, r1_fastq, r2_fastq]
}