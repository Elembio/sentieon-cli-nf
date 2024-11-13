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
 fasta: ${params.fasta}
 fasta_fai: ${params.fasta_fai}
 bwa: ${params.bwa}
 sentieon_ml_modle: ${params.sentieon_ml_model}
 known_sites: ${params.known_sites}
 igenomes_base: ${params.igenomes_base}
 igenomes_ignore: ${params.igenomes_ignore} 
 ======================================================
 """

// Check mandatory parameters
def checkPathParamList = [ params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check input path parameters to see if they exist
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

include { SENTIEON_CLI } from './modules/local/sentieon/sentieon-cli'

def model_file = params.sentieon_ml_model ? file(params.sentieon_ml_model, checkIfExists: true) : [] 
def bed = params.bed ? file(params.bed, checkIfExists: true) : [] 

workflow {  

    ch_versions = Channel.empty()

    ch_genome = [params.fasta, params.fasta_fai]
    
    Channel.value(ch_input)
        .splitCsv ( header:true, sep:',' )
        .set { sheet }

    ch_fastq = sheet.map { row -> [[row.sample], row] }
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
    SENTIEON_CLI (
        ch_grouped_fastq,
        params.bwa,
        ch_genome,
        model_file,
        params.assay,
        bed
    )

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