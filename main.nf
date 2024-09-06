#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Pull in igenomes
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
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
 bwa: ${params.bwa}
 known_sites: ${params.known_sites}
 igenomes_base: ${params.igenomes_base}
 igenomes_ignore: ${params.igenomes_ignore}
 use_elembio_igenomes: ${params.use_elembio_igenomes}
 ======================================================
 """

// Check mandatory parameters
def checkPathParamList = [ params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check input path parameters to see if they exist
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

include { SENTIEON_CLI } from './modules/local/sentieon/sentieon-cli'

def model_file = params.model_file ? file(params.model_file, checkIfExists: true) : [] 
def interval_bed = params.interval_bed ? file(params.model_file, checkIfExists: true) : [] 

workflow {  

    ch_versions = Channel.empty()

    ch_genome = [params.fasta]
    ch_bwa = params.fasta ? params.bwa ? Channel.fromPath(params.bwa).collect() : [] : []

    // parse input samplesheet
    Channel.value(ch_input)
            .splitCsv ( header:true, sep:',' )
            .set { sheet }

    // create ch_fastq
    ch_fastq = sheet.map { row -> [[row.sample], row] }
                .groupTuple()
                .map { meta, rows ->
                    [rows, rows.size()]
                }
                .transpose()
                .map { row, numLanes ->
                    create_fastq_channel(row + [num_lanes:numLanes])
                }

    // convert known sites to ch
    ch_known_sites = Channel.of(params.known_sites)
    
    // fastq -> bam (fq2bam)
    SENTIEON_CLI (
        ch_fastq,
        params.bwa,
        ch_genome,
        ch_known_sites
    )
    ch_versions = ch_versions.mix(SENTIEON_CLI.out.versions.first())

}

def create_fastq_channel(LinkedHashMap row) {

    def fields = [
        'r1_fastq': ['meta': [:], 'read_num': 'R1'],
        'r2_fastq': ['meta': [:], 'read_num': 'R2']
    ]

    def meta = [
        id: row.sample,
        sample: row.sample,
        prefix: row.sample + "__" + row.read_group,
        read_group: row.read_group,
        platform: row.platform,
        gender: row.gender,
        num_lanes: row.num_lanes,
        single_end: false
    ]

    // Add paths of the fastq files to the meta map
    def fastq_files = []

    fields.each { key, value ->
        if (row[key]) {
            def file_path = file(row[key])
            if (!file_path.exists()) {
                error("ERROR: Please check input samplesheet -> ${value.read_num} FastQ file does not exist!\n${row[key]}")
            }
            fastq_files << file_path
        }
    }

    // Determine if the read is single-ended
    if (!row.r2_fastq) {
        meta.single_end = true
    }

    // Return the meta and fastq files list
    return [meta, fastq_files]
}