nextflow.enable.dsl=2

include { RUN_FASTQC } from './workflow/modules/qc/fastqc.nf'

workflow {

    // Expand all FASTQ files individually
    samples_ch = Channel
        .fromPath(params.reads)
        .map { file_path ->
            def sample_id = file_path.getBaseName().replaceFirst(/_\d+$/, '')
            tuple(sample_id, file_path)
        }

    // Run FastQC per file
    RUN_FASTQC(samples_ch)
}


