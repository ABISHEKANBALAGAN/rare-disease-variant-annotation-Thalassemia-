process RUN_FASTQC {
    tag "$sample_id"
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${reads.join(' ')}
    """
}
