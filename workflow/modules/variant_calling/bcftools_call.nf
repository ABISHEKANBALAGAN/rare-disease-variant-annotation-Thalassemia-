process BCFTOOLS_CALL {
    tag "$sample_id"
    publishDir "results/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.vcf")

    script:
    """
    echo "BAM: $sorted_bam"
    echo "REF: $reference"
    ls -l
    bcftools mpileup -Ou -f $reference $sorted_bam | \
    bcftools call -mv -Ov -o ${sample_id}.vcf
    ls -l
    """
}
