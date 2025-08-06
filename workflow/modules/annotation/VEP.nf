process VEP_ANNOTATION {
    input:
        path vcf_file
    output:
        path "*.vep.vcf", emit: vep_out

    publishDir "results/annotation/", mode: 'copy'

    container = 'ensemblorg/ensembl-vep'

    script:
    """
    vep -i ${vcf_file} --assembly GRCh38 --vcf --database -o ${vcf_file.simpleName}.vep.vcf
    """
}
