nextflow.enable.dsl=2

// Load the VEP module
include { VEP_ANNOTATION } from './workflow/modules/annotation/VEP.nf'

// Define input
workflow {
    // Channel to load your filtered VCF file
    vcf_ch = Channel.fromPath("results/variants/SRR32809192.filtered.vcf")

    // Call the VEP annotation process
    VEP_ANNOTATION(vcf_ch)
}
