// todo: should we use trgt merge instead of bcftools merge?
process merge_vcf {
    label "lowMem"

    input:
    tuple(val(chrom), path(vcf_files), path(csi_files))

    output:
    tuple(val(chrom),
          path("${chrom}.merged.vcf.gz"),
          path("${chrom}.merged.vcf.gz.csi"))

    script:
    """
    bcftools merge -Oz -o ${chrom}.merged.vcf.gz $vcf_files
    bcftools index ${chrom}.merged.vcf.gz
    """
}

process count_read_alleles {
    label "lowMem"

    input:
    tuple(val(chrom),
          val(meta), path(trgt_bam), path(trgt_bai),
          path(ref_fa), path(ref_faidx), path(ref_dict),
          path(variants_vcf), path(variants_csi))

    output:
    tuple(val(meta), path("${meta.sample}.allele_read_count.${chrom}.tsv"))

    script:
    """
    allele_read_count.py \
        --bam_file $trgt_bam \
        --genome_fasta $ref_fa \
        --variants_vcf $variants_vcf \
        --sample_name ${meta.sample} \
        --out_count_path ${meta.sample}.allele_read_count.${chrom}.tsv
    """
}

process rbind_allele_read_count {
    label "midMem"
    publishDir params.outdir, mode: 'copy'

    input:
    path(tsv_files)
    val(out_prefix)

    output:
    path("${out_prefix}.combined_allele_read_count.tsv")

    script:
    """
    rbind_df_tsvs.py ${out_prefix}.combined_allele_read_count.tsv ${tsv_files}
    """
}

workflow merge_and_count {
    take:
    vcf_csi_ch
    bam_bai_ch
    gatk_ref_dict_ch
    out_prefix

    main:
    grouped_vcf = vcf_csi_ch.map{
        meta, vcf, csi ->
        tuple(meta.chrom, vcf)
    }.groupTuple()

    grouped_csi = vcf_csi_ch.map{
        meta, vcf, csi ->
        tuple(meta.chrom, csi)
    }.groupTuple()

    //// TODO: skip this step if only 1 sample
    merge_vcf(
        grouped_vcf.combine(
            grouped_csi,
            by: 0
        )
    )

    count_read_alleles(
        bam_bai_ch.map{
            meta, bam, bai ->
            tuple(meta.chrom, meta, bam, bai)
        }.combine(
            gatk_ref_dict_ch,
            by: 0
        ).combine(
            merge_vcf.out,
            by: 0
        )
    )

    rbind_allele_read_count(
        count_read_alleles.out.map{
            meta, tsv -> tsv
        }.collect(),
        out_prefix
    )
}
