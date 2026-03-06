
process filter_vcf {
    label "lowMem"

    input:
    val(mode)
    tuple(val(meta),
          path(in_vcf), path(in_tbi_or_csi))

    output:
    tuple(val(meta),
          path("${in_vcf.baseName}.filtered.gz"),
          path("${in_vcf.baseName}.filtered.gz.csi"))

    script:
    """
    vcf_filter.py \
        --in_vcf $in_vcf \
        --out_vcf tmp.vcf \
        --mode $mode

    bcftools view tmp.vcf -Oz -o ${in_vcf.baseName}.filtered.gz
    bcftools index ${in_vcf.baseName}.filtered.gz
    """
}
