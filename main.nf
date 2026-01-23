include { sort_index as sort_index_orig } from './modules/samtools_sort_index.nf'

process split_bam_by_chrom {
    label "shortJob"
    label "tinyMem"
    
    input:
    tuple val(meta), path(in_bam), path(in_bam_bai)
    each chrom

    output:
    tuple(val(meta), val(chrom),
          path("${meta.sample}.${chrom}.bam"))

    script:
    """
    samtools view \
        -o "${meta.sample}.${chrom}.bam" \
        ${in_bam} \
        ${chrom}
    """
}

process bam_for_cv {
    label "lowMem"
    label "parallel"

    input:
    tuple val(meta), path(in_bam)
    path ref_fa

    output:
    tuple(val(meta), path("${in_bam.baseName}.tagged.bam"))

    script:
    """
    bam_for_cv.py \
        --sample_name ${meta.sample} \
        --bam $in_bam \
        --genome $ref_fa  \
        --outbam ${in_bam.baseName}.tagged.bam \
        --thread $task.cpus
    """
}

process bam_to_fq {
    label "tinyMem"
    label "shortJob"

    input:
    tuple(val(meta), path(in_bam))

    output:
    tuple(val(meta), path("${in_bam.baseName}.fq"))

    script:
    """
    samtools bam2fq $in_bam > ${in_bam.baseName}.fq
    """
}

process index_orig_fasta {
    label "tinyMem"
    label "shortJob"

    input:
    path(ref_fa)

    output:
    path("${ref_fa}.fai")

    script:
    """
    samtools faidx ${ref_fa}
    """
}

process split_fasta {
    label "tinyMem"
    label "shortJob"

    input:
    path(ref_fa)
    path(ref_fa_fai)
    each chrom

    output:
    tuple(val(chrom),
          path("${ref_fa.baseName}.${chrom}.fasta"),
          path("${ref_fa.baseName}.${chrom}.fasta.fai"))

    script:
    """
    samtools faidx ${ref_fa} ${chrom} > ${ref_fa.baseName}.${chrom}.fasta
    samtools faidx ${ref_fa.baseName}.${chrom}.fasta
    """
}

process split_bed {
    label "tinyMem"
    label "shortJob"

    input:
    path(repeats_bed)
    each chrom

    output:
    tuple(val(chrom),
          path("${repeats_bed.baseName}.${chrom}.bed"))

    script:
    """
    grep '^${chrom}[[:space:]]' ${repeats_bed} > ${repeats_bed.baseName}.${chrom}.bed
    """
}

process run_minimap {
    label "highMem"
    label "veryParallel"

    input:
    tuple(val(chrom), val(meta), path(in_fq), path(ref_fa), path(ref_fa_fai))

    output:
    tuple(val(meta), path("${in_fq.baseName}.remapped.bam"))

    script:
    """
    minimap2 \
        -ax splice -uf -C5 \
        -t $task.cpus \
        --secondary=no \
        $ref_fa \
        $in_fq \
        > tmp.sam

    samtools view -bSh -F 2308 tmp.sam > ${in_fq.baseName}.remapped.bam
    """
}

include { sort_index as sort_index_remapped } from './modules/samtools_sort_index.nf'

process gatk_ref_dict {
    label "lowMem"

    input:
    tuple(val(chrom), path(in_fa), path(in_fa_fai))

    output:
    tuple(val(chrom), path(in_fa), path(in_fa_fai), path("${in_fa.baseName}.dict"))

    script:
    """
    gatk CreateSequenceDictionary -R $in_fa
    """
}

// TODO: is gatk still needed now that minimap2 has preset options for
// scisoseq?

process gatk_ncigar {
    label "veryParallel"
    label "highMem"

    input:
    tuple(
        val(chrom),
        val(meta), path(sorted_bam), path(sorted_bam_bai),
        path(ref_fa), path(ref_faidx), path(ref_dict)
    )

    output:
    tuple(val(meta), path(sorted_bam), path(sorted_bam_bai),
          path("${sorted_bam.baseName}.sncr.bam"),
          path("${sorted_bam.baseName}.sncr.bai"))

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}G -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}" \
        SplitNCigarReads \
        -R $ref_fa \
        -I $sorted_bam \
        -O ${sorted_bam.baseName}.sncr.bam
    """
}

process flag_correction {
    label "veryParallel"
    label "highMem"
    
    input:
    tuple(val(meta),
          path(orig_bam), path(orig_bam_bai),
          path(sncr_bam), path(sncr_bam_bai))

    output:
    tuple(val(meta),
          path("${sncr_bam.baseName}.fc.bam"),
          path("${sncr_bam.baseName}.fc.bam.bai"))

    script:
    """
    flagCorrection.py $orig_bam $sncr_bam ${sncr_bam.baseName}.fc.bam $task.cpus
    samtools index -@ $task.cpus ${sncr_bam.baseName}.fc.bam
    """
}

process trgt {
    label "veryParallel"
    label "highMem"
    label "longJob"
    
    input:
    tuple(val(chrom),
          val(meta), path(fc_bam), path(fc_bam_bai),
          path(ref_fa), path(ref_faidx), path(ref_dict),
          path(repeats_bed))

    output:
    tuple(val(meta),
          path("${fc_bam.baseName}.trgt.spanning.sorted.bam"),
          path("${fc_bam.baseName}.trgt.spanning.sorted.bam.bai"),
          path("${fc_bam.baseName}.trgt.sorted.vcf.gz"),
          path("${fc_bam.baseName}.trgt.sorted.vcf.gz.csi"))

    script:
    """
    trgt genotype \
        --genome $ref_fa \
        --repeats $repeats_bed \
        --reads $fc_bam \
        --output-prefix trgt

    bcftools sort -Ob -o ${fc_bam.baseName}.trgt.sorted.vcf.gz trgt.vcf.gz
    bcftools index ${fc_bam.baseName}.trgt.sorted.vcf.gz

    samtools sort -o ${fc_bam.baseName}.trgt.spanning.sorted.bam trgt.spanning.bam
    samtools index ${fc_bam.baseName}.trgt.spanning.sorted.bam
    """
}

process filter_vcf {
    label "lowMem"

    input:
    tuple(val(meta),
          path(trgt_bam), path(trgt_bai),
          path(trgt_vcf), path(trgt_csi))

    output:
    tuple(val(meta),
          path("${trgt_vcf.baseName}.filtered.gz"),
          path("${trgt_vcf.baseName}.filtered.gz.csi"))

    script:
    """
    vcf_filter.py \
        --in_vcf $trgt_vcf \
        --out_vcf tmp.vcf \
        --mode trgt

    bcftools view tmp.vcf -Oz -o ${trgt_vcf.baseName}.filtered.gz
    bcftools index ${trgt_vcf.baseName}.filtered.gz
    """
}

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
    label = "midMem"
    publishDir params.outdir, mode: 'copy'

    input:
    path(tsv_files)

    output:
    path("combined_allele_read_count.tsv")

    script:
    """
    rbind_df_tsvs.py combined_allele_read_count.tsv ${tsv_files}
    """
}

workflow {
    chroms = params.chroms.split(',').collect { "chr$it" }
    
    samples_split = Channel.fromPath(params.sample_sheet).splitCsv(header: true)
    in_bam_ch = samples_split.map{
        row -> tuple([sample: row.sample_name], file(row.path))
    }

    sort_index_orig(in_bam_ch)
    
    split_bam_by_chrom(sort_index_orig.out, chroms)

    // FIXME is adding the sample to the qname really necessary?
    // Worried about downstream issues if we iteratively rerun this on
    // clustered cells (and keep appending the new sample names)
    bam_for_cv(split_bam_by_chrom.out.map{
        meta, chrom, bam ->
        tuple(meta + [chrom: chrom], bam)
    }, params.ref_fa)

    bam_to_fq(bam_for_cv.out)

    index_orig_fasta(params.ref_fa)
    split_fasta(params.ref_fa, index_orig_fasta.out, chroms)

    run_minimap(
        bam_to_fq.out.map{
            meta, fq ->
            tuple(meta.chrom, meta, fq)
        }.combine(
            split_fasta.out,
            by: 0
        )
    )

    sort_index_remapped(run_minimap.out)

    gatk_ref_dict(split_fasta.out)

    gatk_ncigar(
        sort_index_remapped.out.map{
            meta, bam, bai ->
            tuple(meta.chrom, meta, bam, bai)
        }.combine(
            gatk_ref_dict.out,
            by: 0
        )
    )

    flag_correction(gatk_ncigar.out)

    split_bed(params.repeats_bed, chroms)

    trgt(
        flag_correction.out.map{
            meta, bam, bai ->
            tuple(meta.chrom, meta, bam, bai)
        }.combine(
            gatk_ref_dict.out,
            by: 0
        ).combine(
            split_bed.out,
            by: 0
        )
    )

    filter_vcf(trgt.out)

    grouped_vcf = filter_vcf.out.map{
        meta, vcf, csi ->
        tuple(meta.chrom, vcf)
    }.groupTuple()

    grouped_csi = filter_vcf.out.map{
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
        trgt.out.map{
            meta, bam, bai, vcf, csi ->
            tuple(meta.chrom, meta, bam, bai)
        }.combine(
            gatk_ref_dict.out,
            by: 0
        ).combine(
            merge_vcf.out,
            by: 0
        )
    )

    rbind_allele_read_count(
        count_read_alleles.out.map{
            meta, tsv -> tsv
        }.collect()
    )
}
