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
    tuple(val(meta), path("${in_bam.baseName}.fq.gz"))

    script:
    """
    samtools bam2fq $in_bam > ${in_bam.baseName}.fq
    gzip ${in_bam.baseName}.fq
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
        | samtools view -bSh -F 2308 \
        > ${in_fq.baseName}.remapped.bam
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

process samtools_merge_flag_correction {
    label "lowMem"
    label "parallel"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple(val(sample),
          path(bam_files),
          path(bai_files))

    output:
    path("${sample}.flag_corrected.bam")
    path("${sample}.flag_corrected.bam.bai")
         
    
    script:
    """
    samtools merge -@ $task.cpus -o tmp.bam ${bam_files}
    samtools sort -@ $task.cpus -o ${sample}.flag_corrected.bam tmp.bam
    samtools index -@ $task.cpus ${sample}.flag_corrected.bam
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

process deepv {
    //NOTE memory usage grows with num threads
    //label "veryParallel"
    label "parallel"
    label "highMem"
    label "longJob"

    container "docker://google/deepvariant:1.10.0"
    
    input:
    tuple(val(chrom),
          val(meta), path(fc_bam), path(fc_bam_bai),
          path(ref_fa), path(ref_faidx), path(ref_dict),
          path(repeats_bed))

    output:
    tuple(val(meta),
          path("${fc_bam.baseName}.deepv.vcf.gz"),
          path("${fc_bam.baseName}.deepv.vcf.gz.tbi"))

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref "${ref_fa}" \
        --reads "${fc_bam}" \
        --output_vcf "${fc_bam.baseName}.deepv.vcf.gz" \
        --sample_name ${meta.sample} \
        --num_shards "$task.cpus" # Use allocated cores for sharding
    """
}

include { filter_vcf as filter_vcf_trgt } from './modules/filter_vcf.nf'
include { filter_vcf as filter_vcf_deepv } from './modules/filter_vcf.nf'

include {merge_and_count as merge_and_count_trgt} from "./workflows/merge_and_count.nf"

include {merge_and_count as merge_and_count_deepv} from "./workflows/merge_and_count.nf"

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

    flag_correction_bam = flag_correction.out.map{
        meta, bam, bai ->
        tuple(meta.sample, bam)
    }.groupTuple()

    flag_correction_bai = flag_correction.out.map{
        meta, bam, bai ->
        tuple(meta.sample, bai)
    }.groupTuple()

    samtools_merge_flag_correction(
        flag_correction_bam.combine(
            flag_correction_bai,
            by: 0
        )
    )

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

    filter_vcf_trgt("trgt", trgt.out.map{
        meta, bam, bai, vcf, csi ->
        tuple(meta, vcf, csi)
    })

    merge_and_count_trgt(
        filter_vcf_trgt.out,
        trgt.out.map{
            meta, bam, bai, vcf, csi ->
            tuple(meta, bam, bai)
        },
        gatk_ref_dict.out,
        "trgt"
    )

    deepv(
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

    filter_vcf_deepv("dv", deepv.out)

    merge_and_count_deepv(
        filter_vcf_deepv.out,
        flag_correction.out,
        gatk_ref_dict.out,
        "deepv"
    )
}
