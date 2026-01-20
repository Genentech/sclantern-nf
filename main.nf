process bam_for_cv {
    input:
    tuple val(sample_name), path(in_bam)
    path ref_fa

    output:
    tuple(val(sample_name), path("${sample_name}.tagged.bam"))

    script:
    """
    bam_for_cv.py \
        --sample_name $sample_name \
        --bam $in_bam \
        --genome $ref_fa  \
        --outbam ${sample_name}.tagged.bam \
        --thread $task.cpus
    """
}

process bam_to_fq {
    cpus 1

    input:
    tuple(val(sample_name), path(in_bam))

    output:
    tuple(val(sample_name), path("${sample_name}.bam2fq.fq"))

    script:
    """
    samtools bam2fq $in_bam > ${sample_name}.bam2fq.fq
    """
}

process run_minimap {
    cpus params.bigJobCpusMult * params.defaultParallelCpus
    memory params.bigJobMemMult * params.defaultMemory

    input:
    tuple(val(sample_name), path(in_fq))
    path ref_fa

    output:
    tuple(val(sample_name), path("${sample_name}.mapped.bam"))

    script:
    """
    minimap2 \
        -ax splice -uf -C5 \
        -t $task.cpus \
        --secondary=no \
        $ref_fa \
        $in_fq \
        > tmp.sam

    samtools view -bSh -F 2308 tmp.sam > ${sample_name}.mapped.bam
    """
}

process sort_index_remapped {
    input:
    tuple(val(sample_name), path(mapped_bam))

    output:
    tuple(val(sample_name),
          path("${sample_name}.mapped.sorted.bam"),
          path("${sample_name}.mapped.sorted.bam.bai"))

    script:
    """
    samtools sort -@ $task.cpus -o ${sample_name}.mapped.sorted.bam $mapped_bam
    samtools index -@ $task.cpus ${sample_name}.mapped.sorted.bam
    """
}

process gatk_ref_dict {
    cpus 1

    input:
    path in_fa

    output:
    tuple(path("ref.fa"), path("ref.fa.fai"), path("ref.dict"))

    script:
    """
    cp $in_fa ref.fa
    samtools faidx ref.fa
    gatk CreateSequenceDictionary -R ref.fa
    """
}

// TODO: is gatk still needed now that minimap2 has preset options for
// scisoseq?

process gatk_ncigar {
    cpus params.bigJobCpusMult * params.defaultParallelCpus
    memory params.bigJobMemMult * params.defaultMemory

    input:
    tuple(val(sample_name), path(sorted_bam), path(sorted_bam_bai))
    tuple(path(ref_fa), path(ref_faidx), path(ref_dict))

    output:
    tuple(val(sample_name),
          path(sorted_bam), path(sorted_bam_bai),
          path("${sample_name}.sncr.bam"), path("${sample_name}.sncr.bai"))

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}G -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}" \
        SplitNCigarReads \
        -R $ref_fa \
        -I $sorted_bam \
        -O ${sample_name}.sncr.bam
    """
}

process flag_correction {
    cpus params.bigJobCpusMult * params.defaultParallelCpus
    memory params.bigJobMemMult * params.defaultMemory
    
    input:
    tuple(val(sample_name),
          path(orig_bam), path(orig_bam_bai),
          path(sncr_bam), path(sncr_bam_bai))

    output:
    tuple(val(sample_name),
          path("${sample_name}.fc.bam"),path("${sample_name}.fc.bam.bai"))

    script:
    """
    flagCorrection.py $orig_bam $sncr_bam ${sample_name}.fc.bam $task.cpus
    samtools index -@ $task.cpus ${sample_name}.fc.bam
    """
}

process trgt {
    cpus params.bigJobCpusMult * params.defaultParallelCpus
    memory params.bigJobMemMult * params.defaultMemory
    
    input:
    tuple(val(sample_name), path(fc_bam), path(fc_bam_bai))
    tuple(path(ref_fa), path(ref_faidx), path(ref_dict))
    path(repeats_bed)

    output:
    tuple(val(sample_name),
          path("${sample_name}.trgt.spanning.sorted.bam"),
          path("${sample_name}.trgt.spanning.sorted.bam.bai"),
          path("${sample_name}.trgt.sorted.vcf.gz"),
          path("${sample_name}.trgt.sorted.vcf.gz.csi"))

    script:
    """
    trgt genotype \
        --genome $ref_fa \
        --repeats $repeats_bed \
        --reads $fc_bam \
        --output-prefix trgt

    bcftools sort -Ob -o ${sample_name}.trgt.sorted.vcf.gz trgt.vcf.gz
    bcftools index ${sample_name}.trgt.sorted.vcf.gz

    samtools sort -o ${sample_name}.trgt.spanning.sorted.bam trgt.spanning.bam
    samtools index ${sample_name}.trgt.spanning.sorted.bam
    """
}

process filter_vcf {
    cpus 1

    input:
    tuple(val(sample_name),
          path(trgt_bam), path(trgt_bai),
          path(trgt_vcf), path(trgt_csi))

    output:
    path("${sample_name}.filtered.vcf.gz"), emit: vcf
    path("${sample_name}.filtered.vcf.gz.tbi"), emit: tbi

    script:
    """
    vcf_filter.py \
        --in_vcf $trgt_vcf \
        --out_vcf filtered.vcf \
        --mode trgt

    bcftools view filtered.vcf -Oz -o ${sample_name}.filtered.vcf.gz
    tabix ${sample_name}.filtered.vcf.gz
    """
}

// todo: should we use trgt merge instead of bcftools merge?
process merge_vcf {
    input:
    path vcf_files
    path tbi_files

    output:
    path("merged.vcf.gz")

    script:
    """
    bcftools merge -Oz -o merged.vcf.gz $vcf_files
    """
}

process count_read_alleles {
    publishDir params.outdir, mode: 'copy'

    input:
    tuple(val(sample_name),
          path(trgt_bam), path(trgt_bai),
          path(trgt_vcf), path(trgt_csi))
    tuple(path(ref_fa), path(ref_faidx), path(ref_dict))
    path variants_vcf

    output:
    path "${sample_name}.allele_read_counts.csv"

    script:
    """
    allele_read_count.py \
        --bam_file $trgt_bam \
        --genome_fasta $ref_fa \
        --variants_vcf $variants_vcf \
        --out_count_path ${sample_name}.allele_read_counts.csv
    """
}

workflow {
    samples_split = Channel.fromPath(params.sample_sheet).splitCsv(header: true)
    in_bam_ch = samples_split.map{ row -> tuple(row.sample_name, file(row.path))}

    bam_for_cv(in_bam_ch, params.ref_fa)
    bam_to_fq(bam_for_cv.out)
    run_minimap(bam_to_fq.out, params.ref_fa)
    sort_index_remapped(run_minimap.out)

    gatk_ref_dict(params.ref_fa)

    gatk_ncigar(sort_index_remapped.out, gatk_ref_dict.out)

    flag_correction(gatk_ncigar.out)
    trgt(flag_correction.out, gatk_ref_dict.out, params.repeats_bed)

    filter_vcf(trgt.out)

    // TODO: skip this step if only 1 sample
    merge_vcf(filter_vcf.out.vcf.collect(), filter_vcf.out.tbi.collect())

    count_read_alleles(trgt.out, gatk_ref_dict.out, merge_vcf.out)

    // FIXME: need to add the sample name to the cells in the final output table? Or possibly, don't append the sample name to the cell barcode in the bam, I don't think it's needed anymore and may complicate things if the pipeline is run iteratively (clustering and then rerun with the clusters as samples)
}
