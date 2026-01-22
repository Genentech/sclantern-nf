
process sort_index {
    label "lowMem"
    label "parallel"

    input:
    tuple(val(meta), path(in_bam))

    output:
    tuple(val(meta),
          path("${in_bam.baseName}.sorted.bam"),
          path("${in_bam.baseName}.sorted.bam.bai"))

    script:
    """
    samtools sort -@ $task.cpus -o ${in_bam.baseName}.sorted.bam $in_bam
    samtools index -@ $task.cpus ${in_bam.baseName}.sorted.bam
    """
}
