process REGENIE_STEP_2 {

    tag "${phenotype}_${chrom}"

    conda "${moduleDir}/envs/regenie.yaml"

    input:
    tuple val(phenotype), path(phenos), path(covs), val(header), path(bed), path(bim), path(fam),
          val(chrom), path(bgen_file), path(sample_file), path(pred_file)

    output:
    path("${chrom}_${phenotype}_regenie_step_2.*"), emit: s2

    script:
    """
    regenie \
        --step 2 \
        --bgen ${bgen_file} \
        --sample ${sample_file} \
        --ref-first \
        --af-cc \
        --phenoFile ${phenos} \
        --covarFile ${covs} \
        --pred ${pred_file} \
        --bsize 400 --qt --firth --approx --firth-se --pThresh 0.999 --minMAC 5 \
        --test additive \
        --verbose \
        --threads $task.cpus \
        --out ${phenotype}/${chrom}_${phenotype}
    """
}
