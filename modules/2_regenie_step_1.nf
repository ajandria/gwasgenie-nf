process REGENIE_STEP_1 {

    tag "${phenotype}"

    conda "${moduleDir}/envs/regenie.yaml"

    input:
    tuple val(phenotype), path(phenos), path(covs), val(header), path(bed), path(bim), path(fam)

    output:
    tuple val(phenotype), path("${phenotype}_regenie-step_1_pred.list"), emit: s1

    script:
    """
    regenie \
        --step 2 \
        --bgen ${header[0]} \
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