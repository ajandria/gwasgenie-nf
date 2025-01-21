process REGENIE_STEP_2 {

    tag "${phenotype}_${chrom}"

    conda "${moduleDir}/envs/regenie.yaml"

    input:
    tuple val(phenotype), path(phenos), path(covs), val(header), path(bed), path(bim), path(fam), path(pred_s1), val(chrom), path(bgen), path(sample_bgen)

    output:
    path("${phenotype}/${chrom}_${phenotype}_regenie_step_2*"), emit: s2

    script:
    """
    mkdir ${phenotype}
    regenie \
        --step 2 \
        --bgen ${bgen} \
        --sample ${sample_bgen} \
        --ref-first \
        --af-cc \
        --phenoFile ${phenos} \
        --covarFile ${covs} \
        --pred ${pred_s1} \
        --bsize 400 --qt --firth --approx --firth-se --pThresh 0.999 --minMAC 5 \
        --test additive \
        --verbose \
        --threads $task.cpus \
        --out ${phenotype}/${chrom}_${phenotype}_regenie_step_2
    """
}
