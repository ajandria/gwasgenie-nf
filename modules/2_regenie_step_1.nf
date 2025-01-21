process REGENIE_STEP_1 {

    tag "${phenotype}"

    conda "${moduleDir}/envs/regenie.yaml"

    input:
    tuple val(phenotype), path(phenos), path(covs), val(header), path(bed), path(bim), path(fam)

    output:
    path("*"), emit: s1

    script:
    """
    plink --bfile ${header[0]} --mac 100 --write-snplist --out snps_pass

    regenie \
    --step 1 --force-step1 \
    --bed ${header} \
    --phenoFile ${phenos} \
    --covarFile ${covs} \
    --extract snps_pass.snplist \
    --bsize 1000 \
    --qt --lowmem \
    --lowmem-prefix tmp_${phenotype}_regenie-step_1 \
    --out ${phenotype}_regenie-step_1 \
    --threads $task.cpus \
    --verbose
    """
}