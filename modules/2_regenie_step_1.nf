process REGENIE_STEP_1 {

    tag "${phenotype}"

    conda "${moduleDir}/envs/regenie.yaml"

    input:
    tuple val(phenotype), path(phenos), path(covs)
    path(qced_genotypes)

    output:
    path("*"), emit: s1

    script:
    """
    plink --bfile ${qced_genotypes} --mac 100 --write-snplist --out snps_pass

    regenie \
    --step 1 --force-step1 \
    --bed ${qced_genotypes} \
    --phenoFile ${phenos} \
    --covarFile ${covs} \
    --extract snps_pass.snplist \
    --covarColList Gender,Age,BMI,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --phenoColList Total_30_normed \
    --bsize 1000 \
    --qt --lowmem \
    --lowmem-prefix tmp_${phenotype}_regenie-step_1 \
    --out ${phenotype}_regenie-step_1 \
    --threads $task.cpus \
    --verbose
    """
}