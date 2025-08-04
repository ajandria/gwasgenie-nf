process REGENIE_STEP_1 {

    tag "${phenotype}"

    conda "${moduleDir}/envs/regenie.yaml"

    input:
    tuple val(phenotype), path(phenos), path(covs), val(header), path(bed), path(bim), path(fam)

    output:
    tuple val(phenotype), path("${phenotype}_regenie-step_1_pred.list"), emit: s1

    script:
    """
    plink --pfile ${header} --mac 100 --write-snplist --out snps_pass

    awk 'BEGIN {OFS="\t"} NR==1 {print "#FID","IID","SEX"; next} {print 0, \$1, \$2}' ${header}.psam > ${header}.sample
    mv ${header}.sample ${header}.psam

    regenie \
    --step 1 --force-step1 \
    --pgen ${header} \
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
