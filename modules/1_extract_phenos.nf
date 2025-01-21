process EXTRACT_PHENOS {
    tag "Gwas sheet: ${gwas_sheet} | Phenos: ${phenos}"

    conda "${moduleDir}/envs/extract_phenos.yaml"

    input:
    path(gwas_sheet)
    path(phenos)

    output:
    path("*_pheno.txt")     , emit: pheno
    path("*_covariates.txt"), emit: covs

    script:
    template "1_extract_phenos.R"
}