process EXTRACT_PHENOS {
    tag "Gwas sheet: ${gwas_sheet} | Phenos: ${phenos}"

    conda "${moduleDir}/envs/extract_phenos.yaml"

    input:
    path(gwas_sheet)
    path(phenos)
    path(pcs_path)

    output:
    path("*_pheno.txt"), emit: pheno
    path("*_covs.txt") , emit: covs

    script:
    template "1_extract_phenos.R"
}