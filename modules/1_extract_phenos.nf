process EXTRACT_PHENOS {
    tag "${gwas_sheet}"

    conda 'envs/gwasgenie.yaml'

    input:
    path(gwas_sheet)

    output:

    script:
    """
    echo "Placeholder"
    """
}