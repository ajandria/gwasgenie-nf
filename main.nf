#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ajandria/gwasgenie-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ajandria/gwasgenie-nf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_PHENOS } from './modules/1_extract_phenos'
include { REGENIE_STEP_1 } from './modules/2_regenie_step_1'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GWASGENIE {

    take:
    gwas_sheet // channel: samplesheet read in from --gwas_sheet
    phenos
    qced_genotypes

    main:

    //
    // WORKFLOW: Run pipeline
    //
    EXTRACT_PHENOS (
        gwas_sheet,
        phenos
    )

    pheno_covs = EXTRACT_PHENOS.out.pheno
        .flatMap { files -> // Unpack the list
            files.collect { file -> 
                def prefix = file.getBaseName().replace('_pheno', '')
                [prefix, file]
            }
        }
        .join(
            EXTRACT_PHENOS.out.covs.flatMap { files -> // Unpack the list
                files.collect { file -> 
                    def prefix = file.getBaseName().replace('_covariates', '')
                    [prefix, file]
                }
            }
        )
        .map { prefix, phenoFile, covFile ->
            [prefix, phenoFile, covFile]
        }
        .groupTuple()

    pheno_covs.view { prefix, pheno, covariates ->
        log.info """
        \u001B[1m\u001B[34mMerged Result: ${prefix}\u001B[0m
        \u001B[32m------------------------------------------------------------\u001B[0m
        \u001B[33mPhenotype:  \u001B[0m ${pheno}
        \u001B[33mCovariates: \u001B[0m ${covariates}
        \u001B[32m------------------------------------------------------------\u001B[0m
        """
    }

    pheno_covs.view()
    REGENIE_STEP_1 (
        pheno_covs,
        qced_genotypes
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // WORKFLOW: Run main workflow
    //
    GWASGENIE (
        params.gwas_sheet,
        params.phenos,
        params.qced_genotypes
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/