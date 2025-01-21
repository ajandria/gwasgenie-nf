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
include { REGENIE_STEP_2 } from './modules/3_regenie_step_2'

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
    gwas_sheet
    phenos

    main:

    EXTRACT_PHENOS (
        gwas_sheet,
        phenos
    )

    pheno_covs = EXTRACT_PHENOS.out.pheno
        .flatMap { files -> // Unpack phenotypes
            files.collect { file ->
                def prefix = file.name.replace('_pheno.txt', '')
                [prefix, file]
            }
        }
        .join(
            EXTRACT_PHENOS.out.covs.flatMap { files -> // Unpack covariates
                files.collect { file -> 
                    def prefix = file.name.replace('_covs.txt', '')
                    [prefix, file]
                }
            }
        )
        .map { prefix, phenoFile, covFile ->
            def genotypesBase = params.qced_genotypes
            def header = genotypesBase.tokenize('/').last()
            def bedFile = file("${genotypesBase}.bed")
            def bimFile = file("${genotypesBase}.bim")
            def famFile = file("${genotypesBase}.fam")
            [prefix, phenoFile, covFile, header, bedFile, bimFile, famFile]
        }

    REGENIE_STEP_1 (
        pheno_covs
    )

    pheno_chrom_bgen = pheno_covs
        .join(REGENIE_STEP_1.out.s1) // Join phenotypes with REGENIE_STEP_1 output
        .flatMap { pheno_data, pred_file_data ->
            // Unpack phenotype and REGENIE_STEP_1 data
            def (prefix, phenoFile, covFile, header, bedFile, bimFile, famFile) = pheno_data
            def (phenotype, pred_file) = pred_file_data

            // Generate 23 entries for each chromosome
            (1..23).collect { chrom ->
                def chrom_with_x = chrom == 23 ? 'X' : chrom
                def bgen_file = file("${params.imputed_bgen_chrs_path}/chr${chrom_with_x}_imputed_s2m.bgen")
                def sample_file = file("${params.bgen_sample_file}")

                // Return a tuple with all required information
                [
                    phenotype, phenoFile, covFile, header, bedFile, bimFile, famFile,
                    chrom_with_x, bgen_file, sample_file, pred_file
                ]
            }
        }

    // // Step 5: Run REGENIE Step 2
    // REGENIE_STEP_2 (
    //     pheno_chrom_pred
    // )
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
        params.phenos
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/