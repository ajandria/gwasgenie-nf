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

    // Step 1: Extract phenotypes and covariates
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

    // Step 2: Generate chromosome-wise BGEN files
    chromosome_bgen_files = Channel
        .from(1..22, 'X') // Chromosomes 1-22 and X
        .map { chrom ->
            def bgen_file = file("${params.imputed_bgen_chrs_path}/chr${chrom}_imputed_s2m.bgen")
            def sample_file = file("${params.bgen_sample_file}")
            return [chrom, bgen_file, sample_file]
        }
        .filter { chrom, bgen_file, sample_file ->
            bgen_file.exists() && sample_file.exists()
        }
        .flatMap { chrom, bgen_file, sample_file ->
            // Ensure the channel emits one tuple per chromosome
            [[chrom, bgen_file, sample_file]]
        }
    chromosome_bgen_files.view()

    // Step 3: Combine phenotypes (pheno_covs) with chromosomes and their files
    pheno_chrom_bgen = pheno_covs
        .mix(chromosome_bgen_files)
        .map { pheno_data, chrom_data ->
            def (prefix, phenoFile, covFile, header, bedFile, bimFile, famFile) = pheno_data
            def (chrom, bgen_file, sample_file) = chrom_data
            return [prefix, phenoFile, covFile, header, bedFile, bimFile, famFile, chrom, bgen_file, sample_file]
        }
    pheno_chrom_bgen.view()

    // // Step 4: Join with REGENIE_STEP_1 output
    // pheno_chrom_pred = pheno_chrom_bgen
    //     .join(REGENIE_STEP_1.out.s1) { pheno_bgen, step1_output ->
    //         def (phenotype, phenoFile, covFile, header, bedFile, bimFile, famFile, chrom, bgen_file, sample_file) = pheno_bgen
    //         def (phenotype_step1, pred_file) = step1_output
    //         if (phenotype == phenotype_step1) {
    //             return [phenotype, phenoFile, covFile, header, bedFile, bimFile, famFile, chrom, bgen_file, sample_file, pred_file]
    //         }
    //         return null
    //     }
    //     .filter { it != null }

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