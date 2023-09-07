#!/usr/bin/env python
import os
import pandas as pd


def get_scrnaseq_dirs(config, samples_df):
    # Extract variables from configuration file.
    input_dict = config["inputs"]
    scrnaseq_dir = input_dict["sc_rnaseq_dir"]
    individual_list_dir = input_dict["individual_list_dir"]
    samplesheet_filepath = input_dict["samplesheet_filepath"]

    ### Check that all the directories exist and can be accessed
    if not os.path.exists(scrnaseq_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(scrnaseq_dir))
    if not os.path.exists(individual_list_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(individual_list_dir))
    if not os.path.exists(samplesheet_filepath):
        raise Exception("File {} does not exist or you have not mounted a parent directory for the singularity bucket".format(samplesheet_filepath))

    scrnaseq_libs_paths = []
    for pool in samples_df["Pool"]:
        # Check if the pool directory exists.
        scrnaseq_pool_dir = os.path.join(scrnaseq_dir, pool)
        if not os.path.exists(scrnaseq_pool_dir):
            print("Could not find a scRNA-seq directory for pool '{}' in your pool list. Please check that are spelled correctly and you do not have any additional pool names that are not in '{}'.".format(pool, scrnaseq_dir))
            exit()

        # Get the barcode file.
        barcode_file = os.path.join(scrnaseq_dir, pool, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")

        # Get the bam file.
        bam_file = os.path.join(scrnaseq_dir, pool, "outs", "possorted_genome_bam.bam")
        if not os.path.exists(bam_file):
            print("Could not find a bam file for pool '{}' in your pool list. Please check that 'outs/possorted_genome_bam.bam' exist in your pool.".format(pool))
            exit()

        # Get the count h5 file.
        count_h5_file = os.path.join(scrnaseq_dir, pool, "outs", "filtered_feature_bc_matrix.h5")
        if not os.path.exists(count_h5_file):
            print("Could not find a h5 count file for pool '{}' in your pool list. Please check that 'filtered_feature_bc_matrix.h5' exist in your pool.".format(pool))
            exit()

        # Get the individual file.
        individual_file = os.path.join(individual_list_dir, "{}.txt".format(pool))
        if not os.path.exists(individual_file):
            print("Could not find a individual list file for pool '{}' in your pool list. Please check that '{}.txt' exist in your individual_list_dir.".format(pool, pool))
            exit()

        scrnaseq_libs_paths.append([scrnaseq_pool_dir, barcode_file, bam_file, count_h5_file, individual_file])

    scrnaseq_libs_df = pd.DataFrame(scrnaseq_libs_paths,
                                    index=samples_df["Pool"],
                                    columns=["scRNAseqDirectory", "BarcodeFile", "BamFile", "CountH5File", "IndividualFile"])
    return scrnaseq_libs_df
