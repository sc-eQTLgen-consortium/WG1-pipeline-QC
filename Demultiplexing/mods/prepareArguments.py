#!/usr/bin/env python
import os
import pandas as pd
from glob import glob
import re
import subprocess

def matchFolders(x, scrnaseq_dir, dir_list = None):
    for folder in dir_list:
        if os.path.isdir(os.path.join(scrnaseq_dir,folder)):
            if re.search(r'^' + x + "$", folder):
                return(folder)
            elif re.search(r'^' + x + '\D', folder):
                return(folder)
            elif re.search(x + "$", folder):
                return(folder)
            elif re.search(x + "\D", folder):
                return(folder)

def get_barcodes_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if re.search("barcodes.tsv", f)]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath, filename))

def get_bam_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if f.endswith(".bam")]:
            return(os.path.join(dirpath, filename))

def get_matrix_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if re.search("matrix.mtx", f)]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath, filename))
def get_matrix_dirs(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if re.search("matrix.mtx", f)]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath))

def get_individual_files(x, individual_dir = None):
    for filename in individual_dir:
        if re.search(r'^' + x + "$", filename):
            return(filename)
        elif re.search(r'^' + x + '\D', filename):
            return(filename)
        elif re.search(x + "$", filename):
            return(filename)
        elif re.search(x + '\D', filename):
            return(filename)


def get_scrnaseq_dirs(config):
    # Extract variables from configuration file for use within the rest of the pipeline
    input_dict = config["inputs"]
    output_dict = config["outputs"]
    ref_dict = config["refs"]

    # General variables used by the rest of the pipeline
    scrnaseq_dir = input_dict["scRNAseq_dir"]
    individual_list_dir = input_dict["individual_list_dir"]

    ### Check that all the directories exist and can be accessed
    if not os.path.exists(scrnaseq_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(scrnaseq_dir))
    if not os.path.exists(individual_list_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(individual_list_dir))

    # Read in samplesheet from the configfile specified from the command line
    samples = pd.read_csv(input_dict["samplesheet_filepath"], sep = "\t")

    # Expect first colunn to be pools
    pools = samples.iloc[:, 0]

    # Match pools to scrna seq directories to make a list of each scRNA-seq dir
    scrna_seq_dirlist = os.listdir(scrnaseq_dir)
    try:
        scrnaseq_filelist = [os.path.join(scrnaseq_dir, matchFolders(pool, dir_list = scrna_seq_dirlist, scrnaseq_dir = scrnaseq_dir)) for pool in pools]
    except TypeError:
        print("Could not find a scRNA-seq directory for all of the pools in your pool list. Please check that they are spelled correctly and you do not have any additional pool names that are not in {}  ".format(individual_list_dir))
    scrnaseq_filedict = dict(zip(pools, scrnaseq_filelist))
    scrnaseq_libs = pd.Series(scrnaseq_filedict, name="scRNAseq_Directories")

    ### Get the barcode files for each pool
    try:
        barcode_filelist = [get_barcodes_files(pool) for pool in scrnaseq_filelist]
    except:
        print("Could not find a barcode file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'barcodes.tsv' within the name.")
    barcode_filedict = dict(zip(pools, barcode_filelist))
    barcode_libs = pd.Series(barcode_filedict, name="Barcode_Files")

    ### Get the bam files for each pool
    try:
        bam_filelist = [get_bam_files(pool) for pool in scrnaseq_filelist]
    except:
        print("Could not find a bam file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain '.bam' within the name.")
    bam_filedict = dict(zip(pools, bam_filelist))
    bamlibs = pd.Series(bam_filedict, name="Bam_Files")

    ### Get the matrix files for each pool
    try:
        matrix_filelist = [get_matrix_files(pool) for pool in scrnaseq_filelist]
    except:
        print("Could not find a matrix file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'matrix.mtx' within the name.")
    matrix_filedict = dict(zip(pools, matrix_filelist))
    matrix_libs = pd.Series(matrix_filedict, name="Matrix_Files")

    ### Get the directories of all the matrix EXCLUDING the matrix filenames themselves
    try:
        matrix_dirlist = [get_matrix_dirs(pool) for pool in scrnaseq_filelist]
    except:
        print("Could not find a matrix file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'matrix.mtx' within the name.")
    matrix_dirdict = dict(zip(pools, matrix_dirlist))
    matrix_dir_libs = pd.Series(matrix_dirdict, name="Matrix_Directories")

    ### Get the matrix individual list file for each pool
    individual_dirlist = os.listdir(individual_list_dir)
    try:
        individual_filelist = [os.path.join(individual_list_dir, get_individual_files(pool, individual_dir = individual_dirlist)) for pool in pools]
    except:
        print("Could not find a files of individuals in {}. Please check that they exist somewhere in this directory and contain the pool names within the name of the file.".format(individual_list_dir))
    individual_filedict = dict(zip(pools, individual_filelist))
    individual_libs = pd.Series(individual_filedict, name="Individuals_Files")

    dataframe=pd.concat([scrnaseq_libs, barcode_libs, bamlibs, matrix_libs, matrix_dir_libs, individual_libs], axis=1)
    return(dataframe)

