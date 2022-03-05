#!/usr/bin/env python
import os
import pandas as pd
from glob import glob
import re
import subprocess

def parsePaths(config):
    # Strip trailing slashes from paths 
    # scrnaseq_dir, individual_list_dir, output_dir
    config["outputs"]["output_dir"] = (config["outputs"]["output_dir"]).rstrip("/")
    return(config) 


def getPGEN(plink_dir, basename):
	return(plink_dir + "/" + basename + ".pgen")


def getPVAR(plink_dir, basename):
	return(plink_dir + "/" + basename + ".pvar")

def getPSAM(plink_dir, basename):
	return(plink_dir + "/" + basename + ".psam")


def getVCFdir(ref_dir):
	for dirpath, dirnames, filenames in os.walk(ref_dir):
		for filename in [f for f in filenames if re.search("30x-GRCh38_NoSamplesSorted.vcf.gz.tbi", f)]:
			return(dirpath)

def getFASTA(ref_dir):
	for dirpath, dirnames, filenames in os.walk(ref_dir):
		for filename in [f for f in filenames if re.search("Homo_sapiens.GRCh38.dna.primary_assembly.fa$", f)]:
				return(os.path.join(dirpath, filename))

def getMAP(ref_dir):
	for dirpath, dirnames, filenames in os.walk(ref_dir):
		for filename in [f for f in filenames if re.search("genetic_map_hg38_withX.txt.gz", f)]:
			return(os.path.join(dirpath, filename))

def getPHASINGdir(ref_dir):
	for dirpath, dirnames, filenames in os.walk(ref_dir):
		for filename in [f for f in filenames if re.search("chr10.bcf", f)]:
			return(dirpath)

def getIMPUTATIONdir(ref_dir):
	for dirpath, dirnames, filenames in os.walk(ref_dir):
		for filename in [f for f in filenames if re.search("chr10.m3vcf.gz", f)]:
			return(dirpath)