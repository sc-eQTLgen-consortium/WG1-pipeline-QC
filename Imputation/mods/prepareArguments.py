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