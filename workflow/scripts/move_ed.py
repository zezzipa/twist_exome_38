#!/usr/bin/env python3

import glob
import fileinput
import math
import os
import pandas
import pathlib
import shutil
import sys


workdir=os.path.abspath(os.getcwd())

bamfile=snakemake.input[0]

print(bamfile)

for line in open("samples_sv.tsv","r"):
	sample = line.split("\t")[0]
	file = bamfile.split("/")[2]
	name = file.split("_")[0]
	if sample == name:
		if "F" in line:
			shutil.copy("parabricks/pbrun_fq2bam/"+sample+"_N.bam","cnv_sv/exomedepth_F/"+sample+"_N.bam")
			shutil.copy("parabricks/pbrun_fq2bam/"+sample+"_N.bam.bai","cnv_sv/exomedepth_F/"+sample+"_N.bam.bai")
		if "M" in line:
			shutil.copy("parabricks/pbrun_fq2bam/"+sample+"_N.bam","cnv_sv/exomedepth_M/"+sample+"_N.bam")
			shutil.copy("parabricks/pbrun_fq2bam/"+sample+"_N.bam.bai","cnv_sv/exomedepth_M/"+sample+"_N.bam.bai")
