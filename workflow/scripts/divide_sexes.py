#!/usr/bin/env python3

import glob
import fileinput
import math
import os
import pandas
import pathlib
import shutil
import sys


chry=int(1)
chrx=int(1)
ratio=str()
sex=float()

workdir=os.path.abspath(os.getcwd())

with open(workdir+"/qc/samtools_idxstats/idxstats_ratio.txt", "a") as f:
#	f.write("Sample\tRatio\tSex\n")
	for line in open("samples.tsv","r"):
		sample = line.split("\t")[0]
		if sample != "sample":
			for line in open(workdir+"/qc/samtools_idxstats/"+sample+"_N.samtools-idxstats.txt"):
				if "chrX" in line:
					chrx = int(line.split("\t")[2])
				if "chrY" in line:
					chry = int(line.split("\t")[2])
			ratio = str(chry/(chrx+chry))
			sex = float(chry/(chrx+chry))
			if sex < 0.03:
				f.write(sample+"\t"+ratio+"\tF\n")
				shutil.copy("conifer/RPKM/"+sample+".rpkm","/beegfs-storage/projects/wp3/nobackup/Workspace/CoNIFER/hg38/RPKM_F/"+sample+".rpkm")
			elif sex > 0.09:
				f.write(sample+"\t"+ratio+"\tM\n")
				shutil.copy("conifer/RPKM/"+sample+".rpkm","/beegfs-storage/projects/wp3/nobackup/Workspace/CoNIFER/hg38/RPKM_M/"+sample+".rpkm")
			else:
				f.write(sample+"\t"+ratio+"\tUnsolved\n")
				shutil.copy("conifer/RPKM/"+sample+".rpkm","/beegfs-storage/projects/wp3/nobackup/Workspace/CoNIFER/hg38/RPKM_U/"+sample+".rpkm")
				with open(workdir+"/qc/samtools_idxstats/unresolved_sex.txt", "a") as u:
					u.write(sample+"\t"+ratio+"\tCan't resolve\n")

df=pandas.read_csv('samples.tsv', sep='\t')
dr=pandas.read_csv(workdir+'/qc/samtools_idxstats/idxstats_ratio.txt', sep='\t', names=['sample','ratio','sex'])
last_column = dr.iloc[: , -1:]
df['sex']=last_column
df.to_csv('samples_sv.tsv', sep='\t', index=False)

df.to_csv(workdir+'/conifer/RPKM/Done.txt', sep='\t', index=False)
