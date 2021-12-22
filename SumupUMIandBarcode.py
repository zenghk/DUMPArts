import pandas as pd
import numpy as np
import csv, os, sys
import glob
import multiprocessing as mp
from collections import Counter


def ReadUMI(infile):
	titlelist = ['SeqId', 'Sample_id', 'Out5_mis', 'Out5_id', 'Out5_seq', 'Out5_pos',
       'Out3_mis', 'Out3_id', 'Out3_seq', 'Out3_pos', 'Pri5_mis', 'Pri5_id',
       'Pri5_seq', 'Pri5_pos', 'Pri3_mis', 'Pri3_id', 'Pri3_seq', 'Pri3_pos',
       'UMI5', 'UMI3']
	typelist = ["category"] * len(titlelist)
	typedic = dict(zip(titlelist, typelist))
	df = pd.read_csv(infile, sep = "\t", dtype = typedic)
	df = df[df["Sample_id"] != "NO"]
	df = df[~((df["UMI5"].str.contains("N")) & (df["UMI3"].str.contains("N")))]
	return df


def main(kunum):
	infiles = glob.glob("temp_ku%s/file_*/Summary*"%kunum)
	pool = mp.Pool()
	res = pool.map(ReadUMI, infiles)
	totaldf = pd.concat(res)
	totaldf.drop_duplicates(subset = ["SeqId"], inplace = True)
	groups = totaldf.groupby("Sample_id")
	for sample, group in groups:
		group["UMIpair"] = group["UMI5"] + "|" + group["UMI3"]
		Countdic = Counter(group["UMIpair"])
		group["UMInum"] = list(map(lambda x:Countdic[x], group["UMIpair"]))
		del group["UMIpair"]
		group.to_csv("ku%s/Ig.%s/UMIinfo.txt"%(kunum, sample), sep = "\t", index = False)

if __name__ == "__main__":
	kunum = sys.argv[1]
	main(kunum)
