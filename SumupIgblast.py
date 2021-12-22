import pandas as pd 
import numpy as np
import csv, os, sys
import glob
import multiprocessing as mp



def Readfile(infile):
	df = pd.read_csv(infile, sep = "\t")
	df = df[~df["CDR3nt"].isna()]
	return df

def main(output):
	path = os.getcwd()
	infiles = glob.glob("%s/subfile_fasta/seq*"%path)
	pool = mp.Pool()
	res = pool.map(Readfile, infiles)
	totaldf = pd.concat(res)
	totaldf.drop_duplicates(subset = ["SeqId"], inplace = True)
	totaldf.to_csv(output, sep = "\t", index = False)


if __name__ == "__main__":
	output = sys.argv[1]
	main(output)
