import pandas as pd
import numpy as np
import csv, os, sys

def Readseq(Seqfile):
	df = pd.read_csv(Seqfile, sep = "\t", usecols = ["CDR3nt", "SeqId", "Variable", "VHits", "JHits", "SHM", "Vinfo"])
	df["Vgene"] = df["VHits"].str.split("*", expand = True)[0]	
	df["Jgene"] = df["JHits"].str.split("*", expand = True)[0]
	df["Vlength"] = df["Vinfo"].str.split(";", expand = True)[1].astype(int)
	df = df[df["Vlength"] >= 200]
	return df[["CDR3nt", "SeqId", "Variable", "Vgene", "Jgene", "SHM"]]

def ReadUMI(UMIfile):
	df = pd.read_csv(UMIfile, sep = "\t", usecols = ["SeqId", "UMI5", "UMI3", "UMInum"])
	df = df[~(df["UMI5"].str.contains("N")) & ~(df["UMI3"].str.contains("N"))]
	return df

def main(seqfile, UMIfile, output):
	seqdf, UMIdf = Readseq(seqfile), ReadUMI(UMIfile)
	final = pd.merge(seqdf, UMIdf, on = "SeqId")
	final.drop_duplicates(subset = ["SeqId"], inplace = True)
	final.to_csv(output, sep = "\t", index = False)

if __name__ == "__main__":
	seqfile, UMIfile, output = sys.argv[1], sys.argv[2], sys.argv[3]
	main(seqfile, UMIfile, output)
