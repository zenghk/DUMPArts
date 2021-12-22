import pandas as pd
import numpy as np
import csv, os, sys

def Assemble(infile, clonefile):
	df = pd.read_csv(infile, sep = "\t", usecols = ["VHits", "JHits", "CDR3nt", "SHM", "Vinfo"])
	df["Vgene"] = df["VHits"].str.split("\*", expand = True)[0]
	df["Jgene"] = df["JHits"].str.split("\*", expand = True)[0]
	df = df[~df["CDR3nt"].isnull()]
	df["Vlength"] = df["Vinfo"].str.split(";", expand = True)[1].astype(int)
	df = df[df["Vlength"] >= 200]
	groups = df.groupby(["Vgene", "Jgene", "CDR3nt"])
	result = groups.size().to_frame().reset_index()
	SHM = groups["SHM"].apply(np.mean).to_frame().reset_index()
	SHM.columns = ["Vgene", "Jgene", "CDR3nt", "AverageSHM"]
	result.columns = ["Vgene", "Jgene", "CDR3nt", "Size"]
	result["Fre"] = result["Size"]/result["Size"].sum()
	final = pd.merge(result, SHM, on = ["Vgene", "Jgene", "CDR3nt"])
	final.sort_values(by = "Fre", ascending  = False, inplace = True)
	final.to_csv(clonefile, sep = "\t", index = False)

if __name__ == "__main__":
	infile, clonefile = sys.argv[1], sys.argv[2]
	Assemble(infile, clonefile)
