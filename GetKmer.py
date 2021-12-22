import pandas as pd
import numpy as np
import sys, csv, os

def main(infile, output):
	df = pd.read_csv(infile, sep = "\t", usecols = ["UMInum"])
	result = np.unique(df["UMInum"].values, return_counts = True)
	numberArr = result[1]
	posArr = np.ediff1d(numberArr)
	Threshold = np.where(posArr > 0)[0][0]
	out = csv.writer(open(output, 'w'), delimiter = "\t")
	out.writerow(["Threshold", Threshold])
	

if __name__ == "__main__":
	infile, output = sys.argv[1], sys.argv[2]
	main(infile, output)
