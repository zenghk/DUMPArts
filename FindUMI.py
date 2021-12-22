import csv, os, sys
import re
from Bio import SeqIO
import Levenshtein
import glob, argparse
from itertools import chain

'''
Usage: python find_UMI_primer_barcode.py -f1 sub_R1_6.fastq -f2 sub_R2_6.fastq -d summary -o summary_6 -O5S 0 -O5E 10 -O3S 0 -O3E 10 -P5S 25 -P5E 35 -P3S 25 -P3E 35 -B5S 40 -B5E 90 -B3S 39 -B3E 60
'''

def find_UMI(seq, start_5, start_3): #Find the UMI according to the position of outer primer and inner primer
	UMI = seq[start_5: start_3]
	return UMI

#def find_barcode2(barcode_seq, sequence, start, end):
		
def find_barcode(barcode_dict, sequence, start, end): #Find barcode/primer according to the identity.
	label, Tbar, mismatch, best_pos = "-", "NNNNNNNNNN", 100, "-"
	for barcode in barcode_dict.keys():
		barcode_len = len(barcode)
		candidate_bar, candidate_mis, candidate_pos = "NNNNNNNNNN", 100, "-"
		for i in range(start, end):
			misnum = Levenshtein.distance(barcode, sequence[i:i+barcode_len])
			if misnum < candidate_mis:
				candidate_bar, candidate_mis, candidate_pos = sequence[i: i + barcode_len], misnum, i
			else:
				continue
		if candidate_mis <= 5:
			ide1, seq1, pos1 = candidate_mis, candidate_bar, candidate_pos
		else:
			ide1, seq1, pos1 = 100, "NNNNNNNNNN", "-"
		if ide1 < mismatch:
			label, Tbar, mismatch, best_pos = barcode_dict[barcode], seq1, ide1, pos1
		else:
			continue
	return [mismatch, label, Tbar, best_pos]

def get_primer(file): #Get the primer seq as dict, key: primer sequence, values: primer id
	primer_dict = {}
	for line in csv.reader(open(file, 'r'), delimiter = "\t"):
		if line[0].startswith("#"):
			continue
		else:
			primer_dict.setdefault(line[1], line[0])	
	return primer_dict

def get_barcode(file):  #Get the barcode seq as dict, key: barcode sequence, values: barcode id
	reader = csv.reader(open(file, 'r'), delimiter = "\t")
	dict_5, dict_3 = {}, {}
	for line in reader:
		if line[0].startswith("#"):
			continue
		elif line[0].startswith("P5"):
			dict_5.setdefault(line[1], line[0])
		elif line[0].startswith("P3"):
			dict_3.setdefault(line[1], line[0])
	return dict_5, dict_3

def judge_chimera(bar1, bar2): #Judge whether the sequence is chimera or true sequence by the combination of the barcode
	if bar1[1] == bar2[1] == "5":
		label = "5-5"
	elif bar1[1] == bar2[1] == "3":
		label = "3-3"
	else:
		if bar1[1] == "5" and bar2[1] == "3":
			if bar1[2:] == bar2[2:]:
				label = "true_seq"
			else:
				label = "5-3"
	return label

def remove_N(result): #Remove the possible N in the beginning of the sequence
	N_num = result[2].count("N")
	if N_num > 0:
		result[0] = result[0] - N_num
	return result

def get_sample(sample_id):
	sample_dict = {}
	for line in csv.reader(open(sample_id, 'r'), delimiter = "\t"):
		if line[0].startswith("#"):
			continue
		else:
			sample_dict.setdefault("_".join(line[:2]), line[2])
	return sample_dict
	
def main(primer_path):
	outer_barcode = os.path.join(os.getcwd(), 'conf_barcode.seq')
	conf_sample = os.path.join(os.getcwd(), 'conf_sample.txt')
	out_dict = get_primer(outer_barcode)
	primer_dict = get_primer(primer_path)
	sample_dict = get_sample(conf_sample)
	fq1, fq2 = SeqIO.parse(fastq1, 'fastq'), SeqIO.parse(fastq2, 'fastq')
	os.system("rm %s/*_R1.fastq"%outdir)
	os.system("rm %s/*_R2.fastq"%outdir)
	out = csv.writer(open("%s/%s"%(outdir, outfile), 'w'), delimiter = "\t")
	out.writerow(["SeqId","Sample_id", "Out5_mis", "Out5_id", "Out5_seq", "Out5_pos", "Out3_mis", "Out3_id", "Out3_seq", "Out3_pos", "Pri5_mis", "Pri5_id", "Pri5_seq", "Pri5_pos", "Pri3_mis", "Pri3_id", "Pri3_seq", "Pri3_pos", "UMI5", "UMI3"]) 
	try:
		while True:
			rec1, rec2 = fq1.next(), fq2.next()
			#print out_dict, str(rec1.seq), str(rec2.seq)
			out_5 = remove_N(find_barcode(out_dict, str(rec1.seq), start_out_5, end_out_5))
			out_3 = remove_N(find_barcode(out_dict, str(rec2.seq), start_out_3, end_out_3))
			if "P5" in out_5[1]: 
				Rec1, Rec2 = rec1, rec2
			else:
				Rec1, Rec2 = rec2, rec1
				out_5, out_3 = out_3, out_5
			#print out_5, out_3
			primer_5 = find_barcode(primer_dict, str(Rec1.seq), start_primer_5, end_primer_5)
			primer_3 = find_barcode(primer_dict, str(Rec2.seq), start_primer_3, end_primer_3)
			if out_5[0] <= 1 and out_3[0] <= 1:
				bar_pair = "_".join([out_5[1], out_3[1]])
				sample_id = sample_dict.get(bar_pair, "NO")
			else:
				sample_id = "NO"
			if out_5[0] <= 1 and primer_5[0] <= 1 and sample_id.split("-")[-1] != "No-UMI":
				UMI5 = find_UMI(str(Rec1.seq), len(out_5[2]) + out_5[3], primer_5[3])
			else:
				UMI5 = "NNNNNNNNNN"
			if out_3[0] <= 1 and primer_3[0] <= 1 and sample_id.split("-")[-1] != "No-UMI":
				UMI3 = find_UMI(str(Rec2.seq), len(out_3[2]) + out_3[3], primer_3[3])
			else:
				UMI3 = "NNNNNNNNNN"
			out.writerow([rec1.id] + [sample_id] + out_5 + out_3 + primer_5 + primer_3 + [UMI5] + [UMI3])
			if out_5[0] <= 1 and out_3[0] <= 1:
				out1 = open("%s/%s_R1.fastq"%(outdir, sample_id), 'a')
				out2 = open("%s/%s_R2.fastq"%(outdir, sample_id), 'a')
				SeqIO.write(Rec1, out1, 'fastq')
				SeqIO.write(Rec2, out2, 'fastq')
			else:
				out1 = open("%s/No_R1.fastq"%(outdir), 'a')
				out2 = open("%s/No_R2.fastq"%(outdir), 'a')
				SeqIO.write(Rec1, out1, 'fastq')
				SeqIO.write(Rec2, out2, 'fastq')
			#print(rec1.id, sample_id, primer_5, primer_3, out_5, out_3, UMI5, UMI3)
	except StopIteration:
		print("Done")
			
if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='python find_UMI_primer_barcode.py')
	parser.add_argument("-f1", '--fastq1', help = "Input fastq1")
	parser.add_argument("-f2", '--fastq2', help = "Input fastq2")
	parser.add_argument("-d", '--outdir', help = "Output directory")
	parser.add_argument("-o", '--outfile', help = "Output file, tab-split")
	parser.add_argument("-Prim", '--primer_path', help = "Input primer file")
	parser.add_argument('-O5S', '--start_out_5', type = int,  help = "The start point of outer 5 barcode")
	parser.add_argument('-O5E', '--end_out_5', type = int, help = "The end point of outer 5 barcode")
	parser.add_argument('-O3S', '--start_out_3', type = int, help = "The start point of outer 3 barcode")
	parser.add_argument('-O3E', '--end_out_3', type = int, help = "The end point of outer 3 barcode")
	parser.add_argument('-P5S', '--start_primer_5', type = int, help = "The start point of primer 5")
	parser.add_argument('-P5E', '--end_primer_5', type = int, help = "The end point of priemr 5")
	parser.add_argument('-P3S', '--start_primer_3', type = int, help = "The start point of primer 3")
	parser.add_argument('-P3E', '--end_primer_3', type = int, help = "The end point of primer 3")
	args = parser.parse_args()
	fastq1, fastq2 = args.fastq1, args.fastq2
	outdir, outfile = args.outdir, args.outfile
	start_out_5, end_out_5 = args.start_out_5, args.end_out_5
	start_out_3, end_out_3 = args.start_out_3, args.end_out_3
	start_primer_5, end_primer_5 = args.start_primer_5, args.end_primer_5
	start_primer_3, end_primer_3 = args.start_primer_3, args.end_primer_3
	primer_path = args.primer_path
	os.system("mkdir -p %s"%outdir)
	main(primer_path)	
