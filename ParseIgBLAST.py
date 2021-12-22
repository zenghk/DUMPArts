import os,sys,csv
from Bio import SeqIO
import argparse
import re

'''Adding the step to cut the beginning of variable region so that the same primer amplified by differetn primers will be the same'''

def main():
	fasta_dict = SeqIO.index(infasta,'fasta')
	inf1 = open(infile,'rU')
	out = csv.writer(open('%s/%s'%(outdir,outname),'wb'),delimiter='\t')
	out.writerow(["CDR3nt", "CDR3aa", "SeqId", "VHits", "DHits", "JHits", "ChainType", "StopCodon", "Frame", "Productive", "Strand", "SHM", "Vinfo", "Dinfo", "Jinfo", "Variable"])
	igblast_dict = {}
	myid = ""
	for rec in inf1:
		if rec.startswith("# Query:"):
			if len(myid)!=0:
				if len(Vinfo)>0 and len(Jinfo)>0:
					variable_seq = str(myseq1.seq[int(Vinfo[8])-1:int(Jinfo[9])])
					if len(cdr3)>0:
						start_nt=myseq1.seq[int(cdr3[3])-4:int(cdr3[3])-1]
						end_nt=myseq1.seq[int(cdr3[4]):int(cdr3[4])+3]
						whole_nt=start_nt+cdr3[1]+end_nt
						whole_aa=start_nt.translate()+cdr3[2]+end_nt.translate()
						out.writerow([whole_nt,whole_aa,myid]+recom+[mutnum, ";".join(Vinfo[3:12]),";".join(Dinfo[3:12]),";".join(Jinfo[3:12]),variable_seq])
					else:
						out.writerow(["N/A","N/A",myid]+recom+[mutnum, ";".join(Vinfo[3:12]),";".join(Dinfo[3:12]),";".join(Jinfo[3:12]),variable_seq])
			myid,recom,cdr3,Vinfo,Dinfo,Jinfo = "",[],[],[],[],[]
			myid = rec.strip().split(" ")[2]
			vnum, dnum, jnum = 0,0,0
		elif rec.startswith("IG"):
			recominfo = rec.strip().split("\t")
			if len(recominfo)==7:
				recom = [rec.strip().split("\t")[0],"N/A"]+rec.strip().split("\t")[1:]
			elif len(recominfo)==8:
				recom = rec.strip().split("\t")
		elif rec.startswith("CDR3\t"):
			cdr3 = rec.strip().split("\t")
		elif rec.startswith("V\t") and vnum==0:
			Vinfo = rec.strip().split("\t")
			if Vinfo[24]:
				SHMinfo = Vinfo[24]
			else:
				SHMinfo = ""
			mutinfo = re.findall('[0-9]{1,3}[ACGT][AGTC]', SHMinfo)
			mutnum = len(mutinfo)
			if Vinfo[1].startswith("reversed"):
				myseq1 = fasta_dict[myid].reverse_complement()
			else:
				myseq1 = fasta_dict[myid]
			Subjectstart = int(Vinfo[10])
			Vnewmax = max(Subjectstart, 25)
			if Vnewmax == 25:
				Vnewstart = int(Vinfo[8]) + (25 - Subjectstart)
				Vinfo[8] = str(Vnewstart)
			vnum = 1
		elif rec.startswith("D\t") and dnum==0:
			Dinfo = rec.strip().split("\t")
			dnum = 1
		elif rec.startswith("J\t") and jnum==0:
			Jinfo = rec.strip().split("\t")
			Subjectend = int(Jinfo[11])
			Subjectlen = int(Jinfo[15])
			Jnewmin = min(Subjectend, Subjectlen - 8)
			if Jnewmin != Subjectend:
				Jnewend = int(Jinfo[9]) - (Subjectend - (Subjectlen - 8))
				Jinfo[9] = str(Jnewend)
			jnum = 1
	if len(Vinfo)>0 and len(Jinfo)>0:
		variable_seq = str(myseq1.seq[int(Vinfo[8])-1:int(Jinfo[9])])
		if len(cdr3)>0:
			start_nt=myseq1.seq[int(cdr3[3])-4:int(cdr3[3])-1]
			end_nt=myseq1.seq[int(cdr3[4]):int(cdr3[4])+3]
			whole_nt=start_nt+cdr3[1]+end_nt
			whole_aa=start_nt.translate()+cdr3[2]+end_nt.translate()
			out.writerow([whole_nt,whole_aa,myid]+recom+[mutnum, ";".join(Vinfo[3:12]),";".join(Dinfo[3:12]),";".join(Jinfo[3:12]),variable_seq])
		else:
			out.writerow(["N/A","N/A",myid]+recom+[mutnum, ";".join(Vinfo[3:12]),";".join(Dinfo[3:12]),";".join(Jinfo[3:12]),variable_seq])
	
if __name__=='__main__':
	parser = argparse.ArgumentParser(prog='python ParseIgBLAST.py',usage='%(prog)s -f fasta -i igblast -o outname -d outdir')
	parser.add_argument('-i','--igblast',help='Input igblast m7 format result')
	parser.add_argument('-d','--outdir',default=".",help='Output file directory')
	parser.add_argument('-o','--outfile',default="result.txt",help='Output tab format file, which record search results.')
	args = parser.parse_args()
	infasta = args.fasta
	infile = args.igblast
	outdir = args.outdir
	outname = args.outfile
	os.system("mkdir -p %s"%outdir)
	main()
