igblastn \
			 -germline_db_V IGHV.fasta \
			 -germline_db_J IGHJ.fasta \
			 -germline_db_D IGHD.fasta \
			 -organism human \
			 -domain_system imgt \
			 -query $1 \
			 -auxiliary_data human_gl.aux \
			 -num_threads 8 \
			 -outfmt '7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq score frames qframe sframe positive ppos btop staxids stitle sstrand qcovs qcovhsp' \
			 -ig_seqtype Ig \
			 -num_alignments_V 5 \
			 -num_alignments_D 5 \
			 -num_alignments_J 5 \
			 -out IgBLAST.$1.m7.txt
