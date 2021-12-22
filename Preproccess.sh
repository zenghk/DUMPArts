#Split the large raw reads into small files
for i in {1..2};do cd temp_ku${i};zcat *R1*|split -a 3 -d -l 100000  && zcat *2.fq.gz |split -a 4 -d -l 100000 ;cd ../;done

#Find out barcode and UMI in each small files
for j in {1..2};do cd temp_ku${j};for i in {000..999};do python ../find_UMI_primer_barcode.py -f1 x$i -f2 x0$i -d file_x$i 
-o Summary_x$i -Prim conf_primer.seq -O5S 0 -O5E 10 -O3S 0 -O3E 10 -P5S 20 -P5E 30 -P3S 20 -P3E 30;done;cd ../;done

#Merge the pair-end sequences files with pear
for j in {1..2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};
sample=`echo $i|sed "s/Ig.//g"`
pear -j 4 -q 20 -f ${sample}_R1.fastq -r ${sample}_R2.fastq -o merge;cd ../;done;cd ../;done

#Transfer the fastq file into fasta file
for j in {1..2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};seqkit fq2fa merge.assembled.fastq -o merge.assembled.fasta ;cd ../;done;cd ../;done

#Split the fasta files into small fasta file for the Igblast
for j in {1..2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};mkdir -p subfile_fasta; 
python3 split_to_subfiles.py merge.assembled.fasta subfile_fasta;cd ../;done;cd ../;done

#Igblast
for j in {1..2};do cd ku${j};for i in `echo Ig.*`;do cd ${i}/subfile_fasta/;for m in `echo merge*.fasta`;do 
sh IgBLAST4HumanBCR.sh $m;done;cd ../../;done;cd ../;done

$Parse the output of Igblast
for j in {1,2};do cd ku${j};for i in `echo Ig.*`;do cd ${i}/subfile_fasta/;for m in `echo merge*.fasta`;
do "python new_ParseIgBLAST_20200824.py -f ${m} -i IgBLAST.${m}.m7.txt -o seq_${m}.txt";done;cd ../../;done;cd ../;done

#Sum up all result of Igblast
for j in {1,2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};python3 SumupIgblast.py Seqinfo.txt;cd ../;done;cd ../;done

#Merge all information of UMI and barcode into one file
for j in {1,2};do python3 SumupUMIandBarcode.py ${j};done

#Merge barcode information, UMI information, Igblast result into one file
for j in {1,2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};python3 MergeVRandUMI.py Seqinfo.txt UMIinfo.txt SeqUMI.txt;cd ../;done;cd ../;done

#Assemble the clones according to the result of VJ+CDR3nt
for j in {1,2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};python3 AssembleClone.py Seqinfo.txt Clones2.txt;cd ../;done;cd ../;done

#Calculate the K mer in each files
for j in {1,2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};python3 GetKmer.py Seqinfo.txt Kmer.txt;cd ../;done;cd ../;done

#Get the consensus sequences
for j in {1,2};do cd ku${j};for i in `echo Ig.*`;do cd ${i};python3 GetConsensus_V2.py Seqinfo.txt Kmer.txt Consensusfile.txt;cd ../;done;cd ../;done
