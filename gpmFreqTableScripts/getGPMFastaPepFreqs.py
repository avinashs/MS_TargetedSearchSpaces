#Avinash Shanmugam, Aug 6th 2014
#Script will take a list of peptides (from GPMDB) and check if they are present
#in a given fasta file. It will also generate a freq in GPMDB table for these peptides

import os;
import re;

#Fn will read in a fasta file and remove the linebreaks
#in the sequence parts and output a single line fasta.
#Is helpful when quickly looking for seq matches in the file
#using grep. Otherwise seqs broken up by newlines wont be recognized

def singleLineSeqFasta(fastaFileName){

	fastaFileName = "Homo_sapiens.GRCh37.75.pep.all.fa";
	 
	outfileName = fastaFileName.replace(".fa",".singleLineSeq.fa");

	fastaFile = open(fastaFileName);
	outfile = open(outfileName,"w");

	seqNo =0;

	for line in fastaFile:

		line = line.rstrip();

		if line.startswith(">"):

			if seqNo !=0:

				#After the first title, each time we hit a title
				#line, it mean we have finished concatenating seq
				#So writeprev title and seq to outfile
				outfile.write(title+"\n");
				outfile.write(seq+"\n");
				
				print "\rWritten out seqNo: "+str(seqNo),;
				
			title = line;
			seq = "";
			seqNo +=1;
		
		#If not title line, concat to seq
		else:

			seq = seq+line;

	#When loop exits, last seq is still not written
	#So do that.
	outfile.write(title+"\n");
	outfile.write(seq+"\n");

	outfile.close()
	fastaFile.close()

	return outfileName;
}


################# Main ####################

pepFileName = "gpm.peptides.uniq.6.tsv";

fastaFileName = "fasta/Homo_sapiens.GRCh38.76.pep.all.fa";

outfileName = "gpm.peptides.FAmatch.freq.tsv";

#Convert fasta file to a singleLineSeq file. Use this new file
#as the fasta source file henceforth
fastaFileName = singleLineSeqFasta(fastaFileName);

#Open and read in the peptides file
pepFile = open(pepFileName);
outfile = open(outfileName,"w");

#Read and discard title line
pepFile.readline();

#Counter
pepNo = 0;

#Define dictionary to store peptide counts
pepFreq = {};

for pep in pepFile:

	#Get rid of trailing newlines
	pep = pep.rstrip();

	pepNo +=1;

#Commented out block from prev version which would also get the protids that matched
#to each peptide. If uncommenting, change outfile format accordingly (peptide\tprotid string)

#	#Assemble grep command
#	grepCommand = "grep -B1 '"+pep+"' "+fastaFileName+" | grep -v '"+pep+"'";
#	
#	#Use grep command to find occurances of pep in fastaFile and return the title Lines
#	grepReturn = os.popen(grepCommand).read();
#
#	grepReturn = grepReturn.rstrip();
#
#	if grepReturn == "":
#
#		#Write blank out to pepids and type,
#		#and move to next peptide
#		print "No matches to pep:"+ pep+"\n";
#		
#		continue;
#
#	#Split the returned titleLines on newlines
#	titleLines = grepReturn.split("\n");
#
#
#	#Blank string to concat pepids
#	pepidString = "";
#
#	for title in titleLines:
#
#		if not title.startswith(">"):
#
#			continue;
#
#		fields = title.split(" ");
#		
#		if len(fields) < 3:
#
#			print pep;
#			print fields;
#
#		id = fields[0];
#
#		#Format the extracted fields
#		id = id.replace(">","");
#
#		pepidString = pepidString+id+";";
#
#	#Remove trailing semicolon
#	pepidString = pepidString[:-1];
#	
#	outstring = "\t".join((pep,pepidString));
#
#	outfile.write(outstring+"\n");

	if pep in pepFreq:
	
		#If pep is in pepFreq dict, it has already been checked and found present
		#in fasta file. So only increment its freq
		pepFreq[pep] = pepFreq[pep] + 1;
	
	else:
		#If not present in dict, it has not been checked, so check using grep

		#Assemble grep command
		grepCommand = "grep -c '"+pep+"' "+fastaFileName;

		#Use grep command to find no. of occurances of pep in fasta file
		grepReturn = os.popen(grepCommand).read();
		grepReturn = grepReturn.rstrip();

		#If peptide has at least 1 match, include in dictionary
		if grepReturn != 0:
		
			pepFreq[pep] = 1;

		else:
			
			print "No match for peptide"+str(pepNo)+":"+pep+"\n";
	
	if pepNo % 1000000 == 0:
		
		print "Completed checking for peptide no%d" % pepNo;


pepFile.close();

#Begin writing freq table to outfile
outfile.write("peptide\tfreq\n");

for pep in pepFreq:

	outline = "\t".join((pep,str(pepFreq[pep])));

	outfile.write(outline+"\n");

outfile.close();

print "Freq table written to"+outfile;

