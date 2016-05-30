#Avinash Shanmugam, April 16th 2014
#Script will take a list of peptides and check if they are tryptic peptides
#by grepping them against the fasta file

import os;
import re;

pepFileName = "test.freq.tsv";

fastaFileName = "/tpp.data/avinashs/gpmdb/aug102014/fasta/Homo_sapiens.GRCh38.76.pep.all.singleLineSeq.fa";

trypOutName = "gpm.hpeptides.expect0.001filt.istryp.tsv";
nontrypOutName = "gpm.hpeptides.expect0.001filt.nontryp.tsv";

#Open and read in the peptides file
pepFile = open(pepFileName);
trypOut = open(trypOutName,"w");
nontrypOut = open(nontrypOutName,"w");

#Read and discard title line
pepFile.readline();

#Write title line for outfile
trypOut.write("peptide\tfreq\tnmatches\tistryp\n");
nontrypOut.write("peptide\tfreq\tnmatches\tistryp\n");

#Counter
pepNo = 0;

for line in pepFile:

	#Get rid of trailing newlines
	line = line.rstrip();

	#Split line into peptide and freq fields
	fields = line.split();

	pep = fields[0];
	freq = fields[1];

	pepNo +=1;

	if pep.startswith("P"):

		#If peptide begins with Proline, it can only be tryptic if it
		#occurs at the beginning of the protein sequence. So only look
		#for such occurances in grep
		grepCommand = "grep '^"+pep+"' "+fastaFileName;

	elif not pep.endswith(("K","R")):
		#If peptide doesn't end K or R, it can only be tryptic if it
		#occurs at the end of the protein sequence. So only look for
		#such occurances in grep
		grepCommand = "grep '"+pep+"$' "+fastaFileName;

	else:
		#If neither of the above, grep for occurances of pep anywhere
		#in the protein sequence
		grepCommand = "grep '"+pep+"' "+fastaFileName;

	#Use grep command to find occurances of pep in fastaFile
	grepReturn = os.popen(grepCommand).read();

	#Find number of occurances (count number of newlines in return line)
	nMatches = grepReturn.count("\n");

	#Find the A.A. before and after the peptide in each match
	flankingAA = re.findall("(.)?"+pep+"(.)",grepReturn,re.DOTALL);

	#Extract prev and next AAfrom flanking AA list of tuples
	#Element 0 in each tuple is prevAA and element 1 is nextAA
	prevAA = [x[0] for x in flankingAA];
	nextAA = [x[1] for x in flankingAA];

	#Find whether any of the matches have a trytic NTerminal
	if( sum( aa == "K" or aa == "R" or aa == "\n" or aa == '' for aa in prevAA) > 0):
		ntermTryptic =1;

	else:
		ntermTryptic =0;
	
	#Find whether any of the matches have a tryptic CTerminal
	#(Since we already enforce that peptides not ending in K/R must occur
	#at the end of seqs, here we only have to check that nextAA is not a proline)
	if( sum( aa != "P" for aa in nextAA) > 0):

		ctermTryptic =1;

	else:
		ctermTryptic =0;
	

	#Based on the nterminal and cterminal characteristics, decide of the peptide
	#is tryptic or not and write to appropriate outfile
	if ntermTryptic ==1 & ctermTryptic ==1:
		isTryp =1;
		trypOut.write(pep+"\t"+freq+"\t"+str(nMatches)+"\t"+str(isTryp)+"\n");

	else:
		isTryp =0;
		nontrypOut.write(pep+"\t"+freq+"\t"+str(nMatches)+"\t"+str(isTryp)+"\n");

	print "\rFinished searching for peptide no. "+str(pepNo),;

print "\n\nOutfiles written to "+trypOutName+" & "+nontrypOutName;

trypOut.close();
nontrypOut.close();
pepFile.close();



