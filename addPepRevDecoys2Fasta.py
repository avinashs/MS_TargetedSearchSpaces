#Avinash Shanmugam, April 15th 2014
#Script will read in a fasta file, split the sequence at tryptic
#sites, reverse the sequence between the tryptic sites and concat
#again to create decoys.

import re;
import os;
import sys;


def getPepRevDecoy( seq, trypSiteRe = re.compile('(K|R)(?!P)')):

	##Processing the seq

	#Split seq at tryptic sites and store the parts
	seqParts = trypSiteRe.split(seq);

	#For some reason the splitfunction includes a trailing
	#empty space. So removing that last element
	#del seqParts[-1];

	#Copy into decoySeqParts array
	decoySeqParts = seqParts;

	#The even elements (counting from 0) in decoySeqParts are the
	#sequences we want to reverse while the odd elements are the tryptic
	#site itself (K|R). So we go through and reverse the even elements.

	#One additional refinement: if the last character of the seq to be
	#reversed is a proline, making it the first character of the new 
	#reversed seq, we flip the 1st and 2nd characters of the new 
	#reverseq seq making proline the 2nd character. This prevents P 
	#coming after a K|R and destroying a tryptic site that was found in 
	#the forward sequence.

	decoySeqParts[::2] = [ x[-2:]+x[-3::-1] if x.endswith("P") else x[::-1] for x in decoySeqParts[::2]];

	#Join the parts to form back the decoy sequence
	decoySeq = "".join(decoySeqParts);

	return decoySeq;

###### Main ########

if len(sys.argv) != 2:

	print "USAGE python addPepRevDecoys2Fasta.py <fastaFile>";
	sys.exit();

#Get infile, tempfile and outfile names
fastaFileName = sys.argv[1];

decoyFileName = fastaFileName.replace(".fa","_decoySeq.fa");
outfileName = fastaFileName.replace(".fa","_plusPepREV.fa");

#Open files
fastaFile = open(fastaFileName);

decoyFile = open(decoyFileName,'w');


#Blank variable to store seq
seq = "";
seqNo = 0;

for line in fastaFile:

	#Remove trailing newline
	line = line.rstrip();

	if line.startswith(">"):

		#When we hit second title line, process and write out 
		#the previous title and reversed seq
		if seqNo != 0:
			
			#Append rev_ to title
			decoyTitle = title.replace(">",">rev_",1);

			decoySeq = getPepRevDecoy(seq);

			decoyFile.write(decoyTitle+"\n");
			decoyFile.write(decoySeq+"\n");

		#Read in next title line
		title = line;
		seqNo +=1;
		seq = ""

		print "\rProcessing seq no. "+str(seqNo),

	else:

		seq = seq + line;
######
###When the loop exits, all the seqs have been read in, but the last seq
###hasn't been written out (cos of the way the loop is written). So do that.
######

#Append rev_ to title
decoyTitle = title.replace(">",">rev_",1);

#Join the parts to form back the decoy sequence
decoySeq = getPepRevDecoy(seq);

decoyFile.write(decoyTitle+"\n");
decoyFile.write(decoySeq+"\n");

###Last seq written out

#Close files
fastaFile.close()
decoyFile.close()

#The decoy seqs and the original fwd seqs are in different files
#Create a cat command to append them to each other
catCommand = "cat "+fastaFileName+" "+decoyFileName+" > "+outfileName;

os.system(catCommand);

#Separate decoy Seqs file is no longer needed. So delete the file
os.system("rm "+decoyFileName);

print "\nOutput file with appended decoys written to "+outfileName
