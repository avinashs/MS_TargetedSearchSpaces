#Avinash SHanmugam, Oct 27th 2014
#Script will read in a freq.tsv file of tryptic peptides
#and write an outfile only containing those peptides having
#number of missed cleavages less than or equal to the given
#threshold.

import re;

infileName = "gpm.label.humanpeptides.istryp.freq.tsv";

nmissTh = 2;
outfileName = infileName.replace("istryp","istryp.nmiss2")

infile = open(infileName);
outfile = open(outfileName, "w");

trypRe = re.compile('(R|K)(?!P|$)');

for line in infile:

	line = line.rstrip();

	#Split line into fields
	(pep, freq) = line.split("\t");

	#Get number of missed cleavages
	nmiss = len(trypRe.findall(pep));

	#Write out line of nmiss less than threshold
	if nmiss <= nmissTh:

		outfile.write(line+"\n");
	

infile.close();
outfile.close();
