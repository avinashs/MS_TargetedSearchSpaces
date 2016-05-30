#Dec 8th 2014, Avinash Shanmugam
#Script will parse a pepxml file and output as a tsv file

import sys
from pyteomics import pepxml, auxiliary

if len(sys.argv) != 2:

	print "USAGE: python pepxmlparse.py <pepxmlfile>";
	sys.exit();

pepxfile = sys.argv[1];
outfile = pepxfile.replace(".xml","xml.parse.tsv");

#Create pepxml reading iterator
pepxitr = pepxml.read(pepxfile, read_schema=False);

#Open outfile
out = open(outfile,"w");

#Create and write title line to outfile
titleLine = "\t".join(["peptide","modpeptide","specid","mass","charge", "iniprob","pprophprob","isfwd","protidString"]);
out.write(titleLine+"\n");

lineNo = 0;

#Iterate through the pepxml recoprds 

for pepxrec in pepxitr:

	#Extract needed vals from the dict returned
	peptide = pepxrec['search_hit'][0]['peptide'];

	modpeptide = pepxrec['search_hit'][0]['modified_peptide'];

	specid = pepxrec['spectrum'];

	mass = pepxrec['precursor_neutral_mass'];

	charge = pepxrec['assumed_charge'];

	pprophprob = pepxrec['search_hit'][0]['analysis_result'][0]['peptideprophet_result']['probability'];
	
	iniprob = pepxrec['search_hit'][0]['analysis_result'][1]['interprophet_result']['probability'];

	#Read through the protein field list and get all protids that mapped to this peptide
	protids = [];
	isfwd = 0;

	protlist = pepxrec['search_hit'][0]['proteins'];

	for p in protlist:

		protids.append(p['protein']);

		#If none of the protids are fwd then isfwd ==0 (default)
		#But if even one of the mapped protids is fwd
		#the isfwd is set to 1.
		if(not p['protein'].startswith("rev_")):

			isfwd =1;
	
	protidString = ";".join(protids);

	#Convert a list of all extracted fields to str type and join into a tab sep string
	outstring = "\t".join(map(str,[peptide,modpeptide, specid, mass, charge, iniprob, pprophprob, isfwd, protidString]));

	out.write(outstring+"\n");

	lineNo +=1;

	if lineNo % 10000 ==0:

		print "Finished processed line no:"+str(lineNo)+"\r";


print "\n\nAll done!\n";


