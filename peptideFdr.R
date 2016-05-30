#Avinash Shanmugam, Dec 19th 2013
#Script will take a parsed pepXML (tppXMLparser run on pep.xml file)
#with protid and isfwd columns added (using addIsfwd2pepxml.R).
#It will sort this file based on iniProb value, compute FDR at each 
#decreasing iniprob value, and provide the iniProb threshold at which
#a given FDR value is reached. Note: The FDR will be computed on the basis
#of number of unique peptides from fwd and decoys, not number of PSMs

library(data.table)

args = commandArgs(trailingOnly=T);

if( length(args) < 2)
{
	stop("USAGE: peptideFdr.R <.pepxmlparse.tsv> <targetfdr>");
}

pepxmlfile = args[1];
targetFdr = as.numeric(args[2]);

#Create string for outfile name
outfile = sub("pepxml.parse.tsv",paste("pepatfdr.",targetFdr*100,"pc.tsv",sep=""),pepxmlfile);
#outfile = sub("pep.parse.pr.tsv","pepatfdr.1pc.tsv",pepxmlfile);

pepxml = read.table(pepxmlfile,sep="\t",header=T,as.is=T);

names(pepxml) = tolower(names(pepxml));

print(nrow(pepxml));
#In some situations peptides are not mapped to proteins in the fasta file
#Find no. unmapped and subset to only peptides that are mapped to proteins
#nUnmapped = length(unique(pepxml$modpeptide[is.na(pepxml$protid)]));

#pepxml = pepxml[!is.na(pepxml$protid),];

#In order to calculate FDR at the unique peptide level,
#we need to collapse the number of PSMs from the same peptide
#into one row. The maximum iniProb value out of all PSMs is
#stored as maxIniProb.

psmdt = data.table(pepxml);
setkey(psmdt,modpeptide);

#Collapse to unique peptides, (modified peptides treated as separate peptides)
#Also get the maximum iniprob value recorded for each peptide. Additonally
#if a peptide maps to both fwd and decoy sequences, it is treated as a fwd
pepdt = psmdt[,list(maxiniprob = max(iniprob), isfwd=max(isfwd)),by="modpeptide,peptide"];
#pepdt = psmdt[,list(maxiniprob = max(iniprob), isfwd=min(isfwd)),by="modpeptide,peptide"];

#sort the pepdt in decreasing order of maxIniprob
pepdt = pepdt[order(-maxiniprob)];


## Computing R-factor correction for pepdt
rfactor = nrow(pepdt[pepdt$maxiniprob ==0 & pepdt$isfwd ==1,]) / nrow(pepdt[pepdt$maxiniprob ==0 & pepdt$isfwd ==0,]);

print(paste("R-factor=",rfactor));


#From pepdt create miprobGroups data.table, in which rows with same maxiniprob are grouped
#together and the no. of fwd and decoy seqs in each group is also determined

miprobGroups = pepdt[,list( nf =sum(isfwd),nr = (length(isfwd) - sum(isfwd))), by=maxiniprob];


#To the above table add columns storing the cumulative sum of no. forwards and
#no. reverse peptides upto each maxiniprob level
miprobGroups$nfCml = cumsum(miprobGroups$nf);
miprobGroups$nrCml = cumsum(miprobGroups$nr);

#Compute FDR at each maxiniprob level
miprobGroups$fdr = (rfactor * miprobGroups$nrCml) / miprobGroups$nfCml;

#Get the minimum maxiniprob value which has FDR greater than the targetFdr
#i.e. the maxiniprob threshold above which fdr exceeds targetFdr
miprobThreshold = min(miprobGroups$maxiniprob[miprobGroups$fdr <= targetFdr]);

#Using the maxiniprob threshold, get the subset datatable of only peptides
#above the maxiniprob threshold
pepAboveFdr = pepdt[pepdt$maxiniprob >= miprobThreshold,];

nForwards = sum(pepAboveFdr$isfwd ==1);
nDecoys = sum(pepAboveFdr$isfwd == 0);

print(paste("maxiniprob threshold=",miprobThreshold));
print(paste("Nforwards =",nForwards));
print(paste("NDecoys =",nDecoys));

#Write the pepAboveFdr data.table to outfile
write.table(pepAboveFdr,file=outfile,sep="\t",row.names=F,quote=F);

print(paste("Peptides above target fdr written to",outfile));
