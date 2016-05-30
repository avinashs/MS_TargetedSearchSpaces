#Avinash Shanmugam, Dec 19th 2013
#Script will take a parsed pepXML (tppXMLparser run on pep.xml file)
#with protid and isfwd columns added (using addIsfwd2pepxml.R).
#It will sort this file based on iniProb value, compute FDR at each 
#decreasing iniprob value, and provide the iniProb threshold at which
#a given FDR value is reached. It will also output the no. of unique
#peptide seqs at that threshold.


args = commandArgs(trailingOnly=T);

if( length(args) < 2)
{
	stop("USAGE: peptideFdr.R <.pepxmlparse.pr.tsv> <targetfdr>");
}

pepxmlfile = args[1];
targetFdr = as.numeric(args[2]);

#Create string for outfile name
outfile = sub("pepxml.parse.pr.tsv","psmatfdr.1pc.tsv",pepxmlfile);

pepxml = read.table(pepxmlfile,sep="\t",header=T,as.is=T);

names(pepxml) = tolower(names(pepxml));

#Sort df by iniprob
pepxml = pepxml[with(pepxml,order(-iniprob)),];

#Get unique values of iniprob. Since the df is sorted, the unique values
#will also be in descending order
iprobVals = unique(pepxml$iniprob);

#pepxml$fdr = 0;
endIdx=1;

for(iprob in iprobVals)
{
	#Get the isfwd values for PSMs with iniprob value above iprob
	isfwd = pepxml$isfwd[pepxml$iniprob >= iprob];

	f = sum(isfwd ==1);
	r = sum(isfwd ==0);

	fdr = r/f;

#	print(paste("iprob=",iprob,"f=",f,"r=",r,"fdr=",fdr,"idx=",length(isfwd)));
	
	#As long as the fdr stays below the trarget, keep updating
	#the nrows for that iprob as the current endIdx. When it crosses
	#the targetFdr, break and leave the loop.

	if(fdr <= targetFdr)
	{
#		pepxml$fdr[endIdx:length(isfwd)] = fdr;
		endIdx = length(isfwd);
	}
	else
	{
		print(paste("fdr=",fdr,"is greater than target,",targetFdr));
		break;
	}
}

#Once the loop has broken, take the last endIdx that was stored
#and take all PSMs till that row
psmAtFdr = pepxml[1:endIdx,];

iprobThreshold = min(psmAtFdr$iniprob);

#Remove all decoys from psmAtFdr df
psmAtFdr = psmAtFdr[psmAtFdr$isfwd ==1,];

npep = length(unique(psmAtFdr$peptide[psmAtFdr$isfwd ==1]));
nmodpep = length(unique(psmAtFdr$modpeptide[psmAtFdr$isfwd ==1]));

write.table(psmAtFdr,file=outfile,sep="\t",row.names=F,quote=F);

print(paste("Npeptides =",npep));
print(paste("NModpeptides =",nmodpep));
print(paste("Iniprob threshold =",iprobThreshold));











