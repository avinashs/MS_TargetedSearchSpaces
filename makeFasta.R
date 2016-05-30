#Avinash Shanmugam, Oct 8th 2014
#Script takes a peptide freq file and creates creates subsets of
#the dataframe corresponding to defined quantile thresholds
#and then outputs these subset peptide sequences as fasta files.

makeFastafrmVector =function(peptides, outfile)
{
	#Create fastaDf which has 2 identical cols of peptides
	fastaDf = data.frame(peptides,peptides);

	names(fastaDf) = c("title","peptide");

	#Add '>' to title col
	fastaDf$title = paste(">",fastaDf$title,sep="");

	#Write to outfile
	write.table(fastaDf,file=outfile, sep="\n", col.names=F,row.names=F,quote=F);

	return(0);
}


args = commandArgs(trailingOnly = TRUE);

if( length(args) != 2)
{
	stop("USAGE: makeFasta.R <input pepfreqsFile> <outfastaFile prefix>");
}

freqsfile = args[1];
fastaPrefix = args[2];

#Definied quantile thresholds to subset
qntlThresholds =c(0.92,0.95,0.98);

d = read.table(freqsfile,sep="\t",as.is=T);

names(d) = c("peptide","freq");

#Get the frequencies corresponding to the quantile
#thresholds
freqThresholds = quantile(d$freq, qntlThresholds);

for( f in freqThresholds)
{	
	#Subset according to freq Threshold
	dSub = d[d$freq >= f,];
	
	outfile = paste(fastaPrefix,".freq",f,".fa",sep="");

	#Send peptides for fasta file creation
	makeFastafrmVector(dSub$peptide, outfile);

	print(paste("Finished printing fasta file:",outfile));

}

	

	
