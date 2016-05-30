#Avinash Shanmugam, April 10th 2014
#Script will take protxml.parse.tsv, process it as required
#find FDR thresholds and o/p the no of fwds above FDR threshold

library(data.table);
library(plyr);

args = commandArgs(trailingOnly=T);

if( length(args) < 2)
{
	stop("USAGE: proteinFDR.R <protxml.parse.tsv> <targetFdr>");
}

infile = args[1];
targetFdr = args[2];

outfile = sub("prot.parse.tsv","protAtFdr.1pc.tsv",infile);

#Read in dataframe
d = read.table(infile,sep="\t",header=T,as.is=T);

#Convert all names to lowercase for ease of use
names(d) = tolower(names(d));

#Thresholds to filter before considering PSMs
localpwTh = 0.5;
iniprobTh = 0.5;

#Filter d by thresholds
d = d[d$localpw >= localpwTh & d$iniprob >= iniprobTh,];

dt = data.table(d);
setkey(dt, protid);

#Collapse to protid
protdt = dt[,list(maxiniprob = max(iniprob),groupid,siblinggroup,localpw),by="protid"];

protdt = unique(protdt);

#Convert back to data frame
protdf = data.frame(protdt);

#Group and get the protids with the maximum maxiniprob
#andmaximum localpw values from each siblinggroup
topProtDf = ddply(protdf, c("groupid","siblinggroup"),function(df) df[df$maxiniprob == max(df$maxiniprob) & df$localpw == max(df$localpw),]);

#If there is more than one protid with same maxiniprob and
#localpw in the same siblinggroup, just pick the first one
#alphabetically
singleProtDf = ddply(topProtDf, c("groupid","siblinggroup"),function(df) df[df$protid == min(df$protid),]);

#Add isfwd col to dataframe
singleProtDf$isfwd= 1

singleProtDf$isfwd[grep("rev_",singleProtDf$isfwd,fixed=TRUE)] =0;

#Convert back to data.table
singleProtDt = data.table(singleProtDf);
setkey(singleProtDt,maxiniprob);

#Create a data.table grouped by the maxiniprob
miprobGroups = singleProtDt[,list(nf =sum(isfwd), nr=sum(isfwd ==0)),by="maxiniprob"];

#Get cumulative sums of nf and nr
miprobGroups$nfCml = cumsum(miprobGroups$nf);
miprobGroups$nrCml = cumsum(miprobGroups$nr);

#Compute FDR at each maxiniprob level
miprobGroups$fdr = miprobGroups$nrCml / miprobGroups$nfCml;

#Get the minimum maxiniprob val with FDR great than targetFdr
miprobThreshold = min(miprobGroups$maxiniprob[miprobGroups$fdr <= targetFdr]);

#Using the maxiniprob threshold, get the subset datatable of
#only proteins above the maxiniprob threshold
protAboveFdr = singleProtDt[singleProtDt$maxiniprob >= miprobThreshold,];

nForwards = sum(protAboveFdr$isfwd);
nDecoys = sum(protAboveFdr$isfwd ==0);

print(paste("maxiniprob threshold=",miprobThreshold));
print(paste("Nforwards=",nForwards));
print(paste("NDecoys=",nDecoys));

#Write the protAboveFdr data.table to outfile
write.table(protAboveFdr,file=outfile,sep="\t",row.names=F,quote=F);

print(paste("Peptides above target fdr written to",outfile));







