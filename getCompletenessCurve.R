#Avinash Shanmugam, May 20th 2014
#Script will read in a pepatFDR file with gpepfreq values added
#It will find the maximum gpepfreq value, and then loop from
#zero to that value finding the number of peptides at each gpfreq value

library(data.table);

args = commandArgs(trailingOnly = TRUE);

if ( length(args) != 3)
{
	stop("USAGE: getCompletenessCurve.R <pepatfdrFile> <gpm_freqfile> <outfile>");
}

pepfile = args[1];
gpfile = args[2];
outfile = args[3];

d = read.table(pepfile, sep="\t",header=T,as.is=T);
g = read.table(gpfile, sep="\t",header=T,as.is=T);

#Keep only the peptide, modpeptide and gpfreq cols in d
d = d[,c("modpeptide","peptide","gpepfreq")];

#Convert to data table
dd = data.table(d);
gg = data.table(g);


#dd has repeated rows. So get only the unique
#modpeptide rows in dd
setkey(dd,modpeptide);
dd = unique(dd);

#Setkey to the frequencies
setkey(dd, gpepfreq);
setkey(gg, freq);

#Get maximum gpepfreq value
maxgpfreq = max(d$gpepfreq);

print(paste("Max gpfreq val =",maxgpfreq));

#gpfreq vals are closer together in the lower ranges and
#are more spread out more at larger vals. So creating gpfvals
#vector to reflect that.

gpfvals1 = seq(0,100,by=1);
gpfvals2 = seq(110,1000,by=10);
gpfvals3 = seq(1100,10000,by=100);
gpfvals4 = seq(11000,50000,by=1000);

gpfvals = c(gpfvals1, gpfvals2, gpfvals3, gpfvals4);

#Blank numeric vector to store the no. of peptides
npeps = numeric(0);
ngpeps = numeric(0);

for( i in gpfvals)
{

	npep = nrow(subset(dd, gpepfreq >= i));
	ngpep = nrow(subset(gg, freq >= i));

	npeps = append(npeps,npep);
	ngpeps = append(ngpeps,ngpep);

	cat(paste("Finished for val",i,"\r"))
}

npepdf = data.frame(gpfvals,ngpeps,npeps);

write.table(npepdf, file=outfile, sep="\t",row.names=F,quote=F);


