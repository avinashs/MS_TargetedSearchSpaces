#Targeted Search spaces for improving MS/MS peptide identification

Repository contains scripts implementing the methods discussed in the paper
"Effective leveraging of Targeted search spaces for improving peptide identification in Tandem mass spectrometry based proteomics"

Link: http://www.ncbi.nlm.nih.gov/pubmed/26569054

##Usage

#####gpm.label.humanpeptides.istryp.nmiss1.freq.tsv
Sample file of GPMDB peptides present in the Ensembl human proteome file and their frequency in the GPMDB repository. Retrieved and processed on Aug 10th 2014.

#####makeFasta.R
Script takes prev freq file and a frequency percentile threshold to use and creates a fasta file containing only peptides with frequency above given threshold. See Fig. 7 in above paper for details about selecting frequency percentile thresholds.

#####addPepRevDecoys2Fasta.py
Script creates decoy sequences by keeping the tryptic sites on the sequences fixed but reversing the sequences between the tryptic sites. This method of creating decoys keeps the decoy peptides consistent between the targeted and full search DBs and essential to be able to combine the searches later. Also this method of decoy generation has the added benefit of keeping the mass distribution of peptides exactly the same in the target and decoy spaces, which is useful particularly with high mass accuracy data.

#####getCompletenessCurve.R
If you have results from a full DB search on the data, the getCompletenessCurve.R script can be used to help select frequency percentile thresholds.

####Use the targeted fasta files created by the makeFasta script as the database to run searches on your MS/MS data and process using TPP till PeptideProphet.
If using the combined searches workflow (recommended), combine the pep.xml file from the targeted search with a pep.xml file from a previous full DB search using iProphet. The result can then be further processed with ProteinProphet or any other Protein inference tool. 
If using the basic targeted targeted workflow, the pep.xml file from the targeted search can be directly used for Protein inference. 

#####pepxml*Parse.py
Parsing scripts for .pep.xml files resulting from iProphet, MSGF+ w/TPP or X!Tandem w/TPP. Requires the pyteomics package for python. This parsed output is expected for use with FDR calculation scripts.

#####*.FDR.R
Scripts that can determine the False discovery rates at PSM or Peptide levels. Only necessary if you are interested in the intermediate FDR levels. Not needed to be run for performing protein inference and other further analysis. 


###GPM frequency file creation scripts

#####gpmFreqTableScripts/getGPMFastaPepFreqs.py
Main script; takes a text file containing all peptides from GPMDB (Most upto date version can be obtained by downloading and extracting from the GPMDB MySQL dump Link: ftp://ftp.thegpm.org/gpmdb/tables/), compares and filters to only peptides present in a given fasta file and gets the frequency identification of those peptides in GPMDB. Uses system calls to Grep; so only guaranteed to work on Linux. 
Note: A faster alternative to this script might be using some post processing of the exonerate tool from EBI. But haven't played around with it yet.  

#####gpmFreqTableScripts/findIfTryptic.py
Filters a given peptide list to identify only those that occur as tryptic peptides in a given fasta file. 

#####gpmFreqTableScripts/filterTrypticbyNmissed.py
Filters a given list of tryptic peptides to only those having less than a given no. of missed cleavage sites. 

TO-DO: Merge last 2 scripts into the main getGPMFastaPepFreqs.py script. 

###RNAseq custom database creation scripts

A straightforward way to create a custom DB Bioconductor package http://bioconductor.org/packages/release/bioc/html/customProDB.html
But alternative scripts are provided here that can also be used to create a RNAseq based custom database. 

#####RNAseqCustomDBscripts/calcrpkm.R
Script takes a .bam file from Tophat and provides RPKM values for transcripts using the Genomic Features package from Bioconductor. 

#####RNAseqCustomDBscripts/gtf2txdb.R
Optional script to parse and create a txdb from a GTF file if the annotations planned to be used is not directly available online. 

#####RNAseqCustomDBscripts/enst2ensp.R
Scripts maps enst transcript ids to corresponding ensp protein ids in the output of the calcrpkm.R script. 

#####RNAseqCustomDBscripts/filterFastaByRpkm.pl
Script takes the output of the enst2ensp.R script and a fasta file (with protIds corresponding to those in the RPKM file) and filters the fasta file to only have the proteins passing a given RPKM threshold. 


