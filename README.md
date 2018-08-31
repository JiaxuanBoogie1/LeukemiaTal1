# LeukemiaTal1
The combination of ChipPComp and mT1 from alexjgriffith Github

## Refer to alexjgriffith ChipPComp: https://github.com/alexjgriffith/ChipPComp
## Refer to alexjgriffith mT1: https://github.com/alexjgriffith/mT1

## Introduction

Transcription factors are proteins that have the ability to bind to special locations on the DNA based on nucleotide sequence and chromatin conformation. TAL1 is an example of a TF with a critical role in several lineages, including the hematopoietic and the endothelial. It has different functions inside different cellular conditions. By conducting the Chip-seq with four different cellular environments leukemia, hemopoietic stem cell (HSC), erythroid and endothelia colony forming cell (ECFC), peak calling is operated through MACS. The whole programme is composed by two parts: ChipPComp and mT1. 

ChipPComp is the section to identify unique peak regions of specific cellular condition after peak calling. To define regions, alternative feature spaces (AFS) are found as the combination of the peaks into a single unified set, after that density matrix (UDM) is applied to read count pile ups. UDM is the method to normalize the width of peaks and provide a read count density. Finally, principal component analysis (PCA) is applied to the UDM. PCA is used to normalize the density matrix and then preform quantile normalization and return the eigenvectors in a list with the normalized data. By plotting PC graphs, unique peak region of specific cellular environment can be defined, and the peak regions data were added to the peak database to produce the bed file and fasta file for each cellular environment. 

mT1 part is run after the bed files and fasta files about defined peak regions for each cellular environment are produced from ChipPComp. EBox CANNTG and Gata GATAA are two common motif binding sites for Tal1 transcription factor. There are two common motifs binding mechanisms for Tal1. One way is Tal 1 will dimerize with E-protein and directly bind to EBox DNA with bHLH DNA binding domain. Alternatively, Tal1 may interact indirectly through recruitment by other TFs as part of a complex such as GATA. It is important to find the distance between two motifs at different cellular environments to see how distances effect motif functions. In general, the distances between two motifs are 12bp to 15bp, but the accurate distance depends on cells. Besides Ebox-GATA composite motif, the distance between other composite motifs can also be calculated. The preferred distances for Ebox, GATA and a list of motifs (based on personal research interests) are compared as well to find differences among cellular environments. 

## Installation:
Prerequisites:
R (>= 3.4.0) to run LeukemiaTal1

Install from source:
Install devtools at first which contains the command to directly install package from github. To install devtools:
```{r}
>install.packages(“devtools”)
Source the devtools package in R every time before use functions inside devtools:
```{r}
>library(devtools)
```
Install the Leukemia_Tal1 package into R with the help of the function inside devtools:
```{r}
>install_github(“JiaxuanBoogie1/LeukemiaTal1”)
>library(LeukemiaTal1)
```

Install other required packages:

Source Bioconductor which provides tools for the analysis and comprehension of high-throughput genomic data
```{r}
>source(https://bioconductor.org/biocLite.R)
```

Install full genome sequences for Homo sapiens:
```
>biocLite(“biostrings”)
>biocLite(“BSgenome.Hsapiens.USCS.hg19”)
```

To use BSgenome.Hsapiens.USCS.hg19, source it every time before using it:
```
>library(BSgenome.Hsapiens.USCS.hg19)
```

Install GenomicRanges to produce fasta file from defined peak regions:
```
>biocLite(“GenomicRanges”)
To use GenomicRanges, source it every time before use it:
>library(GenomicRanges)
```

## Usage of LeukemiaTal1:
Example for forming the ccca data frame:
>makeCCCA (-r unique_nodups.bed, -p peaks.xls, -c, -mco)

Example for adding regions:
>addRegion (-d, -n, -l) 

Example for finding the distance between two motifs:
>motifDistance (-f .fasta, “-m1”, “-m2”)

Example for finding preferred distances among motifs under peaks:
>mT1( -f fasta file, -m)

## Explanation for important functions:
makeCCCA:	Form the data frame with AFS, UDM, PCA, regions, and categories results from all cell lines. 

plotPrincipalComponents:	Plot request 2-dimensional principal component graphs  

normalizePRC:	Normalize the principal components

addRegion:	Add the defined peak region for each cellular environment by observing principal component graphs

motifDistance:	Find the distances between the composite motifs within the set of genomic range

combHeights:	Determine the number of occurrences for each member along a specific sequence 

mT1:	Investigation of novel preferred distance between a list of motifs. The input motifs are based on personal research interests (at least three motifs). It also includes probability density function results for each motif. 


## makeCCCA:
This is the main function inside Leukemia_Tal1. This function contains sub functions makeAFS (find the alternative feature space), makeUDM (unified all AFS to make unified density matrix), and makePRC (find the principal component analysis) 

Option: 
-r -- rawData
This is the REQUIRED parameter for makeCCCA. rawData is a list of raw read data files inside c(…) and files must be in the name.bed format. 

-p – peaklist
This is the REQUIRED parameter for makeCCCA. Files must be in the .xlsx format which contains a list of all raw peak files c(…) for analysis (from MACS). All .bed files inside rawData list must have a corresponding .xlsx peak file. The list of .xlsx peak must have the same order as the list of .bed raw data. 

-c – category
The is the REQUIRED parameter for makeCCCA. The list of names for all cell lines inside c(…). The order of names must be corresponding to the raw data .bed files and peak .xlsx files.  

-mco – macs cut off (pvalue)
The is the REQUIRED parameter for makeCCCA. Due to peak bias, the macs cut off is used to determine enriched peak locations and cut off bias, the cut offs used must be stringent. The standardized macs cut off value is 20, which is equal to the p value 10-20.


## addRegion:
To define and separate the specific peak regions for each cellular environment based on principal component analysis and all PC comparison plots.

Option:
-d -- database 
The is the REQUIRED parameter for addRegion. The database must be the name of the result database from makeCCCA.  

-n -- tag name
The is the REQUIRED parameter for addRegion. The tag / name of each cellular environment defined region inside the database. The tag will be added to database$reg 

-l -- logic
The is the REQUIRED parameter for addRegion. The logic command to define regions for each cellular environment. The separation logic depends on the normalized principal component (PC) number and standard deviation (SD) number of the PC. The values of PC and SD must be determined through 2-dimensional PC plots.

## motifDistance:
To find the distance between the first and the second motifs inside the genomic ranges. The distances are related to the start of each motif.

Option:
-f -- fasta file 
This is the REQUIRED parameter for motifDistance. The fasta file is extracted from Homo sapiens genome and each defined peak region from ChipPComp.

-m1 -- first motif
This is the REQUIRED parameter for motifDistance. Inside the quotation mark, the exact nucleotide composite of the first motif must be written. For leukemia, the first motif is Ebox “CANNTG”. 

-m2 -- second motif
This is the REQUIRED parameter for motifDistance. Inside the quotation mark, the exact nucleotide composite of the second motif must be showed. For leukemia, the second motif is GATA “GATAA”. 


## mT1:
This is one of the main functions in Leukemia_Tal. mT1 is used to investigate novel preferred motif distances between composite motifs. The input motif list must contain at least three motifs. 

Options:
-f -- fasta file 
This is the REQUIRED parameter for mT1. The fasta file is extracted from Homo sapiens genome and each defined peak region from ChipPComp. -f in mT1 is same as -f in motifDistance.

-m -- motifs under peaks
This is the REQUIRED parameter for mT1. It is a string of motifs to be compared within the fasta sequencing file. The string of motifs contains the first motif Ebox “CANNTG”, the second motif GATA “GATAA” and a list of motifs for your personal research interests. The novel preferred distances among all these motifs are compared with each other.

-verbose 
Optional parameter for mT1. Verbose = FALSE is the default parameter. If it is true, the flag will be printed as motifs completed.

-cl -- cluster
Optional parameter for mT1. cl is a cluster from makeForkCluster. cl=NULL is the default parameter. 


## Output:
1.	The output database contains information about all peaks. Information includes:
  afs (alternative feature space)
  udm (unified density matrix)
  prc$normData (normalized data from principal component analysis)
  prc$eigenVectors (eigenvectors after principal component analysis)
  fasta (the exact homo sapiens nucleotide sequencing for defined peaks in each cellular environment)
  reg (newly defined regions for each cellular environment)
  categories (the list of cellular categories in order)

2.	The output Tag.bed file is a bed file which contains information about peaks in defined peak regions for each cellular environment. You can open it in excel and sort/filter using excel functions. Information includes:
chromosome name
start position of peak
end position of peak

3.	The output name.fasta file is a fasta file which contains information about actual homo sapiens nucleotide sequence of all peaks inside defined peak regions for each cellular environment. To get the fasta file, use getSeq function inside GenomicRanges package. Remember to convert the Tag.bed file into GeneRangelist before running getSeq function. 
Information includes:
DNAStringSet: width of one group of nucleotide sequence and actual nucleotide sequences. 
 

4.	CANNTG-GATAA composite motif preferred distance versus frequency plot. The linear plot with x-axis Ebox CANNTG – GATA GATAA motif distance and y-axis frequency is produced after running motifDistance. The plot helps to find the preferred distance (bp) between composite motifs in different cellular environments. The input composite motifs can be changed based on personal research interest. 


5.	Six composite motifs distance plots. The first plot is the Ebox CANNTG versus probability density function at the specific genomic location and the second one is the GATA GATAA versus probability density. The third plot is the convolution graphs for both motifs. The fourth plot is the frequency of occurrences in specific CANNTG-GATAA motif distance. The fifth one is the peak heights versus frequency plot and the last one is test p-value plot for each motif distance. All these plots represent preferred distances and relationships among Ebox-GATA composite motifs. Differences among all plots with different cellular environments can be observed. 





