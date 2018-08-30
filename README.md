# Leukemia_Tal1
The combination of ChipPComp and mT1 from alexjgriffith Github

## Refer to alexjgriffith ChipPComp: https://github.com/alexjgriffith/ChipPComp
## Refer to alexjgriffith mT1" https://github.com/alexjgriffith/mT1

## Introduction

Transcription factors are proteins that have the ability to bind to special locations on the DNA based on nucleotide sequence and chromatin conformation. TAL1 is an example of a TF with a critical role in several lineages, including the hematopoietic and the endothelial. It has different functions inside different cellular conditions. By conducting the Chip-seq with four different cellular environments leukemia, hemopoietic stem cell (HSC), erythroid and endothelia colony forming cell (ECFC), peak calling is operated through MACS. The whole programme is composed by two parts: ChipPComp and mT1. 

ChipPComp is the section to identify unique peak regions of specific cellular condition after peak calling. To define regions, alternative feature spaces (AFS) are found as the combination of the peaks into a single unified set, after that density matrix (UDM) is applied to read count pile ups. UDM is the method to normalize the width of peaks and provide a read count density. Finally, principal component analysis (PCA) is applied to the UDM. PCA is used to normalize the density matrix and then preform quantile normalization and return the eigenvectors in a list with the normalized data. By plotting PC graphs, unique peak region of specific cellular environment can be defined, and the peak regions data were added to the peak database to produce the bed file and fasta file for each cellular environment. 

mT1 part is run after the bed files and fasta files about defined peak regions for each cellular environment are produced from ChipPComp. EBox CANNTG and Gata GATAA are two common motif binding sites for Tal1 transcription factor. There are two common motifs binding mechanisms for Tal1. One way is Tal 1 will dimerize with E-protein and directly bind to EBox DNA with bHLH DNA binding domain. Alternatively, Tal1 may interact indirectly through recruitment by other TFs as part of a complex such as GATA. It is important to find the distance between two motifs at different cellular environment to see how distances effect motif functions. In general, distance between two motifs is 12bp to 15bp, but the accurate distance depends on cells. In addition, Jaspar database contains a list of other preference composite motifs for Tal1. The preferred distances for Ebox, GATA and all motifs inside Jaspar are compared as well to find differences among cellular environments. 

## Installation:
Prerequisites:
R (>= 3.4.0) to run Leukemia_Tal1

Install from source:
Install devtools at first which contains the command to directly install package from github. To install devtools:
```{r}
>install.packages(“devtools”)
Source the devtools package in R every time before use functions inside devtools:
```{r}
>library(devtools)
```
Install the ChipPComp package into R with the help of the function inside devtools:
```{r}
>install_github(“JiaxuanBoogie1/Leukemia_Tal1”)
>library(ChipPComp)
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

To use BSgenome.Hsapiens.USCS.hg19, source it every time before use it:
```
>library(BSgenome.Hsapiens.USCS.hg19)
```

Install GeomicRanges to produce fasta file from defined peak regions:
```
>biocLite(“GenomicRanges”)
To use GenomicRanges, source it every time before use it:
>library(GenomicRanges)
```





