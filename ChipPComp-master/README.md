# ChipPComp

## Refer to alexjgriffith github: https://github.com/alexjgriffith/ChipPComp.

These examples use the data sets that will be in ChipPCompData, change the root directory to wherever you have this package saved.

## Installation

Install devtools and ChipPComp from github. ChipPComp package is used to define specifc PCA peak regions for each celluar environment:
```{r}
>install.packages("devtools")
>library(devtools)

>install_github("alexjgriffith/ChipPComp")
>library(ChipPComp)

>source("https://bioconductor.org/biocLite.R")

>biocLite("Biostrings")
>library(Biostrings)
```

## Example

List the work directory as root and list all bed files, peak files and categories for all 15 cell lines (Notice: all file names inside peaklist, rawData and categories should corresponding to others and keep in the same files order):
```{r
>root <- "/home/work/Desktop/Data"

>peaklist <- c("00-000-A_Prima2_Tal1_peaks.xls", "00-000-A_Prima5_Tal1_peaks.xls", "00-000-B_Jurkat_Tal1_peaks.xls", "00-000-C_Erythroid_Tal1_peaks.xls", "00-000-C_Jurkat_Tal1_peaks.xls", "00-000-D_RPMI_Tal1_peaks.xls", "00-000-E_CCRF-CEM_Tal1_peaks.xls", "00-000-F_ECFC_Tal1_peaks.xls", "00-000-G_HSPC_Tal1_peaks.xls", "00-000-H_Meka_Tal1_peaks.xls", "00-000_HSPC_Tal1_peaks.xls", "00-000-I_CD34_Tal1_peaks.xls", "00-000-J_Erythroid-Fetal_peaks.xls", "00-000-K_Erythroid-Adult_peaks.xls", "00-000-L_K562_Tal1_peaks.xls")

>rawData <- c("00-000-A_Prima2_Tal1_Rep1-2_unique_nodups.bed", "00-000-A_Prima5_Tal1_Rep1-2_unique_nodups.bed", "00-000-B_Jurkat_Tal1_Rep1-2_unique_nodups.bed", "00-000-C_Erythroid_Tal1_Rep1-2_unique_nodups.bed", "00-000-C_Jurkat_Tal1_Rep1-2_unique_nodups.bed", "00-000-D_RPMI_Tal1_Rep1-2_unique_nodups.bed", "00-000-E_CCRF-CEM_Tal1_Rep1-2_unique_nodups.bed", "00-000-F_ECFC_Tal1_Rep1_unique_nodups.bed", "00-000-G_HSPC_Tal1_Rep1-2_unique_nodups.bed", "00-000-H_Meka_Tal1_Rep1_unique_nodups.bed", "00-000_HSPC_Tal1_Rep1_unique_nodups.bed", "00-000-I_CD34_Tal1_Rep1_unique_nodups.bed", "00-000-J_Erythroid-Fetal_Tal1_Rep1_unique_nodups.bed", "00-000-K_Erythroid-Adult_Tal1_Rep1_unique_nodups.bed", "00-000-L_K562_Tal1_Rep1-2_unique_nodups.bed")

>categories <- c("Prima2", "Prima5", "Jurkat_Rep1", "Erythroid", "Jurkat_Rep1-2", "PRMI", "CEM", "ECFC", "CD133", "Meka", "HSPC_Rep1", "CD34", "Erythroid_fetal", "Erythroid_adult", "K562")
```


Find and analyze the alternative feature space, unified density matrix and principal component analysis of all 15 cell lines and show the result inside the ccca database:
```{r}
>ccca <- makeCCCA(rawData, peaklist, categories, 20)
```


Plot all principle component analysis graphs for the ccca database:
```{r}
>ChipPComp::plotPrincipalComponents(ccca, c(1,2))
>ChipPComp::plotPrincipalComponents(ccca, c(1,3))
>ChipPComp::plotPrincipalComponents(ccca, c(1,4))
```


Setup categories for cell lines, and the column numbers are based on the order of column on ccca$afs:
```{r
>head(ccca$afs)

>Leukemia.columns <- c(4,5,6,8,9,10)
>Erythroid.columns <- c(7, 16, 17, 18)
>HSC.columns <- c(12:15)
>ECFC.columns <- 11
```


Create vector to indicate each peak only present in one cell type. rowSums is used to count the number of cell lines in peaks, then convert it into boolean to check if the result is greater than zero:
```{r}
>Peaks.in.Leukemia <- rowSums(ccca$afs[, Leukemia.columns]) >0
>Peaks.in.Erythroid <- rowSums(ccca$afs[, Erythroid.columns]) >0
>Peaks.in.HSC <- rowSums(ccca$afs[, HSC.columns]) >0
>Peaks.in.ECFC <- ccca$afs[, ECFC.columns]>0
```


Create the boolean to check if each peak is only present in one cell type. Then use the cbind to join all peaks from four cellular environments into one matrix:
```{r}
>Peaks.in.one.type <- rowSums(cbind(Peaks.in.Leukemia, Peaks.in.Erythroid, Peaks.in.HSC, Peaks.in.ECFC))==1
```


Check PC plots without coloring:
```{r}
>plot(ccca$prc$eigenVectors[Peaks.in.one.type, c(1,2)])
>plot(ccca$prc$eigenVectors[Peaks.in.one.type, c(1,4)])
```


Setup an empty vector for coloring:
```{r}
>my.colours <- NULL


Define coloring rule for each cellular environment on PC plots:
>my.colours[Peaks.in.one.type & Peaks.in.Leukemia] <- "blue"
>my.colours[Peaks.in.one.type & Peaks.in.Erythroid] <- "red"
>my.colours[Peaks.in.one.type & Peaks.in.HSC] <- "yellow"
>my.colours[Peaks.in.one.type & Peaks.in.ECFC] <- "green"
```


Plot PC plots for peaks and coloring peaks based on categories:
```{r}
>plot(ccca$prc$eigenVectors[Peaks.in.one.type, c(1,2)], 
     col=my.colours[Peaks.in.one.type])
>plot(ccca$prc$eigenVectors[Peaks.in.one.type, c(1,4)], 
     col=my.colours[Peaks.in.one.type])
```



Normalize all principal component analysis data and then add specific defined region of each cellular environment into the ccca database. The specific values of PC and sd are based on PC plots above:
```{r}
>nm <- normalizePRC(ccca$prc)

>ccca <- addRegion(ccca, "Leukemia", 
                  nm[,1]<(mean(nm[,1])-1*sd(nm[,1])))
>ccca <- addRegion(ccca, "Erythroid", 
                  nm[,1]>(mean(nm[,1])+1*sd(nm[,1])))
>ccca <- addRegion(ccca, "HSC", 
                  nm[,4]>(mean(nm[,4])+1*sd(nm[,4])))
>ccca <- addRegion(ccca, "ECFC", 
                  nm[,4]<(mean(nm[,4])-1*sd(nm[,4])))
```


Name all defined region according to the order on ccca$reg (same as the above addRegion order):
```{r}
>Leukemia <- data.frame(chr=ccca$afs$chr[ccca$reg[,1]],
                       start=ccca$afs$start[ccca$reg[,1]], 
                       end = ccca$afs$end[ccca$reg[,1]])
>Erythroid <- data.frame(chr=ccca$afs$chr[ccca$reg[,2]],
                        start=ccca$afs$start[ccca$reg[,2]], 
                        end = ccca$afs$end[ccca$reg[,2]])
>HSC <- data.frame(chr=ccca$afs$chr[ccca$reg[,3]],
                  start=ccca$afs$start[ccca$reg[,3]], 
                  end = ccca$afs$end[ccca$reg[,3]])
>ECFC <- data.frame(chr=ccca$afs$chr[ccca$reg[,4]],
                   start=ccca$afs$start[ccca$reg[,4]], 
                   end = ccca$afs$end[ccca$reg[,4]])


>Write out defined peak regions .bed file for each cellular environment:
>write.table(ECFC, "ECFC.bed", sep = "\t", quote = FALSE, row.names = FALSE)
>write.table(Erythroid, "Erythroid.bed", sep = "\t", quote = FALSE, row.names = FALSE)
>write.table(Leukemia, "Leukamia.bed", sep = "\t", quote = FALSE, row.names = FALSE)
>write.table(HSC, "HSC.bed", sep = "\t", quote = FALSE, row.names = FALSE)


Add peak name inside data frame for fasta file motif identifiers:
>Leukemia$name <- paste0("Leukemeia",Leukemia$chr,":",Leukemia$start,"-",Leukemia$end)
>Erythroid$name <- paste0("Erythroid",Erythroid$chr,":",Erythroid$start,"-",Erythroid$end)
>HSC$name <- paste0("HSC",HSC$chr,":",HSC$start,"-",HSC$end)
>ECFC$name <- paste0("ECFC",ECFC$chr,":",ECFC$start,"-",ECFC$end)
```


Intall GenomicRanges and full Homo sapiens genomic seuqneces:
```{r
>source("https://bioconductor.org/biocLite.R")

>biocLite("GenomicRanges")
>library(GenomicRanges)

>biocLite("BSgenome.Hsapiens.UCSC.hg19")
>library(BSgenome.Hsapiens.UCSC.hg19)
```


Firslty convert all defined region .bed files into genomic range lists. Then convert and write out all genomic range lists into the fasta files for motif analysis:
```{r}
>Leukemia.granges <- GRanges(seqnames = Leukemia$chr, 
                            ranges = IRanges(start = Leukemia$start, 
                                             end =Leukemia$end))
>Erythroid.granges <- GRanges(seqnames = Erythroid$chr, 
                             ranges = IRanges(start = Erythroid$start, 
                                              end = Erythroid$end))
>HSC.granges <- GRanges(seqnames = HSC$chr, 
                       ranges = IRanges(start = HSC$start, 
                                        end = HSC$end))
>ECFC.granges <- GRanges (seqnames = ECFC$chr, 
                         ranges = IRanges (start = ECFC$start, 
                                           end=ECFC$end))

>fastaLeukemia <- getSeq(BSgenome.Hsapiens.UCSC.hg19, Leukemia.granges)
>fastaErythroid <- getSeq(BSgenome.Hsapiens.UCSC.hg19, Erythroid.granges)
>fastaHSC <- getSeq(BSgenome.Hsapiens.UCSC.hg19, HSC.granges)
>fastaECFC <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ECFC.granges)


>names(fastaLeukemia) <- Leukemia$name
>names(fastaErythroid) <- Erythroid$name
>names(fastaHSC) <- HSC$name
>names(fastaECFC) <- ECFC$name

>writeXStringSet(fastaLeukemia, "Leukemia.fasta")
>writeXStringSet(fastaErythroid, "Erythroid.fasta")
>writeXStringSet(fastaHSC, "HSC.fasta")
>writeXStringSet (fastaECFC, "ECFC.fasta")
```

