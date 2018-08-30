# mT1
## Motifs in Tandem With One Another


**mT1** identifies the distance between DNA motifs in subsets of the genome. First it builds the empirical PDF from individual motifs, then it generates the expectations for the distances between the two motifs. **mT1** was designed to be applied to ChIP-Seq peaks.
___
### Installation
Currently only the source is available, compiled versions for Windows and Mac will be made available once mT1 is submitted to Bioconductor. There are two easy ways to install the package, via command line, and within R itself. The `devtools` package is required to install from within R.


From the command line:

```sh
## Clone the repository
git clone https://github.com/alexjgriffith/mT1.git
cd mT1

## Install 
R CMD INSTALL .
```
From within R:


```R
## Check if devtools is installed
if(!require(devtools)){
	install.packages("devtools")
	library(devtools)
}
## Install mT1
devtools::install_github("alexjgriffith/mT1")

## Install mT1 from the develop branch
## devtools::install_github("alexjgriffith/mT1",ref="develop")
```

The main branch is guaranteed to pass `R CMD check .` with no warnings once mT1 is build using `R CMD build --resave-data .`. The merges to the develop branch should always pass, however it is not a guarantee.


Note that mT1 also requires the Bioconductor package `Biostrings`.

To install `Biostrings`:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```

___
### Examples
Compare the distances between two motifs in a set of genome subsets.

Find distances between two motifs CANNTG Ebox and GATAA GATA on each string they share at four different cellular environments. The distances are relative to the start of each motif:
```{r}
distanceLeukeima<- motifDistance(fastaLeukemia, "CANNTG", "GATAA")
distanceErythroid <- motifDistance(fastaErythroid, "CANNTG", "GATAA")
distanceHSC <- motifDistance(fastaHSC, "CANNTG", "GATAA")
distancesECFC <- motifDistance(fastaECFC, "CANNTG", "GATAA")
```

Plot the motif-motif distance CANNTG-GATAA frequency results between -20 and 20 bp:
```{r}
x <- seq(-20, 20)

yLeukemia <- combHeights(x, distanceLeukeima[,2])[[1]]
yErythroid <- combHeights(x, distanceErythroid[,2])[[1]]
yHSC <- combHeights(x, distanceHSC[,2])[[1]]
yECFC <- combHeights(x,distancesECFC[,2])[[1]]

plot(x, yLeukemia, 
     main = "CANNTG-GATAA", 
     xlab = "distance(bp)", ylab = "Frequency", 
     type = "l")
plot(x, yErythroid, 
     main = "CANNTG-GATAA", 
     xlab = "distance(bp)", ylab = "Frequency", 
     type = "l")
plot(x, yHSC, 
     main = "CANNTG-GATAA", 
     xlab = "distance(bp)", ylab = "Frequency", 
     type = "l")
plot(x,yECFC,
     main="CANNTG-GATAA",
     xlab="distance(bp)", ylab="Frequency",
     type="l")
```


Write the list of motifs (The list of motifs can be written based on personal interests)
```{r}
motifs<-c("CANNTG","GATAA", "MCWTNT", "HGATAA", "GGGGA", "TAAWNG", "TGCCC", "TGACA", "TTGGC", "GCCAA")
```


Find preferred distaces among all motifs and identify composite motifs. Every motif is compared with all motifs inside the motifs list:
```{r}
objMT1_Leukemia <-mT1(fastaLeukemia,motifs)
objMT1_Erythroid <- mT1(fastaErythroid, motifs)
objMT1_HSC <- mT1(fastaHSC, motifs)
objMT1_ECFC<-mT1(fastaECFC,motifs)
```

Plots all motif-motif distance relationship graphs based on the string:
```{r}
plot(objMT1_Leukemia, "CANNTG", "GATAA")
plot(objMT1_Erythroid, "CANNTG", "GATAA")
plot(objMT1_HSC, "CANNTG", "GATAA")
plot(objMT1_ECFC,"CANNTG","GATAA")
```



