% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.r
\name{eMP}
\alias{eMP}
\title{Expected Motif Probability}
\usage{
eMP(a, b, width, nb)
}
\arguments{
\item{a}{numerical vector PDF of motif a}

\item{b}{numerical vector PDF of motif b}

\item{width}{width of fasta indicies}

\item{nb}{Number of Basepairs seperating the start of Motif A
From motif B}
}
\value{
expectation for each motif distance
}
\description{
Determine the expected probability of finding two motifs at a specific
distance based on the empirical PDFs of the individual motifs
}
