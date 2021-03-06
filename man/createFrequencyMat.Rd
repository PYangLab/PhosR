% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMotifFreq.R
\name{createFrequencyMat}
\alias{createFrequencyMat}
\title{Create frequency matrix}
\usage{
createFrequencyMat(substrates.seq)
}
\arguments{
\item{substrates.seq}{A substrate sequence}
}
\value{
A frequency matrix of amino acid from substrates.seq.
}
\description{
Create frequency matrix
}
\examples{

data("phospho_L6_ratio_pe")

# We will create a frequency matrix of Tfg S198 phosphosite.
idx = which(grepl("TFG\\\\;S198\\\\;", rownames(phospho.L6.ratio.pe)))
substrate.seq = Sequence(phospho.L6.ratio.pe)[idx]
freq.mat = createFrequencyMat(substrate.seq)

}
