% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_cvm_density}
\alias{plot_cvm_density}
\title{Plot distributions and CVM score for a gene.}
\usage{
plot_cvm_density(cvmobj, gene_name)
}
\arguments{
\item{cvmobj}{An \code{\link{CVMomics}} object, typically returned via a call
to \code{\link{calculate_cvm}}.}

\item{gene_name}{The gene to visualize. The name should be defined as a row
name in \code{cvmobj$cvm}.}
}
\value{
A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
not assigned, a plot will be drawn.
}
\description{
The data for the specified gene is retrieved from
\code{cvmobj$data}. \code{outcomes} is used to divide the data into
distributions for each group, which are then visualized as
density distributions. The calculated CVM score for the specified gene is
displayed in the plot title.
}
\examples{
# 100 genes, 100 samples
dat <- matrix(rnorm(10000), nrow=100, ncol=100)
rownames(dat) <- paste("gene", 1:100, sep="")
colnames(dat) <- paste("sample", 1:100, sep="")

# "A": first 50 samples; "B": next 30 samples; "C": final 20 samples
outcomes <- c(rep("A",50), rep("B",30), rep("C",20))
names(outcomes) <- colnames(dat)

results <- calculate_cvm(dat, outcomes, nperm=10, parallel=FALSE)
plot_cvm_density(results, "gene5")
}
\seealso{
\code{\link{calculate_cvm}} \code{\link[ggplot2]{ggplot}}
}

