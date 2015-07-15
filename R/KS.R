#' @export
#' @title Calculate the Komolgorov-Smirnov test statistic and q-values.
#' 
#' @description This is the main user interface to the \pkg{EMDomics} package, and
#' will usually the only function needed when conducting an analysis using the Komolgorov-Smirnov
#' algorithm. Analyses can also be conducted with the EMD algorithm using
#' \code{calculate_emd}.
#'
#' The algorithm is used to compare genomics data between any number of groups. 
#' Usually the data will be gene expression
#' values from array-based or sequence-based experiments, but data from other
#' types of experiments can also be analyzed (e.g. copy number variation).
#'
#' Traditional methods like Significance Analysis of Microarrays (SAM) and Linear
#' Models for Microarray Data (LIMMA) use significance tests based on summary
#' statistics (mean and standard deviation) of the two distributions. This
#' approach tends to give non-significant results if the two distributions are
#' highly heterogeneous, which can be the case in many biological circumstances
#' (e.g sensitive vs. resistant tumor samples).
#'
#' Komolgorov-Smirnov instead calculates a test statistic that is the maximum distance between
#' two cumulative distribution functions (CDFs). Unlike the EMD score, the KS test statistic
#' summarizes only the maximum difference (while EMD considers quantity and distance between
#' differences).
#'
#' The KS algorithm implemented in \pkg{EMDomics} has two main steps.
#' First, a matrix (e.g. of expression data) is divided into data for each of the groups.
#' Every possible pairwise KS score is then computed and stored in a table. The KS score
#' for a single gene is calculated by averaging all of the pairwise KS scores. The p-values
#' from the KS test are adjusted using the Benjamini-Hochberg method.
#' Next, the labels for each of the groups are randomly
#' permuted a specified number of times, and an EMD score for each permutation is
#' calculated. The median of the permuted scores for each gene is used as
#' the null distribution, and the False Discovery Rate (FDR) is computed for
#' a range of permissive to restrictive significance thresholds. The threshold
#' that minimizes the FDR is defined as the q-value, and is used to interpret
#' the significance of the EMD score analogously to a p-value (e.g. q-value
#' < 0.05 = significant). The q-values returned by the KS test (and adjusted for multiple 
#' significance testing) can be compared to the permuted q-values.
#'
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param outcomes A vector containing group labels for each of the samples provided
#' in the \code{data} matrix. The names should be the sample identifiers provided in \code{data}.
#' @param nperm An integer specifying the number of randomly permuted EMD
#' scores to be computed. Defaults to 100.
#' @param verbose Boolean specifying whether to display progress messages.
#' @param parallel Boolean specifying whether to use parallel processing via
#' the \pkg{BiocParallel} package. Defaults to \code{TRUE}.
#' @return The function returns an \code{\link{EMDomics}} object.
#' 
#' @seealso \code{\link{EMDomics}} \code{\link{ks.test}}
calculate_ks <- function(data, outcomes, nperm=100, 
                         verbose=TRUE, parallel=TRUE) {
  bpparam <- BiocParallel::bpparam()
  
  if (parallel == FALSE)
    bpparam <- BiocParallel::SerialParam()
  
  # transpose and coerce to df (for bplapply)
  data.df <- as.data.frame(t(data))
  sample_names <- rownames(data.df)
  
  # ---------- pairwise ks table: test statistic and p-value -----------
  
  # generate pairwise emd table for each gene
  if (verbose)
    message("Calculating pairwise ks scores and p-values...", appendLF=FALSE)
  
  # all possible pairwise comparisons
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  names <- apply(pairs,2,function(x){paste(x[1],'vs',x[2])})
  
  ks.tab <- BiocParallel::bplapply(data.df, .ks_pairwise_table, sample_names, 
                                    outcomes, pairs, verbose,
                                    BPPARAM = bpparam)
  
  matrix(unlist(ks.tab), nrow=ncol(data.df), ncol=ncol(pairs)*2, byrow=TRUE)
}

# table of pairwise K-S results
.ks_pairwise_table <- function(geneData, sample_names, outcomes, 
                               pairs, verbose) {
  names(geneData) <- sample_names
  
  KS <- matrix(NA, nrow=1, ncol=ncol(pairs)*2)

  for (p in 1:ncol(pairs))
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    KS <- .ks_gene_pairwise(geneData,src.lab,sink.lab)
    KS[1,p] <- KS$p.value
    KS[1,ncol(pairs)+p] <- unname(KS$statistic)
  }
  
  KS <- as.numeric(KS)
  
  KS
}

# pairwise K-S for a single gene
.ks_gene_pairwise <- function(vec, idxA, idxB) {
  dataA <- vec[idxA]
  dataB <- vec[idxB]
  
  dataA <- as.numeric(dataA)
  dataB <- as.numeric(dataB)
  
  ks.test(dataA,dataB)
}

