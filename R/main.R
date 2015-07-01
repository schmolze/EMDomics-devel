#' Earth Mover's Distance algorithm for differential analysis of genomics data.
#'
#' \code{\link{calculate_emd}} will usually be the only function needed.
#'
#'
#' @import emdist
#' @import BiocParallel
#' @import matrixStats
#' @import ggplot2
#' @name emdomics-package
#' @docType package
NULL


#' @export
#' @title Earth Mover's Distance for differential analysis of genomics data
#' @description This is the main user interface to the \pkg{EMDomics} package, and
#' will usually the only function needed.
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
#' The Earth Mover's Distance algorithm instead computes the "work" needed
#' to transform one distribution into another, thus capturing possibly
#' valuable information relating to the overall difference in shape between
#' two heterogeneous distributions.
#'
#' The EMD-based algorithm implemented in \pkg{EMDomics} has two main steps.
#' First, a matrix (e.g. of expression data) is divided into data for each of the groups.
#' Every possible pairwise EMD score is then computed and stored in a table. The EMD score
#' for a single gene is calculated by summing all of the pairwise EMD scores.
#' Next, the labels for each of the groups are randomly
#' permuted a specified number of times, and an EMD score for each permutation is
#' calculated. The median of the permuted scores for each gene is used as
#' the null distribution, and the False Discovery Rate (FDR) is computed for
#' a range of permissive to restrictive significance thresholds. The threshold
#' that minimizes the FDR is defined as the q-value, and is used to interpret
#' the significance of the EMD score analogously to a p-value (e.g. q-value
#' < 0.05 = significant.)
#'
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param outcomes A vector containing group labels for each of the samples provided
#' in the \code{data} matrix. The names should be the sample identifiers provided in \code{data}.
#' @param binSize The bin size to be used when generating histograms of
#' the data for each group. Defaults to 0.2.
#' @param nperm An integer specifying the number of randomly permuted EMD
#' scores to be computed. Defaults to 100.
#' @param verbose Boolean specifying whether to display progress messages.
#' @param parallel Boolean specifying whether to use parallel processing via
#' the \pkg{BiocParallel} package. Defaults to \code{TRUE}.
#' @return The function returns an \code{\link{EMDomics}} object.
#' @examples
#' # 100 genes, 100 samples
#' dat <- matrix(rnorm(10000), nrow=100, ncol=100)
#' rownames(dat) <- paste("gene", 1:100, sep="")
#' colnames(dat) <- paste("sample", 1:100, sep="")
#'
#' # "A": first 50 samples; "B": next 30 samples; "C": final 20 samples
#' outcomes <- c(rep("A",50), rep("B",30), rep("C",20))
#' names(outcomes) <- colnames(dat)
#' 
#' results <- calculate_emd(dat, outcomes, nperm=10, parallel=FALSE)
#' head(results$emd)
#' 
#' @seealso \code{\link{EMDomics}} \code{\link[emdist]{emd2d}}
calculate_emd <- function(data, outcomes, binSize=0.2,
                            nperm=100, verbose=TRUE, parallel=TRUE) {

  bpparam <- BiocParallel::bpparam()

  if (parallel == FALSE)
    bpparam <- BiocParallel::SerialParam()

  # transpose and coerce to df (for bplapply)
  data.df <- as.data.frame(t(data))
  sample_names <- rownames(data.df)
  
  # ---------- pairwise emd table -----------
  
  # generate pairwise emd table for each gene
  if (verbose)
    message("Calculating pairwise emd scores...", appendLF=FALSE)
  
  # all possible pairwise comparisons
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  names <- apply(pairs,2,function(x){paste(x[1],'vs',x[2])})
  
  emd.tab <- BiocParallel::bplapply(data.df, .emd_pairwise_table, sample_names, 
                                           outcomes, pairs,
                                           binSize, verbose,
                                           BPPARAM = bpparam)
  
  emd.tab <- matrix(unlist(emd.tab), nrow=dim(data.df)[2], byrow=TRUE)
  rownames(emd.tab) <- colnames(data.df)
  colnames(emd.tab) <- names
  
  if (verbose)
    message("done.")
  
  # ---------- emd ------------

  # calculate emd for each gene
  if (verbose)
    message("Calculating emd...", appendLF=FALSE)
  
  emd <- apply(emd.tab, 1, function(x){mean(as.numeric(x))})

  emd <- as.matrix(emd)
  colnames(emd) <- "emd"

  if (verbose)
    message("done.")

  # --------------- permuted emd scores ----------------

  sample_count <- length(outcomes)

  # matrix to hold permuted emd values
  emd.perm <- matrix(nrow=ncol(data.df), ncol=nperm)
  rownames(emd.perm) <- colnames(data.df)
  colnames(emd.perm) <- as.character(1:nperm)

  for (i in 1:nperm) {

    msg <- paste("Calculating permuted emd #", i, " of ",
                 nperm, "...", sep="")

    if (verbose)
      message(msg, appendLF=FALSE)

    # permute samples
    idx.perm <- sample(1:sample_count, replace=FALSE)
    data.perm <- data.df[idx.perm, ]

    # calculate emd for permuted samples
    perm.val <- BiocParallel::bplapply(data.perm, calculate_emd_gene,
                                                   outcomes, rownames(data.perm),
                                                   binSize,
                                                   BPPARAM = bpparam)
  
    emd.perm[,i] <- unlist(perm.val)
    
    if (verbose)
      message("done.")

  }

  # ------------------ q-values --------------------

  if (verbose)
    message("Calculating q-values...", appendLF=FALSE)

  perm.medians <- apply(emd.perm,1,function(x){median(x)})

  # generate thresholds and qval matrix
  thr_upper <- ceiling(max(emd))
  thr <- seq(thr_upper, 0, by = -0.001)
  qvals <- matrix(1, nrow=nrow(emd), ncol=length(thr))

  colnames(qvals) <- thr
  rownames(qvals) <- rownames(emd)

  # calculate fdr at each threshold
  j <- 0
  for (d in thr) {

    j <- j+1

    # calculate true discoveries at this threshold
    idx <- which(emd > d)
    n.signif <- length(idx)
    genes.signif <- rownames(emd)[idx]

    # calculate false discoveries at this threshold
    idx <- which(perm.medians > d)
    n.fd <- length(idx)

    fdr <- n.fd/n.signif
    qvals[genes.signif, j] <- fdr

  }

  # final q-value = smallest fdr
  emd.qval <- apply(qvals, 1, min)

  emd.qval <- as.matrix(emd.qval)
  colnames(emd.qval) <- "q-value"

  # leave q-values of 0 as-is

  if (verbose)
    message("done.")

  emd <- cbind(emd, emd.qval)

  EMDomics(data, outcomes, emd, emd.perm, emd.tab)

}


#' @export
#' @title Calculate EMD score for a single gene
#' @details All possible combinations of the classes are used as pairwise comparisons.
#' The data in \code{vec} is divided based on class labels based on the \code{outcomes}
#' identifiers given. For each pairwise computation, the \code{\link{hist}} function is
#' used to generate histograms for the two groups. The densities are then retrieved
#' and passed to  \code{\link[emdist]{emd2d}} to compute the pairwise EMD score. The 
#' total EMD score for the given data is the sum of the pairwise EMD scores.
#' 
#' @param vec A named vector containing data (e.g. expression data) for a single
#' gene.
#' @param outcomes A vector of group labels for the samples. The names must correspond
#' to the names of \code{vec}.
#' @param sample_names A character vector with the names of the samples in \code{vec}.
#' @param binSize The bin size to be used when generating histograms for each of the groups.
#' @return The emd score is returned.
#' 
#' @examples
#' # 100 genes, 100 samples
#' dat <- matrix(rnorm(10000), nrow=100, ncol=100)
#' rownames(dat) <- paste("gene", 1:100, sep="")
#' colnames(dat) <- paste("sample", 1:100, sep="")
#'
#' # "A": first 50 samples; "B": next 30 samples; "C": final 20 samples
#' outcomes <- c(rep("A",50), rep("B",30), rep("C",20))
#' names(outcomes) <- colnames(dat)
#'
#' calculate_emd_gene(dat[1,], outcomes, colnames(dat))
#' 
#' @seealso \code{\link[emdist]{emd2d}}
calculate_emd_gene <- function(vec, outcomes, sample_names, binSize=0.2) {
  
  names(vec) <- sample_names
  
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  
  EMD.tab <- matrix(NA, nrow=1, ncol=dim(pairs)[2])
  colnames<-list()
  
  for (p in 1:dim(pairs)[2])
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    EMD <- .emd_gene_pairwise(vec,src.lab,sink.lab,binSize)
    EMD.tab[1,p] <- EMD
  }
  
  EMD.tab <- as.numeric(EMD.tab)
  names(EMD.tab) <- colnames
  
  mean(EMD.tab)
}


#' @export
#' @title Create an EMDomics object
#' @description This is the constructor for objects of class 'EMDomics'. It
#' is used in \code{\link{calculate_emd}} to construct the return value.
#' 
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param outcomes A vector of group labels for each of the sample identifiers. The
#' names of this vector must correspond to the column names of \code{data}.
#' @param emd A matrix containing a row for each gene in \code{data}, and with
#' the following columns:
#' \itemize{
#' \item \code{emd} The calculated emd score.
#' \item \code{q-value} The calculated q-value.
#' }
#' The row names should specify the gene identifiers for each row.
#' @param emd.perm A matrix containing a row for each gene in \code{data}, and
#' with a column containing emd scores for each random permutation calculated
#' via \code{\link{calculate_emd}}.
#' @param pairwise.emd.table A table containing the EMD scores for each pairwise
#' comparison for each gene. For a two-class problem, there should be only one column
#' comparing class 1 and class 2. The row names should be gene identifiers. The column
#' names should be in the format "<class 1> vs <class 2>" (e.g. "1 vs 2" or "A vs B").
#' 
#' @return The function combines its arguments in a list, which is assigned class
#' 'EMDomics'. The resulting object is returned.
#' 
#' @seealso \code{\link{calculate_emd}}
EMDomics <- function(data, outcomes, emd, emd.perm, pairwise.emd.table) {

  structure(list("data"=data, "outcomes"=outcomes,
                 "emd"=emd, "emd.perm"=emd.perm, "pairwise.emd.table"=pairwise.emd.table),
            class = "EMDomics")

}

# Creates a table of all the pairwise EMD scores for one gene
.emd_pairwise_table <- function(geneData, sample_names, outcomes, pairs, binSize, verbose) {
  
  names(geneData) <- sample_names
  
  EMD.tab <- matrix(NA, nrow=1, ncol=dim(pairs)[2])
  colnames<-list()
  
  for (p in 1:dim(pairs)[2])
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    col.name <- as.character(paste(src,'vs',sink))
    colnames<-c(colnames,col.name)
    
    EMD <- .emd_gene_pairwise(geneData,src.lab,sink.lab,binSize)
    EMD.tab[1,p] <- EMD
  }

  EMD.tab <- as.numeric(EMD.tab)
  names(EMD.tab) <- colnames
  
  EMD.tab
}

# computes pairwise EMD
.emd_gene_pairwise <- function(vec, idxA, idxB, binSize=0.2) {
  dataA <- vec[idxA]
  dataB <- vec[idxB]
  
  dataA <- as.numeric(dataA)
  dataB <- as.numeric(dataB)
  
  bins <- seq(floor(min(c(dataA, dataB))),
              ceiling(max(c(dataA, dataB))),
              by=binSize )
  
  histA <- hist(dataA, breaks=bins, plot=FALSE)
  histB <- hist(dataB, breaks=bins, plot=FALSE)
  
  densA <- as.matrix(histA$density)
  densB <- as.matrix(histB$density)
  
  emdist::emd2d(densA, densB)
}
