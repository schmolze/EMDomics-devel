#' @export
#' @title Calculate the Komolgorov-Smirnov test statistic and q-values for differential gene expression
#' analysis.
#' 
#' @description This is only function needed when conducting an analysis using the Komolgorov-Smirnov
#' algorithm. Analyses can also be conducted with the EMD algorithm using
#' \code{calculate_emd} or the Cramer Von Mises (CVM) algorithm using \code{calculate_cvm}.
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
#' summarizes only the maximum difference (while EMD considers quantity and distance between all
#' differences).
#'
#' The KS algorithm implemented in \pkg{EMDomics} has two main steps.
#' First, a matrix (e.g. of expression data) is divided into data for each of the groups.
#' Every possible pairwise KS score is then computed and stored in a table. The KS score
#' for a single gene is calculated by averaging all of the pairwise KS scores. If the user
#' sets \code{pairwise.p} to true, then the p-values
#' from the KS test are adjusted using the Benjamini-Hochberg method and stored in a table.
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
#' @param pairwise.p Boolean specifying whether the user wants the pairwise p-values. Pairwise
#' p-values returned by \code{\link{ks.test}} are adjusted within pairwise comparison using the
#' Benjamini-Hochberg (BH) method. Defaults to \code{FALSE}.
#' @param seq Boolean specifying if the given data is RNA Sequencing data and ought to be
#' normalized. Set to \code{TRUE}, if passing transcripts per million (TPM) data or raw
#' data that is not scaled. If \code{TRUE}, data will be normalized by first multiplying by 1E6, then adding
#' 1, then taking the log base 2. If \code{FALSE}, the data will be handled as is (unless 
#' \code{quantile.norm} is \code{TRUE}). Note that as a distribution comparison function, K-S will
#' compute faster with scaled data. Defaults to \code{FALSE}.
#' @param quantile.norm Boolean specifying is data should be normalized by quantiles. If
#' \code{TRUE}, then the \code{\link[preprocessCore]{normalize.quantiles}} function is used.
#' Defaults to \code{FALSE}.
#' @param verbose Boolean specifying whether to display progress messages.
#' @param parallel Boolean specifying whether to use parallel processing via
#' the \pkg{BiocParallel} package. Defaults to \code{TRUE}.
#' @return The function returns an \code{\link{KSomics}} object.
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
#' results <- calculate_ks(dat, outcomes, nperm=10, parallel=FALSE)
#' head(results$ks)
#' 
#' @seealso \code{\link{EMDomics}} \code{\link{ks.test}}
calculate_ks <- function(data, outcomes, nperm=100, 
                         pairwise.p=FALSE, seq=FALSE,
                         quantile.norm=FALSE,
                         verbose=TRUE, parallel=TRUE) {
  
  bpparam <- BiocParallel::bpparam()
  
  if (parallel == FALSE)
    bpparam <- BiocParallel::SerialParam()
  
  
  if (seq)
  {
    data<-data*1E6
    data<-log2(data+1)
  }
  
  if (quantile.norm)
  {
    preprocessCore::normalize.quantiles(data)
  }
  
  # transpose and coerce to df (for bplapply)
  data.df <- as.data.frame(t(data))
  sample_names <- rownames(data.df)
  
  # ---------- pairwise ks table: test statistic and p-value -----------
  
  # generate pairwise ks table for each gene
  if (verbose)
    message("Calculating pairwise KS scores...", appendLF=FALSE)
  
  # all possible pairwise comparisons
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  names <- apply(pairs,2,function(x){paste(x[1],'vs',x[2])})
  
  ks.tab <- BiocParallel::bplapply(data.df, .ks_pairwise_table, sample_names, 
                                    outcomes, pairs, pairwise.p, verbose,
                                    BPPARAM = bpparam)
  
  if (pairwise.p) {
    ks.tab <- matrix(unlist(ks.tab), nrow=ncol(data.df), ncol=ncol(pairs)*2, byrow=TRUE)
  
    ks.p <- ks.tab[,seq(1,ncol(pairs)*2,2)]
    ks.stat <- ks.tab[,seq(2,ncol(pairs)*2,2)]
    
    # adjust p-values using BH/FDR by each pairwise comparison
    if (verbose)
      message("Adjusting p-values with B-H method by pairwise comparison...", appendLF=FALSE)
    
    ks.p <- apply(ks.p,2,function(x){ p.adjust(x, method='fdr') })
    rownames(ks.p) <- colnames(data.df)
    colnames(ks.p) <- names
    
    if (verbose)
      message("done.")
  } else {
    ks.tab <- matrix(unlist(ks.tab), nrow=ncol(data.df), ncol=ncol(pairs), byrow=TRUE)
    ks.stat <- ks.tab
    ks.p <- NULL
  }
  
  if (verbose)
    message("done.")
  
  rownames(ks.stat) <- colnames(data.df)
  colnames(ks.stat) <- names
  
  # ------------------ aggregate ks score for each gene (mean) -----------------------
  if (verbose)
    message("Calculating KS...", appendLF=FALSE)
  
  ks <- apply(ks.stat, 1, function(x){mean(x)})
  
  ks <- as.matrix(ks)
  rownames(ks) <- colnames(data.df)
  colnames(ks) <- 'ks'
  
  if (verbose)
    message("done.")
  
  # ------------------ permuted ks value ------------------------
  sample_count <- length(outcomes)
  
  # matrix to hold permuted ks values
  ks.perm <- matrix(nrow=ncol(data.df), ncol=nperm)
  rownames(ks.perm) <- colnames(data.df)
  colnames(ks.perm) <- as.character(1:nperm)
  
  for (i in 1:nperm) {
    
    msg <- paste("Calculating permuted KS #", i, " of ",
                 nperm, "...", sep="")
    
    if (verbose)
      message(msg, appendLF=FALSE)
    
    # permute samples
    idx.perm <- sample(1:sample_count, replace=FALSE)
    sample.id <- names(outcomes)
    outcomes.perm <- outcomes[idx.perm]
    names(outcomes.perm) <- sample.id
    
    # calculate ks for permuted samples
    perm.val <- BiocParallel::bplapply(data.df, calculate_ks_gene,
                                       outcomes.perm, rownames(data.df),
                                       BPPARAM = bpparam)
    
    ks.perm[,i] <- as.numeric(unlist(sapply(perm.val,"[",1)))
    
    if (verbose)
      message("done.")
  }
  
  # -------------- calculate permutation-based q-values ------------------
  if (verbose)
    message("Calculating q-values...", appendLF=FALSE)
  
  perm.medians <- apply(ks.perm,1,function(x){median(x)})
  
  # generate thresholds and qval matrix
  thr_upper <- ceiling(max(ks))
  thr <- seq(thr_upper, 0, by = -0.001)
  qvals <- matrix(1, nrow=nrow(ks), ncol=length(thr))
  
  colnames(qvals) <- thr
  rownames(qvals) <- rownames(ks)
  
  # calculate fdr at each threshold
  j <- 0
  for (d in thr) {
    
    j <- j+1
    
    # calculate true discoveries at this threshold
    idx <- which(ks > d)
    n.signif <- length(idx)
    genes.signif <- rownames(ks)[idx]
    
    # calculate false discoveries at this threshold
    idx <- which(perm.medians > d)
    n.fd <- length(idx)
    
    fdr <- n.fd/n.signif
    qvals[genes.signif, j] <- fdr
    
  }
  
  # final q-value = smallest fdr
  ks.qval <- apply(qvals, 1, min)
  
  ks.qval <- as.matrix(ks.qval)
  colnames(ks.qval) <- "permutation q-value"
  
  # leave q-values of 0 as-is
  
  if (verbose)
    message("done.")
  
  ks <- cbind(ks, ks.qval)
  
  KSomics(data, outcomes, ks, ks.perm, ks.stat, ks.p)
}

#' @export
#' @title Calculate KS score for a single gene
#' @details All possible combinations of the classes are used as pairwise comparisons.
#' The data in \code{vec} is divided based on class labels based on the \code{outcomes}
#' identifiers given. For each pairwise computation, \code{\link{ks.test}} is used to compute 
#' the pairwise KS scores. The 
#' total KS score for the given data is the average of the pairwise KS scores.
#' 
#' @param vec A named vector containing data (e.g. expression data) for a single
#' gene.
#' @param outcomes A vector of group labels for the samples. The names must correspond
#' to the names of \code{vec}.
#' @param sample_names A character vector with the names of the samples in \code{vec}.
#' @return The KS score is returned.
#' 
#' @examples
#' # 100 genes, 100 samples
#' dat <- matrix(rnorm(100000), nrow=100, ncol=1000)
#' rownames(dat) <- paste("gene", 1:100, sep="")
#' colnames(dat) <- paste("sample", 1:1000, sep="")
#'
#' # assign outcomes
#' outcomes <- c(rep(1,500), rep(2,300), rep(3,200))
#' names(outcomes) <- colnames(dat)
#'
#' calculate_ks_gene(dat[1,], outcomes, colnames(dat))
#' 
#' @seealso \code{\link{ks.test}}
calculate_ks_gene <- function(vec, outcomes, sample_names) {
  
  names(vec) <- sample_names
  
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  
  KS.tab <- matrix(NA, nrow=1, ncol=ncol(pairs))
  
  for (p in 1:ncol(pairs))
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]

    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    KS <- ks.test(vec[src.lab],vec[sink.lab])
    KS.tab[1,p] <- unname(KS$statistic)
  }
  
  KS.tab <- as.numeric(KS.tab)
  
  mean(KS.tab,na.rm=T)
}

#' @export
#' @title Create an KSomics object
#' @description This is the constructor for objects of class 'KSomics'. It
#' is used in \code{\link{calculate_ks}} to construct the return value.
#' 
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param outcomes A vector of group labels for each of the sample identifiers. The
#' names of this vector must correspond to the column names of \code{data}.
#' @param ks A matrix containing a row for each gene in \code{data}, and with
#' the following columns:
#' \itemize{
#' \item \code{ks} The calculated KS score.
#' \item \code{q-value} The calculated q-value (by permutation analysis).
#' }
#' The row names should specify the gene identifiers for each row.
#' @param ks.perm A matrix containing a row for each gene in \code{data}, and
#' with a column containing KS scores for each random permutation calculated
#' via \code{\link{calculate_ks}}.
#' @param pairwise.ks.score A table containing the KS scores for each pairwise
#' comparison for each gene. For a two-class problem, there should be only one column
#' comparing class 1 and class 2. The row names should be gene identifiers. The column
#' names should be in the format "<class 1> vs <class 2>" (e.g. "1 vs 2" or "A vs B").
#' @param pairwise.ks.q A table of the same dimensions as \code{pairwise.ks.score} with
#' the q-values for the pairwise comparisons. Q-values are computed by adjusting the p-value
#' using the Benjamini-Hochberg method within each pairwise comparison.
#' 
#' @return The function combines its arguments in a list, which is assigned class
#' 'KSomics'. The resulting object is returned.
#' 
#' @seealso \code{\link{calculate_ks}}
KSomics <- function(data, outcomes, ks, ks.perm, 
                        pairwise.ks.score, pairwise.ks.q=NULL) {
  structure(list("data"=data, "outcomes"=outcomes,
                 "ks"=ks, "ks.perm"=ks.perm, "pairwise.ks.score"=pairwise.ks.score, 
                 "pairwise.ks.q"=pairwise.ks.q),
            class = "KSomics")
}

# table of pairwise K-S results
.ks_pairwise_table <- function(geneData, sample_names, outcomes, 
                               pairs, pairwise.p, verbose) {
  names(geneData) <- sample_names
  
  if (pairwise.p) {
    KS.tab <- matrix(NA, nrow=1, ncol=ncol(pairs)*2)
  } else {
    KS.tab <- matrix(NA, nrow=1, ncol=ncol(pairs))
  }
  
  for (p in 1:ncol(pairs))
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    KS <- ks.test(geneData[src.lab],geneData[sink.lab])
    
    if (pairwise.p) {
      KS.tab[1,p] <- KS$p.value
      KS.tab[1,ncol(pairs)+p] <- unname(KS$statistic)
    } else {
      KS.tab[1,p] <- KS$statistic
    }
  }
  
  KS.tab <- as.numeric(KS.tab)
  
  KS.tab
}
