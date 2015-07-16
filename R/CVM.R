#' @export
#' @title Earth Mover's Distance for differential analysis of genomics data
#' @description This is the main user interface to the \pkg{EMDomics} package, and
#' will usually the only function needed when conducting an analysis using the CVM
#' algorithm. Analyses can also be conducted with the Komolgorov-Smirnov Test using
#' \code{calculate_ks} or the Cramer Von Mises algorithm using \code{calculate_cvm}.
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
#' The Cramer von Mises (CVM) algorithm generates a test statistic that is the sum of 
#' the squared values of the differences between two cumulative distribution functions (CDFs). 
#' As a result, the test statistic tends to overestimate the similarity between two distributions and
#' cannot effectively handle partial matching like EMD does.
#' 
#' The CVM-based algorithm implemented in \pkg{EMDomics} has two main steps.
#' First, a matrix (e.g. of expression data) is divided into data for each of the groups.
#' Every possible pairwise CVM score is then computed and stored in a table. The CVM score
#' for a single gene is calculated by averaging all of the pairwise CVM scores.
#' Next, the labels for each of the groups are randomly
#' permuted a specified number of times, and an CVM score for each permutation is
#' calculated. The median of the permuted scores for each gene is used as
#' the null distribution, and the False Discovery Rate (FDR) is computed for
#' a range of permissive to restrictive significance thresholds. The threshold
#' that minimizes the FDR is defined as the q-value, and is used to interpret
#' the significance of the CVM score analogously to a p-value (e.g. q-value
#' < 0.05 is significant.)
#'
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param outcomes A vector containing group labels for each of the samples provided
#' in the \code{data} matrix. The names should be the sample identifiers provided in \code{data}.
#' @param nperm An integer specifying the number of randomly permuted CVM
#' scores to be computed. Defaults to 100.
#' @param pairwise.p Boolean specifying whether the permutation-based q-values should
#' be computed for each pairwise comparison. Defaults to \code{FALSE}.
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
#' @return The function returns an \code{\link{CVMomics}} object.
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
#' results <- calculate_cvm(dat, outcomes, nperm=10, parallel=FALSE)
#' head(results$cvm)
#' 
#' @seealso \code{\link{CVMomics}} \code{\link[CDFt]{CramerVonMisesTwoSamples}}
calculate_cvm <- function(data, outcomes,
                          nperm=100, pairwise.p=FALSE, seq=FALSE, 
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
  
  # ---------- pairwise cvm table -----------
  
  # generate pairwise cvm table for each gene
  if (verbose)
    message("Calculating pairwise cvm scores...", appendLF=FALSE)
  
  # all possible pairwise comparisons
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  names <- apply(pairs,2,function(x){paste(x[1],'vs',x[2])})
  
  cvm.tab <- BiocParallel::bplapply(data.df, .cvm_pairwise_table, sample_names, 
                                    outcomes, pairs, verbose,
                                    BPPARAM = bpparam)
  
  cvm.tab <- matrix(unlist(cvm.tab), nrow=ncol(data.df), ncol=ncol(pairs), byrow=TRUE)
  rownames(cvm.tab) <- colnames(data.df)
  colnames(cvm.tab) <- names
  
  if (verbose)
    message("done.")
  
  # ---------- cvm ------------
  
  # calculate cvm for each gene
  if (verbose)
    message("Calculating cvm...", appendLF=FALSE)
  
  cvm <- apply(cvm.tab, 1, function(x){mean(as.numeric(x))})
  
  cvm <- as.matrix(cvm)
  colnames(cvm) <- "cvm"
  
  if (verbose)
    message("done.")
  
  # pairwise q-value computation if specified
  if (pairwise.p) {
    cvm.pairwise.q <- .cvm_pairwise_q(data.df, cvm, outcomes, nperm, verbose, bpparam)
  } else {
    cvm.pairwise.q <- NULL
  }
  
  # --------------- permuted cvm scores ----------------
  message("Calculating for overall q-values...")
  
  sample_count <- length(outcomes)
  
  # matrix to hold permuted cvm values
  cvm.perm <- matrix(nrow=ncol(data.df), ncol=nperm)
  rownames(cvm.perm) <- colnames(data.df)
  colnames(cvm.perm) <- as.character(1:nperm)
  
  for (i in 1:nperm) {
    
    msg <- paste("Calculating permuted cvm #", i, " of ",
                 nperm, "...", sep="")
    
    if (verbose)
      message(msg, appendLF=FALSE)
    
    # permute samples
    idx.perm <- sample(1:sample_count, replace=FALSE)
    sample.id <- names(outcomes)
    outcomes.perm <- outcomes[idx.perm]
    names(outcomes.perm) <- sample.id
    
    # calculate cvm for permuted samples
    perm.val <- BiocParallel::bplapply(data.df, calculate_cvm_gene,
                                       outcomes.perm, rownames(data.df),
                                       BPPARAM = bpparam)
    
    cvm.perm[,i] <- unlist(sapply(perm.val,"[",1))
    
    if (verbose)
      message("done.")
    
  }
  
  # ------------------ q-values --------------------
  
  if (verbose)
    message("Calculating q-values...", appendLF=FALSE)
  
  perm.medians <- apply(cvm.perm,1,function(x){median(x)})
  
  # generate thresholds and qval matrix
  thr_upper <- ceiling(max(cvm))
  thr <- seq(thr_upper, 0, by = -0.001)
  qvals <- matrix(1, nrow=nrow(cvm), ncol=length(thr))
  
  colnames(qvals) <- thr
  rownames(qvals) <- rownames(cvm)
  
  # calculate fdr at each threshold
  j <- 0
  for (d in thr) {
    
    j <- j+1
    
    # calculate true discoveries at this threshold
    idx <- which(cvm > d)
    n.signif <- length(idx)
    genes.signif <- rownames(cvm)[idx]
    
    # calculate false discoveries at this threshold
    idx <- which(perm.medians > d)
    n.fd <- length(idx)
    
    fdr <- n.fd/n.signif
    qvals[genes.signif, j] <- fdr
    
  }
  
  # final q-value = smallest fdr
  cvm.qval <- apply(qvals, 1, min)
  
  cvm.qval <- as.matrix(cvm.qval)
  colnames(cvm.qval) <- "q-value"
  
  # leave q-values of 0 as-is
  
  if (verbose)
    message("done.")
  
  cvm <- cbind(cvm, cvm.qval)
  
  CVMomics(data, outcomes, cvm, cvm.perm, cvm.tab, cvm.pairwise.q)
  
}


#' @export
#' @title Calculate CVM score for a single gene
#' @details All possible combinations of the classes are used as pairwise comparisons.
#' The data in \code{vec} is divided based on class labels based on the \code{outcomes}
#' identifiers given. For each pairwise computation, the \code{\link{hist}} function is
#' used to generate histograms for the two groups. The densities are then retrieved
#' and passed to  \code{\link[CDFt]{CramerVonMisesTwoSamples}} to compute the pairwise CVM score. The 
#' total CVM score for the given data is the average of the pairwise CVM scores.
#' 
#' @param vec A named vector containing data (e.g. expression data) for a single
#' gene.
#' @param outcomes A vector of group labels for the samples. The names must correspond
#' to the names of \code{vec}.
#' @param sample_names A character vector with the names of the samples in \code{vec}.
#' @return The cvm score is returned.
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
#' calculate_cvm_gene(dat[1,], outcomes, colnames(dat))
#' 
#' @seealso \code{\link[CDFt]{CramerVonMisesTwoSamples}}
calculate_cvm_gene <- function(vec, outcomes, sample_names) {
  
  names(vec) <- sample_names
  
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  
  CVM.tab <- matrix(NA, nrow=1, ncol=dim(pairs)[2])
  colnames<-list()
  
  for (p in 1:dim(pairs)[2])
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    CVM <- CDFt::CramerVonMisesTwoSamples(vec[src.lab],vec[sink.lab])
    CVM.tab[1,p] <- CVM
  }
  
  CVM.tab <- as.numeric(CVM.tab)
  
  mean(CVM.tab)
}


#' @export
#' @title Create an CVMomics object
#' @description This is the constructor for objects of class 'CVMomics'. It
#' is used in \code{\link{calculate_cvm}} to construct the return value.
#' 
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param outcomes A vector of group labels for each of the sample identifiers. The
#' names of this vector must correspond to the column names of \code{data}.
#' @param cvm A matrix containing a row for each gene in \code{data}, and with
#' the following columns:
#' \itemize{
#' \item \code{cvm} The calculated cvm score.
#' \item \code{q-value} The calculated q-value.
#' }
#' The row names should specify the gene identifiers for each row.
#' @param cvm.perm A matrix containing a row for each gene in \code{data}, and
#' with a column containing cvm scores for each random permutation calculated
#' via \code{\link{calculate_cvm}}.
#' @param pairwise.cvm.table A table containing the CVM scores for each pairwise
#' comparison for each gene. For a two-class problem, there should be only one column
#' comparing class 1 and class 2. The row names should be gene identifiers. The column
#' names should be in the format "<class 1> vs <class 2>" (e.g. "1 vs 2" or "A vs B").
#' 
#' @return The function combines its arguments in a list, which is assigned class
#' 'CVMomics'. The resulting object is returned.
#' 
#' @seealso \code{\link{calculate_cvm}}
CVMomics <- function(data, outcomes, cvm, cvm.perm, pairwise.cvm.table, pairwise.q.table) {
  
  structure(list("data"=data, "outcomes"=outcomes,
                 "cvm"=cvm, "cvm.perm"=cvm.perm, "pairwise.cvm.table"=pairwise.cvm.table,
                 "pairwise.q.table"=pairwise.q.table),
            class = "CVMomics")
  
}

# Creates a table of all the pairwise q-values for all genes
.cvm_pairwise_q <- function(data.df, cvm, outcomes, nperm, verbose, bpparam) {
  
  classes <- unique(outcomes)
  pairs <- combn(classes,2)
  names <- apply(pairs,2,function(x){paste(x[1],'vs',x[2])})
  
  q.tab <- matrix(NA, nrow=ncol(data.df), ncol=ncol(pairs))
  colnames(q.tab) <- names
  rownames(q.tab) <- colnames(data.df)
  
  for (p in 1:ncol(pairs)) {
    
    # ------------------ calculate permuted cvm scores ---------------------
    sample_count <- length(outcomes)
    
    # matrix to hold permuted cvm values
    cvm.perm <- matrix(nrow=ncol(data.df), ncol=nperm)
    rownames(cvm.perm) <- colnames(data.df)
    colnames(cvm.perm) <- as.character(1:nperm)
    
    msg <- paste("Beginning pairwise q-value computation for",names[p])
    if (verbose)
      message(msg)
    
    for (i in 1:nperm) {
      
      msg <- paste("Calculating permuted cvm #", i, " of ",
                   nperm, "...", sep="")
      
      if (verbose)
        message(msg, appendLF=FALSE)
      
      # permute samples
      idx.perm <- sample(1:sample_count, replace=FALSE)
      sample.id <- names(outcomes)
      outcomes.perm <- outcomes[idx.perm]
      names(outcomes.perm) <- sample.id
      
      # calculate cvm for permuted samples
      perm.val <- BiocParallel::bplapply(data.df, calculate_cvm_gene,
                                         outcomes.perm, rownames(data.df),
                                         BPPARAM = bpparam)
      
      cvm.perm[,i] <- unlist(sapply(perm.val,"[",1))
      
      if (verbose)
        message("done.")
    }
    
    # ------------------ q-values --------------------
    
    if (verbose)
      message("Calculating pairwise q-values...", appendLF=FALSE)
    
    perm.medians <- apply(cvm.perm,1,function(x){median(x)})
    
    # generate thresholds and qval matrix
    thr_upper <- ceiling(max(cvm))
    thr <- seq(thr_upper, 0, by = -0.001)
    qvals <- matrix(1, nrow=nrow(cvm), ncol=length(thr))
    
    colnames(qvals) <- thr
    rownames(qvals) <- rownames(cvm)
    
    # calculate fdr at each threshold
    j <- 0
    for (d in thr) {
      
      j <- j+1
      
      # calculate true discoveries at this threshold
      idx <- which(cvm > d)
      n.signif <- length(idx)
      genes.signif <- rownames(cvm)[idx]
      
      # calculate false discoveries at this threshold
      idx <- which(perm.medians > d)
      n.fd <- length(idx)
      
      fdr <- n.fd/n.signif
      qvals[genes.signif, j] <- fdr
      
    }
    
    # final q-value = smallest fdr
    cvm.qval <- apply(qvals, 1, min)
    
    q.tab[,p] <- cvm.qval
    
    if (verbose)
      message("done.")
  }
  
  q.tab
}

# Creates a table of all the pairwise CVM scores for one gene
.cvm_pairwise_table <- function(geneData, sample_names, outcomes, pairs, verbose) {
  
  names(geneData) <- sample_names
  
  CVM.tab <- matrix(NA, nrow=1, ncol=dim(pairs)[2])
  colnames<-list()
  
  for (p in 1:dim(pairs)[2])
  {
    inds <- pairs[,p]
    src <- inds[1]
    sink <- inds[2]
    
    src.lab <- names(outcomes[outcomes==src])
    sink.lab <- names(outcomes[outcomes==sink])
    
    CVM <- CDFt::CramerVonMisesTwoSamples(geneData[src.lab],geneData[sink.lab])
    CVM.tab[1,p] <- CVM
  }
  
  CVM.tab <- as.numeric(CVM.tab)
  
  CVM.tab
}
