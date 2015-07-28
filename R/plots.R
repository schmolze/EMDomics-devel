#' @export
#' @title Plot null distribution of permuted EMD scores vs. calculated EMD
#' scores.
#' @description The median of the randomly permuted EMD scores (i.e. the null
#' distribution) is plotted on the x-axis, vs. the observed EMD scores on the
#' y-axis. The line \code{y=x} is superimposed.
#' @param emdobj An \code{\link{EMDomics}} object, typically returned via a call
#' to \code{\link{calculate_emd}}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
#' @examples
#' # 100 genes, 100 samples
#' dat <- matrix(rnorm(10000), nrow=100, ncol=100)
#' rownames(dat) <- paste("gene", 1:100, sep="")
#' colnames(dat) <- paste("sample", 1:100, sep="")
#'
#' # "group A" = first 50, "group B" = second 50
#' groups <- c(rep("A",50),rep("B",50))
#' names(groups) <- colnames(dat)
#'
#' results <- calculate_emd(dat, groups, nperm=10, parallel=FALSE)
#' plot_emdnull(results)
#' @seealso \code{\link{calculate_emd}} \code{\link[ggplot2]{ggplot}}
plot_emdnull <- function(emdobj) {

  emd <- emdobj$emd
  emd.perm <- emdobj$emd.perm
  rms <- rowMedians(emd.perm)

  data <- as.data.frame(cbind(emd[,"emd",drop=FALSE], rms))

  title <- "Null distribution vs. observed emd scores"

  ggplot(data, aes(rms, emd)) + geom_point(alpha=0.3) +
    geom_segment(x=0, y=0, xend=10, yend=10, colour="red") +
    xlab("median of permuted emd scores")  +
    ylab("observed emd scores") +
    ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24))

}

#' @export
#' @title Plot null distribution of permuted ks scores vs. calculated ks
#' scores.
#' @description The median of the randomly permuted KS scores (i.e. the null
#' distribution) is plotted on the x-axis, vs. the observed KS scores on the
#' y-axis. The line \code{y=x} is superimposed.
#' @param ksobj An \code{\link{KSomics}} object, typically returned via a call
#' to \code{\link{calculate_ks}}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
#' @examples
#' # 100 genes, 100 samples
#' dat <- matrix(rnorm(10000), nrow=100, ncol=100)
#' rownames(dat) <- paste("gene", 1:100, sep="")
#' colnames(dat) <- paste("sample", 1:100, sep="")
#'
#' # "group A" = first 50, "group B" = second 50
#' groups <- c(rep("A",50),rep("B",50))
#' names(groups) <- colnames(dat)
#'
#' results <- calculate_ks(dat, groups, nperm=10, parallel=FALSE)
#' plot_ksnull(results)
#' @seealso \code{\link{calculate_ks}} \code{\link[ggplot2]{ggplot}}
plot_ksnull <- function(ksobj) {
  
  ks <- ksobj$ks
  ks.perm <- ksobj$ks.perm
  rms <- rowMedians(ks.perm)
  
  data <- as.data.frame(cbind(ks[,"ks",drop=FALSE], rms))
  
  title <- "Null distribution vs. observed KS scores"
  
  ggplot(data, aes(rms, ks)) + geom_point(alpha=0.3) +
    geom_segment(x=0, y=0, xend=10, yend=10, colour="red") +
    xlab("median of permuted ks scores")  +
    ylab("observed ks scores") +
    ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24))
  
}

#' @export
#' @title Plot null distribution of permuted cvm scores vs. calculated cvm
#' scores.
#' @description The median of the randomly permuted CVM scores (i.e. the null
#' distribution) is plotted on the x-axis, vs. the observed CVM scores on the
#' y-axis. The line \code{y=x} is superimposed.
#' @param cvmobj An \code{\link{CVMomics}} object, typically returned via a call
#' to \code{\link{calculate_cvm}}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
#' @examples
#' # 100 genes, 100 samples
#' dat <- matrix(rnorm(10000), nrow=100, ncol=100)
#' rownames(dat) <- paste("gene", 1:100, sep="")
#' colnames(dat) <- paste("sample", 1:100, sep="")
#'
#' # "group A" = first 50, "group B" = second 50
#' groups <- c(rep("A",50),rep("B",50))
#' names(groups) <- colnames(dat)
#'
#' results <- calculate_cvm(dat, groups, nperm=10, parallel=FALSE)
#' plot_cvmnull(results)
#' @seealso \code{\link{calculate_cvm}} \code{\link[ggplot2]{ggplot}}
plot_cvmnull <- function(cvmobj) {
  
  cvm <- cvmobj$cvm
  cvm.perm <- cvmobj$cvm.perm
  rms <- rowMedians(cvm.perm)
  
  data <- as.data.frame(cbind(cvm[,"cvm",drop=FALSE], rms))
  
  title <- "Null distribution vs. observed CVM scores"
  
  ggplot(data, aes(rms, cvm)) + geom_point(alpha=0.3) +
    geom_segment(x=0, y=0, xend=10, yend=10, colour="red") +
    xlab("median of permuted cvm scores")  +
    ylab("observed cvm scores") +
    ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24))
  
}

#' @export
#' @title Plot histogram of EMD scores calculated via random permutation.
#' @description The permuted EMD scores stored in \code{emdobj$emd.perm} are
#' plotted as a histogram.
#' @param emdobj An \code{\link{EMDomics}} object, typically returned via a call
#' to \code{\link{calculate_emd}}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
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
#' plot_emdperms(results)
#' @seealso \code{\link{calculate_emd}} \code{\link[ggplot2]{ggplot}}
plot_emdperms <- function(emdobj) {

  emd.perm <- as.data.frame(emdobj$emd.perm)

  # to appease CRAN
  x <- NULL

  colnames(emd.perm) <- "x"

  title <- "Histogram of permuted emd scores"

  ggplot(emd.perm, aes(x)) + geom_histogram(alpha=0.7) +
    xlab("emd score")  + ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24))

}

#' @export
#' @title Plot histogram of KS scores calculated via random permutation.
#' @description The permuted KS scores stored in \code{ksobj$ks.perm} are
#' plotted as a histogram.
#' @param ksobj An \code{\link{KSomics}} object, typically returned via a call
#' to \code{\link{calculate_ks}}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
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
#' plot_ksperms(results)
#' @seealso \code{\link{calculate_ks}} \code{\link[ggplot2]{ggplot}}
plot_ksperms <- function(ksobj) {
  
  ks.perm <- as.data.frame(ksobj$ks.perm)
  
  # to appease CRAN
  x <- NULL
  
  colnames(ks.perm) <- "x"
  
  title <- "Histogram of permuted KS scores"
  
  ggplot(ks.perm, aes(x)) + geom_histogram(alpha=0.7) +
    xlab("ks score")  + ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24))
  
}

#' @export
#' @title Plot histogram of CVM scores calculated via random permutation.
#' @description The permuted CVM scores stored in \code{cvmobj$cvm.perm} are
#' plotted as a histogram.
#' @param cvmobj An \code{\link{CVMomics}} object, typically returned via a call
#' to \code{\link{calculate_cvm}}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
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
#' plot_cvmperms(results)
#' @seealso \code{\link{calculate_cvm}} \code{\link[ggplot2]{ggplot}}
plot_cvmperms <- function(cvmobj) {
  
  cvm.perm <- as.data.frame(cvmobj$cvm.perm)
  
  # to appease CRAN
  x <- NULL
  
  colnames(cvm.perm) <- "x"
  
  title <- "Histogram of permuted CVM scores"
  
  ggplot(cvm.perm, aes(x)) + geom_histogram(alpha=0.7) +
    xlab("cvm score")  + ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24))
  
}


#' @export
#' @title Plot distributions and EMD score for a gene.
#' @description The data for the specified gene is retrieved from
#' \code{emdobj$data}. \code{outcomes} is used to divide the data into 
#' distributions for each group, which are then visualized as
#' density distributions. The calculated EMD score for the specified gene is
#' displayed in the plot title.
#' @param emdobj An \code{\link{EMDomics}} object, typically returned via a call
#' to \code{\link{calculate_emd}}.
#' @param gene_name The gene to visualize. The name should be defined as a row
#' name in \code{emdobj$emd}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
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
#' plot_emd_density(results, "gene5")
#' 
#' @seealso \code{\link{calculate_emd}} \code{\link[ggplot2]{ggplot}}
plot_emd_density <- function(emdobj, gene_name) {

  data <- emdobj$data
  outcomes <- emdobj$outcomes
  classes <- unique(outcomes)
  
  emd_score <- emdobj$emd[gene_name, "emd"]
  gene.data <- data[gene_name,]
  
  # to appease CRAN
  group <- NULL
  exp <- NULL
  
  df<-data.frame(row.names=colnames(data), group=outcomes, exp=gene.data)

  title <- paste(gene_name, "\n", "(emd score = ",
                 round(emd_score, 2), ")", sep="")

  ggplot(df, aes(exp, fill=group)) + geom_density(alpha=0.5) +
    xlab("data")  + ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24),
          legend.text = element_text(size = 24),
          legend.title = element_text(size=24))
}

#' @export
#' @title Plot distributions and KS score for a gene.
#' @description The data for the specified gene is retrieved from
#' \code{ksobj$data}. \code{outcomes} is used to divide the data into 
#' distributions for each group, which are then visualized as
#' density distributions. The calculated KS score for the specified gene is
#' displayed in the plot title.
#' @param ksobj An \code{\link{KSomics}} object, typically returned via a call
#' to \code{\link{calculate_ks}}.
#' @param gene_name The gene to visualize. The name should be defined as a row
#' name in \code{ksobj$ks}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
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
#' plot_ks_density(results, "gene5")
#' 
#' @seealso \code{\link{calculate_ks}} \code{\link[ggplot2]{ggplot}}
plot_ks_density <- function(ksobj, gene_name) {
  
  data <- ksobj$data
  outcomes <- ksobj$outcomes
  classes <- unique(outcomes)
  
  ks_score <- ksobj$ks[gene_name, "ks"]
  gene.data <- data[gene_name,]
  
  # to appease CRAN
  group <- NULL
  exp <- NULL
  
  df<-data.frame(row.names=colnames(data), group=outcomes, exp=gene.data)
  
  title <- paste(gene_name, "\n", "(ks score = ",
                 round(ks_score, 2), ")", sep="")
  
  ggplot(df, aes(exp, fill=group)) + geom_density(alpha=0.5) +
    xlab("data")  + ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24),
          legend.text = element_text(size = 24),
          legend.title = element_text(size=24))
}

#' @export
#' @title Plot distributions and CVM score for a gene.
#' @description The data for the specified gene is retrieved from
#' \code{cvmobj$data}. \code{outcomes} is used to divide the data into 
#' distributions for each group, which are then visualized as
#' density distributions. The calculated CVM score for the specified gene is
#' displayed in the plot title.
#' @param cvmobj An \code{\link{CVMomics}} object, typically returned via a call
#' to \code{\link{calculate_cvm}}.
#' @param gene_name The gene to visualize. The name should be defined as a row
#' name in \code{cvmobj$cvm}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the value is
#' not assigned, a plot will be drawn.
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
#' plot_cvm_density(results, "gene5")
#' 
#' @seealso \code{\link{calculate_cvm}} \code{\link[ggplot2]{ggplot}}
plot_cvm_density <- function(cvmobj, gene_name) {
  
  data <- cvmobj$data
  outcomes <- cvmobj$outcomes
  classes <- unique(outcomes)
  
  cvm_score <- cvmobj$cvm[gene_name, "cvm"]
  gene.data <- data[gene_name,]
  
  # to appease CRAN
  group <- NULL
  exp <- NULL
  
  df<-data.frame(row.names=colnames(data), group=outcomes, exp=gene.data)
  
  title <- paste(gene_name, "\n", "(cvm score = ",
                 round(cvm_score, 2), ")", sep="")
  
  ggplot(df, aes(exp, fill=group)) + geom_density(alpha=0.5) +
    xlab("data")  + ggtitle(title) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24),
          legend.text = element_text(size = 24),
          legend.title = element_text(size=24))
}
