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
#'
#' results <- calculate_emd(dat, groups, nperm=10)
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
#' 
#' results <- calculate_emd(dat, outcomes, nperm=10, parallel=FALSE)
#' plot_perms(results)
#' @seealso \code{\link{calculate_emd}} \code{\link[ggplot2]{ggplot}}
plot_perms <- function(emdobj) {

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
#' 
#' results <- calculate_emd(dat, outcomes, nperm=10, parallel=FALSE)
#' plot_density(results, "gene5")
#' 
#' @seealso \code{\link{calculate_emd}} \code{\link[ggplot2]{ggplot}}
plot_density <- function(emdobj, gene_name) {

  data <- emdobj$data
  outcomes <- emdobj$outcomes
  classes <- unique(outcomes)
  
  emd_score <- emdobj$emd[gene_name, "emd"]
  gene.data <- data[gene_name,]
  
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
