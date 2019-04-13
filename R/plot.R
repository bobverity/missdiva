
#------------------------------------------------
#' @title Plot matrix of coverage for all samples and all loci
#'
#' @description Plot matrix of coverage for all samples and all loci.
#'
#' @param x matrix of coverage.
#'
#' @import reshape2
#' @import ggplot2
#' @export

plot_coverage_matrix <- function(x = NULL) {
  
  # check inputs
  assert_matrix(x)
  
  # produce basic plot
  plot1 <- ggplot() + theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.y = element_blank(),
                            axis.ticks.y = element_blank())
  
  # add raster
  df_plot <- reshape2::melt(x)
  plot1 <- plot1 + geom_raster(aes_(x = ~V1, y = ~V2, fill = log(~value)), data = df_plot)
  
  # legends and scales
  plot1 <- plot1 + scale_fill_viridis_c(option = "plasma", name = "log(coverage)")
  plot1 <- plot1 + xlab("sample") + ylab("locus")
  
  return(plot1)
}

#------------------------------------------------
#' @title PCoA plot
#'
#' @description Plot the first two PCoA axes.
#'
#' @param pcoa result of running PCoA, e.g. through the function
#'   \code{get_pcoa()}.
#' @param title plot title.
#'
#' @import ggplot2
#' @export

plot_pcoa <- function(pcoa = NULL, title = NULL) {
  
  # get data into ggplot format
  df_plot <- as.data.frame(pcoa$points)
  
  # get variance explained by first two components
  var1 <- sprintf("PC1 = %s%%", round(pcoa$var_explained[1]*100, digits = 2))
  var2 <- sprintf("PC2 = %s%%", round(pcoa$var_explained[2]*100, digits = 2))
  
  # plot PCoA first two axes
  plot1 <- ggplot() + theme_bw()
  plot1 <- plot1 + geom_point(aes_(x = ~V1, y = ~V2), data = df_plot)
  plot1 <- plot1 + ggtitle(title) + xlab(var1) + ylab(var2)
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot quality of each locus
#'
#' @description Plot quality of each locus.
#'
#' @param missing_matrix matrix of missingness encoded as a binary variable.
#'
#' @import ggplot2
#' @export

plot_locus_quality <- function(missing_matrix) {
  
  # check inputs
  assert_matrix(missing_matrix)
  
  # get quality per locus
  q <- 1 - colMeans(missing_matrix, na.rm = TRUE)
  df_plot <- data.frame(pos = 1:length(q), quality = q)
  
  # basic plot
  plot1 <- ggplot() + theme_bw() + theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank())
  plot1 <- plot1 + geom_line(aes_(x = ~pos, y = ~quality), data = df_plot)
  plot1 <- plot1 + xlab("pos (relative)") + ylab("quality")
  plot1 <- plot1 + ylim(c(0,1))
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot quality of each sample
#'
#' @description Plot quality of each sample.
#'
#' @param missing_matrix matrix of missingness encoded as a binary variable.
#'
#' @import ggplot2
#' @export

plot_sample_quality <- function(missing_matrix) {
  
  # check inputs
  assert_matrix(missing_matrix)
  
  # get quality per sample
  q <- 1 - rowMeans(missing_matrix, na.rm = TRUE)
  df_plot <- data.frame(samp = 1:length(q), quality = q)
  
  # basic plot
  plot1 <- ggplot() + theme_bw() + theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank())
  plot1 <- plot1 + geom_line(aes_(x = ~samp, y = ~quality), data = df_plot)
  plot1 <- plot1 + xlab("sample") + ylab("quality")
  plot1 <- plot1 + ylim(c(0,1))
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot histogram of quality of each locus
#'
#' @description Plot histogram of quality of each locus.
#'
#' @param missing_matrix matrix of missingness encoded as a binary variable.
#'
#' @import ggplot2
#' @export

plot_locus_quality_hist <- function(missing_matrix) {
  
  # check inputs
  assert_matrix(missing_matrix)
  
  # get quality per locus
  q <- 1 - colMeans(missing_matrix, na.rm = TRUE)
  df_plot <- data.frame(x = q)
  
  # basic plot
  plot1 <- ggplot() + theme_bw()
  plot1 <- plot1 + geom_histogram(aes_(x = ~x), breaks = seq(0,1,0.02), data = df_plot)
  plot1 <- plot1 + xlab("quality") + ylab("frequency")
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot histogram of quality of each sample
#'
#' @description Plot histogram of quality of each sample.
#'
#' @param missing_matrix matrix of missingness encoded as a binary variable.
#'
#' @import ggplot2
#' @export

plot_sample_quality_hist <- function(missing_matrix) {
  
  # check inputs
  assert_matrix(missing_matrix)
  
  # get quality per sample
  q <- 1 - rowMeans(missing_matrix, na.rm = TRUE)
  df_plot <- data.frame(x = q)
  
  # basic plot
  plot1 <- ggplot() + theme_bw()
  plot1 <- plot1 + geom_histogram(aes_(x = ~x), breaks = seq(0,1,0.02), data = df_plot)
  plot1 <- plot1 + xlab("quality") + ylab("frequency")
  
  return(plot1)
}
