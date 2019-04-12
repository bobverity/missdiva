
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
  plot1 <- plot1 + geom_raster(aes(x = Var1, y = Var2, fill = log(value)), data = reshape2::melt(x))
  
  # legends and scales
  plot1 <- plot1 + scale_fill_viridis_c(option = "plasma", name = "log(coverage)")
  plot1 <- plot1 + xlab("sample") + ylab("locus")
  
  return(plot1)
}





