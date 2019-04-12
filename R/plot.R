
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
  plot_df <- as.data.frame(pcoa$points)
  
  # get variance explained by first two components
  var1 <- sprintf("PC1 = %s%%", round(pcoa$var_explained[1]*100, digits = 2))
  var2 <- sprintf("PC2 = %s%%", round(pcoa$var_explained[2]*100, digits = 2))
  
  # plot PCoA first two axes
  plot1 <- ggplot() + theme_bw()
  plot1 <- plot1 + geom_point(aes(x = V1, y = V2), data = plot_df)
  plot1 <- plot1 + ggtitle(title) + xlab(var1) + ylab(var2)
  
  return(plot1)
}

#------------------------------------------------
#' @title Make exploratory plots of missingness by position and sample 
#'
#' @description Make exploratory plots of missingness by position and sample 
#'
#' @param coverage coverage matrix.
#' @param threshold below this threshold is coded as missing.
#'
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom stats dist
#'
#' @export

coverage_plots <- function(coverage, threshold) {
  
  # first converting the coverage matrix into a data frame 
  coverage_df <- as.data.frame(coverage)
  coverage_df$sample_id <- rownames(coverage_df)
  
  # need to turn data into a long format
  coverage_long <- coverage_df %>%
    gather(key = "POS", value = "coverage", -sample_id)
  
  # now filtering 
  coverage_long$filt_pass <- NA
  coverage_long$filt_pass[coverage_long$coverage >= threshold] <- 1
  coverage_long$filt_pass[coverage_long$coverage < threshold] <- 0
  coverage_long$filt_pass[is.na(coverage_long$coverage)] <- 0
  cov_table <- table(coverage_long$filt_pass)
  
  # now need to plot 
  coverage_pos <- coverage_long %>%
    group_by(POS) %>%
    summarise(pass_pct = mean(filt_pass, na.rm = T))
  
  
  # plotting the chromosome position vs. the proportion of reads that were > 20
  plot1 <- ggplot(coverage_pos, aes(x =POS, y = pass_pct, group = 1)) + geom_line() + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    ggtitle("% of samples that pass filtering at each position")
  
  # now we want to look at the proportion of sites that passed filtering for each sample 
  coverage_samp <- coverage_long %>%
    group_by(sample_id) %>%
    summarise(pass_pct = mean(filt_pass, na.rm = T))
  
  # we need to order the bars by the level of coverage 
  coverage_samp$sample_id <- factor(coverage_samp$sample_id, levels = coverage_samp$sample_id[order(-coverage_samp$pass_pct)])
  
  plot2 <- ggplot(coverage_samp, aes(x =sample_id, y = pass_pct))+ geom_bar(stat = "identity") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    ggtitle("% of sites that pass filtering for each sample")
  
  
  # histogram of coverage passing 
  plot3 <- ggplot(coverage_samp, aes(x = pass_pct)) + geom_histogram(breaks = seq(0,1, by =0.02)) + 
    ggtitle("Histogram of sample quality") + 
    theme_bw()
  
  return(list(cov_table, plot1, plot2, plot3))
}



