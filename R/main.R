
#------------------------------------------------
#' @useDynLib missdiva, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
#' @title Check that missdiva package has loaded successfully
#'
#' @description Simple function to check that missdiva package has loaded
#'   successfully. Prints "missdiva loaded successfully!" if so.
#'
#' @export

check_missdiva_loaded <- function() {
  message("missdiva loaded successfully!")
}

#------------------------------------------------
#' @title Dummy function
#'
#' @description Simple test function that demonstrates some of the features of
#'   this package.
#'
#' @details Takes a vector of values, returns the square.
#'
#' @param x vector of values
#'
#' @export
#' @examples
#' # Find square of first 100 values
#' dummy1(1:100)

dummy1 <- function(x = 1:5) {
  
  # print message to console
  message("running R dummy1 function")
  
  # get arguments in list form
  args <- list(x = x)
  
  # run C++ function with these arguments
  output_raw <- dummy1_cpp(args)
  
  # some optional processing of output
  message("processing output")
  ret <- output_raw$x_squared
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Crude method for stitching together vcfs
#'
#' @description Crude method for stitching together vcfs.
#'
#' @param input_paths vector of file paths to input vcf files.
#' @param output_path single file path specifying thr output vcf file.
#' @param trim_lines how many lines at the top of the vcf file to drop when
#'   stitching.
#' @param check_overwrite Boolean, whether to prompt user before overwriting
#'   existing files.
#'
#' @export

vcf_merge_crude <- function(input_paths, output_path, trim_lines = 5, check_overwrite = TRUE) {
  
  # check inputs
  assert_vector(input_paths)
  assert_single(output_path)
  assert_single_pos_int(trim_lines, zero_allowed = TRUE)
  assert_single_logical(check_overwrite)
  
  # check that all input files exist
  if (!all(file.exists(input_paths))) {
    w <- which(!file.exists(input_paths))[1]
    stop(sprintf("input file %s could not be found", input_paths[w]))
  }
  
  # check whether to overwrite existing
  if (file.exists(output_path) & check_overwrite) {
    if (!user_yes_no("output file already exists. Overwrite? (Y/N)")) {
      message("returning without overwriting")
      return()
    }
    message("overwriting existing")
  }
  
  # loop through input files
  nf <- length(input_paths)
  ret <- NULL
  for (i in 1:nf) {
    
    # read input file
    con <- file(input_paths[i], "r")
    dat <- readLines(con)
    close(con)
    
    # trim header lines
    if (i > 1 & trim_lines > 0) {
      dat <- dat[-(1:trim_lines)]
    }
    
    # stitch together
    ret <- c(ret, dat)
  }
  
  # write output to file
  con <- file(output_path, "w")
  writeLines(ret, con)
  close(con)
}

#------------------------------------------------
#' @title Read in vcf object
#'
#' @description Read in vcf object.
#'
#' @param file path to vcf file.
#' @param verbose if reading from file, whether to read in verbose manner.
#'
#' @export

read_vcf <- function(file, verbose = TRUE) {
  
  # check inputs
  assert_single_logical(verbose)
  
  # read in vcf
  vcf <- vcfR::read.vcfR(file = file, verbose = verbose)
  
  return(vcf)
}

#------------------------------------------------
#' @title Extract coverage matrix from vcf object
#'
#' @description Extract coverage matrix from vcf object.
#'
#' @param file path to vcf file.
#' @param vcf object of class \code{vcfR}.
#' @param verbose if reading from file, whether to read in verbose manner.
#'
#' @import ggplot2
#' @importFrom stats prcomp
#' @export

get_coverage_matrix <- function(file = NULL, vcf = NULL, verbose = TRUE) {
  
  # check inputs
  if (!xor(!is.null(file), !is.null(vcf))) {
    stop("Must specify one input: either a raw vcf file path or a vcfR object")
  }
  
  # get vcf object
  if (!is.null(vcf)){
    assert_custom_class(vcf, "vcfR")
  } else {
    assert_file_exists(file)
    vcf <- vcfR::read.vcfR(file = file, verbose = verbose)
  }
  
  # extract coverage matrix
  coverage <- t(vcfR::extract.gt(vcf, element = "DP", as.numeric = T))
  
  return(coverage)
}


#------------------------------------------------
#' @title Make a matrix of missingness using coverage matrix
#'
#' @description Use coverage values extracted by the coverage function and using a variable threshold make a matrix of missingness, with a binary variable.
#'
#' @param coverage coverage matrix, should be numeric values, can contain NAs.
#' @param threshold numeric object saying what the minimum coverage threshold is.
#'
#' @export

make_missing_matrix <- function(coverage=NULL, threshold=NULL) {
  
  # convert matrix to binary using threshold 0 is non-missing 1 is missing
  missing_matrix<- ifelse(is.na(coverage)==TRUE,1,ifelse(coverage<=threshold,1,0))
  return(missing_matrix)
}

#------------------------------------------------
#' @title Get pairwise distances based on concordance
#'
#' @description Get pairwise distances based on concordance.
#'
#' @param missing_matrix matrix of missingness encoded as a binary variable
#'
#' @export

get_missing_distance <- function(missing_matrix = NULL) {
  
  # compute concordance matrix
  ret <- as.matrix(dist(missing_matrix, method = "manhattan") / ncol(missing_matrix))
  
  return(ret)
}

#------------------------------------------------
#' @title TODO
#'
#' @description TODO
#'
#' @param distance_matrix matrix of distance
#'
#' @importFrom stats prcomp
#' @export

missing_pca <- function(distance_matrix = NULL) {
  
  # perform PCA
  pca <- prcomp(distance_matrix, scale = FALSE)
  
  return(pca)
}
  
  
#------------------------------------------------
#' @title Do a PCA plot of missingness 
#'
#' @description plot the results of PCA using missingness distance matrix
#'
#' @param pca_distance PCA on distance matrix
#'
#' @import ggplot2
#' @export

pca_missingness_plot <- function(pca_dist=NULL,title=NULL) {
  
  #get percent variance explained by first 2 pr comps
  eigs<- pca_dist$sdev^2
  
  pct_var_exp_1<- as.character(eigs[1]/sum(eigs))
  pct_var_exp_2<- as.character(eigs[2]/sim(eigs))
  
  pca_dist$x<- as.data.frame(pca_dist$x)
  #plot PCA
  pca_missingness_plot<-ggplot(data=pca_on_dist_matrix$x,aes(x=PC1,y=PC2))
  pca_missingness_plot+geom_point()+ggtitle(title)+xlab(pct_var_exp_1) + ylab(pct_var_exp_2)
  
  return(pca_missingness_plot)
}


#------------------------------------------------
#' @title Making exploratory plots of missingness by position and sample 
#'
#' @description from the binary missingness matrix compute and plot missingness 
#'
#' @param missing_matrix matrix of missingness encoded as a binary variable
#'
#'
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stats prcomp
#' @importFrom stats dist
#'
#' @export


## 
coverage_plots <- function(coverage, threshold) { 
  # first convering the coverage matrix into a data frame 
  coverage_df <- as.data.frame(coverage)
  coverage_df$sample_id <- rownames(coverage_df)
  
  ### need to turn data into a long format
  coverage_long <- coverage_df %>%
    gather(key = "POS", value = "coverage", -sample_id)
  
  ## now filtering 
  coverage_long$filt_pass <- NA
  coverage_long$filt_pass[coverage_long$coverage >= threshold] <- 1
  coverage_long$filt_pass[coverage_long$coverage < threshold] <- 0
  coverage_long$filt_pass[is.na(coverage_long$coverage)] <- 0
  cov_table <- table(coverage_long$filt_pass)
  
  ## now need to plot 
  coverage_pos <- coverage_long %>%
    group_by(POS) %>%
    summarise(pass_pct = mean(filt_pass, na.rm = T))
  
  
  ### plotting the chromosome position vs. the proportion of reads that were > 20
  plot1 <- ggplot(coverage_pos, aes(x =POS, y = pass_pct, group = 1)) + geom_line() + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    ggtitle("% of samples that pass filtering at each position")
  
  ### now we want to look at the proportion of sites that passed filtering for each sample 
  coverage_samp <- coverage_long %>%
    group_by(sample_id) %>%
    summarise(pass_pct = mean(filt_pass, na.rm = T))
  
  # we need to order the bars by the level of coverage 
  coverage_samp$sample_id <- factor(coverage_samp$sample_id, levels = coverage_samp$sample_id[order(-coverage_samp$pass_pct)])
  
  plot2 <- ggplot(coverage_samp, aes(x =sample_id, y = pass_pct))+ geom_bar(stat = "identity") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    #theme_bw() +
    ggtitle("% of sites that pass filtering for each sample")
    
  
  # histogram of coverage passing 
  plot3 <- ggplot(coverage_samp, aes(x = pass_pct)) + geom_histogram(breaks = seq(0,1, by =0.02)) + 
    ggtitle("Histogram of sample quality") + 
    theme_bw()
  
  return(list(cov_table, plot1, plot2, plot3))
}

