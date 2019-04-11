
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


