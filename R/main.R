
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
#' @title Crude method for concatenating vcfs
#'
#' @description A vector of input vcf file paths is specified. Each file is read
#'   in, and any lines starting with the "#" symbol are taken to be header
#'   lines. Header lines are retained for the first input file, but dropped
#'   thereafter when concatenating files. The final concatenated data is printed
#'   to file.
#'
#' @param input_paths vector of file paths to input vcf files.
#' @param output_path single file path specifying thr output vcf file.
#' @param check_overwrite Boolean, whether to prompt user before overwriting
#'   existing files.
#'
#' @export

vcf_concat_crude <- function(input_paths, output_path, check_overwrite = TRUE) {
  
  # check inputs
  assert_vector(input_paths)
  assert_single(output_path)
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
    
    # trim lines that start with # symbol
    if (i > 1) {
      w <- which(as.vector(mapply(function(x) substr(x,1,1) == "#", dat)))
      dat <- dat[-w]
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
#' @description Read in vcf object using the \code{vcfR} function \code{read.vcf}.
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
#' @description Extract coverage matrix, i.e. the "DP" element, from vcf object.
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
#' @description Use coverage values extracted by the function
#'   \code{get_coverage_matrix()} and, using a variable threshold, make a binary
#'   matrix of missingness.
#'
#' @param coverage coverage matrix, should be numeric values, can contain NAs.
#' @param threshold coverage values below this threshold are considered missing
#'   data.
#'
#' @export

make_missing_matrix <- function(coverage = NULL, threshold = NULL) {
  
  # convert matrix to binary using threshold 0 is non-missing 1 is missing
  missing_matrix <- ifelse(is.na(coverage) == TRUE, 1, ifelse(coverage <= threshold, 1, 0))
  
  return(missing_matrix)
}

#------------------------------------------------
#' @title Get pairwise missingness distances based on discordance
#'
#' @description Discordance here is defined as the proportion of loci that do
#'   not match in their missingness status, i.e. one sample is missing at that
#'   locus and the other is not. This is equivalent to the Manhattan distance
#'   between rows of the missingness matrix.
#'
#' @param missing_matrix matrix of missingness encoded as a binary variable.
#'
#' @importFrom stats dist
#' @export

get_missing_distance <- function(missing_matrix) {
  
  # compute concordance matrix
  ret <- as.matrix(dist(missing_matrix, method = "manhattan") / ncol(missing_matrix))
  
  return(ret)
}

#------------------------------------------------
#' @title Principal coordinates analysis (PCoA) on distance matrix
#'
#' @description Carry out principal coordinates analysis (PCoA) on a distance
#'   matrix.
#'
#' @param distance_matrix matrix of distances.
#'
#' @importFrom stats cmdscale
#' @export

get_pcoa <- function(distance_matrix = NULL) {
  
  # perform pcoa
  ret <- cmdscale(distance_matrix, eig = TRUE)
  
  # add variance explained to ret object
  ret$var_explained <- abs(ret$eig)/sum(abs(ret$eig))
  
  return(ret)
}
