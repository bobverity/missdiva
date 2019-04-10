
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
#' @title Simulate biallelic data
#'
#' @description Simulate biallelic data from a simple statistical model. Inputs
#'   include the complexity of infection (COI), population-level allele
#'   frequencies (PLAF) and some parameters dicating skew and error
#'   distributions. Outputs include the phased haplotypes and the un-phased read
#'   count and coverage data.
#'
#' @details Simulated data are drawn from a simple statistical model:
#'   \enumerate{
#'     \item Strain proportions are drawn from a symmetric Dirichlet
#'     distribution with shape parameter \code{alpha}.
#'     \item Phased haplotypes are drawn at every locus, one for each
#'     \code{COI}. The allele at each locus is drawn from a Bernoulli
#'     distribution with probability given by the \code{PLAF}.
#'     \item The "true" within-sample allele frequency at every locus is
#'     obtained by multiplying haplotypes by their strain proportions, and
#'     summing over haplotypes. Errors are introduced through the equation
#'     \deqn{wsaf_error = wsaf*(1-e) + (1-wsaf)*e}where \eqn{wsaf} is the WSAF
#'     without error and \eqn{e} is the error parameter \code{epsilon}.
#'     \item Final read counts are drawn from a beta-binomial distribution with
#'     expectation \eqn{w_error}. The raw number of draws is given by the
#'     \code{coverage}, and the skew of the distribution is given by the
#'     \code{overdispersion} parameter. If \code{overdispersion = 0} then the
#'     distribution is binomial, rather than beta-binomial.
#'   }
#'
#' @param COI complexity of infection.
#' @param PLAF vector of population-level allele frequencies at each locus.
#' @param coverage coverage at each locus. If a single value then the same
#'   coverage is applied over all loci.
#' @param prob_missing the probabilty of each locus failing. If a locus fails
#'   then coverage is set to \code{NA}. Can be a vector of values over all loci,
#'   or a single value in which case the same value applies over all loci.
#' @param alpha shape parameter of the symmetric Dirichlet prior on strain
#'   proportions.
#' @param overdispersion the extent to which counts are over-dispersed relative
#'   to the binomial distribution. Counts are Beta-binomially distributed, with
#'   the beta distribution having shape parameters \code{p/overdispersion} and
#'   \code{(1-p)/overdispersion}.
#' @param epsilon the probability of a single read being mis-called as the other
#'   allele. Applies in both directions.
#'
#' @export

sim_biallelic <- function(COI = 3,
                          PLAF = runif(10,0,0.5),
                          coverage = 100,
                          prob_missing = 0.1,
                          alpha = 1,
                          overdispersion = 0,
                          epsilon = 0) {
  
  # check inputs
  assert_single_pos_int(COI)
  assert_vector(PLAF)
  assert_bounded(PLAF)
  L <- length(PLAF)
  if (length(coverage) == 1) {
    coverage <- rep(coverage, L)
  }
  assert_vector(coverage)
  assert_pos_int(coverage, zero_allowed = TRUE)
  assert_bounded(prob_missing)
  if (length(prob_missing) == 1) {
    prob_missing <- rep(prob_missing, L)
  }
  assert_same_length_multiple(PLAF, coverage, prob_missing)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion, zero_allowed = TRUE)
  assert_single_pos(epsilon, zero_allowed = TRUE)
  assert_bounded(epsilon)
  
  # generate strain proportions
  w <- rdirichlet(rep(alpha, COI))
  
  # generate true WSAF levels by summing binomial draws over strain proportions
  m <- mapply(function(x) rbinom(COI, 1, x), x = PLAF)
  p_levels <- colSums(sweep(m, 1, w, "*"))
  
  # add in genotyping error
  p_error <- p_levels*(1-epsilon) + (1-p_levels)*epsilon
  
  # draw read counts, taking into account overdispersion
  if (overdispersion == 0) {
    counts <- rbinom(L, size = coverage, prob = p_error)
  } else {
    counts <- rbetabinom(L, k = coverage, alpha = p_error/overdispersion, beta = (1-p_error)/overdispersion)
  }
  
  # calculate within-sample allele frequencies
  wsaf <- counts/coverage
  
  # make some loci missing
  missing_draw <- which(rbinom(L, 1, prob_missing) == 1)
  coverage[missing_draw] <- 0
  counts[missing_draw] <- 0
  wsaf[missing_draw] <- NA
  
  # return list
  ret <- list(COI = COI,
              strain_proportions = w,
              phased = m,
              data = data.frame(PLAF = PLAF,
                                coverage = coverage,
                                counts = counts,
                                wsaf = wsaf))
  return(ret)
}

