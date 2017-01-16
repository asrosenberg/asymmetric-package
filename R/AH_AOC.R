#' Calculate Area Outside the Curve
#'
#' \code{AH_AOC} calculates the Area Outside of the Curve for asymmetric boundaries. It is implemented for a variety of techniques and can return bootstrapped confidence intervals.
#'
#' @param dat A Dataset
#' @param input Independent variable
#' @param output Dependent variable
#' @param method What boundary method should be used? One of: SFA, QR, Polynomial, Kernel.
#' @param print_style "ascii" or "r"?
#' @param CI Should confidence intervals be returned?
#' @param nboots Number of bootstraps to generate the confidence intervals
#' @import data.table
#' @examples
#' DATA <- read.csv("dat.csv")
#' AH_AOC(DATA, input = DATA$x, output = DATA$y, method = "QR",
#'    print_style = "ascii", CI = TRUE, nboots = 5000)
#' @export
AH_AOC <- function(dat, input, output, method, print_style = "ascii",
  CI = FALSE, nboots)
{
  AOC_estimates <- lapply(method, which_technique_to_use, dat = dat,
    input = input, output = output)
  estimates <- round(do.call(rbind, AOC_estimates), digits = 3)
  rownames(estimates) <- method
  colnames(estimates) <- "AOC"
  if(CI == TRUE)
  {
    samples <- lapply(method, replicates, nboots = nboots, dat = dat,
      input = input, output = output)
    CIs <- lapply(samples, quantile, probs = c(0.025, 0.5, 0.975))
    CIs <- round(do.call(rbind, CIs), digits = 3)
    CI_estimates <- cbind(estimates, CIs)
    CI_estimates <- cbind(method, CI_estimates)
    estimates_df <- as.data.frame(CI_estimates)
  }
  if(CI == FALSE)
  {
    add_rows_names <- cbind(method, estimates)
    estimates_df <- as.data.frame(add_rows_names)
  }
  if(print_style == "ascii")
  {
    asciify(estimates_df)
  }
  if(print_style == "R")
  {
    return(estimates_df)
  }
}
