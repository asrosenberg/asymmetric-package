#' Calculate Area Outside the Curve
#'
#' @param dat A Dataset
#' @param input Independent variable
#' @param output Dependent variable
#' @param method What boundary method should be used?
#' @param print_style "ascii" or "r"?
#' @param CI Should confidence intervals be returned?
#' @param nboots Number of bootstraps to generate the confidence intervals
#' @import data.table
#' @import quantreg
#' @import frontier
#' @import npbr
#' @export
estimate_boundary <- function(dat, input, output, sufficient = FALSE,
  method, AOC = FALSE)
{
  dat$x <- input
  dat$y <- output
  id <- order(dat$x)
  if(method == "Polynomial")
  {
    odeg <- opt_degree(x, y, x, prange = 0:20)
    boundary <- poly_estimate(x, y, x, deg = odeg)
  }
  if(method == "Kernel")
  {
    bw <- kern_smooth_bw(x, y)
    boundary <- kernel_smoothing(x, y, x, h = bw)
  }
  if(method == "SFA")
  {
    options(warn = -1)
    ka_sfa <- sfa(y ~ x | x, ineffDecrease = TRUE, data = dat)
    boundary <- ka_sfa$mleParam[1] + ka_sfa$mleParam[2] * input
    options(warn = 0)
  }
  if(method == "QR")
  {
    ka_qr <- rq(output ~ input, tau = 0.95, data = dat)
    boundary <- ka_qr$coefficients[1] + ka_qr$coefficients[2] * input
  }
  if(AOC == TRUE)
  {
    AOC <- sum(diff(sort(x)) * rollmean(sort(boundary), 2))
    AOC_lower <- sum(diff(dat$x[id]) * ymin)
    AOC_difference <- AOC - AOC_lower
    AOC_percent <- AOC_difference/divide_by_this
    AOC_result <- 1 - AOC_percent
  }
  return(boundary)
}
