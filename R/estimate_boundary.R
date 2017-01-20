#' Estimate Boundary of Asymmetric Hypothesis
#'
#' This function uses a variety of techniques to estimate data boundaries for asymmetric hypotheses.
#'
#' @param dat A Dataset
#' @param input Independent variable
#' @param output Dependent variable
#' @param method What boundary technique? One of: SFA, QR, Polynomial, Kernel.
#' @param AOC Should AOC be returned? Default is FALSE.
#' @param sufficient Is the data sufficient? Default is FALSE.
#' @import data.table
#' @import quantreg
#' @import frontier
#' @import npbr
#' @examples
#' DATA <- read.csv("dat.csv")
#' estimate_boundary(dat = DATA, input = DATA$x, output = DATA$y, sufficient = FALSE, merthod = "QR", AOC = FALSE)
#' @export
#'
estimate_boundary <- function(dat, input, output, sufficient = FALSE,
  method, AOC = FALSE)
{
  dat$x <- input
  dat$y <- output
  id <- order(dat$x)
  if(method == "Polynomial")
  {
    if(sufficient == TRUE)
    {
      dat$x <- 1 - dat$x
      dat$y <- 1 - dat$y
      odeg <- opt_degree(dat$x, dat$y, dat$x, prange = 0:20)
      boundary <- poly_estimate(dat$x, dat$y, dat$x, deg = odeg)
      boundary < 1 - boundary
    }
    if(sufficient == FALSE)
    {
      odeg <- opt_degree(dat$x, dat$y, dat$x, prange = 0:20)
      boundary <- poly_estimate(dat$x, dat$y, dat$x, deg = odeg)
    }
  }
  if(method == "Kernel")
  {
    if(sufficient == TRUE)
    {
      dat$x <- 1 - dat$x
      dat$y <- 1 - dat$y
      bw <- kern_smooth_bw(dat$x, dat$y)
      boundary <- kernel_smoothing(dat$x, dat$y, dat$x, h = bw)
      boundary <- 1 - boundary
    }
    if(sufficient == FALSE)
    {
      bw <- kern_smooth_bw(dat$x, dat$y)
      boundary <- kernel_smoothing(dat$x, dat$y, dat$x, h = bw)
    }

  }
  if(method == "SFA")
  {
    if(sufficient == TRUE)
    {
      dat$x <- 1 - dat$x
      dat$y <- 1 - dat$y
      ka_sfa <- suppressWarnings(sfa(dat$y ~ dat$x | dat$x, 
        ineffDecrease = TRUE, data = dat))
      flip_boundary <- ka_sfa$mleParam[1] + ka_sfa$mleParam[2] * dat$x
      boundary <- 1 - flip_boundary
    }
    if(sufficient == FALSE)
    {
      ka_sfa <- suppressWarnings(sfa(dat$y ~ dat$x | dat$x, 
        ineffDecrease = TRUE, data = dat))
      boundary <- ka_sfa$mleParam[1] + ka_sfa$mleParam[2] * input
    }
  }
  if(method == "QR")
  {
    if(sufficient == TRUE)
    {
      dat$x <- 1 - dat$x
      dat$y <- 1 - dat$y
      ka_qr <- rq(dat$y ~ dat$x, tau = 0.95, data = dat)
      flip_boundary <- ka_qr$coefficients[1] + ka_qr$coefficients[2] * dat$x
      boundary <- 1 - flip_boundary
    }
    if(sufficient == FALSE)
    {
      ka_qr <- rq(output ~ input, tau = 0.95, data = dat)
      boundary <- ka_qr$coefficients[1] + ka_qr$coefficients[2] * input
    }
  }
  if(AOC == TRUE)
  {
    id <- order(dat$x)
    divide_by_this <- (max(dat$x) - min(dat$x)) * (max(dat$y) - min(dat$y))
    AOC <- sum(diff(sort(dat$x)) * rollmean(sort(boundary), 2))
    AOC_lower <- sum(diff(dat$x[id]) * min(dat$y))
    AOC_difference <- AOC - AOC_lower
    AOC_percent <- AOC_difference/divide_by_this
    AOC_result <- 1 - AOC_percent
    results_list <- list(boundary, AOC_result)
    names(results_list) <- c("boundary", "AOC")
    return(results_list)
  }
  if(AOC == FALSE)
  {
    return(boundary)
  }
}
