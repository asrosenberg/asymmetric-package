boundary_routine <- function(dat, input = input, output = output, method)
{
  dat$x <- input
  dat$y <- output
  id <- order(dat$x)
  divide_by_this <- (max(dat$x) - min(dat$x)) * (max(dat$y) - min(dat$y))
  ymin <- min(dat$y)
  idx <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
  sample.dta <- dat[idx,]
  y <- sample.dta$y
  x <- sample.dta$x
  if(method == "Polynomial")
  {
    odeg <- opt_degree(x, y, x, prange=0:20)
    polfront <- poly_estimate(x, y, x, deg = odeg)
    AUC_poly <- sum(diff(sort(x)) * rollmean(sort(polfront), 2))
    AUC_lower <- sum(diff(dat$x[id]) * ymin)
    AUC_difference <- AUC_poly - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    result <- 1 - AUC_percent
  }
  if(method == "Kernel")
  {
    bw <- kern_smooth_bw(x, y)
    kernsmooth <- kernel_smoothing(x, y, x, h = bw)
    AUC_kern <- sum(diff(sort(x)) * rollmean(sort(kernsmooth), 2))
    AUC_lower <- sum(diff(dat$x[id]) * ymin)
    AUC_difference <- AUC_kern - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    result <- 1 - AUC_percent
  }
  if(method == "SFA")
  {
    ka_sfa <- sfa(y ~ x | x, ineffDecrease = TRUE, data = sample.dta)
    y_sfa <- ka_sfa$mleParam[1] + ka_sfa$mleParam[2] * x
    AUC_sfa <- sum(diff(sort(dat$x)) * rollmean(sort(y_sfa), 2))
    AUC_lower <- sum(diff(dat$x[id]) * ymin)
    AUC_difference <- AUC_sfa - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    result <- 1 - AUC_percent
  }
  if(method == "QR")
  {
    ka_qr <- rq(output ~ input, tau = 0.95, data = sample.dta)
    y_qr <- ka_qr$coefficients[1] + ka_qr$coefficients[2] * input
    AUC_qr <- sum(diff(sort(input)) * rollmean(sort(y_qr), 2))
    AUC_lower <- sum(diff(input[id]) * ymin)
    AUC_difference <- AUC_qr - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    result <- 1 - AUC_percent
  }
  result
}

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

which_technique_to_use <- function(method, dat, input = input, output = output)
{
  id <- order(input)
  divide_by_this <- (max(input) - min(input)) * (max(output) - min(output))
  ymin <- min(dat$y)
  if(method == "QR"){
    ka_qr <- rq(output ~ input, tau = 0.95, data = dat)
    y_qr <- ka_qr$coefficients[1] + ka_qr$coefficients[2] * input
    AUC_qr <- sum(diff(input[id]) * rollmean(sort(y_qr), 2))
    AUC_lower <- sum(diff(input[id]) * ymin)
    AUC_difference <- AUC_qr - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    qr_result <- 1 - AUC_percent
    return(qr_result)
  }
  if(method == "SFA"){
    ka_sfa <- sfa(output ~ input | input, ineffDecrease = TRUE, data = dat)
    y_sfa <- ka_sfa$mleParam[1] + ka_sfa$mleParam[2] * input
    AUC_sfa <- sum(diff(input[id]) * rollmean(sort(y_sfa), 2))
    AUC_lower <- sum(diff(input[id]) * ymin)
    AUC_difference <- AUC_sfa - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    sfa_result <- 1 - AUC_percent
    return(sfa_result)
  }

  if(method == "Polynomial"){
    odeg <- opt_degree(xtab = input, ytab = output, x = input, prange = 1:20)
    polfront <- poly_estimate(input, output, x = input, deg = odeg)
    AUC_poly <- sum(diff(input[id]) * rollmean(sort(polfront), 2))
    AUC_lower <- sum(diff(input[id]) * ymin)
    AUC_difference <- AUC_poly - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    poly_result <- 1 - AUC_percent
    return(poly_result)
  }

  if(method == "Kernel"){
    bw <- kern_smooth_bw(input, output, method = "u", technique = "noh",
      bw_method="bic")
    kernsmooth <- kernel_smoothing(input, output, input, h = bw)
    AUC_kern <- sum(diff(input[id]) * rollmean(sort(kernsmooth), 2))
    AUC_lower <- sum(diff(input[id]) * ymin)
    AUC_difference <- AUC_kern - AUC_lower
    AUC_percent <- AUC_difference/divide_by_this
    kern_result <- 1 - AUC_percent
    return(kern_result)
  }
}

replicates <- function(method, nboots, dat, input, output)
{
  samples <- replicate(nboots, boundary_routine(dat, method = method,
    input = input, output = output))
}



asciify <- function(df, pad = 1, ...) {
  ## error checking
  stopifnot(is.data.frame(df))
  ## internal functions
  SepLine <- function(n, pad = 1) {
    tmp <- lapply(n, function(x, pad) paste(rep("-", x + (2* pad)),
                                            collapse = ""),
                  pad = pad)
    paste0("+", paste(tmp, collapse = "+"), "+")
  }
  Row <- function(x, n, pad = 1) {
    foo <- function(i, x, n) {
      fmt <- paste0("%", n[i], "s")
      sprintf(fmt, as.character(x[i]))
    }
    rowc <- sapply(seq_along(x), foo, x = x, n = n)
    paste0("|", paste(paste0(rep(" ", pad), rowc, rep(" ", pad)),
                      collapse = "|"),
           "|")
  }
  ## convert everything to characters
  df <- as.matrix(df)
  ## nchar in data
  mdf <- apply(df, 2, function(x) max(nchar(x)))
  ## nchar in names
  cnames <- nchar(colnames(df))
  ## max nchar of name+data per elements
  M <- pmax(mdf, cnames)
  ## write the header
  sep <- SepLine(M, pad = pad)
  writeLines(sep)
  writeLines(Row(colnames(df), M, pad = pad))
  writeLines(sep)
  ## write the rows
  for(i in seq_len(nrow(df))) {
    ## write a row
    writeLines(Row(df[i,], M, pad = pad))
    ## write separator
    writeLines(sep)
  }
  invisible(df)
}


desMthing <- function(i, xtab = xtab)
{
  xtab^i
}

gridMthing <- function(i, x = x, xtab = xtab)
{
  x^i
}
poly_estimate <- function(xtab, ytab, x, deg)
{
  opt_coef <- (1/(1:(deg + 1))) * (max(xtab)^(1:(deg + 1)) -
    min(xtab)^(1:(deg + 1)))
  desM <- sapply(0:deg, desMthing, xtab = xtab)
  dir <- c(rep(">=", length(xtab)))
  bounds <- list(lower = list(ind = 1:(deg + 1), val = rep(-Inf,
    deg + 1)), upper = list(ind = 1:(deg + 1), val = rep(Inf, deg + 1)))
  Sol <- Rglpk_solve_LP(obj = opt_coef, mat = desM,
    dir = c(rep(">=", length(xtab))), rhs = ytab,
    bounds, types = NULL, max = FALSE)
  OPT <- Sol$solution
  gridM <- sapply(0:deg, gridMthing, xtab = xtab, x = x)
  fitt <- gridM %*% OPT
  return(fitt)
}

difs <- function(i, xs, ys, x = x)
{
  log(sum(abs(ys - poly_estimate(xtab = xs,
    ytab = ys, x = x, i)))) + (i + 1)/length(xs)
}

opt_degree <- function(xtab, ytab, x, prange = 0:20)
{
  criteria <- cbind(prange, sapply(prange, difs, xs = xtab, ys = ytab, x = x))
  return(prange[which.min(criteria[, 2])])
}

bootstrap_degrees <- function(dta)
{
  opt_degree(xtab = dta$gdp, ytab = dta$life)
}

bs.sample <- function(dta, ...)
{
  idx <- sample(1:nrow(dta), nrow(dta), replace=TRUE)
  dta <- dta[idx,]
  rownames(dta) <- NULL
  dta
}

bootstrap_estimates <- function(i)
{
  dta <- bootstrapped_data[[i]]
  deg <- degrees[[i]]
  poly_estimate(dta[, gdp], dta[, life], x = mydat$gdp, deg = deg)
}

bs.routine <- function(dta)
{
  idx <- sample(1:nrow(dta), nrow(dta), replace=TRUE)
  sample.dta <- dta[idx,]
  y <- sample.dta$y
  x <- sample.dta$x
  odeg <- opt_degree(x, y, analysis_dat$x, prange=0:20)
  polfront <- poly_estimate(x, y, analysis_dat$x, deg = odeg)
  #ys <- 1 - polfront
  #xs <- 1 - x
  #list(cbind(xs, ys))
  #list(cbind(1 - x, 1 - polfront))
  c(1 - polfront)
}

bs.routine.beta <- function(dta)
{
  idx <- sample(1:nrow(dta), nrow(dta), replace=TRUE)
  sample.dta <- dta[idx,]
  y <- sample.dta$y
  x <- sample.dta$x
  bw <- kern_smooth_bw(x, y, method="u", technique="noh", bw_method="bic")
  kernsmooth <- kernel_smoothing(x, y, analysis_dat$x, h = bw)
  c(1 - kernsmooth)
}

bs.routine_general <- function(xtab, dta, method)
{
  idx <- sample(1:nrow(dta), nrow(dta), replace=TRUE)
  sample.dta <- dta[idx,]
  y <- sample.dta$y
  x <- sample.dta$x
  if(method == "Kernel")
  {
    bw <- kern_smooth_bw(x, y, method="u", technique="noh", bw_method="bic")
    kernsmooth <- kernel_smoothing(x, y, x, h = bw)
    return(kernsmooth)
  }
  if(method == "Poly")
  {
    odeg <- opt_degree(x, y, xtab, prange=0:20)
    polfront <- poly_estimate(x, y, xtab, deg = odeg)
    return(polfront)
  }
  if(method == "QR")
  {
    qr <- rq(y ~ x, tau = 0.95, data = sample.dta)
    return(qr$fitted.values)
  }
}

comb_dat <- function(dat, ...)
{
  dat[, gdp]
}

put_wealth_health_together <- function(i)
{
  cbind(wealth[[i]], hw_bs_est[[i]])
}

kern_smooth_bw<-function(xtab, ytab, method="u", technique="noh", bw_method="bic")
{
  n<-length(xtab)
  ndata<-length(xtab) # number of data points

  # sorting step to use the Priestly-Chao estimator
  oind<-order(xtab)
  xtab<-xtab[oind]
  ytab<-ytab[oind]

  if (technique=="noh" && bw_method=="bic")
  {
    h_min<-2*max(diff(sort(xtab)))
    h_max<-(max(xtab)-min(xtab))
    h_grid<-seq(h_min,h_max,length=184)

    BIC<-NULL

    for ( h in h_grid)
    {
      Phi<-kernel_smoothing(xtab,ytab,xtab,h,method,technique="noh")

      # calculating model complexity
      xtab2<-xtab[-1]
      ytab2<-ytab[-1]
      DX<-diff(xtab)

      A<-dnorm(outer(xtab,xtab2,'-')/h)/h
      B<-rep(1,n) %*% t(ytab2*DX)
      S<-A * rep(1,n) %*% t(DX)
      SS<-S[-1,]
      DF<-sum(diag(SS))

      BIC<-c(BIC,log(sum(abs(ytab-Phi)))+log(length(ytab))*DF/(2*length(ytab)))
    }

    if (which.min(BIC)==1) {mind<-which.min(BIC[-1])+1}
    if (which.min(BIC)>1) {mind<-which.min(BIC)}
    hopt<- h_grid[mind]
  }

  if (technique=="noh" && bw_method=="cv")
  {
    bw<-npregbw(xdat=xtab, ydat=ytab, regtype="lc", ckertype="gaussian")
    hopt<-max(max(diff(sort(xtab))),bw$bw)
  }


  if (technique=="pr" && bw_method=="bic")
  {
    h_min<-2*max(diff(sort(xtab)))
    h_max<-(max(xtab)-min(xtab))
    h_grid<-seq(h_min,h_max,length=30)

    BIC<-NULL

    for ( h in h_grid)
    {
      Phi<-kernel_smoothing(xtab,ytab,xtab,h,method,technique="pr")

      # calculating model complexity
      xtab2<-xtab[-1]
      ytab2<-ytab[-1]
      DX<-diff(xtab)

      A<-dnorm(outer(xtab,xtab2,'-')/h)/h
      B<-rep(1,n) %*% t(ytab2*DX)
      S<-A * rep(1,n) %*% t(DX)
      SS<-S[-1,]
      DF<-sum(diag(SS))

      BIC<-c(BIC,log(sum(abs(ytab-Phi)))+log(length(ytab))*DF/(2*length(ytab)))
    }
    if (which.min(BIC)==1) {mind<-which.min(BIC[-1])+1}
    if (which.min(BIC)>1) {mind<-which.min(BIC)}
    hopt<- h_grid[mind]

  }

  if (technique=="pr" && bw_method=="cv")
  {
    bw<-npregbw(xdat=xtab, ydat=ytab, regtype="lc", ckertype="gaussian")
    hopt<-max(max(diff(sort(xtab))),bw$bw)

  }
  return(hopt)

}

kernel_smoothing <- function (xtab, ytab, x, h, method = "u", technique = "noh")
{
  stopifnot(method %in% c("u", "m", "mc"))
  stopifnot(technique %in% c("pr", "noh"))
  consX <- unique(sort(c(xtab, seq(min(xtab), max(xtab), length = 201))))
  Cstable = 10000
  n <- length(xtab)
  ndata <- length(xtab)
  ncons <- length(consX)
  oind <- order(xtab)
  xtab <- xtab[oind]
  ytab <- ytab[oind]
  if (technique == "noh") {
    oind <- order(xtab)
    xtab <- xtab[oind]
    ytab <- ytab[oind]
    xtab2 <- xtab[-1]
    ytab2 <- ytab[-1]
    DX <- diff(xtab)
    r_end <- max(xtab)
    l_end <- min(xtab)
    opt_coef <- (pnorm((r_end - xtab2)/h) - pnorm((l_end -
      xtab2)/h)) * DX * ytab2
    C0 <- .fitMat(xtab, ytab, xtab, h, type = 0)
    C1 <- .fitMat(xtab, ytab, consX, h, type = 1)
    C2 <- .fitMat(xtab, ytab, consX, h, type = 2)
    obj <- c(opt_coef, -opt_coef)
    if (method == "u") {
      mat_temp <- C0
      mat <- rbind(cbind(mat_temp, -mat_temp), rep(1, length(obj)))
      rhs <- c(ytab, Cstable)
      dir <- c(rep(">=", ndata), "<=")
    }
    if (method == "m") {
      mat_temp <- rbind(C0, C1)
      mat <- rbind(cbind(mat_temp, -mat_temp), rep(1, length(obj)))
      rhs <- c(ytab, rep(0, ncons), Cstable)
      dir <- c(rep(">=", ndata + ncons), "<=")
    }
    if (method == "mc") {
      mat_temp <- rbind(C0, C1, C2)
      mat <- rbind(cbind(mat_temp, -mat_temp), rep(1, length(obj)))
      rhs <- c(ytab, rep(0, 2 * ncons), Cstable)
      dir <- c(rep(">=", ndata + ncons), rep("<=", ncons),
        "<=")
    }
    Sol <- Rglpk_solve_LP(obj, mat, dir, rhs, types = NULL,
        max = FALSE)
    OPT_temp <- Sol$solution
    OPT <- OPT_temp[1:length(opt_coef)] - OPT_temp[(length(opt_coef) +
        1):(length(opt_coef) * 2)]
    fitt <- .fitMat(xtab, ytab, x, h, 0) %*% OPT
  }
  if (technique == "pr") {
    A <- .fitMat_nw(xtab, ytab, xtab, h, type = 0)
    A.deriv.1 <- .fitMat_nw(xtab, ytab, consX, h, type = 1)
    A.deriv.2 <- .fitMat_nw(xtab, ytab, consX, h, type = 2)
    if (method == "u") {
      p.hat <- solve.QP(Dmat = diag(n), dvec = rep(1, n),
        Amat = t(A), bvec = ytab, meq = 0)$solution
    }
    if (method == "m") {
      p.hat <- solve.QP(Dmat = diag(n), dvec = rep(1, n),
        Amat = cbind(t(A), t(A.deriv.1)), bvec = c(ytab,
        rep(0, ncons)), meq = 0)$solution
    }
    if (method == "mc") {
      p.hat <- solve.QP(Dmat = diag(n), dvec = rep(1, n),
        Amat = cbind(t(A), t(A.deriv.1), -t(A.deriv.2)),
        bvec = c(ytab, rep(0, ncons), rep(0, ncons)),
        meq = 0)$solution
    }
    fitt <- .fitMat_nw(xtab, ytab, x, h, type = 0) %*% p.hat
  }
  return(as.vector(fitt))
}

solve.QP <- function (Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE)
{
  n <- nrow(Dmat)
  q <- ncol(Amat)
  if (missing(bvec))
    bvec <- rep(0, q)
  if (n != ncol(Dmat))
    stop("Dmat is not symmetric!")
  if (n != length(dvec))
    stop("Dmat and dvec are incompatible!")
  if (n != nrow(Amat))
    stop("Amat and dvec are incompatible!")
  if (q != length(bvec))
    stop("Amat and bvec are incompatible!")
  if ((meq > q) || (meq < 0))
    stop("Value of meq is invalid!")
  iact <- rep(0, q)
  nact <- 0
  r <- min(n, q)
  sol <- rep(0, n)
  lagr <- rep(0, q)
  crval <- 0
  work <- rep(0, 2 * n + r * (r + 5)/2 + 2 * q + 1)
  iter <- rep(0, 2)
  res1 <- .Fortran(.QP_qpgen2, as.double(Dmat), dvec = as.double(dvec),
    as.integer(n), as.integer(n), sol = as.double(sol), lagr = as.double(lagr),
    crval = as.double(crval), as.double(Amat), as.double(bvec),
    as.integer(n), as.integer(q), as.integer(meq), iact = as.integer(iact),
    nact = as.integer(nact), iter = as.integer(iter), work = as.double(work),
    ierr = as.integer(factorized))
  if (res1$ierr == 1)
    stop("constraints are inconsistent, no solution!")
  else if (res1$ierr == 2)
    stop("matrix D in quadratic function is not positive definite!")
  list(solution = res1$sol, value = res1$crval, unconstrained.solution = res1$dvec,
       iterations = res1$iter, Lagrangian = res1$lagr, iact = res1$iact[1:res1$nact])
}

.fitMat <- function (X, Y, consX, h, type = 0)
{
  n <- length(X)
  stopifnot(sum((order(X) - 1:n)^2) == 0)
  X2 <- X[-1]
  Y2 <- Y[-1]
  DX <- diff(X)
  A <- dnorm(outer(consX, X2, "-")/h)/h
  B <- rep(1, length(consX)) %*% t(Y2 * DX)
  if (type == 0) {
    return(A * B)
  }
  if (type == 1) {
    A2 <- .dderi(outer(consX, X2, "-")/h)/h^2
    return(A2 * B)
  }
  if (type == 2) {
    A3 <- .d2deri(outer(consX, X2, "-")/h)/h^3
    return(A3 * B)
  }
}

.dderi <- function (x)
{
  return(-x * dnorm(x))
}

.d2deri <- function (x)
{
  return((x^2 - 1) * dnorm(x))
}

