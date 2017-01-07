library(asymmetric)

colors <- brewer.pal(7, "Set1")

# Make fake asymmetric data
x <- runif(75, min = 0, max = 1)
y <- 1/2*x^(1/2)*runif(75, min = 0, max = 2)
analysis_dat <- data.frame(x = x, y = y)
plot(x, y)

# QR, SFA, poly
qr_fake <- rq(y ~ x, tau = 0.95)
sfa_fake <- sfa(y ~ x | x, ineffDecrease = TRUE)
ols_fake <- lm(y ~ x)
bw <- kern_smooth_bw(x, y, method="u", technique="noh", bw_method="bic")
kernsmooth <- kernel_smoothing(x, y, y, h = bw)
par(mar=c(4,3,3,0), pty = "s")
plot(x, y, main = "Comparison of Each Technique on Simulated Asymmetric Data",
  xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,1), type = "p", pch = 19,
  col="gray", xaxp = c(0, 1, 1), yaxp = c(0, 1, 1), cex.lab = 0.8,
  cex.axis = 0.8, cex = 1.5)
abline(a = qr_fake$coefficients[1] , b = qr_fake$coefficients[2], col = colors[4],
  lty = "dotted", lwd = 2)
abline(a = sfa_fake$mleParam[1], b = sfa_fake$mleParam[2], col = colors[5],
  lty = "dotdash", lwd = 2)
points(x = sort(x), y = sort(kernsmooth), type = "l", lwd = 2,
  col = colors[1], lty = "twodash")
abline(a = ols_fake$coefficients[1], b = ols_fake$coefficients[2], col = colors[7],
  lty = 1, lwd = 2)
legend.txt <- c("Kernel", "QR (.95)", "SFA", "OLS")
legend('topleft', legend.txt,
  lty=c("twodash", "dotted", "dotdash", "solid"),
  col=c(colors[1], colors[4],colors[5], colors[7]),
  bty='n', cex=.75, lwd = 2)

# AOC on simulated data (Table 2)
AOC_fd <- AH_AOC(dat = analysis_dat, input = analysis_dat$x,
  output = analysis_dat$y, method = c("Kernel", "QR", "SFA"), print_style = "R",
  CI = TRUE, nboots = 5000)
toLatex(AOC_fd)
