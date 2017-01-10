library(asymmetric)
library(countrycode)
library(RColorBrewer)
library(memisc)
dat <- load_dataset("KA_data")
countryname <- dat$X
y <- dat$Supranationalist.government
x <- dat$Strong.regional.governance * dat$High.policy.conformity
y <- 1 - y
x <- 1 - x
analysis_dat <- data.frame(x, y, countryname)
bw <- kern_smooth_bw(analysis_dat$x, analysis_dat$y, method="u",
  technique="noh", bw_method="bic")
kernsmooth <- kernel_smoothing(analysis_dat$x, analysis_dat$y, analysis_dat$x,
  h = bw)

nboots <- 10000

bs.est <- replicate(nboots, bs.routine(dta = analysis_dat))
ymatrix <- apply(bs.est, 2, sort)

interval_lines_y <- apply(ymatrix, 1, quantile, probs = c(0.1, 0.9))
mean_y <- apply(ymatrix, 1, quantile, probs = 0.5)

xs <- with(dat, Strong.regional.governance*High.policy.conformity)
dat$abbrv <- countrycode(as.character(dat$X), "country.name", "iso3c")
par(mar=c(4,3,3,0), pty = "s")
set.seed(2)
plot(dat$Strong.regional.governance * dat$High.policy.conformity,
  dat$Supranationalist.government,
  main = "Kernel Boundary Estimate for Koenig-Archibugi (2004)",
  xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,1), type = "p", pch = 19,
  col="#00000055", xaxp = c(0, 1, 1), yaxp = c(0, 1, 1), cex.lab = 0.8,
  cex.axis = 0.8, cex = 1.5)
with(dat, text(jitter(Supranationalist.government, factor = 1) ~ xs,
  labels = dat$abbrv, pos = 3, cex = 0.75))
mtext("Support for Supranational CFSP", side=2, line=1)
mtext("Combination of Regional Governance and Policy Conformity",
  side=1, line=1)
polygon(x = c(sort(1 - analysis_dat$x), rev(sort(1 - analysis_dat$x))),
  c(sort(interval_lines_y[1, ]), rev(sort(interval_lines_y[2, ]))),
  col = "grey87", border = NA)
points(x = sort(1 - analysis_dat$x), y = sort(1 - kernsmooth), type = "l",
  lwd = 2, col = colors[1], lty = "longdash")
points(x = sort(1 - analysis_dat$x), y = sort(mean_y), type = "l",
  lwd = 2, col = colors[2], lty = "solid")
legend.txt <- c("Kernel", "Kernel Boot (0.5)")
legend('bottomright', legend.txt, lty=c("longdash","solid"),
       col=c(colors[1], colors[2]), bty='n', cex=.75, lwd = 2)

# AOC for Koenig-Archibugi example (Table 3)
KA_AOC <- AH_AOC(dat = analysis_dat, input = analysis_dat$x,
  output = analysis_dat$y, method = "Kernel", print_style = "R", CI = TRUE,
  nboots = 5000)
toLatex(KA_AOC)
