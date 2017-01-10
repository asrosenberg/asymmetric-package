library(asymmetric)
library(RColorBrewer)
library(zoo)
library(pracma)
library(plot3D)
options(stringsAsFactors = FALSE)
colors <- brewer.pal(7, "Set1")
load("inst/extdata/Coxappend.RData")

# Make Figure 6 and Figure 7
CGG.test <- sfa(x$enps ~ x$eneth * x$lnml | x$eneth * x$lnml, ineffDecrease=TRUE)
yhat <- CGG.test$mleParam[1] + x$eneth*CGG.test$mleParam[2] +
  x$lnml*CGG.test$mleParam[3] + x$eneth*x$lnml*CGG.test$mleParam[4]
xgrid <- seq(min(x$eneth), max(x$eneth), length=10)
wgrid <- seq(min(x$lnml), max(x$lnml), length=10)
frontier.pred <- function(x, w){
  CGG.test$mleParam[1] + x*CGG.test$mleParam[2] + w*CGG.test$mleParam[3] +
  x*w*CGG.test$mleParam[4]
}
ygrid <- outer(xgrid, wgrid, frontier.pred)
par(mfrow=c(1,2))
par(mar=c(2,2,2,2))
for(r in c(240,300)){
  p <- persp(xgrid, wgrid, ygrid, col = "#FFFFFF22", border = "#00000022",
  phi = 10, theta = r, xlab = "Electoral system permissiveness",
  ylab = "Social heterogeneity", zlab = "Party system size", cex.lab = 0.7)
  obs <- trans3d(x$eneth, x$lnml, x$enps, p)
  pred <- trans3d(x$eneth, x$lnml, yhat, p)
  points(obs, col = "black", pch = 19, cex = 0.2)
  segments(obs$x, obs$y, pred$x, pred$y, col = "#00000022")
}

CGG.test.reduced <- sfa(x$enps ~ x$eneth * x$lnml, ineffDecrease=TRUE)
lrtest <- -2*CGG.test.reduced$mleLogl + 2*CGG.test$mleLogl
pchisq(lrtest, df=4, lower.tail=FALSE)

resid <- x$enps - yhat
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
    sep="", collapse=" ")
}

x$var12[45] <- "St. Kitts and Nevis"
x$var12[46] <- "St. Lucia"
x$var12[47] <- "St. Vincent & Grenadines"
x$var12[51] <- "United Kingdom"
x$var12[52] <- "United States"

##############
# Figure 7
##############

dotchart(sort(resid), labels=sapply(tolower(x$var12[order(resid)]), simpleCap),
  cex=0.4, pch=19, lcolor="gray90", main="Party System Size 'Inefficiency'")
abline(v=0)

#############
# Do AOC Test
#############
CGG.test <- sfa(x$enps ~ x$eneth * x$lnml | x$eneth * x$lnml, ineffDecrease=TRUE)
yhat <- CGG.test$mleParam[1] + x$eneth*CGG.test$mleParam[2] +
  x$lnml*CGG.test$mleParam[3] + x$eneth*x$lnml*CGG.test$mleParam[4]

cgg_function <- function(x, y)
{
  CGG.test$mleParam[1] + x * CGG.test$mleParam[2] + y * CGG.test$mleParam[3] +
    x * y * CGG.test$mleParam[4]
}

plane_function <- function(x, y)
{
  CGG.test$mleParam[1]
}

AOC_SFA_3D <- function(x, y, fun)
{
  integral <- integral2(fun = cgg_function, xmin = min(x), xmax = max(x),
    ymin = min(y), ymax = max(y))
  AUC_lower <- integral2(fun = plane_function, xmin = min(x),
    xmax = max(x), ymin = min(y), ymax = max(y))
  AUC_dif <- integral$Q - AUC_lower$Q
  AUC <- (AUC_dif)/((max(x) - min(x)) * (max(y) -
    min(y)) * (max(cgg_function(x, y)) - min(cgg_function(x, y))))
  1 - AUC
}

#################
# AOC for Table 6
#################
AOC_SFA_3D(x = x$eneth, y = x$lnml, fun = cgg_function)

boot_CGG <- function(dat)
{
  idx <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
  sample.dta <- dat[idx,]
  x <- sample.dta$eneth
  y <- sample.dta$lnml
  CGG.test <- sfa(sample.dta$enps ~ x * y | x * y, ineffDecrease=TRUE)
  cgg_function <- function(x, y)
  {
    CGG.test$mleParam[1] + x * CGG.test$mleParam[2] + y * CGG.test$mleParam[3] +
      x * y * CGG.test$mleParam[4]
  }
  plane_function <- function(x, y)
  {
    CGG.test$mleParam[1]
  }
  integral <- integral2(fun = cgg_function, xmin = min(x), xmax = max(x),
    ymin = min(y), ymax = max(y))
  AUC_lower <- integral2(fun = plane_function, xmin = min(x),
    xmax = max(x), ymin = min(y), ymax = max(y))
  AUC_dif <- integral$Q - AUC_lower$Q
  AUC <- (AUC_dif)/((max(x) - min(x)) * (max(y) -
    min(y)) * (max(cgg_function(x, y)) - min(cgg_function(x, y))))
  1 - AUC
}
dat <- x
rm(x)
CGG_reps <- replicate(10000, boot_CGG(dat = dat))
quantile(CGG_reps, c(0.025, 0.975))


# generate CI with for loop
nboots <- 5000
boots <- c()
for(i in 1:nboots){
idx <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
sample.dta <- dat[idx,]
x <- sample.dta$eneth
y <- sample.dta$lnml
CGG.test <- sfa(sample.dta$enps ~ x * y | x * y, ineffDecrease=TRUE)
cgg_function <- function(x, y)
{
  CGG.test$mleParam[1] + x * CGG.test$mleParam[2] + y * CGG.test$mleParam[3] +
    x * y * CGG.test$mleParam[4]
}
plane_function <- function(x, y)
{
  CGG.test$mleParam[1]
}
integral <- integral2(fun = cgg_function, xmin = min(x), xmax = max(x),
  ymin = min(y), ymax = max(y))
AUC_lower <- integral2(fun = plane_function, xmin = min(x),
  xmax = max(x), ymin = min(y), ymax = max(y))
AUC_dif <- integral$Q - AUC_lower$Q
AUC <- (AUC_dif)/((max(x) - min(x)) * (max(y) -
  min(y)) * (max(cgg_function(x, y)) - min(cgg_function(x, y))))
boots[i] <- 1 - AUC
}

#################
# CIs for Table 6
#################
quantile(boots, c(0.025, 0.5, 0.975))
