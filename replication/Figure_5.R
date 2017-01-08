library(asymmetric)
library(data.table)
# Do QR on HW
dat_gdp <- load_dataset("gdp_ppp")
dat_life <- load_dataset("life_expect")
# log GDP
dat_gdp[, 2:217] <- log(dat_gdp[, 2:217])
columnsgdp <- c("state", paste("gdp",seq(from = 1800, to = 2015), sep = "_"))
columnslife <- c("state", paste("life",seq(from = 1800, to = 2015), sep = "_"))
colnames(dat_gdp) <- columnsgdp
colnames(dat_life) <- columnslife
dat_life <- dat_life[,-1]
dat <- data.frame(cbind(dat_gdp, dat_life))
dat <- dat[complete.cases(dat),]
gdps <- dat[2:217]
dat <- melt(dat)
dat$variable <- as.character(dat$variable)
variable <- do.call(rbind, strsplit(dat$variable, split =  "_"))[,1]
date <- do.call(rbind, strsplit(dat$variable, split =  "_"))[,2]
dat$variable <- variable
dat$date <- as.numeric(date)
dat <- dcast(dat, state + date ~ variable)
setDT(dat)
setnames(dat, c("state", "year", "gdp", "life"))


qr_funct <- function(date, DT) {
  sub_dat <- DT[year == date]
  qr <- rq(life ~ gdp, tau = 0.95, data = sub_dat)
  qr$fitted
}
get_gdps <- function(date)
{
  sub_dat <- dat[year == date]
  sub_dat$gdp
}
dates <- seq(from = 1810, to = 2015, by = 20)
yrs20 <- seq(1, 215, by = 20)
gdps <- as.data.table(sapply(dates, get_gdps))
results <- as.data.table(sapply(dates, qr_funct, DT = dat))
setnames(gdps, as.character(yrs20))
setnames(results, as.character(yrs20))
gdps[["1810"]]
points(gdps[["1810"]], results[["1810"]])

# Get Bootstraps for every 20 years
qr_20_years_bootstrap <- function(dat, date, nboots)
{
  boots <- function(dat, date){
    dta <- dat[year == date]
    idx <- sample(1:nrow(dta), nrow(dta), replace=TRUE)
    sample.dta <- dta[idx,]
    qr <- rq(life ~ gdp, tau = 0.95, data = sample.dta)
    qr$coefficients[1] + dta[, gdp] * qr$coefficients[2]
  }
  reps <- apply(replicate(nboots, boots(dat, date)), 2, sort)
  quantiles <- apply(reps, 1, quantile, probs = c(0.025, 0.975))
  quantiles
}

dates <- seq(from = 1810, to = 2010, by = 20)
bootstraps <- lapply(dates, qr_20_years_bootstrap, dat = dat, nboots = 5000)
names(bootstraps) <- seq(1, 215, by = 20)

# Make Figure 5
yrs20 <- seq(1, 215, by = 20)
par(mar=c(2,2,0,0))
par(mfrow=c(2,6), oma = c(0, 0, 2, 0))
for(i in seq(from = 1, to = 215, by = 20)){
  plot(0, xlim=c(5,12), ylim=c(15,90), type = "n", xlab = " ", ylab = " ",
    main = " ", axes=FALSE)
  axis(2, labels=FALSE, col="#88888888")
  if(i==1 | i==121){axis(2, labels=TRUE, col="#88888888", col.axis="#88888888")}
  if(i>120){axis(1, labels=c("150","1,500","15,000","150,000"),
    at=c(5.01,7.31,9.61,11.91), col="#88888888", col.axis="#88888888")}
  points(dat[year == 1799 + i, gdp], dat[year == 1799 + i, life], pch=".",
    col="#88888888", cex = 2)
  points(eval(parse(text = paste0("gdps[['", i, "']]"))),
    eval(parse(text = paste0("results[['", i, "']]"))), type = "l",
    col = "black", lwd = 1)
  ylow <- eval(parse(text = paste0("bootstraps$`", i, "`")))[1,]
  yhigh <- eval(parse(text = paste0("bootstraps$`", i, "`")))[2,]
  polygon(x = c(sort(eval(parse(text = paste0("gdps[['", i, "']]")))),
    rev(sort(eval(parse(text = paste0("gdps[['", i, "']]")))))), c(sort(ylow),
    rev(sort(yhigh))), col = alpha("grey30", 0.25), border = NA)
  text(6,75, as.character(1809 + i))
}
mtext("Life Expectancy vs. GDP Across Generations", outer = TRUE, cex = 1.5)
