library(asymmetric)
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
m <- as.matrix(dat[2:217])
dat <- melt(dat)
dates <- seq(from = 1810, to = 2015, by = 20)
dat$variable <- as.character(dat$variable)
variable <- do.call(rbind, strsplit(dat$variable, split =  "_"))[,1]
date <- do.call(rbind, strsplit(dat$variable, split =  "_"))[,2]
dat$variable <- variable
dat$date <- as.numeric(date)
dat <- dcast(dat, state + date ~ variable)
setDT(dat)
setnames(dat, c("state", "year", "x", "y"))
subset_hw_data <- function(date, dat)
{
  dta <- dat[year == date]
  dta[, state := NULL]
  dta[, year := NULL]
}
years_data <- lapply(dates, subset_hw_data, dat = dat)
analysis_dat <- years_data[[11]]

get_time_series_aocs <- function(DAT)
{
  analysis_dat <- DAT
  AH_AOC(dat = analysis_dat, method = c("QR"), input = analysis_dat$x,
    output = analysis_dat$y, print_style = "R", CI = TRUE, nboots = 5000)
}

ts_aoc_hw <- lapply(years_data, get_time_series_aocs)
ts_aoc_hw_table <- do.call(rbind, ts_aoc_hw)
ts_aoc_hw_table[,1] <- dates
colnames(ts_aoc_hw_table) <- c("year", "AOC", "2.5%", "50%", "97.5%")
ts_aoc_hw_table <- apply(ts_aoc_hw_table, 2, as.numeric)
xtable(ts_aoc_hw_table)

# let's get AOC over time
all_years <- seq(1800, 2010, 1)
all_years_data <- lapply(all_years, subset_hw_data, dat = dat)
all_years_aoc <- lapply(all_years_data, get_time_series_aocs)
aoc_data_all_years <- do.call(rbind, all_years_aoc)
aoc_data_all_years$year <- all_years
aoc_data_all_years$AOC <- as.numeric(aoc_data_all_years$AOC)
aoc_data_all_years$year <- as.numeric(aoc_data_all_years$year)
aoc_data_all_years$above_2 <- ifelse(aoc_data_all_years$AOC >= 0.2, 1, 0)
aoc_data_all_years$pre1950 <- ifelse(aoc_data_all_years$year < 1950, 1, 0)
setDT(aoc_data_all_years)
aoc_data_all_years[, sum(above_2)/184]
aoc_data_all_years[, mean(AOC), by = pre1950]
prop_table <- t(aoc_data_all_years[, sum(above_2)/184, by = pre1950])
prop_table[1,] <- c("Before 1950", "After 1950")
xtable(prop_table)
ggplot(aoc_data_all_years, aes(x = year, y = AOC)) + geom_point() +
  geom_smooth() + theme_bw() + geom_hline(yintercept = 0.2) +
  ggtitle("How Asymmetric is Health and Wealth from 1800 to 2010?")
