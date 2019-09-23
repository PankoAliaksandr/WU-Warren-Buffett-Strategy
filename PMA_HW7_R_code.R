# Libraries
library(quantmod)
library(PerformanceAnalytics)
library(xtable)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(psych)
library(lubridate)
library(dplyr)
#library(tseries)


#Set working directory
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Downloading data and packages
load("PerfMBerkshireData.RData")

# Evaluating data structure and values


# Check missing values
nacols <- function(df) {
  colnames(df)[unlist(lapply(df, function(x) any(is.na(x))))]
} #Extracting the columns with missing data
nacols(DT)

# Calculating excess return for BRK.ret
DT$BRK.ret <- DT$BRK.ret - DT$RF
colnames(DT)[ncol(DT)] <- "BRK.ret_e"

#Plot each column of our data:
#df <- melt(DT,  id.vars = 'date', variable.name = 'series')

# plot on same grid, each series colored differently -- 
# good if the series have same scale

# or plot on different plots
#Plot <- ggplot(df, aes(date,value)) + geom_line() + facet_grid(series ~ .)

#Calculating the sd of each column
#DT_sd <- as.data.frame(colSds(as.matrix(subset(DT, select = -date))))
#row.names(DT_sd) <- colnames(subset(DT, select = -date))
#colnames(DT_sd) <- c("sd_values")
#DT_sd[order(DT_sd$sd_values),c(1,1)]

# Descriptive statistics

# Find the starting date of a period
x <- as.Date(range(DT$date)[1])

# Split the whole period into 10 years blocks
DT_10 <- DT[DT$date >= x & DT$date <= x %m+% years(10),]
DT_20 <- DT[DT$date >= x %m+% years(10) & DT$date <= x %m+% years(20),]
DT_30 <- DT[DT$date >= x %m+% years(20) & DT$date <= x %m+% years(30),]

# Convert to timeseries
DT_10 <- xts(DT_10[,-8], order.by=as.Date(DT_10[,8], "%m/%d/%Y"))
DT_20 <- xts(DT_20[,-8], order.by=as.Date(DT_20[,8], "%m/%d/%Y"))
DT_30 <- xts(DT_30[,-8], order.by=as.Date(DT_30[,8], "%m/%d/%Y"))
DT <- xts(DT[,-8], order.by=as.Date(DT[,8], "%m/%d/%Y"))

# Get base statistics for each period
des_10 <- describe(subset(DT_10, select = BRK.ret_e));
des_10 <- subset(des_10, select=-c(trimmed,mad))

des_15 <- describe(subset(DT_15, select = BRK.ret_e));
des_15 <- subset(des_15, select=-c(trimmed,mad))

des_20 <- describe(subset(DT_20, select = BRK.ret_e));
des_20 <- subset(des_20, select=-c(trimmed,mad)) 

des_25 <- describe(subset(DT_25, select = BRK.ret_e));
des_25 <- subset(des_25, select=-c(trimmed,mad))

des_30 <- describe(subset(DT_30, select = BRK.ret_e));
des_30 <- subset(des_30, select=-c(trimmed,mad))

des_full <- describe(subset(DT, select = BRK.ret_e));
des_full <- subset(des_full, select=-c(trimmed,mad)) 

# Descriptive statistics table
table_BRK <- rbind(des_10, des_20, des_30, des_full)
row.names(table_BRK) <- c("1988-1998", "1999-2008","2009-2018", "full")
print(table_BRK, digits = 3)
xtable(table_BRK, digits = 3)

# Sharpe Ratio
SR_10 <- SharpeRatio.annualized(DT_10$BRK.ret_e, Rf = 0, scale = 12, geometric = TRUE)
SR_20 <- SharpeRatio.annualized(DT_20$BRK.ret_e, Rf = 0, scale = 12, geometric = TRUE)
SR_30 <- SharpeRatio.annualized(DT_30$BRK.ret_e, Rf = 0, scale = 12, geometric = TRUE)
SR_full <- SharpeRatio.annualized(DT$BRK.ret_e, Rf = 0, scale = 12, geometric = TRUE)
SR_all <- c(SR_10, SR_20, SR_30, SR_full)

# Treynor Ratio
TR_10 <- TreynorRatio(DT_10$BRK.ret_e, DT_10$Mkt, Rf = 0, scale = 12, modified = FALSE)
TR_20 <- TreynorRatio(DT_20$BRK.ret_e, DT_20$Mkt, Rf = 0, scale = 12, modified = FALSE)
TR_30 <- TreynorRatio(DT_30$BRK.ret_e, DT_30$Mkt, Rf = 0, scale = 12, modified = FALSE)
TR_full <- TreynorRatio(DT$BRK.ret_e, DT$Mkt, Rf = 0, scale = 12, modified = FALSE)
TR_all <- c(TR_10, TR_20, TR_30, TR_full)

# incormation ratio
IR_10 <- InformationRatio(DT_10$BRK.ret_e, DT_10$Mkt, scale = 12)
IR_20 <- InformationRatio(DT_20$BRK.ret_e, DT_20$Mkt, scale = 12)
IR_30 <- InformationRatio(DT_30$BRK.ret_e, DT_30$Mkt, scale = 12)
IR_full <- InformationRatio(DT$BRK.ret_e, DT$Mkt, scale = 12)
IR_all <- c(IR_10, IR_20, IR_30, IR_full)

#CAPM
JA_10 <- CAPM.jensenAlpha(DT_10$BRK.ret_e, DT_10$Mkt, Rf = 0)
JA_20 <- CAPM.jensenAlpha(DT_20$BRK.ret_e, DT_20$Mkt, Rf = 0)
JA_30 <- CAPM.jensenAlpha(DT_30$BRK.ret_e, DT_30$Mkt, Rf = 0)
JA_full <- CAPM.jensenAlpha(DT$BRK.ret_e, DT$Mkt, Rf = 0)
JA_all <- c(JA_10, JA_20, JA_30, JA_full)

# Multifactor 5 factors

model_10 <- lm(DT_10$BRK.ret_e ~ DT_10$Mkt + DT_10$SMB + DT_10$HML + DT_10$RMW + DT_10$CMA)
coeff_10 <- summary(model_10)$coefficients
tvalues_10 <- coeff_10[,3]
# annualized
alpha_10 <- coeff_10[1][1] * 12

model_20 <- lm(DT_20$BRK.ret_e ~ DT_20$Mkt+ DT_20$SMB+ DT_20$HML + DT_20$RMW + DT_20$CMA)
coeff_20<-summary(model_20)$coefficients
tvalues_20 <- coeff_20[,3]
# annualized
alpha_20 <- coeff_20[1][1] * 12

model_30 <- lm(DT_30$BRK.ret_e~DT_30$Mkt+ DT_30$SMB+ DT_30$HML + DT_30$RMW + DT_30$CMA)
coeff_30<-summary(model_30)$coefficients
tvalues_30 <- coeff_30[,3]
# annualized
alpha_30 <- coeff_30[1][1] * 12

model_full <- lm(DT$BRK.ret_e ~ DT$Mkt+ DT$SMB+ DT$HML + DT$RMW + DT$CMA)
coeff_full <- summary(model_full)$coefficients
tvalues_full <- coeff_full[,3]
# annualized
alpha_full <- coeff_full[1][1] * 12

alphas_all <- c(alpha_10, alpha_20,alpha_30, alpha_full)

summary_table <- cbind(SR_all, TR_all, IR_all, JA_all, alphas_all)
colnames(summary_table) <- c("SR", "TR", "IR", "JensenAlpha","MultifAlpha")
row.names(summary_table) <- c("1988-1998", "1999-2008","2009-2018", "full")
xtable(summary_table, digits = 3)

coeffs_all <- rbind(coeff_10, coeff_20, coeff_30, coeff_full)
tvalues_all <- rbind(tvalues_10, tvalues_20, tvalues_30, tvalues_full)
row.names(tvalues_all) <- c("1988-1998", "1999-2008","2009-2018", "full")
xtable(tvalues_all, digits = 3)


