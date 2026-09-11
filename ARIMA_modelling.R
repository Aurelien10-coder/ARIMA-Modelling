# ===============================================================
# Time Series Project - ARIMA Modelling
# ===============================================================

# ===========================
# Load Required Libraries
# ===========================

require(zoo)
library(ggplot2)
library(dplyr)
library(fUnitRoots)
require(tseries)
library(tibble)
library(ellipse)
# ===========================
# Load and Explore the Dataset
# ===========================

# Read the CSV file (skip first 3 rows, use ";" as separator)
data <- read.csv("data/valeurs_mensuelles.csv",
                 sep = ";",
                 stringsAsFactors = FALSE,
                 skip = 3)
colnames(data) <- c("Date", "X", "Code")
data <- data[, 1:2]  # Keep only relevant columns

# Clean and transform the dataset
clean_data <- data %>%
  mutate(
    X = as.numeric(X),
    Date = as.yearmon(Date, "%Y-%m"),
    Date_num = as.numeric(format(Date, "%Y")) + (as.numeric(format(Date, "%m")) - 1) / 12
  )

# Plot the original time series
ggplot(clean_data, aes(x = Date_num, y = X)) +
  geom_line(color = "steelblue") +
  geom_point(color = "darkblue") +
  labs(title = "Original Series: X", x = "Date", y = "X") +
  theme_minimal()


# ===========================
# Stationarity Analysis - Differencing
# ===========================

# Ljung-Box test helper
Qtests <- function(series, k, fitdf = 0) {
  pvals <- sapply(1:k, function(l) {
    if (l <= fitdf) return(NA)
    Box.test(series, lag = l, type = "Ljung-Box", fitdf = fitdf)$p.value
  })
  return(cbind(lag = 1:k, pval = pvals))
}

# Augmented Dickey-Fuller test with autocorrelation check
adfTest_valid <- function(series, kmax, adftype) {
  k <- 0
  while (TRUE) {
    cat(sprintf("ADF with %d lags: residuals OK? ", k))
    adf <- adfTest(series, lags = k, type = adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[, "pval"]
    if (all(pvals >= 0.05, na.rm = TRUE)) {
      cat("OK\n")
      break
    } else {
      cat("Not OK\n")
      k <- k + 1
    }
  }
  return(adf)
}

# Phillips-Perron test wrapper
do_pp_test <- function(series) {
  pp.test(series)
}

# Original series tests
cat("\n===\nOriginal Series\n===\n")

model_orig <- lm(X ~ Date_num, data = clean_data)
print(summary(model_orig))

spread_orig <- zoo(clean_data$X, order.by = clean_data$Date)
adf_orig <- adfTest_valid(spread_orig, 40, adftype = "ct")
print(adf_orig)

cat("Phillips-Perron Test:\n")
print(do_pp_test(spread_orig))

# First differencing
data_diff <- clean_data %>%
  arrange(Date) %>%
  mutate(
    X_diff = c(NA, diff(X)),
    Date_num = as.numeric(format(Date, "%Y")) + (as.numeric(format(Date, "%m")) - 1) / 12
  ) %>%
  filter(!is.na(X_diff)) %>%
  dplyr::select(Date_num, X_diff, Date)



# Plot differenced series
ggplot(data_diff, aes(x = Date_num, y = X_diff)) +
  geom_line(color = "steelblue") +
  geom_point(color = "darkblue") +
  labs(title = "Differenced Series: X_diff", x = "Date", y = "X_diff") +
  theme_minimal()

model_diff <- lm(X_diff ~ Date_num, data = data_diff)
print(summary(model_diff))

spread_diff <- zoo(data_diff$X_diff, order.by = data_diff$Date)
adf_diff <- adfTest_valid(spread_diff, 40, adftype = "ct")
print(adf_diff)

cat("Phillips-Perron Test:\n")
print(do_pp_test(spread_diff))

# Check variance stability across thirds
n <- nrow(data_diff)
v1 <- var(data_diff$X_diff[1:(n/3)])
v2 <- var(data_diff$X_diff[(n/3 + 1):(2*n/3)])
v3 <- var(data_diff$X_diff[(2*n/3 + 1):n])
cat("Variance - First Third:", v1, "\n")
cat("Variance - Second Third:", v2, "\n")
cat("Variance - Third Third:", v3, "\n")

# Log + differencing
data_log_diff <- clean_data %>%
  arrange(Date) %>%
  mutate(
    X_log = log(X),
    X_log_diff = c(NA, diff(log(X)))
  ) %>%
  filter(!is.na(X_log_diff)) %>%
  dplyr::select(Date_num, X_log_diff, Date)



ggplot(data_log_diff, aes(x = Date_num, y = X_log_diff)) +
  geom_line(color = "steelblue") +
  geom_point(color = "darkblue") +
  labs(title = "Log-Differenced Series: X_log_diff", x = "Date", y = "X_log_diff") +
  theme_minimal()

model_log_diff <- lm(X_log_diff ~ Date_num, data = data_log_diff)
print(summary(model_log_diff))

spread_log_diff <- zoo(data_log_diff$X_log_diff, order.by = data_log_diff$Date)
adf_log_diff <- adfTest_valid(spread_log_diff, 40, adftype = "ct")
print(adf_log_diff)

cat("Phillips-Perron Test:\n")
print(do_pp_test(spread_log_diff))

# Variance check for log-differenced series
n <- nrow(data_log_diff)
cat("Variance - First Third:", var(data_log_diff$X_log_diff[1:(n/3)]), "\n")
cat("Variance - Second Third:", var(data_log_diff$X_log_diff[(n/3 + 1):(2*n/3)]), "\n")
cat("Variance - Third Third:", var(data_log_diff$X_log_diff[(2*n/3 + 1):n]), "\n")


# ===========================
# Visualization
# ===========================

ggplot(clean_data, aes(x = Date_num, y = X)) +
  geom_line(color = "steelblue") +
  geom_point(color = "darkblue") +
  labs(title = "Original Series: X", x = "Date", y = "X") +
  theme_minimal()

ggplot(data_log_diff, aes(x = Date_num, y = X_log_diff)) +
  geom_line(color = "steelblue") +
  geom_point(color = "darkblue") +
  labs(title = "Log-Differenced Series", x = "Date", y = "X_log_diff") +
  theme_minimal()


# ===========================
# ARMA Modelling
# ===========================

# ACF & PACF plots
y <- data_log_diff$X_log_diff
par(mfrow = c(1, 2))
acf(y, lag.max = 24, main = "ACF - log(diff(X))")
pacf(y, lag.max = 24, main = "PACF - log(diff(X))")

# Model selection via AIC & BIC
pmax <- 3; qmax <- 3
AICs <- BICs <- matrix(NA, nrow = pmax + 1, ncol = qmax + 1)
rownames(AICs) <- colnames(AICs) <- paste0("p=", 0:pmax)
colnames(BICs) <- paste0("q=", 0:qmax)

for (p in 0:pmax) {
  for (q in 0:qmax) {
    model <- try(arima(y, order = c(p, 0, q)), silent = TRUE)
    if (!inherits(model, "try-error")) {
      AICs[p+1, q+1] <- AIC(model)
      BICs[p+1, q+1] <- BIC(model)
    }
  }
}
print("AIC Table:"); print(AICs)
print("BIC Table:"); print(BICs)
cat("\n est AIC → ARMA(3,3) | Best BIC → ARMA(0,1)\n")

# Fit candidate models
models <- list(
  arma01 = arima(y, order = c(0, 0, 1)),
  arma10 = arima(y, order = c(1, 0, 0)),
  arma11 = arima(y, order = c(1, 0, 1)),
  arma22 = arima(y, order = c(2, 0, 2)),
  arma33 = arima(y, order = c(3, 0, 3))
)

# Model diagnostic function
significance <- function(model) {
  coefs <- model$coef
  se <- sqrt(diag(model$var.coef))
  t_vals <- coefs / se
  pvals <- 2 * (1 - pnorm(abs(t_vals)))
  return(rbind(coef = coefs, std_err = se, p_val = pvals))
}

arima_diagnostics <- function(name, model) {
  cat("\n=== Model:", name, "===\n")
  print(significance(model))
  cat("\n-- Ljung-Box Test (lag=12) --\n")
  print(Box.test(model$residuals, lag = 12, type = "Ljung-Box", fitdf = length(model$coef)))
  cat("\n-- Ljung-Box Test (lags 1 to 24) --\n")
  print(round(Qtests(model$residuals, 24, fitdf = length(model$coef)), 3))
}

# Run diagnostics
for (name in names(models)) {
  arima_diagnostics(name, models[[name]])
}


# ===========================
# Final ARIMA Model
# ===========================

# We select an ARIMA(0,1,1) model on the log-transformed series, equivalent to ARIMA(0,0,1) on the log-differenced series.
# This model, chosen using BIC, is parsimonious with a significant MA(1) coefficient and acceptable residuals.

# Residual diagnostics
arma_model <- arima(y, order = c(0, 0, 1))
resid <- residuals(arma_model)
std_resid <- resid / sqrt(arma_model$sigma2)

par(mfrow = c(2, 3))
plot(std_resid, type = "l", main = "Standardized Residuals")
abline(h = 0, col = "red", lty = 2)
acf(std_resid, main = "ACF of Residuals")
qqnorm(std_resid); qqline(std_resid, col = "red")
hist(std_resid, breaks = 20, col = "lightblue", freq = FALSE, main = "Residual Histogram")
curve(dnorm(x), col = "red", lwd = 2, add = TRUE)
lb_pvalues <- sapply(1:20, function(lag) Box.test(std_resid, lag = lag, type = "Ljung-Box", fitdf = length(arma_model$coef))$p.value)
plot(1:20, lb_pvalues, type = "b", main = "Ljung-Box Test P-values")
abline(h = 0.05, col = "blue", lty = 2)
par(mfrow = c(1, 1))


# ===========================
# Forecasting
# ===========================

h <- 2  # Forecast horizon
z <- qnorm(0.975)  # 95% CI

forecast <- predict(arma_model, n.ahead = h)


# Rebuild log(X)
logX_base <- cumsum(y) + log(clean_data$X[1])
logX_last <- tail(logX_base, 1)

# Forecast in log space
logX_forecast <- cumsum(forecast$pred) + logX_last
logX_upper <- cumsum(forecast$pred + z * forecast$se) + logX_last
logX_lower <- cumsum(forecast$pred - z * forecast$se) + logX_last

mu <- c(logX_forecast[1], logX_forecast[2])
exp(mu)

# Convert back to X
X_pred <- exp(logX_forecast)
X_upper <- exp(logX_upper)
X_lower <- exp(logX_lower)

# Build date vectors
dates_hist <- tail(data_log_diff$Date, 100)
last_date <- tail(dates_hist, 1)
future_dates <- seq(last_date + 1/12, by = 1/12, length.out = h)

# Create dataframes
X_hist <- exp(tail(logX_base, 100))
df_hist <- tibble(Date = dates_hist, X = X_hist, Type = "Historical")
df_forecast <- tibble(Date = future_dates, X = X_pred, Upper = X_upper, Lower = X_lower, Type = "Forecast")
df_forecast <- bind_rows(
  df_hist %>% tail(1) %>% mutate(Type = "Forecast"),
  df_forecast
)
df_all <- bind_rows(df_hist, df_forecast)

# Plot forecasts
ggplot(df_all, aes(x = Date, y = X)) +
  geom_line(data = df_hist, aes(color = Type), size = 1.2) +
  geom_line(data = df_forecast, color = "red", size = 1.2) +
  geom_ribbon(data = df_forecast %>% filter(!is.na(Upper)),
              aes(ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.4) +
  scale_color_manual(values = c("Historical" = "steelblue")) +
  labs(title = "Forecast of Series X with ARMA(0,1)",
       x = "Date", y = "Value X") +
  theme_minimal(base_size = 14)




# Plot of the ellipse
cov_matrix <- matrix(c(forecast$se[1]^2, 0,
                       0, forecast$se[2]^2), 
                     nrow = 2)


ellipse_data <- as.data.frame(ellipse(cov_matrix, centre = mu, level = 0.95))
colnames(ellipse_data) <- c("x", "y")

ellipse_data <- ellipse_data %>%
  mutate(x = exp(x), y = exp(y))

center_point <- data.frame(x = X_pred[1], y = X_pred[2])

ggplot() +
  geom_path(data = ellipse_data, aes(x = x, y = y), color = "black") +
  geom_point(data = center_point, aes(x = x, y = y), shape = 18, size = 3) +
  labs(title = "95% Joint Confidence Ellipse (Forecasts T+1 & T+2)",
       x = "Forecast for March 2025 (T+1)",
       y = "Forecast for April 2025 (T+2)") +
  theme_minimal(base_size = 14)
