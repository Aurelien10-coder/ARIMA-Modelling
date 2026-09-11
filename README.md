# ARIMA Modelling of a French Industrial Production Index

This project studies the monthly Industrial Production Index (IPI) for the pharmaceutical industry in metropolitan France, published by INSEE. The dataset contains more than 400 monthly observations from January 1990 to February 2025.

## Overview

The analysis follows a standard time-series modelling workflow:

- Exploration of the original time series
- Stationarity analysis using Augmented Dickey-Fuller (ADF) and Phillips-Perron (PP) tests
- First differencing and log-differencing
- ACF/PACF analysis for ARMA identification
- ARMA model selection using AIC and BIC
- Residual diagnostics using Ljung-Box tests, ACF and Q-Q plots
- Estimation of an ARIMA(0,1,1) model
- Two-step-ahead forecasting with 95% confidence intervals
- Construction of a joint 95% confidence ellipse for the forecasts

## Methodology

The original series exhibits an upward trend and increasing variance. First differencing improves stationarity, but variance remains unstable. A logarithmic transformation followed by first differencing provides a more stable series.

The ADF and Phillips-Perron tests support the stationarity of the log-differenced series.

Several ARMA specifications are compared using information criteria and residual diagnostics. While AIC favors the more complex ARMA(3,3), BIC selects the more parsimonious ARMA(0,1) specification.

## Final Model

The final model is an ARIMA(0,1,1) applied to the logarithm of the original series:

$$
(1-B)\\log(X_t) = \\varepsilon_t - 0.6991\\varepsilon_{t-1}.
$$

This is equivalent to an ARMA(0,1) model for the log-differenced series.

## Forecasting

The model is used to produce two-step-ahead forecasts for March and April 2025, together with 95% confidence intervals.

A joint 95% confidence ellipse is also constructed for the two forecasts.

## Repository Structure

```text
.
├── README.md
├── ARIMA-Modelling.Rproj
├── ARIMA_modelling.R
├── data/
│   └── valeurs_mensuelles.csv
└── Rapport.pdf
```


## Requirements
The analysis was implemented in R using the following packages:
zoo
ggplot2
dplyr
fUnitRoots
tseries
tibble
ellipse

## Report
[![Report](https://img.shields.io/badge/Read%20Report-PDF-red?style=for-the-badge&logo=adobeacrobatreader)](Rapport.pdf)

The complete analysis, including the mathematical derivations, figures, diagnostics and forecasting results, is available in the link above.

## Author
Aurélien Tarroux
