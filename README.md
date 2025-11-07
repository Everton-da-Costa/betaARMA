# betaARMA

[![Status](https://img.shields.io/badge/Status-In_Development-blue.svg)](https://github.com/Everton-da-Costa/betaARMA)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

An R package for fitting, forecasting, and simulating Beta Autoregressive Moving Average $(\beta\text{ARMA})$ models. This package provides a comprehensive and user-friendly toolkit for modeling time series data bounded on the (0, 1) interval, such as rates, proportions, and indices.

---

## 📚 Table of Contents

- [🎯 Project Motivation](#-project-motivation)
- [✨ Core Features](#-core-features)
- [🗺️ Development Roadmap](#️-development-roadmap)
- [🛠️ Installation](#️-installation)
- [🚀 Getting Started](#-getting-started)
- [📂 Repository Structure](#-repository-structure)
- [🎓 Citation](#-citation)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📬 Contact](#-contact)

---

## 🎯 Project Motivation

The Beta Autoregressive Moving Average $(\beta\text{ARMA})$ model is a powerful tool for analyzing time series data bounded between 0 and 1. While foundational models exist, there is a need for a unified R package that simplifies the entire modeling workflow—from fitting flexible AR, MA, and ARMA structures to performing diagnostics, forecasting, and simulation.

This project aims to create the `betaARMA` package as a go-to resource for researchers and practitioners working with bounded time series data. The focus is on a clean interface, robust implementation, and strong documentation.

---

## ✨ Core Features

* **Unified Model Fitting:** A single core function, `barma()`, now handles $\beta$AR, $\beta$MA, and $\beta$ARMA models through a single, clean interface.
* **Object-Oriented Design:** The package uses a modern S3 class system. The `barma()` output object works directly with standard R generics:
    * `print()` for a concise model overview.
    * `summary()` for a detailed table of coefficients, std. errors, p-values, and information criteria.
    * `coef()` to extract the coefficient vector.
    * `fitted()` to extract the NA-padded fitted values as a `ts` object.
    * `residuals()` to calculate and extract standardized residuals.
* **Forecasting Engine:** A `forecast()` method is implemented to generate dynamic, multi-step-ahead point forecasts from a fitted model.
* **Simulation Tools:** Includes a `simu_barma()` function to generate time series from known $\beta$ARMA processes for testing and validation.

---

## 🗺️ Development Roadmap

This is the development plan for the `betaARMA` package.

### Phase 1: Research, Architecture, and Setup (Deadline: October 31, 2025)
- [ ] **Research:** Analyze reference packages (e.g., `btsr`, `arima2::arima`, `forecast::Arima`, `stats::arima`) to find efficient and stable way of implementation.
- [ ] **Continuous Integration:** Set up GitHub Actions to run `R-CMD-check`. (User will implement this later).
- [X] **Architecture:** Define the use of the S3 object system for model objects (class `"barma"`).
- [X] **Optimization:** Select and test the optimization algorithm (`stats::optim` with the `BFGS` method).
- [X] **Setup:** Create the package skeleton and initialize version control with Git.

### Phase 2: Core Model Implementation (Deadline: November 14, 2025)
- [X] **Main Function:** Develop `barma()` to unify AR, MA, and ARMA model fitting (without regressors).
- [ ] **Regressors:** Implement support for static regressors (`xreg`).
- [X] **S3 Object:** Structure the `barma` class with a standardized list of outputs.
- [X] **Basic Methods:** Create the essential S3 methods: `print()`, `summary()`, `coef()`, and `fitted()`.

### Phase 3: Essential Functionality (Deadline: November 28, 2025)
- [ ] **Add Regressor Support:** Enhance the `barma()` function to support static regressors via an `xreg` argument.
- [X] **Forecasting:** Implement the `forecast.barma()` method.
- [X] **Simulation:** Create the `simu_barma()` function.
- [X] **Residuals:** Implement the `residuals.barma()` method.
- [ ] **Diagnostics:** Develop a `plot.barma()` method for residual analysis.

### Phase 4: Documentation & Polishing (Deadline: December 12, 2025)
- [ ] **Datasets:** Add and document the seasonal and non-seasonal datasets.
- [ ] **Help Pages:** Finalize documentation for all exported functions.
- [ ] **Vignette:** Write a complete tutorial (package vignette) demonstrating a full workflow.
- [ ] **Review:** Conduct a final review of all code and documentation.

---

## 🛠️ Installation
Once the first version is stable, the package will be installable directly from GitHub.

First, ensure you have the `remotes` package:
```R
if (!require("remotes")) {
  install.packages("remotes")
}
```

Then, install the package from GitHub (note: this link will be active once the repository is public):
```R
remotes::install_github("everton-da-costa/betaARMA", 
                        dependencies = TRUE)
```

---

## 🚀 Getting Started

Once installed, the best way to get started will be through the package vignette, which will provide a detailed, narrated code example.

```R
library(betaARMA)

# 1. Simulate some data
set.seed(123)
y <- simu_barma(n = 100, ar = 1, varphi = 0.5, phi = 20)

# 2. Fit a model
fit <- barma(y, ar = 1)

# 3. Get a detailed summary
summary(fit)

# 4. Get 10-step-ahead forecasts
forecast_h10 <- forecast(fit, h = 10)
print(forecast_h10)
```

---

## 📂 Repository Structure

The repository is structured as a standard R package for clarity and reproducibility.

```plaintext
.
├── R/                  # Source code for all R functions.
├── man/                # R package documentation files (generated by roxygen2).
├── validation/         # Scripts for testing and validation.
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace (generated by roxygen2).
├── LICENSE             # MIT License file.
└── README.md           # This file.
```

---

## 🎓 Citation

Once the package is developed, you will be able to get citation information by running the following command in R:
```R
citation("betaARMA")
```

---

## 🤝 Contributing
Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

## 📄 License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## 📬 Contact
For questions, suggestions, or issues related to the code, please contact:

Everton da Costa
📧 everto.cost@gmail.com