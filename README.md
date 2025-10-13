---
title: "betaARMA: A Package for Beta Autoregressive Moving Average Models"
author: "Everton da Costa"
output: github_document
---

# betaARMA

[![Status](https://img.shields.io/badge/Status-In_Development-blue.svg)](https://github.com/Everton-da-Costa/betaARMA)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An R package for fitting, forecasting, and simulating Beta Autoregressive Moving Average $(\beta\text{ARMA})$ models. This package provides a comprehensive and user-friendly toolkit for modeling time series data bounded on the (0, 1) interval, such as rates, proportions, and indices.

---

## 📚 Table of Contents

- [🎯 Project Motivation](#-project-motivation)
- [✨ Planned Features](#-planned-features)
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

This project aims to create the `betaARMA` package as a go-to resource for researchers and practitioners working with bounded time series data. The focus is on a clean interface, and strong documentation.

---

## ✨ Planned Features

* **Unified Model Fitting:** A single core function, `betaARMA()`, for fitting AR, MA, and ARMA models with support for regressors (`xreg`).
* **Object-Oriented Design:** A clean S3 class system, allowing for intuitive use of standard R generics like `predict()`, `summary()`, `plot()`, and `simulate()`.
* **Forecasting Engine:** A powerful `predict()` method to generate multi-step-ahead forecasts.
* **Simulation Tools:** A `simulate()` method to generate sample paths from a fitted model for analysis and testing.
* **Model Diagnostics:** Built-in functions for residual analysis and model validation, including Portmanteau tests.

---

## 🗺️ Development Roadmap

This is the development plan for the `betaARMA` package.

### Phase 1: Research, Architecture, and Setup (Deadline: October 31, 2025)
- [ ] **Research:** Analyze reference packages (e.g., `btsr`, `stats::arima`) to find the most efficient and stable way to implement the **recursion of the dynamic component**.
- [ ] **Architecture:** Define the use of the S3 object system for model objects (class `"betaARMA"`).
- [ ] **Optimization:** Select and test the optimization algorithm (e.g., `stats::optim` with the `L-BFGS-B`, `lbfgs` method).
- [ ] **Setup:** Create the package skeleton and initialize version control with Git.

### Phase 2: Core Model Implementation (Deadline: November 14, 2025)
- [ ] **Main Function:** Develop `betaARMA()` to unify AR, MA, and ARMA model fitting (without regressors initially).
- [ ] **Regressors:** Implement support for static regressors (`xreg`).
- [ ] **S3 Object:** Structure the `betaARMA` class with a standardized list of outputs (coefficients, residuals, vcov, etc.).
- [ ] **Basic Methods:** Create the essential S3 methods: `print.betaARMA()`, `summary.betaARMA()`, `coef.betaARMA()`, and `fitted.betaARMA()`.

### Phase 3: Essential Functionality (Deadline: November 28, 2025)
- [ ] **Add Regressor Support:** Enhance the `betaARMA()` function to support static regressors via an `xreg` argument.
- [ ] **Forecasting:** Implement the `predict.betaARMA()` method with support for `n.ahead` and `newxreg`.
- [ ] **Simulation:** Create the `simulate.betaARMA()` method to generate sample paths from a fitted model.
- [ ] **Diagnostics:** Develop a `check_residuals()` function with options for Portmanteau tests (Ljung-Box, Monti, `Q_4`).

### Phase 4: Documentation & Polishing (Deadline: December 12, 2025)
- [ ] **Datasets:** Add and document the seasonal and non-seasonal datasets.
- [ ] **Help Pages:** Document all exported functions and datasets using `roxygen2`.
- [ ] **Vignette:** Write a tutorial (package vignette) demonstrating a complete workflow.
- [ ] **Review:** Conduct a final review of all code and documentation before tagging a "version 1.0".

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
                        dependencies = TRUE,
                        build_vignettes = TRUE)
```

---

## 🚀 Getting Started

Once installed, the best way to get started will be through the package vignette, which will provide a detailed, narrated code example.

```R
# This command will work once the first vignette is complete
vignette("intro_betaARMA", package = "betaARMA")
```

---

## 📂 Repository Structure

The repository is structured as a standard R package for clarity and reproducibility.

```plaintext
.
├── R/                  # Source code for all R functions.
├── data/               # Processed data included in the package (.rda).
├── data-raw/           # Raw data and scripts used to process it.
├── man/                # R package documentation files for functions.
├── vignettes/          # Detailed tutorial and case study (.Rmd).
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace.
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