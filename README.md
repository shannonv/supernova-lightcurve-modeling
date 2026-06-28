# 🌟 Supernova Light Curve & Radio Emission Modeling

> **Python models for optical and radio emission from Type IIb supernovae, including shock-cooling light-curve modeling, synthetic photometry, Monte Carlo uncertainty analysis, parameter-space grid searches, and synchrotron self-absorption modeling.**

This repository contains a collection of Python scripts developed during my graduate research in astronomy and astrophysics to model the early optical and radio emission from Type IIb supernovae.

The project combines semi-analytic physical models with observational data to constrain progenitor star properties and investigate the physical processes governing supernova explosions. The analysis emphasizes reproducible scientific computing, uncertainty quantification, and physically motivated model fitting.

---

# Project Objectives

This project investigates several scientific questions:

* How do the physical properties of a progenitor star influence the observed optical light curve?
* What combinations of envelope mass and radius best reproduce the observed data?
* How do observational uncertainties propagate into the inferred stellar properties?
* Can synchrotron self-absorption models reproduce observed radio light curves and spectra?

---

# Tools & Technologies

| Category             | Tools                                                                            |
| -------------------- | -------------------------------------------------------------------------------- |
| Programming          | Python                                                                           |
| Scientific Computing | NumPy, SciPy                                                                     |
| Astronomy            | Astropy                                                                          |
| Visualization        | Matplotlib                                                                       |
| Data Formats         | JSON, Pickle                                                                     |
| Methods              | Monte Carlo simulation, Grid Search, Synthetic Photometry, Numerical Integration |

---

# Repository Structure

```text
.
├── optical_LC_modeling.py
├── optical_grid_search_16gkg_MC.py
├── radio_SSA_model.py
├── README.md
└── requirements.txt
```

---

# Analysis Workflow

```text
Observed Supernova Photometry
            │
            ▼
Data Preparation
            │
            ▼
Shock-Cooling Light Curve Model
            │
            ▼
Synthetic Photometry
            │
            ▼
Monte Carlo Sampling
            │
            ▼
Grid Search Parameter Estimation
            │
            ▼
Best-Fit Stellar Parameters
```

---

# Repository Contents

## Optical Light Curve Modeling

Implements a semi-analytic shock-cooling model (Piro 2015) to calculate:

* Luminosity evolution
* Photospheric temperature
* Photospheric radius
* Synthetic optical magnitudes

Model parameters include:

* Core mass
* Envelope mass
* Envelope radius
* Explosion energy
* Opacity

---

## Monte Carlo Parameter Estimation

Observational uncertainties are propagated by repeatedly sampling photometric measurements from their error distributions.

A multidimensional grid search is then performed over:

* Envelope radius
* Envelope mass
* Explosion time offset

to identify physically plausible solutions and estimate parameter uncertainties.

---

## Radio Emission Modeling

Implements synchrotron self-absorption (SSA) and synchrotron self-absorption with free-free absorption (SSA+FFA) models to reproduce observed radio light curves and spectra.

Additional helper functions calculate:

* Stellar radius
* Escape velocity
* Luminosity scaling
* Distance uncertainty propagation

---

# Methods Demonstrated

* Scientific computing
* Numerical modeling
* Monte Carlo uncertainty propagation
* Parameter estimation
* Grid search optimization
* Numerical integration
* Synthetic photometry
* Data visualization
* Computational astrophysics

---

# Scientific Context

These scripts were developed as part of graduate research on stripped-envelope (Type IIb) supernovae.

The methods are representative of computational workflows used in observational astrophysics, where physical models are compared with multi-wavelength observations to infer stellar and explosion properties.

---

# Skills Demonstrated

* Python programming
* Scientific software development
* Numerical simulation
* Model fitting
* Error propagation
* Statistical analysis
* Data visualization
* Reproducible research
* Technical documentation

---

# Future Improvements

Potential future enhancements include:

* Bayesian parameter estimation (MCMC)
* Interactive visualization
* Modular object-oriented implementation
* Automated model comparison
* Jupyter notebook examples

---

# About Me

**Shannon Vanderwoude**

PhD in Astronomy & Astrophysics
University of Toronto

My interests include scientific computing, machine learning, statistical modeling, data visualization, and applying quantitative methods to problems in astronomy, climate science, healthcare, and public-sector analytics.

