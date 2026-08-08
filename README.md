# Some Advances in Progressive Type-II Censoring Scheme Under Weibull Lifetime Model

[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the R scripts, algorithmic implementations, simulation engines, and empirical data applications for the paper **"Some Advances in Progressive Type-II Censoring Scheme Under Weibull Lifetime Model"** by Okechukwu J. Obulezi (Nnamdi Azikiwe University, Awka, Nigeria).

---

## 📌 Overview

Progressive Type-II Censoring (PC-T2) is a fundamental framework in reliability engineering and survival analysis. However, classical progressive procedures under Weibull baseline distributions face challenges when encountering dynamic environments, non-linear degradation, high-dimensional dependent competing risks, and model misspecification. 

This repository provides R implementations across **four novel methodological frontiers**:

1. **DRL-AP2** *(Dynamic Reinforcement Learning Adaptive PC-T2)*: Formulates progressive unit removal choices as a Markov Decision Process (MDP) balancing information gain against operational budgets.
2. **PINN-ProgCens** *(Physics-Informed Neural ODE)*: Integrates physical ODE degradation trajectories within a progressively censored Weibull likelihood.
3. **Prog-Vine** *(Progressive Regularized Vine Copula)*: Decomposes high-dimensional dependent competing risks subject to cause masking using C-Vine/D-Vine copulas and EM/Data Augmentation algorithms.
4. **DRO-PC** *(Distributionally Robust Progressive Optimization)*: Uses Wasserstein ambiguity balls and optimal transport theory to derive robust inference bounds resilient to out-of-distribution contamination.

---

## 📁 Repository Structure

```text
├── R/
│   ├── 01_drl_ap2.R          # MDP state observer, MLE engine, and Metropolis-within-Gibbs sampler
│   ├── 02_pinn_progcens.R    # Continuous Adjoint sensitivity solver & SGLD Bayesian engine
│   ├── 03_prog_vine.R        # EM-Vine MLE solver & Data Augmentation Gibbs sampler for masked risks
│   └── 04_dro_pc.R           # Primal-Dual convex solver & Wasserstein extreme-point sampler
├── data/
│   └── jet_engine_data.csv   # Industrial lifetime / degradation benchmark datasets
├── simulations/
│   ├── run_drl_sims.R        # Monte Carlo simulation pipeline for DRL-AP2
│   ├── run_pinn_sims.R       # Monte Carlo simulation pipeline for PINN-ProgCens
│   ├── run_vine_sims.R       # Monte Carlo simulation pipeline for Prog-Vine
│   └── run_dro_sims.R        # Monte Carlo simulation pipeline for DRO-PC
├── README.md                 # Project description and guide
└── LICENSE                   # MIT License
