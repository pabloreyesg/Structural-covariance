# Structural Covariance Network (SCN) Analysis App

This Shiny application provides a user-friendly interface to compute and visualize **Structural Covariance Networks (SCNs)** from neuroimaging-derived regional data. The tool enables neuroscientists and cognitive researchers to explore inter-regional covariance patterns across subject groups, compare network properties, and perform statistical inference.

## ✨ Features

- Upload custom `.csv` files with morphometric or imaging-derived metrics.
- Compute SCNs based on:
  - **Pearson correlation**
  - **Spearman correlation**
  - **Covariance**
- Adjust for covariates (e.g., age, sex) using linear regression residuals.
- Visualize:
  - SCN heatmaps
  - Graph-theoretical networks (interactive)
  - Network-level metrics (density, clustering, path length)
  - Group comparison statistics (parametric, non-parametric, Bayesian)
- Flexible graph layouts (`circle`, `Fruchterman-Reingold`, `Kamada-Kawai`, `tree`, etc.)

## 📂 Input Format

The app expects a `.csv` file with the following structure:

| ID | Group | Edad | Region1 | Region2 | ... |
|----|-------|------|---------|---------|-----|
| 01 | MCI   | 72   | 2.31    | 3.45    | ... |
| 02 | Control | 66 | 2.12    | 3.51    | ... |

- **Group column**: Indicates group membership (e.g., `Control`, `MCI`, `AD`).
- **Covariate column**: Optional numeric variable (e.g., `Edad`).
- **Regional columns**: Brain regions or other neurobiological measurements.

## 🧠 Network Construction

SCNs are computed within each group by:
1. Extracting residuals (if covariate is selected).
2. Computing inter-regional Pearson/Spearman/covariance matrices.
3. Thresholding absolute values below a user-defined threshold.
4. Generating undirected weighted graphs from adjacency matrices.
5. Visualizing and computing graph-theoretical measures.

## 📊 Statistical Comparisons

You can compare SCNs between any two groups using:

- **t-test**: For normally distributed differences.
- **Wilcoxon**: For non-parametric comparison.
- **Bayesian inference**: Using Bayes Factor (`BayesFactor` package).

Results are shown as group-wise statistics and a heatmap of differences.

## 🖥️ Installation

Ensure the following R packages are installed:

```r
install.packages(c(
  \"shiny\", \"igraph\", \"ggplot2\", \"reshape2\", \"dplyr\",
  \"visNetwork\", \"BayesFactor\"
))
