# Distorted M-Quantiles in R

## Extended Abstract
**Distorted M-Quantiles: Depth, Central Regions, and Multiple Output Regression in R**  
This repository contains a collection of R functions designed to compute distorted $M$-quantiles, central regions, data depth, and regression models for multivariate datasets, as introduced in the article "Data Depth and Multiple Output Regression, the Distorted $M$-Quantiles Approach" by Maicol Ochoa and Ignacio Cascos (Mathematics, 2022). The $M$-quantiles generalize classical quantiles and expectiles by solving asymmetric minimization problems with power loss functions, where the parameter $r$ determines the loss (e.g., $r=1$ for quantiles, $r=2$ for expectiles). The distortion functions (e.g., trim, sigmoid) adjust the weighting of distribution tails, enhancing robustness against outliers. These concepts are leveraged to define multivariate central regions via intersections of halfspaces and associated depth functions, which measure the centrality of points in $\mathbb{R}^d$. The code extends these ideas to multiple output regression, computing conditional and regression regions that capture trends, variability, and dependencies in the response variables given predictors.

The implemented functions include:
- `dist_mquantile`: Computes distorted $M$-quantiles for univariate data with a specified power $r$ and distortion function $g$.
- `rexpectile` and `rquant`: Robust expectile and quantile calculations with trimming options.
- `lsreg`: Fits $M$-quantile regression models, supporting robust multivariate regression.
- `rextreme.points`, `qextreme.points`, `mqextreme.points`, etc.: Generate extreme points for central regions using halfspace intersections.
- `conditional.regions`: Constructs conditional regression regions for multivariate responses, visualizing relationships between predictors and responses.

The theoretical foundation stems from statistical depth literature (e.g., Tukey’s halfspace depth) and $M$-quantile advancements (e.g., Breckling and Chambers, 1988), with distortions inspired by non-additive probabilities and Choquet expectations. This implementation provides practical tools for robust multivariate analysis, as demonstrated with simulated and real datasets in the original article.

## Installation
To use this code, ensure you have R installed along with the following packages:
```R
install.packages(c("MASS", "geometry", "depth"))
```

## Usage
See the example visualizations in the script for generating mregions.pdf (central regions) and conditional_regions.pdf (conditional regression regions).

# Acknowledgments
The code was created by I. Cascos (https://github.com/icascos) and M. Ochoa. The theoretical description of the procedures is available at https://doi.org/10.3390/math10183272.
