# Bayesian Decision-Making Forests (BaD-MF)

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
- [License](./LICENSE)
- [Issues](https://github.com/ebridge2/lol/issues)
- [Citation](#citation)

## Overview

A package for fitting Bayesian Decision-Making Forests (BaD-MFs). At each node of a decision tree, a subset of features is chosen from the target feature space for determination of the optimal split criterion. In a simple classification problem, for instance, a single feature is selected from `d < p` possible choices. These features are selected uniformly and at random,and no data outside of the data seen by the tree may be used for assessing the split criterion.

In a random forest, numerous decision trees are aggregated to construct an ensemble learner. Trees are constructed using 
Bootstrap AGGregation (bagging), through which each tree sees a nonparametric sample from the target dataset upon which 
a decision tree is constructed. Bagging has been shown by numerous investigators to reduce the variance and help avoid overfitting of individual trees.

As each tree only obtains a subset of the dataset, we risk choosing features that may, or may not, accurately represent the
target dataset. While growing the forest may seem an easy solution, larger forests are computationally expensive, and we risk
not seeing the features that really matter for our classification problem at a high frequency representative the quantity of signal they contain.

Through `badmf`, we train an initial random forest using the traditional approach, where features are selected for assessment 
uniformly and at random using a non-informative Dirichlet Prior (each feature sampled with equal probability). Assuming the sampling distribution for the features is Multinomial, we obtain a Dirichlet Posterior by combining information from our non-informative prior with the relative counts of each feature in the trained forest. We finally train a new classifier, where features are instead sampled using the resulting posterior at each split node. Optimally, the resulting forest is small (running less risk of being highly variant and overfit) and efficient (with the optimal features selected more frequently, we obtain similar performance to a much larger and more complex forest using far fewer nodes).

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): package documentation, and usage of the `lolR` package on many real and simulated data examples.
- [man](./man): package manual for help in R session.
- [tests](./tests): `R` unit tests written using the `testthat` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.

# System Requirements

## Hardware Requirements

The `badmf` package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3.3 GHz) and internet of speed 25 Mbps. Note that multiple cores greatly improves execution, particularly with wider forests, as the package is parallelized over trees.

## Software Requirements

### OS Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Ubuntu 16.04  
Mac OSX:  
Windows:  

The CRAN package should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the `badmf` package, users should have `R` version 3.5.3 or higher, and several packages set up from CRAN.

# Installation Guide

## Development Version

### Package dependencies

Users should install the following packages prior to installing `lolR`, from an `R` terminal:

```
install.packages(c('parallel', 'MCMCpack', 'forcats'))
```

### Package Installation

From an `R` session, type:

```
require(devtools)
install_github('ebridge2/badmf', build_vignettes=TRUE, force=TRUE)  # install lol with the vignettes
require(badmf)
vignette("badmf", package="badmf")  # view one of the basic vignettes
```

The package should take approximately 40 seconds to install with vignettes on a recommended computer. 

# Demo

## Functions

For interactive demos of the functions, please check out the vignettes built into the package. They can be accessed as follows:

```
require(badmf)
vignette('badmf')
vignette('rf')
vignette('decision_tree')
```
