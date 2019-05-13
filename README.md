# Bayesian Decision-Making Forests (BaD-MF)

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Issues](https://github.com/ebridge2/badmf/issues)
- [Citation](#citation)

![A BaD-MF](https://media0.giphy.com/media/l2YWxte7sJB2XuE8M/giphy.gif)


## Overview

A package for fitting Bayesian Decision-Making Forests (BaD-MFs) using the `R` programming language. At each node of a decision tree, a subset of features is chosen from the target feature space for determination of the optimal split criterion. In a simple classification problem, for instance, a single feature is selected from `d < p` possible choices. These features are selected uniformly and at random, and no data outside of the data seen by the tree may be used for assessing the split criterion or determining the features to attempt at a given split.

In a random forest, numerous decision trees are aggregated to construct an ensemble learner. Trees are build using 
Bootstrap AGGregation (bagging), through which each tree sees a nonparametric sample from the target dataset upon which 
a decision tree is trained. Bagging has been shown by numerous investigators to reduce the variance and help avoid overfitting of individual trees.

As each tree only obtains a subset of the dataset, we risk choosing features that may, or may not, accurately represent the
target dataset. While growing the forest may seem an easy solution, larger forests are computationally expensive, and we risk
not seeing the features that really matter for our classification problem at a high frequency representative the quantity of signal they contain.

Through `badmf`, we train an initial random forest using a slight augmentation to the traditional approach, where we sample the feature probability vector using a user-provided Dirichlet prior, defaulting to a non-informative prior (`Dir(1,1,...)`). Assuming the sampling distribution for the features is Multinomial, we obtain a Dirichlet Posterior by combining information from our prior with the relative counts of each feature in the trained forest. We finally train a new classifier, where the feature probability vectors are instead sampled using the resulting posterior at each split node. In this fashion, each constructed tree is more likely to attempt splits along the features containing a high portion of signal. Optimally, the resulting forest is small (running less risk of being highly variant and overfit) and efficient (with the optimal features selected more frequently, we obtain similar performance to a much larger and more complex forest using far fewer nodes and greater adherence of feature importance in the resulting forest to the truth).

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): package documentation, and usage of the `badmf` package on simulated examples.
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

Users should install the following packages prior to installing `badmf`, from an `R` terminal:

```
install.packages(c('parallel', 'MCMCpack', 'forcats', 'MASS', 'ggplot2'))
```

### Package Installation

From an `R` session, type:

```
require(devtools)
install_github('ebridge2/badmf', build_vignettes=TRUE, force=TRUE)  # install badmf with the vignettes
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
