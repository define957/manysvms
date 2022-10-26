# many SVMs
![GitHub RCMDCHECK](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)
![GitHub license](https://img.shields.io/github/license/define957/manysvms)
![GitHub ISSUES](https://img.shields.io/github/issues/define957/manysvms)
![GitHub STARS](https://img.shields.io/github/stars/define957/manysvms)
![GitHub STARS](https://img.shields.io/github/forks/define957/manysvms)

MBSVM :
<div align=center><img src = man\figures\MBSVM.png width="60%"></div>
SVR :
<div align=center><img src = man\figures\SVR.png width="60%"></div>


## Introduction

There are many different types of SVMs in this repository. 

## Organization structure

<div align=center><img src = man\figures\manysvms.svg width="100%"></div>

## Why I created this project ?

In order to learn `SVMs` better, I built this repository to implement support vector machines. The package is under active development.

## How to install manysvms

Make sure you have installed `devtools`, if you don't have `devtools` installed, please run the following command first 
```{r}
install.packages("devtools")
```

and run the following command :
```{r}
devtools::install_github("define957/manysvms")
```

## SVMs for classification

+ Twin-SVM for binary classification [Rcpp acceleration is available]
+ Twin-SVM for Multi-classification (One versus Rest strategy) [Rcpp acceleration is available]
+ Multiple Birth SVM for Multi-classification [Rcpp acceleration is available]
+ Twin K SVM for Multi-classification (One versus One versus Rest strategy) [Rcpp acceleration is available]
+ Ramp Twin K SVM for Multi-Classification (One versus One versus Rest strategy) [Rcpp acceleration is available]

## SVMs for regression

+ Epsilon Support Vector Regression [Rcpp acceleration is available]
+ Twin Support Vector Regression [Rcpp acceleration is available]
+ Least Squares Support Vector Regression

## Kernel options

+ Linear kernel
+ RBF kernel
+ Polynomial kernel

## Development environment and dependency

enviroment: R 4.2.1, windows 10 x64 &#x2705;

dependency: 

+ ggplot2 
+ Rcpp
+ RcppArmadillo

## Bug report

If you find bug in this package, please post an issue on the [issue](https://github.com/define957/manysvms/issues) website.

## Contact us

&#x2709; Email : zhangjiaqi957957@outlook.com

## Licenses

GNU GENERAL PUBLIC LICENSE Version 3 (GPL-3.0)
