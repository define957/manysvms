# **manysvms** <img src="man/figures/Logo.png" align="right" width="150" />
             
![GitHub RCMDCHECK](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)
![GitHub license](https://img.shields.io/github/license/define957/manysvms)
![GitHub ISSUES](https://img.shields.io/github/issues/define957/manysvms)
![GitHub STARS](https://img.shields.io/github/stars/define957/manysvms)
![GitHub FORKS](https://img.shields.io/github/forks/define957/manysvms)

***


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
After you installed `devtools`, remember to instal `RTools`. You can find suitable `RTools` on the website https://cran.r-project.org/bin/windows/Rtools/. （on Mac OSX you should install suitable `clang` and `gfortran`, see more: https://cran.r-project.org/bin/macosx/tools/ and https://mac.r-project.org/tools/）
After above steps, please run the following command :
```{r}
devtools::install_github("define957/manysvms")
```
Then you can have `manysvms` package on your PC。
## Important notice
This project is currently being redesigned !!! 🔧
The old version of the function will be deprecated in the new version !!!

## SVMs for classification

+ Support Vector Machine for Multi-classification (One versus Rest strategy) 
+ + [Rcpp acceleration is available]
+ Least Squares Support Vector Machine for Multi-classification (One versus Rest strategy) 
+ + [Rcpp acceleration is available]
+ Twin-SVM for binary classification 
+ + [Rcpp acceleration is available]
+ Twin-SVM for Multi-classification (One versus Rest strategy) 
+ + [Rcpp acceleration is available]
+ + [Parallel execution of grid searches is available]
+ Multiple Birth SVM for Multi-classification 
+ + [Rcpp acceleration is available] 
+ + [Parallel execution of grid searches is available]
+ Twin K SVM for Multi-classification (One versus One versus Rest strategy) 
+ + [Rcpp acceleration is available]
+ Ramp Twin K SVM for Multi-Classification (One versus One versus Rest strategy)
+ + [Rcpp acceleration is available]

## SVMs for regression

+ Epsilon Support Vector Regression [Rcpp acceleration is available]
+ Twin Support Vector Regression [Rcpp acceleration is available]
+ Least Squares Support Vector Regression

## Kernel options

+ Linear kernel
+ RBF kernel
+ Polynomial kernel

## Development environment and dependency

My enviroment: R 4.2.1, windows 10 x64 &#x2705;

Other test environment detail: 
+ Windows 10/11 x64 &#x2705;
+ Mac osx (ARM platform) &#x2705; 
+ Linux : We haven't tested it yet &#x2753;

Dependency: 

+ ggplot2 
+ Rcpp
+ RcppArmadillo
+ foreach
+ doParallel
+ doSNOW

## Bug report

If you find bug in this package, please post an issue on the [issue](https://github.com/define957/manysvms/issues) website.

## Contact us

&#x2709; Email : zhangjiaqi957957@outlook.com

## Licenses

GNU GENERAL PUBLIC LICENSE Version 3 (GPL-3.0)

## Examples

MBSVM :
<div align=center><img src = man\figures\MBSVM.png width="60%"></div>
SVR :
<div align=center><img src = man\figures\SVR.png width="60%"></div>
