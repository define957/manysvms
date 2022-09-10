# many SVMs
![GitHub RCMDCHECK](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)
![GitHub license](https://img.shields.io/github/license/define957/manysvms)
![GitHub ISSUES](https://img.shields.io/github/issues/define957/manysvms)
![GitHub STARS](https://img.shields.io/github/stars/define957/manysvms)
![GitHub STARS](https://img.shields.io/github/forks/define957/manysvms)
<div align=center><img src = man\figures\MBSVM.png width="60%"></div>


## Introduction

There are many different types of SVMs in this repository. 

## Why I created this project ?**

In order to learn `SVMs` better, I built this repository to implement support vector machines. The package is under active development.

## How to install manysvms

Make sure you have installed devtools and run the following command :
```{r}
devtools::install_github('define957/manysvms')
```

## SVMs

+ Twin-SVM for binary classification
+ Twin-SVM for Multi-classification (Ones versus Rest strategy) 
+ Twin-SVM for Multi-classification (Ones versus Rest strategy and K-fold cross validation) 
+ Multiple Birth SVM for Multi-classification
+ Epsilon Support Vector Regression [Rcpp acceleration is available]

## Kernel options

+ linear kernel
+ rbf kernel

## Development environment and dependency

R Version : 4.2.1

## Contact us

Email : zhangjiaqi957957@outlook.com

## Licenses

GNU GENERAL PUBLIC LICENSE Version 3 (GPL-3.0)
