# Bayesian negative binomial generalized linear regression for single-cell RNA-seq

Single-cell RNA-seq can generate highly accurate and sensitive data on the abundance of mRNA in single cells. Using unique molecular identifiers (UMIs), it is now possible to actually count the number of mRNA molecules in single cells, modulo the losses incurred during reverse transcription. 

This package contains the [Stan](http://mc-stan.org/) code for a generalized linear regression model, which can be used to estimate the true underlying expression level of a gene in a population of single cells. Such estimates can then be used to ask questions about differences between cell types, over time or between experiments. 
