# nbregression command-line tool

This package contains a commmand-line tool that implements a Bayesian negative binomial generalized linear regression model for single-cell RNA-seq. The algorithm was described in Zeisel et al. *Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq* **Science** 2015 (PMID: [25700174](http://www.ncbi.nlm.nih.gov/pubmed/25700174), doi: [10.1126/science.aaa1934](http://dx.doi.org/10.1126/science.aaa1934)). Please cite this paper if you use the nbregression algorithm in your work.

The regression model was implemented in [Stan](http://mc-stan.org/) by Sten Linnarsson, with a Python wrapper and command-line tool by Peter Lönnerberg. nbregression takes input in CEF format and produces an annotated CEF file as output. Use [ceftools](https://github.com/linnarsson-group/ceftools) to create and manipulate CEF files.
 
...work in progress...
