# GWAS_meta
This package provids three methods to get optimal ABF of single SNP and a phenotype in Bayesian GWAS meta-analysis including Subset-Exhaustive, MCMC and shotgun stochastic search. 

Install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/GWAS_meta")
```
This package is modified based on the Bayesian meta-analysis package "metabf" which is available at https://github.com/trochet/metabf. We appreciate the authors' work and details about the original edition can be found at [[1]](#rf1).

## Software
An web tool to conduct GWAS meta-analysis using SMetABF is developed based on R Shiny. You can use the tool at https://sunjianle-sjtu.shinyapps.io/analycode/. Please uplaod your file according to the example format. The example data can be found at the folder .

## Reference
<div id="rf1"></div>
[1] Trochet, H., Pirinen, M., Band, G., Jostins, L., McVean, G., & Spencer, C. (2019). Bayesian meta-analysis across genome-wide association studies of diverse phenotypes. Genetic epidemiology, 43(5), 532–547. https://doi.org/10.1002/gepi.22202

## Lisence
This package is available under the GNU General Public License version 3.
