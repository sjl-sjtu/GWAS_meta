# GWAS_meta: A rapid Bayesian methods for subset selection and GWAS meta-analysis.
This package provids three methods to get optimal ABF for the variant-trait association in Bayesian GWAS meta-analysis across numerous studies. The three approaches include Subset-Exhaustive, MCMC and Shotgun Stochastic Search. 

Install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/GWAS_meta")
```
Usage and examples can be found at https://github.com/sjl-sjtu/GWAS_meta/blob/main/GWASmeta_0.1.0.pdf. You can cite the work and find more methodological details about the package at [[1]](#rf1). This package is modified based on the Bayesian meta-analysis package "metabf" which is available at https://github.com/trochet/metabf. We appreciate the authors' work and details about the original edition can be found at [[2]](#rf2).

## Software
An web tool to conduct GWAS meta-analysis using SMetABF is developed based on R Shiny. You can use the tool at https://sunjianle-sjtu.shinyapps.io/analycode/. Please uplaod your file according to the example format. The example data can be found at the folder https://github.com/sjl-sjtu/GWAS_meta/tree/main/data.

## Reference
<div id="rf1"></div>
[1] Sun, J., Lyu, R., Deng, L., Li, Q., Zhao, Y., & Zhang, Y. (2022). SMetABF: A rapid algorithm for Bayesian GWAS meta-analysis with a large number of studies included. <i>PLoS computational biology</i>, 18(3), e1009948. https://doi.org/10.1371/journal.pcbi.1009948
<div id="rf2"></div>
[2] Trochet, H., Pirinen, M., Band, G., Jostins, L., McVean, G., & Spencer, C. (2019). Bayesian meta-analysis across genome-wide association studies of diverse phenotypes. <i>Genetic epidemiology</i>, 43(5), 532–547. https://doi.org/10.1002/gepi.22202

## Lisence
This package is available under the GNU General Public License version 3.

The latest updation for the package is at Sep. 16th, 2023.

If there is any problem, please contact me: Jianle Sun (sjl-2017@sjtu.edu.cn)
