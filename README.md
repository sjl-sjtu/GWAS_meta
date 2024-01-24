# GWAS_meta: A rapid Bayesian methods for subset selection and GWAS meta-analysis.
This package provids three methods to get optimal ABF for the variant-trait association in Bayesian GWAS meta-analysis across numerous studies. The three approaches include Subset-Exhaustive, MCMC and Shotgun Stochastic Search. 

Latest revision: Jan. 24th, 2024.

## Installation
Install the package as the following:
```R
library(devtools)
devtools::install_github("sjl-sjtu/GWAS_meta")
```
Then use the package as
```R
library(GWASmeta)
```
Usage and examples can be found at https://github.com/sjl-sjtu/GWAS_meta/blob/main/GWASmeta_0.1.0.pdf. You can cite the work and find more methodological details about the package at [[1]](#rf1). This package is modified based on the Bayesian meta-analysis package "metabf" which is available at https://github.com/trochet/metabf. We appreciate the authors' work and details about the original edition can be found at [[2]](#rf2).

## Tutorial
### 1. GWAS meta-analysis for single locus
First you need to organize the data set. As shown in the figure below,  with each row corresponding to a study, the organized dataset should contain two columns. The first column is the effect size ($\beta$) between the SNP and the phenotype obtained in each study, which is usually the coefficient of linear regression, or Log OR value in logistic regression. The second column is the standard error of $\beta$ in the corresponding study in the first column.
![image](https://github.com/sjl-sjtu/GWAS_meta/blob/main/software/www/single.png)

Here is an example using the built-in sample dataset
```R
data(single)
betas <- single$betas
ses <- single$ses

re <- exh_abf(betas,ses) # exhausive
re <- mcmc_abf(betas,ses) # MCMC
re <- shotgun_abf(betas,ses) # SSS
```

If you want to obtain the ABF, and at the same time obtain the studies included in the final effect subset, you can use the following code
```R
re <- shotgun_abf_model(betas,ses)
```
The result will return a list, the first element is the calculated ABF, and the second element is a vector composed of 0 and 1, representing whether each study is included in the effect subset.

There are some important parameters that can be specified, including
* `prior.sigma`: the prior variance of the impact of each study
* `prior.cor`, the specification of the relationship between studies in the subset, which can be specified as `indep` (independent effect), `fixed` (fixed effect) or `correlated` (correlated effect). When specifying the model with a correlation effect, you can specify the prior of the correlation coefficient by specifying `prior.rho`. In addition,  ​​
* `log` and `log10`: two Boolean values used to specify whether the calculated ABF should be taken as the natural logarithm or the common logarithm before output.

### 2. GWAS meta-analysis for multiple loci
First, you need to organize the data, as shown in the figure below. Each row represents a SNP (locus). The first column is the name of the SNP, and every two columns after that represent the effect size ($\beta$) of each SNP on the phenotype in a study and the corresponding estimation standard (se), i.e., the second column is the estimated effect size of each SNP on the trait in the first study, the third column is the estimated standard error of the corresponding SNP's effect on the trait in the first study, and the fourth column is the estimated effect value of each SNP on the trait in the second study, and the 5th column is the estimated standard error of the corresponding SNP's effect on the trait in the second study, and so on. Missing values ​​are represented by NA.
![image](https://github.com/sjl-sjtu/GWAS_meta/blob/main/software/www/multiple.png)

Here is an example using the built-in sample dataset
```R
data(multi)
re <- multi_shotgun_abf(multi)
```
The result is a data frame, including the name of each SNP, the calculated ABF, a 01 vector representing whether each study included in the analysis belongs to the effective subset of the SNP (The length of the vector is equal to the number of non-NA studies), the number of studies included in the analysis of the SNP (i.e., all The number of studies excluding studies with NA estimates for the SNP, and a 01 vector indicates whether each study was included in the analysis of the SNP (The length of the vector is equal to the number of all studies in the input dataset).

## Online Software
An web tool to conduct GWAS meta-analysis using SMetABF is developed based on R Shiny. You can use the tool at https://sunjianle-sjtu.shinyapps.io/analycode/. Please uplaod your file according to the example format. The example data can be found at the folder https://github.com/sjl-sjtu/GWAS_meta/tree/main/data. The output data frame is consistent with the output of the above code.

**Note**: Please do not run large datasets on online tools!

## Reference
<div id="rf1"></div>
[1] Sun, J., Lyu, R., Deng, L., Li, Q., Zhao, Y., & Zhang, Y. (2022). SMetABF: A rapid algorithm for Bayesian GWAS meta-analysis with a large number of studies included. <i>PLoS computational biology</i>, 18(3), e1009948. https://doi.org/10.1371/journal.pcbi.1009948
<div id="rf2"></div>
[2] Trochet, H., Pirinen, M., Band, G., Jostins, L., McVean, G., & Spencer, C. (2019). Bayesian meta-analysis across genome-wide association studies of diverse phenotypes. <i>Genetic epidemiology</i>, 43(5), 532–547. https://doi.org/10.1002/gepi.22202

## Lisence
This package is available under the GNU General Public License version 3.

The latest updation for the package is at Sep. 16th, 2023.

If there is any problem, please contact me: Jianle Sun (sjl-2017@sjtu.edu.cn)
