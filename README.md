# PhosR

<!-- badges: start -->
[![R build status](https://github.com/PYangLab/PhosR/workflows/R-CMD-check/badge.svg)](https://github.com/PYangLab/PhosR/actions)
  <!-- badges: end -->

`PhosR` is a package for the comprenhensive analysis of phosphoproteomic data. There are two major components to PhosR: processing and downstream analysis. PhosR consists of various processing tools for phosphoproteomic data including filtering, imputation, normalisaton and batch correction, which enables integration of multiple phosphoproteomic datasets. Downstream analytical tools consists of site- and protein-centric pathway analysis to evaluate activities of kinases and signalling pathways, large-scale kinase-substrate annotation from dynamic phosphoproteomic profiling, and visualisation and construction of signalomes present in the phosphoproteomic data of interest.

## PhosR overview

<img src="https://raw.githubusercontent.com/PYangLab/PhosR/master/inst/graphical_abstract.png" align="center"/>


## Installation

```r
library(devtools)
devtools::install_github("PYangLab/PhosR")
```

## Vignette 

<!--You can find the vignette at our website: https://PYangLab.github.io/PhosR/articles/webOnly/index.html-->
Please note that the vignette generated in the link below is using the full dataset. 
Our vignette prepared for BioC (`PhosR.Rmd`) utilisesthe subset of these datasets and the results may not be identical.


* Processing of phosphoproteomic data 
     * Filtering and imputation (https://PYangLab.github.io/PhosR/articles/webOnly/imputation.html)
     * Batch correction (https://PYangLab.github.io/PhosR/articles/webOnly/batch_correction.html)
* Downstream analysis of phosphoproteomic data
     * Pathway analysis (https://PYangLab.github.io/PhosR/articles/webOnly/pathway_analysis.html)
     * Site- and gene-centric analysis (https://PYangLab.github.io/PhosR/articles/webOnly/site_gene_analysis.html)
     * Kinase-substrate relationship scoring and signalome construction (https://PYangLab.github.io/PhosR/articles/webOnly/signalomes.html)

## Code relating to our preprint

* Data imputation: [https://pyanglab.github.io/PhosR/Analysis-of-FL83B-and-Hepa-1-6.html]
* Identification of stably phosphorylated sites: [https://pyanglab.github.io/PhosR/IdentifySPS_dataIntegration.html]
* 1-3 dimensional enrichmet analysis: [https://pyanglab.github.io/PhosR/1,-2,-3D-Enrichment-Analysis.html]

## Contact us

If you have any enquiries, especially about performing PhosR to analyse your phosphoproteomic data, please contact taiyun.kim@sydney.edu.au or jieun.kim@sydney.edu.au. We are also happy to receive any suggestions and comments.


