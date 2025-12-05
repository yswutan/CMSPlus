# CMSPlus
An R package for molecular subtyping by integrating the inter-tumor and intra-tumor heterogeneity of colorectal cancer

```
devtools::install_github("yswutan/CMSPlus")
library(CMSPlus)
library(GSVA)

## exp2symbol: a dataframe with Gene Expression Profiles data values,
##       samples in columns, genes in rows, rownames corresponding to gene symbols
## parallel.sz: If parallel is loaded and this argument is left with its default value (parallel.sz=0) then it will use all available core processors unless we set this argument with a smaller number.

CMSPlusLabels <- CMSPlus(exp2symbol, plot=TRUE, parallel.sz=0)

```

