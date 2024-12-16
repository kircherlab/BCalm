# BCalm and analyze your MPRA data

BCalm is a package that provides a modification from [the mpralm package](https://github.com/hansenlab/mpra/tree/master), an R package that provides tools for differential analysis in MPRA studies.
BCalm allows users to use individual barcodes as model input.
See the vignette for examples on how to run BCalm.

BCalm requires R >=3.5, <= 4.4.0 and can be installed using `devtools` or `remotes`.

### Installation guide:

#### Using conda
The installation guide assumes conda as package management tool is installed

```bash
conda create -n BCalm_env r-base=4.4.0
```
After activating the environment you can start the R terminal you can use either devtools or remotes

```R
install.packages("remotes", repos='http://cran.us.r-project.org', dependencies=TRUE)
remotes::install_github("kircherlab/bcalm")
```

```R
install.packages("devtools", repos='http://cran.us.r-project.org', dependencies=TRUE)
devtools::install_github("kircherlab/bcalm")
```
### Vignette
You can either follow the prepared vignette directly in the editor of your choice (`vignettes/BCalm.Rmd`) or scroll through it by opening it in your browser (`vignettes/BCalm.html`).


mpralm publication:
Myint, Leslie, Dimitrios G. Avramopoulos, Loyal A. Goff, and Kasper D. Hansen. *Linear models enable powerful differential activity analysis in massively parallel reporter assays.* BMC Genomics 2019, 209. doi: 10.1186/s12864-019-5556-x.
