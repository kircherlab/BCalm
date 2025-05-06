# BCalm and analyze your MPRA data

BCalm is a package that provides a modification from [the mpralm package](https://github.com/hansenlab/mpra/tree/master), an R package that provides tools for differential analysis in MPRA studies.
BCalm allows users to use individual barcodes as model input.
See the [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-025-06065-9) for a detailed description, and the vignette for examples on how to run BCalm.

BCalm requires R >=3.5, <= 4.4.0 and can be installed using `devtools` or `remotes`.

### Installation guide:

#### Using conda
We suggest using conda as a package management tool. Its installation guide can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

`BCalm` is available from GitHub. Here we give installation instructions using either `devtools` or `remotes`.

**Remotes:**
```bash
conda create -n BCalm_env r-base=4.4.0 r-remotes
```

**Devtools:**
```bash
conda create -n BCalm_env r-base=4.4.0 r-devtools
```
After activating the environment (`conda activate BCalm_env`) you can start the R terminal and install `BCalm` using either `devtools` or `remotes`, by loading the chosen package and running:
```R
install_github("kircherlab/BCalm")
```

After installation, you can start using `BCalm` (after loading it with `library(BCalm)`).
For a more extensive user guide, please see the vignette (installation described below).


> **Note:**
> If you have problems installing BCalm dependencies (e.g. similar to [issue #10](https://github.com/kircherlab/BCalm/issues/10)) you can install them via conda. Due to one dependency (`bioconductor-genomeinfodbdata`) we have to use the gcc7 label for the bioconda channel. R base version might be different to `4.4.0` but BCALm should work on all R versions `bioconductor-mpra` is supported.
> ```bash
> conda install -c bioconda/label/gcc7 -c conda-forge bioconductor-mpra r-devtools r-tidyr r-ggplot2
> ```


### Vignette

In the following code snippets, we build the vignette and install all necessary packages for the vignette. 

**Remotes:**
```R
remotes::install_github("kircherlab/BCalm", build_vignette=TRUE, dependencies=TRUE)
```

**Devtools:**
```R
devtools::install_github("kircherlab/BCalm", build_vignette=TRUE, dependencies=TRUE)
```

After this you can open the built vignette by `vignette('BCalm')` (Alternatively, we provide a pre-built vignette in `vignettes/BCalm.html`)

You can either follow the prepared vignette directly in the editor of your choice (`vignettes/BCalm.Rmd`) or scroll through it by opening it in your browser (`vignettes/BCalm.html`).


### How to cite: 
Keukeleire, P., Rosen, J.D., Göbel-Knapp, A. et al. Using individual barcodes to increase quantification power of massively parallel reporter assays. BMC Bioinformatics 26, 52 (2025). https://doi.org/10.1186/s12859-025-06065-9

The paper for the `BCalm` package can be found [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-025-06065-9). When using `BCalm`, please also cite the mpra package (Myint et al. 2019) and the limma-voom framework (Law et al. 2014).

#### References
Myint, Leslie, Dimitrios G. Avramopoulos, Loyal A. Goff, and Kasper D. Hansen. *Linear models enable powerful differential activity analysis in massively parallel reporter assays.* BMC Genomics 2019, 209. doi: 10.1186/s12864-019-5556-x.
 
Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 2014, 15: R29. doi: 10.1186/gb-2014-15-2-r29
