# BCalm and analyze your MPRA data

BCalm is a package that provides a modification from [the mpralm package](https://github.com/hansenlab/mpra/tree/master), an R package that provides tools for differential analysis in MPRA studies.
BCalm allows users to use individual barcodes as model input.
See the vignette for examples on how to run BCalm.

BCalm requires R >=3.5, <= 4.4.0 and can be installed using `devtools` or `remotes`.

### Installation guide:

#### Using conda
The installation guide assumes conda as package management tool is installed. How to install can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

The provided package is available from Github and we will show how to install it with either `devtools` or `remotes`.

**Remotes:**
```bash
conda create -n BCalm_env r-base=4.4.0 r-remotes
```

**Devtools:**
```bash
conda create -n BCalm_env r-base=4.4.0 r-devtools
```
After activating the environment you can start the R terminal and install `BCalm` depending on which of the aforementioned packages you have installed. 

```bash
conda activate BCalm_env
```

If you just want to install `BCalm` without installing packages only necessary for the vignette skip the following two code blocks.  

In the following code snippets, we build the vignette and install all necessary packages for the vignette. 

**Remotes:**
```R
remotes::install_github("kircherlab/BCalm", build_vignette=TRUE, dependencies=TRUE)
```

**Devtools:**
```R
devtools::install_github("kircherlab/BCalm", build_vignette=TRUE, dependencies=TRUE)
```

After this you can open the build vignette by `vignette('BCalm')` (This is not possible if you are working on WSL but the vignette can be found at `vignettes/BCalm.html`)

**NOTE:** All required packages to run `BCalm` are installed if you use the following code as well. This command will not install additional packages only used by the vignette nor build the vignette.
**Remotes:**
```R
remotes::install_github("kircherlab/BCalm")
```

**Devtools:**
```R
devtools::install_github("kircherlab/BCalm")
```
### Vignette
You can either follow the prepared vignette directly in the editor of your choice (`vignettes/BCalm.Rmd`) or scroll through it by opening it in your browser (`vignettes/BCalm.html`).


### How to cite: 
The BCalm package is still unpublished, citing details will be provided later. When using BCalm, please cite the mpra package (Myint et al. 2019) and the limma-voom framework (Law et al. 2014).

#### References
Myint, Leslie, Dimitrios G. Avramopoulos, Loyal A. Goff, and Kasper D. Hansen. *Linear models enable powerful differential activity analysis in massively parallel reporter assays.* BMC Genomics 2019, 209. doi: 10.1186/s12864-019-5556-x.
 
Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 2014, 15: R29. doi: 10.1186/gb-2014-15-2-r29
