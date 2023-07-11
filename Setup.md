# Setup
All of the calculations and analyses are done in the iPython notebook. The code is written in Python and partially in R. The R is used to get data from GEO in a CEL file and process it to RNA-seq type data.

We recommend installing our Python virtual environment to perform the analysis. For detailed instructions, please refer to the “Preparation of the environment” section


## Environment Requirements 
* Python 3.10
  * The packages are in requirements.txt
* R 4.0.0 or higher
  * The packages are in install_R_packages.R
* Jupyter notebook with its R and Python kernels
* WSL (for Windows users)

## WSL installation
If you are a Windows user, install WSL on your computer. If your operating system is Unix-based, skip this step.

To install WSL, follow the instructions provided on the [WSL installation webpage](https://learn.microsoft.com/en-us/windows/wsl/install).

## Preparation of the environment
Please follow the instruction steps in the given order:

***Installation of python using apt source***


    sudo apt-get update
    sudo apt-get install python3.10-venv python3.10-dev python3-pip
    
***Installation of python environment via pip***


    clone https://github.com/BostonGene/MFP
    cd MFP
    bash make_tme_environment.sh

make_tme_environment.sh creates python3.10 environment with all necessary packages and creates ipykernel core for the environment with name tme_env


***Installation of python using conda***


If you want to create a python environment via conda please follow this [link](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/installing-with-conda.html) 


***Install the jupyter kernel for your environment (for conda)***


    python -m ipykernel install --user --name=MFP_env


**If your data is in CEL format you have to download R and all of the required packages.**


***Installation and preparation of R***


    sudo apt install r-base-core 
    Rscript -e "install.packages('IRkernel')"
    Rscript -e "IRkernel::installspec(user = FALSE)"
    Rscript -e 'install.packages("httr", repos="http://cran.rstudio.com/")' 
    Rscript -e 'install.packages("RJSONIO", repos="http://cran.rstudio.com/")' 
    sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev
    Rscript install_R_packages.R


