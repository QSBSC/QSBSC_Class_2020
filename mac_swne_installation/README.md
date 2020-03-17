# Installation of swne R package on Macs
If you tried to install the *LLNM* or *swne* R packages on a Mac and got an error that looked like:  

`clang: error: unsupported option '-fopenmp'`  

`make: *** [RcppExports.o] Error 1`  

it means that you don't have a compiler on your computer that supports OpenMP (the default Apple-provided compiler does not.)

To overcome this problem you should install a compiler that supports OpenMP, such as `gcc` or `llvm`. The instructions below are how to install `llvm` using [Homebrew](https://brew.sh).

### Brew Installation
These steps will all be done in Terminal.
* Make sure R and RStudio are closed.
* Install Homebrew by typing:
`/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`
* Use Homebrew to install the llvm compiler: `brew install llvm`

### Let R know the location of the new compilers
* First make sure you are in your home directory by typing: `cd ~`  
* Check for presence of hidden `.R` directory by typing: `ls -la`  
* If not, make a new (hidden) directory by typing: `mkdir .R`
* Download the file [`Makevar.txt`](https://raw.githubusercontent.com/QSBSC/QSBSC_Class_2020/master/mac_swne_installation/Makevars.txt) (right click and select `Download Linked File As...`) and note the location (e.g. `~/Desktop/Makevars.txt`)
* Move/rename the file into the new hidden folder you made by typing in *Terminal*: `mv <path_to_Makevars.txt> ~/.R/Makevars`

After that you should be able to install the libraries *in R* using:  

`if(!require(devtools)){ install.packages("devtools")}`  
`devtools::install_url("https://cran.r-project.org/src/contrib/Archive/NNLM/NNLM_0.4.3.tar.gz", type="source")`  
`devtools::install_github("yanwu2014/swne")`  

#### Thanks to stackoverflow!
The information on how to do this was obtained from the response to a [question on stackoverflow]  (https://stackoverflow.com/questions/43595457/alternate-compiler-for-installing-r-packages-clang-error-unsupported-option)
