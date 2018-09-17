# MetaMapp
R codes for creating metamapp graphs and files

MetaMapp installation method.

In R, run following code.
```
if (!require("devtools"))
install.packages('devtools', repos="http://cran.rstudio.com/")
if (!require("opencpu"))
install.packages('opencpu', repos="http://cran.rstudio.com/")
if (!require("RCurl"))
install.packages('RCurl', repos="http://cran.rstudio.com/")
library(devtools)
library(RCurl)
source('https://bioconductor.org/biocLite.R')
install.packages("https://github.com/barupal/metamapp/blob/master/MetaMapp2017_2.0.1.tar.gz?raw=true", repos = NULL, type = "source")
library(MetaMapp2017)
library(opencpu)
opencpu$browse('/library/MetaMapp2017/www')
```
## Online version 

 [MetaMapp Online Version](http://metamapp.fiehnlab.ucdavis.edu)

## Citation

[Barupal, D.K., Haldiya, P.K., Wohlgemuth, G., Kind, T., Kothari, S.L., Pinkerton, K.E. and Fiehn, O., 2012. MetaMapp: mapping and visualizing metabolomic data by integrating information from biochemical pathways and chemical and mass spectral similarity. BMC bioinformatics, 13-1, p.1.] (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-99) 

## Docker image 

[MetaMapp Docker image] (https://hub.docker.com/r/barupal/metamapp-docker/)
