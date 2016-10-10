# MetaMapp
R codes for creating metamapp graphs and files

MetaMapp installation method.

In R, run following code.
```
if (!require("devtools"))
install.packages('devtools', repos="http://cran.rstudio.com/")
library(devtools)
library(RCurl)
source('https://bioconductor.org/biocLite.R')
install_github('barupal/metamapp')
library(MetaMapp2016)
library(opencpu)
opencpu$browse('/library/MetaMapp2016/www')
```

## Citation

[Barupal, D.K., Haldiya, P.K., Wohlgemuth, G., Kind, T., Kothari, S.L., Pinkerton, K.E. and Fiehn, O., 2012. MetaMapp: mapping and visualizing metabolomic data by integrating information from biochemical pathways and chemical and mass spectral similarity. BMC bioinformatics, 13-1, p.1.] (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-99) 

