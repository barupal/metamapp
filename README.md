# metamapp
R codes for creating metamapp graphs and files

MetaMapp installation method.

In R, run following code.
```
if (!require("devtools"))
install.packages('devtools', repos="http://cran.rstudio.com/")
library(devtools)
source('https://bioconductor.org/biocLite.R')
install_github('barupal/metamapp')
library(MetaMapp2016)
library(opencpu)
opencpu$browse('/library/MetaMapp2016/www')
```

