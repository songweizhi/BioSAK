# ICoVeR

source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("opencpu")
install.packages(devtools)


library(devtools)
library(opencpu)
install_local(file.path(getwd(), "R.ICoVeR"))

opencpu$stop() # It starts at a random port, which is annoying.
opencpu$start(8000)
browseURL("http://localhost:8000/ocpu/library/ICoVeR/www/", browser = getOption("browser"), encodeIfNeeded = FALSE)
