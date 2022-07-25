
# check.packages function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = 1)}

# install packages if not installed
packages<-c("optparse", "ape", "vegan")
invisible(suppressMessages(check.packages(packages)))

option_list = list(
  make_option(c("-a", "--treeo"),  type="character", default=NULL, help="the first tree"),
  make_option(c("-b", "--treet"),  type="character", default=NULL, help="the second tree"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

TREE1 = read.tree(opt$treeo)
TREE2 = read.tree(opt$treet)

D1 = cophenetic(TREE1)
D1 = D1[order(row.names(D1)),order(row.names(D1))]
D2 = cophenetic(TREE2)
D2 = D2[order(row.names(D2)),order(row.names(D2))]

mantel(xdis = D1, ydis = D2, permutations = 999)

