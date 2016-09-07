## Sequence of commands to check and build the package 
## and generate html page documents

#install.packages(devtools)
#install_github('hadley/staticdocs')

require(devtools)
require(staticdocs)

#package.dir <- normalizePath("../../../Gistic2GE")
package.dir <-"/Users/amyli/Desktop/git_projects/iEDGE"




cat("documenting package...\n")
document(package.dir) #creates help pages

cat("building package..\n")
build(package.dir)

cat("checking package...\n")
check(package.dir) #package checking
#check(package.dir, args = "--no-examples")
cat("loading package...\n")
load_all(package.dir) #building package
cat("library package...\n")
library(iEDGE)

##install package
#cat("installing package...\n")
#install.packages(package.dir, repos = NULL, type = "source")


#PAT <- "04fe676593e46b6bda5a5d09431156e8a500349a"
#install_github("montilab/iEDGE",auth_token = PAT)


#generate html pages
setwd(package.dir)
cat("creating html help pages for package...\n")
build_site(pkg = package.dir, examples = TRUE, launch = TRUE)


