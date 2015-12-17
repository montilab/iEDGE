## Sequence of commands to check and build the package 
## and generate html page documents

#install.packages(devtools)
#install_github('hadley/staticdocs')

require(devtools)
require(staticdocs)

#package.dir <- normalizePath("../../../Gistic2GE")
package.dir <-"/Users/amyli/Desktop/git_projects/Gistic2GE"

cat("documenting package...\n")
document(package.dir) #creates help pages
cat("checking package...\n")
#check(package.dir) #package checking
check(package.dir, args = "--no-examples")
cat("loading package...\n")
load_all(package.dir) #building package
cat("library package...\n")
library(Gistic2GE)

#generate html pages
setwd(package.dir)
cat("creating html help pages for package...\n")
build_site(pkg = package.dir, examples = TRUE, launch = TRUE)

##install package
cat("installing package...\n")
install.packages(package.dir, repos = NULL, type = "source")


