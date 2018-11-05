# Script to render vignettes for source package.
# Not to be included in package build.
library(rmarkdown)

render("polyRADtutorial.Rmd") # html_vignette

# markdown version viewable on GitHub
render("polyRADtutorial.Rmd", 
       output_format = md_document("markdown_github"))

# manually add table of contents
mdlines <- readLines("polyRADtutorial.md")

mdlines <- c("# polyRAD tutorial", "", "## Table of Contents",
             "* [Introduction](#introduction)", 
             "* [Summary of available functions](#functions)",
             "* [Estimating genotype probabilities in a mapping population](#mapping)",
             "* [Estimating genotype probabilities in a diversity panel](#diversity)",
             "* [Considerations for RAM and processing time](#considerations)",
             "* [Citing polyRAD](#citation)", "",
             mdlines)

cat(mdlines, sep = "\n", file = "polyRADtutorial.md")
