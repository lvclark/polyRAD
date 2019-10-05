# Script to render vignettes for source package.
# Not to be included in package build.
library(rmarkdown)

render("polyRADtutorial.Rmd") # html_vignette

# markdown version viewable on GitHub
render("polyRADtutorial.Rmd", 
       output_format = github_document(toc = TRUE, html_preview = FALSE))

mylines <- readLines("polyRADtutorial.md")
mylines <- gsub("$\\\\frac{ploidy - 1}{ploidy}$", "(ploidy - 1)/ploidy", fixed = TRUE,
                mylines)
mylines <- gsub("\\(\\frac{ploidy - 1}{ploidy}\\)", "(ploidy - 1)/ploidy", fixed = TRUE,
                mylines)
mylines <- gsub("\\(H_{ind}/H_E\\)", "*H*<sub>*i**n**d*</sub>/*H*<sub>*E*</sub>",
                fixed = TRUE, mylines)
writeLines(mylines, con = "polyRADtutorial.md")
