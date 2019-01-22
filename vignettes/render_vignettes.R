# Script to render vignettes for source package.
# Not to be included in package build.
library(rmarkdown)

render("polyRADtutorial.Rmd") # html_vignette

# markdown version viewable on GitHub
render("polyRADtutorial.Rmd", 
       output_format = github_document(toc = TRUE, html_preview = FALSE))
