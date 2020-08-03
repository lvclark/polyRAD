## Delete this section

Thank you for contributing to polyRAD!  Below is a template for your pull
request to help both of us make sure that all the needed commits have
been made before it is merged, and to document why the changes were made.
Please edit it as appropriate.  You can always make the pull request,
make a few more commits, and then edit this document further.

## Purpose

Use this section to describe the purpose of your pull request.
Was there a bug that needed to be fixed?  Is there a reason why a new feature
would be helpful or why the documentation needed to be updated?

## Approach

What did you do to address the problem?  What sorts of changes were made
to the code and/or documentation, and why?  Note that you don't need to
put too much (or possibly any) code here since your commits are also
visible.

## Checklist

You can still make the pull request if you haven't completed this checklist,
but you may be asked to make additional commits before the request is merged.
You can delete or strikethough some of these items if you are certain they are
not necessary.  Feel free also to add any items that you didn't already discuss
in the above section.

- [x] Here is an example of a checked box.
- [ ] Code style matches polyRAD in terms of argument names, line length
(try not to go over 80 characters) and whitespace.
- [ ] Documentation in the `man` folder has been updated to reflect changes in
the code, including new functions and arguments.
- [ ] The list of functions in `vignettes/polyRADtutorial.Rmd` has been updated
to add any new functions.
- [ ] If vignettes have been updated, `vignettes/render_vignettes.R` has been
run (with `vignettes` as the working directory).
- [ ] New functions and/or methods have been listed in `NAMESPACE`.
- [ ] If this is the first functionality added since the last release, 
the minor version number (_e.g._ 2.4 -> 2.5) in `DESCRIPTION` has been updated.
- [ ] If this is the first bug fix since the last release, and the minor
version number has not increased, the patch version number (_e.g._ 2.4 -> 2.4.1
or 2.4.16 -> 2.4.17) in `DESCRIPTION` has been updated.
- [ ] `NEWS` has been updated to summarize the changes.  Be sure to put the
changes under the correct version number, according to the above two points.
- [ ] If you are contributing anything more substantial than a typo fix, you
have added yourself as a contributor with the "ctb" role in `Author` and
`Authors@R` in the `DESCRIPTION` file.
- [ ] The date in `DESCRIPTION` has been updated.
- [ ] `README.md` has been updated if necessary, particularly if any new input
or output formats are supported.
- [ ] If any Rcpp functions have been added or their arguments modified,
`Rcpp::compileAttributes` has been run on the package.
- [ ] After running `R CMD build polyRAD` with the current version of R, the package
passes `R CMD check polyRAD_x.x.x.tar.gz` (both of these commands are run from the
terminal, one directory up from the polyRAD directory).  If there are any errors,
warnings, or notes that you can't resolve (or think are not your fault), please
list them in this document.
