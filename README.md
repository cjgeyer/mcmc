The main branch is morph (not master).  We will probably never go back
to master being the main branch.

This is the source tree for the R contributed package mcmc.
The version for users is at CRAN (http://cran.r-project.org/package=mcmc)

This package suggests packages `Iso` and `xtable`.  So if don't have need to
add to R.  In R do

    install.packages(c("Iso", "xtable"))

To check do

    cd package
    rm -f mcmc_*.tar.gz
    R CMD build mcmc
    R CMD check mcmc_*.tar.gz

Since one of the vignettes takes a long time, probably want

    R CMD check mcmc_*.tar.gz --no-vignettes

except for one last check before push.

