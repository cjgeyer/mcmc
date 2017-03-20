
 library(mcmc, lib.loc = "../../package/mcmc.Rcheck")

 bits <- .Machine$sizeof.pointer * 8

 rdas <- list.files(pattern = "\\.rda$")

 myrda <- paste(bits, "bit.rda", sep = "")

 if (! myrda %in% rdas) {
     source("../../package/mcmc/tests/temp-par-witch.R")
     save.image(myrda)
 }

 if (length(rdas) == 2)
     cat("Ready to roll!\n")
