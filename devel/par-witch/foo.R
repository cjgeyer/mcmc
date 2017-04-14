
 library(mcmc, lib.loc = "../../package/mcmc.Rcheck")

 bits <- .Machine$sizeof.pointer * 8

 rdas <- list.files(pattern = "\\.rda$")

 myrda <- paste(bits, "bit.rda", sep = "")

 DEBUG <- TRUE

 source("../../package/mcmc/tests/temp-par-witch.R")

 if (myrda %in% rdas) {
     fred <- new.env(parent = emptyenv())
     load(myrda, fred)
     fredOK <- TRUE
     cat("checking objects saved in", myrda, "\n")
     for (sally in names(fred)) {
         cat("checking", sally, "... ")
         if (exists(sally)) {
             sally1.obj <- get(sally)
             sally2.obj <- fred[[sally]]
             if (inherits(sally1.obj, "tempering")) {
                 sally1.obj$time <- NULL
                 sally2.obj$time <- NULL
             }
             sallyOK <- identical(sally1.obj, sally2.obj)
             if (sallyOK) cat("OK\n") else cat("differs\n")
             fredOK <- fredOK && sallyOK
         } else {
             cat("absent\n")
         }
     }
     if (! fredOK) {
         rm(fred, fredOK, sally, sallyOK, sally1.obj, sally2.obj)
         cat("re-saving", myrda, "\n")
         save.image(myrda)
     }
 } else {
     save.image(myrda)
 }

 if (length(rdas) == 2)
     cat("Ready to roll!\n")
