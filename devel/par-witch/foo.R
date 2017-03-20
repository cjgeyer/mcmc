
 library(mcmc, lib.loc = "../../package/mcmc.Rcheck")

 bits <- .Machine$sizeof.pointer * 8

 rdas <- list.files(pattern = "\\.rda$")

 myrda <- paste(bits, "bit.rda", sep = "")

 if (! myrda %in% rdas) {
     source("../../package/mcmc/tests/temp-par-witch.R")
     save.image(myrda)
 }

 # redo of above, get one just written
 rdas <- list.files(pattern = "\\.rda$")

 if (length(rdas) == 2) {

     envir.rda <- sub("\\.rda$", "", rdas)
     envir.rda <- sub("^", "env", envir.rda)

     for (i in 1:2) {
         assign(envir.rda[i], new.env(parent = emptyenv()))
         load(rdas[i], get(envir.rda[i]))
     }
 }

 e1 <- get(envir.rda[1])
 e2 <- get(envir.rda[2])

 identical(e1$save.initial.seed, e1$out$initial.seed)
 identical(e2$save.initial.seed, e2$out$initial.seed)
 identical(e1$save.initial.seed, e2$save.initial.seed)

 options(digits=4) # avoid rounding differences

 # discrepancy one

 e1$out$acceptx
 e2$out$acceptx

 names(e1$out)
 dim(e1$out$initial)
 dim(e1$out$state)
 dim(e1$out$neighbors)

 # so there are 6 distributions and the state space is 3-dimensional

 # copied from vignettes/debug.Rnw
 #
 #         which: a logical vector of length niter, the type of update
 #                for each iteration: within component (TRUE) or
 #                swap components (FALSE)
 #    unif.which: a numeric vector of length niter, the Uniform(0, 1) random
 #                variate used to decide which type of update is done
 #         state: a numeric niter by ncomp by d array, the state before
 #                iteration i is state[i, , ]
 #      proposal: a numeric niter by d + 1 matrix, the proposal for
 #                iteration i is proposal[i, ] (explanation below)
 #    coproposal: a numeric niter by d + 1 matrix, the coproposal for
 #                iteration i is coproposal[i, ] (explanation below)
 #  log.hastings: a numeric vector of length niter, the logarithm of the
 #                Hastings ratio for each iteration
 # unif.hastings: a numeric vector of length niter, the Uniform(0, 1)
 #                random variate compared to the Hastings ratio for each
 #                iteration or NA if none is needed (when the log Hastings
 #                ratio is nonnegative)
 #       acceptd: a logical vector of length niter, the decision for each
 #                iteration: accept the proposal (TRUE) or reject it (FALSE)
 #          norm: a numeric niter by d matrix, the vector of standard normal
 #                random variates used to generate the proposal for
 #                iteration i is z[i, ] unless none are needed (for swap
 #                updates) when it is NA
 #   unif.choose: a numeric niter by 2 matrix, the vector of Uniform(0, 1)
 #                random variates used to choose the components to update
 #                in iteration i is unif.choose[i, ]; in a swap update two
 #                are used; in a within-component update only one is used
 #                and the second is NA
 #
 # In a within-component update, one component, say j, is chosen for
 # update.  The coproposal is the current value of the state for this
 # component, which is a vector of length d + 1, the first
 # component of which is j and the rest of which is state[i, j, ]
 # if we are in iteration i.
 # The proposal is a similar vector, the first
 # component of which is again j and the rest of which is a multivariate
 # normal random vector centered at state[i, j, ].
 # The coproposal is the current state; the proposal is the possible value
 # (if accepted) of the state at the next time.
 #
 # In a swap update, two components, say j1 and j2, are chosen for
 # update.  Strictly, speaking the coproposal is the pair of vectors
 # c(j1, state[i, j1, ]) and c(j2, state[i, j2, ])
 # and the proposal is these swapped, that is, the pair of vectors
 # c(j2, state[i, j1, ]) and c(j1, state[i, j2, ])
 # if we are in iteration i.
 # Since, however, there is a lot of redundant information here,
 # the vector c(j1, state[i, j1, ]) is output as coproposal[i, ]
 # and the vector c(j2, state[i, j2, ]) is output as proposal[i, ].

 k <- 1
 identical(e1$out$which[k], e2$out$which[k])
 if (e1$out$which[k]) "within" else "swap"
 identical(e1$out$unif.which[k], e2$out$unif.which[k])
 identical(e1$out$state[k, , ], e1$out$initial)
 identical(e2$out$state[k, , ], e2$out$initial)
 e1$out$state[k, , ]
 identical(e1$out$proposal[k, ], e2$out$proposal[k, ])
 identical(e1$out$coproposal[k, ], e2$out$coproposal[k, ])
 identical(e1$out$log.hastings[k], e2$out$log.hastings[k])
 identical(e1$out$unif.hastings[k], e2$out$unif.hastings[k])
 identical(e1$out$acceptd[k], e2$out$acceptd[k])
 e1$out$acceptd[k]
 identical(e1$out$norm[k, ], e2$out$norm[k, ])
 identical(e1$out$unif.choose[k, ], e2$out$unif.choose[k, ])

 # first iteration trivially accepts because states of two distributions
 # are same at start

 k <- 2
 identical(e1$out$state[k, , ], e1$out$state[k - 1, , ]) # because didn't move
 identical(e1$out$which[k], e2$out$which[k])
 if (e1$out$which[k]) "within" else "swap"
 identical(e1$out$unif.which[k], e2$out$unif.which[k])
 e1$out$unif.which[k]
 e2$out$unif.which[k]
 identical(e1$out$state[k, , ], e1$out$initial)
 identical(e2$out$state[k, , ], e2$out$initial)
 e1$out$state[k, , ]
 identical(e1$out$proposal[k, ], e2$out$proposal[k, ])
 identical(e1$out$coproposal[k, ], e2$out$coproposal[k, ])
 identical(e1$out$log.hastings[k], e2$out$log.hastings[k])
 identical(e1$out$unif.hastings[k], e2$out$unif.hastings[k])
 identical(e1$out$acceptd[k], e2$out$acceptd[k])
 e1$out$acceptd[k]
 identical(e1$out$norm[k, ], e2$out$norm[k, ])
 identical(e1$out$unif.choose[k, ], e2$out$unif.choose[k, ])
 e1$out$unif.choose[k, ]
 e2$out$unif.choose[k, ]

