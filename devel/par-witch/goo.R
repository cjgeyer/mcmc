
 e32 <- new.env(parent = emptyenv())
 e64 <- new.env(parent = emptyenv())

 load(envir = e32, file = "32bit.rda")
 load(envir = e64, file = "64bit.rda")

 identical(names(e32), names(e64))

 for (n in names(e32)) {
     cat(n)
     o32 <- e32[[n]]
     o64 <- e64[[n]]
     if (inherits(o32, "tempering")) {
         o32$time <- NULL
         o64$time <- NULL
         if (identical(o32, o64)) cat(" same except for time component\n") else
             cat(" different\n");
     } else {
         if (identical(o32, o64)) cat(" same\n") else cat(" different\n");
     }
 }
 
 o32 <- e32$out
 o64 <- e64$out

 identical(names(o32), names(o64))
 names(o32)
 identical(dim(o32$state), dim(o64$state))
 dim(o32$state)

 for (i in 1:dim(o32$state)[1]) {
     state32 <- o32$state[i, , ]
     state64 <- o64$state[i, , ]
     if (! identical(state32, state64)) break
 }

    ## regular output (if parallel)
    ##
    ##     batch      (if outfun) nbatch x nout matrix
    ##                (if no outfun) nbatch x ncomp x nx array
    ##     acceptx    vector of length ncomp
    ##     accepti    ncomp x ncomp matrix
    ##     initial    copy of initial state
    ##     final      final state
    ##
    ## debug output (if parallel)
    ##
    ##     which         vector of length niter (TRUE if within-component)
    ##     unif_which    vector of length niter (uniform for deciding which)
    ##     state         niter x ncomp x nx array (state before)
    ##     coproposal    niter x (nx + 1) matrix
    ##     proposal      niter x (nx + 1) matrix
    ##     log_hastings  niter vector
    ##     unif_hastings niter vector (uniform for deciding acceptance)
    ##     acceptd       niter vector (TRUE if accept)
    ##     norm          niter x nx matrix (std normals)
    ##     unif_choose   niter x 2 matrix (uniforms for choosing components
    ##                       to update)
    ##
    ##     for within-component move coproposal and proposal have natural
    ##         meaning
    ##     for swap move we have 2 coproposals (i, x_i) and (j, x_j)
    ##         and 2 proposals (i, x_j) and (j, x_i) but since the information
    ##         here is quite redundant we just store (i, x_i) in "coproposal"
    ##         and (j, x_j) in "proposal" -- the checker can figure it out

 for (i in 1:length(o32$which)) {
     if (! identical(o32$which[i], o64$which[i])) break
     if (! identical(o32$unif.which[i], o64$unif.which[i])) break
     if (! identical(o32$state[i, , ], o64$state[i, , ])) break
     if (! identical(o32$coproposal[i, ], o64$coproposal[i, ])) break
     if (! identical(o32$log.hastings[i], o64$log.hastings[i])) break
 }
 i

 identical(o32$state[i, , ], o64$state[i, , ])
 # state agrees at start of iteration i
 identical(o32$unif.which[i], o64$unif.which[i])
 o32$unif.which[i]
 identical(o32$which[i], o64$which[i])
 o64$which[i]
 # which kind of update agrees, will attempt swap move
 identical(o32$proposal[i, ], o64$proposal[i, ])
 identical(o32$coproposal[i, ], o64$coproposal[i, ])
 o32$proposal[i, ]
 o32$coproposal[i, ]
 # we are attempting to swap distributions 6 and 4 of current state
 # swap must succeed because they are same position
 o32$state[i, , ]
 o32$state[i, , ]
 identical(o32$log.hastings[i], o64$log.hastings[i])
 ##### what ?????
 o32$log.hastings[i]
 o64$log.hastings[i]
 # that's crazy.  I though it was only the 32 bit chips that could
 # do stuff like that
 identical(o32$lud, o64$lud)
 o32$lud
 ncomp <- dim(o32$state)[2]
 ncomp
 d <- dim(o32$state)[3]
 d
 # copied from mcmc/tests/temp-par-witch.R
 witch.which <- 1 - (1 / 2)^(1 / d) * (1 / 4)^(seq(0, 5) / d)
 witch.which

 o32$lud(o32$proposal[i, ])
 o32$lud(o32$coproposal[i, ])

 # we've seen enough to see that purely from the vagarities of computer
 # arithmetic the RNG streams get out of phase

 o32$unif.hastings[i]
 o64$unif.hastings[i]

 # This is correct!  When the log hastings ratio is >= 1 we automatically
 # accept the proposal (in this case a swap) and don't need to use a
 # uniform(0, 1) output of the RNG system.    When the log hastings ratio
 # is < 1 we cannot automatically accept but must accept with probability
 # exp(log_hastings) and we do need to use a
 # uniform(0, 1) output of the RNG system for that.
 # Now that we are out of phase, the future of the Markov chain will
 # be quite different.

 # so just leave temp-par-witch.Rout.save out of the package as we
 # did before

