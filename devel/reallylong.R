
library(mcmc, lib.loc = "../package/mcmc.Rcheck")

lud <- function(x) dt(x, df=3, log=TRUE)

set.seed(42)

out <- metrop(lud, 0, blen=1000, nbatch=1)
out$accept
out <- metrop(out, scale=4)
out$accept
out <- metrop(out, scale=6)
out$accept

out <- metrop(out, nbatch = 1000)
out$accept

out <- metrop(out, blen = 1e6)
out$accept
shapiro.test(as.vector(out$batch))

mout <- morph.metrop(lud, 0, blen=1000, nbatch=1, morph=morph(b=1))
mout$accept
mout <- morph.metrop(mout, scale=4)
mout$accept

mout <- morph.metrop(mout, nbatch = 1000)
mout$accept

mout <- morph.metrop(mout, blen = 1e6)
mout$accept
shapiro.test(as.vector(mout$batch))

save(out, mout, file = "reallylong.rda")

