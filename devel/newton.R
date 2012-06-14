
fred <- quote(x + (x - r)^p)
sally <- D(fred, "x")

newpath <- function(y, r, p) {
    stopifnot(is.atomic(y))
    stopifnot(is.numeric(y))
    stopifnot(is.finite(y))
    stopifnot(length(y) == 1)
    stopifnot(is.atomic(r))
    stopifnot(is.numeric(r))
    stopifnot(is.finite(r))
    stopifnot(length(r) == 1)
    stopifnot(is.atomic(p))
    stopifnot(is.numeric(p))
    stopifnot(is.finite(p))
    stopifnot(length(p) == 1)
    stopifnot(r >= 0)
    stopifnot(p > 2)
    stopifnot(y >= r)

    x <- r + y^(1 / p)
    xlist <- NULL
    slist <- NULL
    elist <- NULL
    repeat {
        err <- y - eval(fred)
        step <- err / eval(sally)
        xlist <- c(xlist, x)
        slist <- c(slist, step)
        elist <- c(elist, err)
        x <- x + step
        if (abs(err) < sqrt(.Machine$double.eps)) break
    }
    return(data.frame(x = xlist, error = elist, step = slist))
}

newpath(y = 5, r = 5, p = 2.5)
newpath(y = 10, r = 5, p = 2.5)
newpath(y = 20, r = 5, p = 2.5)
newpath(y = 40, r = 5, p = 2.5)
newpath(y = 80, r = 5, p = 2.5)
newpath(y = 160, r = 5, p = 2.5)
newpath(y = 320, r = 5, p = 2.5)
newpath(y = 640, r = 5, p = 2.5)
newpath(y = 1280, r = 5, p = 2.5)
newpath(y = 2560, r = 5, p = 2.5)
newpath(y = 5120, r = 5, p = 2.5)
newpath(y = 10240, r = 5, p = 2.5)

newpath(y = 5, r = 5, p = 5.5)
newpath(y = 10, r = 5, p = 5.5)
newpath(y = 20, r = 5, p = 5.5)
newpath(y = 40, r = 5, p = 5.5)
newpath(y = 80, r = 5, p = 5.5)
newpath(y = 160, r = 5, p = 5.5)
newpath(y = 320, r = 5, p = 5.5)
newpath(y = 640, r = 5, p = 5.5)
newpath(y = 1280, r = 5, p = 5.5)
newpath(y = 2560, r = 5, p = 5.5)
newpath(y = 5120, r = 5, p = 5.5)
newpath(y = 10240, r = 5, p = 5.5)


