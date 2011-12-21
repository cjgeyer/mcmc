
morph <- function(type = c("none", "exponential", "subexponential"), R, p, b)
{
    type <- match.arg(type)
    if (type == "none")
        return(list(type = "none"))
    if (missing(R) || missing(p))
        stop("must specify arguments R and p for type", type)
    stopifnot(is.numeric(R))
    stopifnot(length(R) == 1)
    stopifnot(R >= 0)
    stopifnot(is.numeric(p))
    stopifnot(length(p) == 1)
    stopifnot(R > 2)
    if (type == "exponential")
        return(list(type = "exponential", R = R, p = p))
    if (missing(b))
        stop("must specify argument b for type subexponential")
    stopifnot(is.numeric(b))
    stopifnot(length(b) == 1)
    stopifnot(b > 0)
    return(list(type = "subexponential", R = R, p = p, b = b))
}

