initseq <- function(x) {
    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))
    .Call(C_initseq, x - mean(x))
}
