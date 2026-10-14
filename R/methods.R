## ---------------------------------------------------------------------------
## Methods for "sandwichrand" objects.
##
##   print()  -> the design and the first rows of the list
##   summary()-> design, arm balance and the verification table
##   plot()   -> the running imbalance, the signature of the sandwich
## ---------------------------------------------------------------------------

#' Methods for \pkg{sandwichrand} objects
#'
#' \code{print} shows the design and the head of the allocation list,
#' \code{summary} adds the arm balance and the verification table,
#' \code{as.data.frame} extracts the allocation list itself, and \code{plot}
#' draws the running imbalance.
#'
#' @param x,object An object of class \code{"sandwichrand"}.
#' @param ... Further arguments passed to or from other methods.
#' @param n Integer. Number of allocation rows echoed by \code{print}.
#' @param row.names Passed to \code{as.data.frame}.
#' @param optional Passed to \code{as.data.frame}.
#' @param y Ignored, present for compatibility with \code{plot}.
#' @param type Character. Which curve to draw: \code{"imbalance"} (the
#'   default) for the running imbalance, or \code{"cumulative"} for the
#'   cumulative deviation of each arm from its target.
#' @param xlab,ylab,main,col,lwd Graphical parameters.
#'
#' @return \code{as.data.frame} returns the allocation data frame; \code{format}
#'   returns a character vector; \code{summary} returns a
#'   \code{"summary.sandwichrand"} object; the print and plot methods return their
#'   input invisibly.
#'
#' @name sandwichrand-methods
#' @aliases print.sandwichrand
#' @aliases summary.sandwichrand
#' @aliases print.summary.sandwichrand
#' @aliases plot.sandwichrand
#' @aliases as.data.frame.sandwichrand
#' @aliases format.sandwichrand
NULL

.smr_header <- function(d) {
    c("SMR allocation sequence (sandwich mixed randomization)",
      sprintf("  N             = %d", d$N),
      sprintf("  arms          = %s   (ratio %s, ratio_sum = %d, %d arms)",
              paste(d$arms, collapse = " : "), d$ratio, d$ratio_sum,
              d$n_arms),
      sprintf("  layout        = %s   (path %d)", .smr_layout(d), d$path),
      if (d$path != 1L)
          sprintf("  first block   = n1 = %d, imbalance %s (target [%s, %s])",
                  d$n1, format(d$imbalance),
                  format(d$imbalance_target[1L]),
                  format(d$imbalance_target[2L])),
      if (d$path != 1L)
          sprintf("  middle        = %d block(s), %d slots, candidate sizes %s",
                  d$n_mid_blocks, d$n_mid,
                  paste(d$block.sizes, collapse = " ")),
      if (isTRUE(d$stream_used))
          "  random        = independent stream supplied by the caller"
      else
          sprintf("  seed          = %d   (%s)", d$seed, d$seed_note),
      sprintf("  block columns = %s",
              if (isTRUE(d$keep.blocks)) "kept"
              else "suppressed in $alloc (ICH E9 / GCP)"))
}

#' @rdname sandwichrand-methods
#' @export
format.sandwichrand <- function(x, ...) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    .smr_header(x$design)
}

#' @rdname sandwichrand-methods
#' @export
print.sandwichrand <- function(x, n = 6L, ...) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    cat(format(x), sep = "\n")
    n <- min(max(as.integer(n), 0L), nrow(x$alloc))
    if (n > 0L) {
        cat(sprintf("\nAllocation list (first %d of %d rows)\n",
                    n, nrow(x$alloc)))
        print(head(x$alloc, n), row.names = FALSE)
    }
    cat("\nArm balance\n")
    print(x$balance, row.names = FALSE)
    if (length(x$warnings)) {
        cat("\nWarnings\n")
        for (w in x$warnings) cat("  - ", w, "\n", sep = "")
    }
    invisible(x)
}

#' @rdname sandwichrand-methods
#' @export
summary.sandwichrand <- function(object, ...) {
    if (!inherits(object, "sandwichrand")) {
        stop("`object` must be an object of class \"sandwichrand\".",
             call. = FALSE)
    }
    structure(list(call = object$call,
                   design = object$design,
                   balance = object$balance,
                   blocks = object$blocks,
                   verification = smr_verify(object)),
              class = "summary.sandwichrand")
}

#' @rdname sandwichrand-methods
#' @export
print.summary.sandwichrand <- function(x, ...) {
    cat(.smr_header(x$design), sep = "\n")
    cat("\nArm balance\n")
    print(x$balance, row.names = FALSE)
    if (!is.null(x$blocks)) {
        cat("\nBlock structure\n")
        print(x$blocks, row.names = FALSE)
    }
    cat("\nVerification\n")
    print(x$verification)
    invisible(x)
}

#' @rdname sandwichrand-methods
#' @export
plot.sandwichrand <- function(x, y = NULL, type = "imbalance",
                           main = "Running imbalance",
                           xlab = "allocation number",
                           ylab = NULL, col = 1, lwd = 2, ...) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    type <- match.arg(type, c("imbalance", "cumulative"))
    a <- x$alloc
    d <- x$design
    N <- d$N
    if (type == "imbalance") {
        imb <- smr_imbalance(x)
        plot(seq_len(N), imb, type = "l", lwd = lwd, col = col,
             xlab = xlab, ylab = if (is.null(ylab)) "spread of deviations"
                 else ylab,
             main = main, ...)
        abline(h = d$imbalance_target, lty = 2, col = "grey50")
        mtext(sprintf("target window [%.0f, %.0f]",
                      d$imbalance_target[1L], d$imbalance_target[2L]),
              side = 3, line = 0.4, col = "grey40", cex = 0.7)
    } else {
        mult <- d$mult
        rs <- d$ratio_sum
        mi <- match(a$arm, d$arms)
        M <- matrix(0L, nrow = N, ncol = length(mult))
        M[cbind(seq_len(N), mi)] <- 1L
        cum <- apply(M, 2L, cumsum)          # N x n_arms cumulative counts
        target <- outer(seq_len(N), mult) / rs
        dev <- cum - target                   # N x n_arms
        matplot(seq_len(N), dev, type = "l", lwd = lwd, col = col,
                xlab = xlab,
                ylab = if (is.null(ylab)) "cumulative deviation" else ylab,
                main = main, ...)
        abline(h = 0, lty = 2, col = "grey50")
        legend("topright", legend = d$arms, col = col, lty = 1, bty = "n")
    }
    invisible(x)
}

#' @rdname sandwichrand-methods
#' @export
as.data.frame.sandwichrand <- function(x, row.names = NULL, optional = FALSE,
                                    ...) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    x$alloc
}
