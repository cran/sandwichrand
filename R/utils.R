## ---------------------------------------------------------------------------
## Internal helpers: independent random streams, shuffling, apportionment,
## imbalance measures.
##
## These are the R counterparts of the temporary-array mechanics of the SAS
## macro SanMR. They are deliberately kept free of any package state so that
## every sequence is a pure, reproducible function of its seed (or of the
## stream that is handed to it).
## ---------------------------------------------------------------------------

## Run `expr` with the RNG state `state` and restore the caller's state on
## exit. This is the R counterpart of CALL STREAMINIT in the SAS macro: every
## module consumes its own sub-stream and never touches the stream of another
## module (or of the user's session).
##
## `state` is an integer vector (the raw .Random.seed) of a L'Ecuyer-CMRG
## stream. Assigning it to .Random.seed fully determines the subsequent draws.
.with_rng_state <- function(state, expr) {
    had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        if (had) {
            assign(".Random.seed", old, envir = .GlobalEnv, inherits = FALSE)
        } else if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }
    }, add = TRUE)
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
    force(expr)
}

## Initialise the session RNG to the L'Ecuyer-CMRG generator (the only
## generator for which parallel::nextRNGStream() guarantees non-overlapping
## sub-streams). Must be called before nextRNGStream is used.
.ensure_lecuyer <- function() {
    kind <- RNGkind()
    if (!identical(kind[1L], "L'Ecuyer-CMRG")) {
        RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    }
    invisible(kind)
}

## Derive a fresh L'Ecuyer-CMRG stream from `seed` (an integer seed).
##
## The returned value is the raw .Random.seed state of the stream. It is the
## first sub-stream, i.e. equivalent to `set.seed(seed); .Random.seed`.
.stream_from_seed <- function(seed) {
    .ensure_lecuyer()
    set.seed(seed)
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

## Derive `n` non-overlapping L'Ecuyer-CMRG sub-streams from `seed`.
##
## This is the simple, efficient way to give each of `n` independent calls its
## own stream: the streams are derived once, up front, and are statistically
## independent (each is roughly 2^127 draws long). It is what makes
## high-frequency calling safe: the k-th stream is completely different from
## the (k-1)-st, so a loop never re-uses the same random numbers even when the
## same base `seed` is used throughout.
.smr_streams <- function(seed, n) {
    n <- as.integer(n)
    .ensure_lecuyer()
    set.seed(seed)
    st <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    out <- vector("list", n)
    out[[1L]] <- st
    for (i in seq_len(n - 1L)) {
        st <- nextRNGStream(st)
        out[[i + 1L]] <- st
    }
    out
}

## A single stream for one call, derived from `seed`.
##
## `key` (>= 1) selects which sub-stream of the `seed`-derived family is used,
## mirroring the SAS macro's `stream_id` convention. key = 1 is the default
## (the stream straight from the seed). The same (seed, key) always yields the
## same stream, so generation is reproducible; different keys yield
## independent streams.
.stream_for_call <- function(seed, key = 1L) {
    key <- as.integer(key)
    .ensure_lecuyer()
    set.seed(seed)
    st <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    for (i in seq_len(pmax(0L, key - 1L))) st <- nextRNGStream(st)
    st
}

## Run `expr` using the RNG `state`. If `state` is NULL, `seed` is used to
## build a state first. Returns the value of `expr` and always restores the
## caller's RNG state.
.with_rng <- function(state, expr, seed = NULL) {
    if (is.null(state)) state <- .stream_for_call(seed)
    .with_rng_state(state, expr)
}

## Fisher-Yates shuffle.
##
## All uniforms are drawn in one runif() call, in exactly the order the
## elementwise loop would consume them (i = n, n-1, ..., 2), so the result is
## identical to the naive loop but roughly an order of magnitude faster for
## long vectors. runif() can in principle return 1, so the index is clamped
## exactly as the SAS code clamps FLOOR(RAND('UNIFORM') * i).
.smr_shuffle <- function(x) {
    n <- length(x)
    if (n > 1L) {
        i <- n:2L
        j <- 1L + floor(runif(length(i)) * i)
        over <- j > i
        if (any(over)) j[over] <- i[over]
        for (k in seq_along(i)) {
            ii <- i[k]
            jj <- j[k]
            tmp <- x[ii]
            x[ii] <- x[jj]
            x[jj] <- tmp
        }
    }
    x
}

## Largest-remainder apportionment of `b` slots over the treatment arms.
##
## Arm g receives floor(b * mult_g / ratio_sum) slots; the slots left over go
## to the largest remainders, ties broken at random. Filling a block by
## cycling through the allocation pattern instead (as earlier versions did)
## piles every leftover slot onto the first arm and leaves both the block and
## the whole sequence further from the target.
.smr_apportion <- function(b, mult, ratio_sum) {
    ng <- length(mult)
    base <- floor(b * mult / ratio_sum)
    left <- b - sum(base)
    if (left > 0L) {
        frac <- b * mult - base * ratio_sum
        key <- runif(ng)
        ord <- order(-frac, -key)
        take <- ord[seq_len(left)]
        base[take] <- base[take] + 1L
    }
    as.integer(base)
}

## Spread of the deviations from the target allocation:
##
##   spread = max_g(count_g - n * mult_g / ratio_sum)
##            - min_g(count_g - n * mult_g / ratio_sum)
##
## For a 1:1 allocation this is the difference of the two arm sizes, i.e.
## exactly the imbalance used in the SMR paper; for unequal ratios it is the
## natural generalisation (spread = 0 means exactly proportional).
.smr_spread <- function(cnt, n, mult, ratio_sum) {
    dev <- cnt - n * mult / ratio_sum
    max(dev) - min(dev)
}

## Is a whole number? (tolerates the floating-point representation)
.smr_is_whole <- function(x, tol = 1e-8) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && abs(x - round(x)) < tol
}

## Whole-number test over a vector.
.smr_is_whole_set <- function(x, tol = 1e-8) {
    is.numeric(x) && all(is.finite(x)) && all(abs(x - round(x)) < tol)
}

## Frequency table of an allocation, aligned on the declared arm order.
.smr_tabulate <- function(x, levels) {
    as.integer(tabulate(match(as.character(x), levels),
                        nbins = length(levels)))
}

## Compact "n1 + blocks + n1" description of a design, used by print()/summary().
.smr_layout <- function(d) {
    if (d$path == 1L) return(sprintf("%d", d$N))
    if (d$path == 3L) return(sprintf("%d + %d", d$n1, d$n1))
    sprintf("%d + %d x %s + %d", d$n1, d$n_mid_blocks,
            paste(d$block.sizes, collapse = "/"), d$n1)
}

## Human-readable allocation ratio, e.g. "2:2:1".
.smr_ratio_label <- function(mult) paste(mult, collapse = ":")
