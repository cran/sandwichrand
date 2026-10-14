## ---------------------------------------------------------------------------
## The user-facing generator.
## ---------------------------------------------------------------------------

#' Allocation sequences by sandwich mixed randomization
#'
#' Draws a reproducible allocation sequence (randomization list) for a
#' randomized controlled trial with the sandwich mixed randomization (SMR)
#' scheme.
#'
#' The SMR sequence is a sandwich. Writing \code{n1} for the size of the first
#' block and \code{N} for the total sample size, the layout is
#'
#' \preformatted{
#'   |<--- n1 --->|<-------- N - 2*n1 -------->|<--- n1 --->|
#'   +------------+----------------------------+------------+
#'   | block 1    | blocks 2 .. K              | block K+1  |
#'   | imbalanced | permuted blocks of         | complement |
#'   |            | varying sizes              | of block 1 |
#'   +------------+----------------------------+------------+
#' }
#'
#' Block 1 and block K+1 are the two halves of one balanced pool of
#' \code{2*n1} slots, so the last block is the exact complement of the first
#' one: overall balance is guaranteed by construction, while the imbalance
#' deliberately introduced in block 1 keeps the next allocation unpredictable.
#'
#' The head-block size \code{n1} is always sampled at random from the
#' admissible multiples of \code{ratio_sum} in \code{[2*ratio_sum, N/2]}, so
#' there is no need to pre-specify it and the block structure is never revealed
#' in advance.
#'
#' Two further paths handle the boundary cases:
#' \describe{
#'   \item{path 1}{\code{N < 4 * ratio_sum}: complete randomization with a
#'     pre-determined \code{N}, in one single block.}
#'   \item{path 3}{\code{N == 2 * n1}: head and tail blocks only, no middle
#'     section.}
#' }
#'
#' @param N Integer. Total sample size.
#' @param arms Specification of the treatment arms and of the allocation ratio.
#'   Either a character vector with one token per allocation-ratio unit --
#'   \code{c("TestTrt", "Placebo")} for 1:1 and \code{c("A", "A", "B", "B",
#'   "C")} for 2:2:1, a single blank-separated string \code{"A A B B C"} is
#'   also accepted -- or a named vector of ratio units such as
#'   \code{c(TestTrt = 2, Placebo = 2, Control = 1)}. Arm identity is
#'   case-insensitive.
#' @param seed Integer. Random seed. \code{NULL} (the default) draws a seed
#'   from the system clock; the seed actually used is stored in the result and
#'   is what has to be archived to reproduce the list.
#' @param stream Numeric vector. The raw state of an independent
#'   L'Ecuyer-CMRG sub-stream (for example one of the values returned by
#'   \code{\link{smr_streams}}). Supplying \code{stream} is the most efficient
#'   way to give one call an independent random stream when the package is
#'   called at high frequency: each stream is non-overlapping with every
#'   other, so a loop never re-uses the same random numbers. At most one of
#'   \code{seed} and \code{stream} may be supplied.
#' @param block.sizes Integer vector. Candidate sizes of the middle blocks.
#'   Defaults to \code{c(2, 3, 4) * ratio_sum}.
#' @param imbalance Numeric vector of length two, \code{c(min, max)}. Admissible
#'   imbalance of the first block (default \code{c(1, 4)}; 4 is the safety
#'   factor proposed in the SMR paper). The imbalance is the \emph{spread} of
#'   the deviations from the target allocation, which for a 1:1 allocation is
#'   the between-group count difference and is zero for an exactly
#'   proportional block.
#' @param max.cand Integer. Maximum number of candidate values kept for
#'   \code{n1} (default 10).
#' @param max.iter Integer. Maximum number of rejection-sampling attempts per
#'   module (default 1000). If no candidate qualifies, the closest one is used,
#'   so a sequence is always produced.
#' @param keep.blocks Logical. \code{TRUE} (default) keeps the \code{block} and
#'   \code{block_size} columns in \code{$alloc}; \code{FALSE} drops them, as
#'   required by ICH E9 and good clinical practice when the list is released to
#'   the trial sites. The block structure is still kept in \code{$blocks},
#'   which is the copy for the statistician and is never written out by
#'   \code{\link{smr_write_csv}} or \code{\link{smr_rtf}} when
#'   \code{keep.blocks = FALSE}.
#' @param verbose Logical. Emit the design notes and warnings as messages
#'   (default \code{FALSE}; they are always stored in the returned object).
#'
#' @return An object of class \code{"sandwichrand"}, a list with the elements
#'   \describe{
#'     \item{\code{call}}{the matched call.}
#'     \item{\code{design}}{list describing the design actually used:
#'       \code{N}, \code{arms}, \code{ratio}, \code{ratio_sum}, \code{n_arms},
#'       \code{path}, \code{n1}, \code{n_mid}, \code{n_mid_blocks},
#'       \code{imbalance_observed}, \code{imbalance_target}, \code{seed}, ...}
#'     \item{\code{alloc}}{data frame with columns \code{id}, \code{block},
#'       \code{block_size} (the latter two only when \code{keep.blocks = TRUE})
#'       and \code{arm}.}
#'     \item{\code{balance}}{data frame with \code{arm}, \code{n},
#'       \code{expected} and \code{deviation}.}
#'     \item{\code{blocks}}{data frame with one row per block: \code{block},
#'       \code{size}, \code{first_id}, \code{last_id} and the arm counts.}
#'     \item{\code{notes}, \code{warnings}}{character vectors.}
#'     \item{\code{spec}}{the normalized specification, used by
#'       \code{\link{smr_reproduce}}.}
#'   }
#'
#' @references
#' Wang B, Wang X, Ni C, Xie J, Wu S, Wu Y, Shen Y, Dong Z, Lv W, He Y (2026).
#' \dQuote{Reducing selection bias risk to enhance RCT validity: sandwich mixed randomization outperforms permuted block design.}
#' \emph{BMC Medical Research Methodology}, \bold{26}(1), 152.
#' \doi{10.1186/s12874-026-02877-2}.
#'
#' @seealso \code{\link{smr_verify}}, \code{\link{smr_balance}},
#'   \code{\link{smr_imbalance}}, \code{\link{smr_blocks}},
#'   \code{\link{smr_rtf}}, \code{\link{smr_write_csv}},
#'   \code{\link{smr_streams}}.
#'
#' @examples
#' ## Two-arm 1:1 trial, N = 48 -- the illustration of the SMR paper
#' x <- smr_alloc(N = 48, arms = c("TestTrt", "Placebo"), seed = 84529052)
#' x
#' head(x$alloc)
#' x$balance
#'
#' ## Three-arm 2:2:1 trial, N = 150, released without the block structure
#' y <- smr_alloc(N = 150, arms = c("A", "A", "B", "B", "C"), seed = 123456,
#'                keep.blocks = FALSE)
#' names(y$alloc)
#' smr_balance(y)
#'
#' ## Small trial: path 1, complete randomization in one block
#' z <- smr_alloc(N = 9, arms = c("Drug", "Placebo"), seed = 777)
#' z$design$path
#'
#' ## High-frequency calling: derive independent streams up front so that no
#' ## two calls ever draw from the same random stream
#' ss <- smr_streams(seed = 42, n = 3)
#' head(smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[1]])$alloc)
#' head(smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[2]])$alloc)
#'
#' @export
smr_alloc <- function(N,
                      arms,
                      seed         = NULL,
                      stream       = NULL,
                      block.sizes  = NULL,
                      imbalance    = c(1, 4),
                      max.cand     = 10L,
                      max.iter     = 1000L,
                      keep.blocks  = TRUE,
                      verbose      = FALSE) {

    cl <- match.call()

    if (missing(arms)) {
        stop("`arms` is required, e.g. c(\"TestTrt\", \"Placebo\") for 1:1 ",
             "or c(\"A\", \"A\", \"B\", \"B\", \"C\") for 2:2:1.", call. = FALSE)
    }
    if (!is.logical(keep.blocks) || length(keep.blocks) != 1L ||
        is.na(keep.blocks)) {
        stop("`keep.blocks` must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
    }

    spec <- .smr_spec(N = N, arms = arms, seed = seed, stream = stream,
                      block.sizes = block.sizes, imbalance = imbalance,
                      max.cand = max.cand, max.iter = max.iter)

    notes <- spec$notes
    warns <- spec$warnings

    ## Preserve the caller's RNG state across the WHOLE call. Stream derivation
    ## below calls set.seed()/RNGkind(), which mutate the global state; this
    ## guard restores it on exit so a call to smr_alloc is invisible to the
    ## session's random numbers.
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

    ## The base state is either the supplied `stream` or one derived from the
    ## `seed`. Four independent sub-streams, one per module, are derived from it
    ## so that the modules of one call are independent of one another and every
    ## call is a deterministic function of its base.
    base <- if (is.null(spec$stream)) .stream_for_call(spec$seed)
            else spec$stream
    st_sel  <- nextRNGStream(base)
    st_m1   <- nextRNGStream(st_sel)
    st_m2   <- nextRNGStream(st_m1)
    st_p1   <- nextRNGStream(st_m2)

    gen <- .with_rng_state(base,
                           .smr_generate(spec, st_sel, st_m1, st_m2, st_p1))
    notes <- c(notes, gen$notes)
    warns <- c(warns, gen$warns)
    out <- list(alloc = gen$alloc, design = gen$design)

    ## =======================================================================
    ## Balance and internal verification
    ## =======================================================================
    alloc_full <- out$alloc
    if (nrow(alloc_full) != spec$N) {
        stop("internal check failed: the sequence has ", nrow(alloc_full),
             " rows but N = ", spec$N, ".", call. = FALSE)
    }

    n_obs <- .smr_tabulate(alloc_full$arm, spec$arms)
    expected <- spec$N * spec$mult / spec$ratio_sum
    balance <- data.frame(arm = spec$arms,
                          n = n_obs,
                          expected = expected,
                          deviation = n_obs - expected,
                          stringsAsFactors = FALSE)
    dev_max <- max(abs(balance$deviation))
    if (dev_max > 1) {
        warns <- c(warns,
                   sprintf(paste("the largest deviation from the target",
                                 "allocation is %s (expected <= 1)."),
                           format(dev_max)))
    }

    blocks <- .smr_block_table(alloc_full)

    design <- c(out$design,
                list(arms = spec$arms,
                     ratio = .smr_ratio_label(spec$mult),
                     ratio_sum = spec$ratio_sum,
                     n_arms = spec$n_arms,
                     mult = spec$mult,
                     tmap = spec$tmap,
                     block.sizes = spec$block.sizes,
                     imbalance_target = c(spec$ndiffmin, spec$ndiffmax),
                     max.cand = spec$max.cand,
                     max.iter = spec$max.iter,
                     keep.blocks = keep.blocks,
                     seed = spec$seed,
                     stream_used = !is.null(spec$stream),
                     seed_note = spec$seed_note,
                     dev_max = dev_max))

    alloc <- if (isTRUE(keep.blocks)) {
        alloc_full
    } else {
        alloc_full[, c("id", "arm"), drop = FALSE]
    }

    res <- structure(list(call = cl,
                          design = design,
                          alloc = alloc,
                          balance = balance,
                          blocks = blocks,
                          notes = notes,
                          warnings = warns,
                          spec = spec),
                     class = "sandwichrand")

    if (isTRUE(verbose)) {
        for (m in notes) message("NOTE: smr_alloc - ", m)
        for (m in warns) message("WARNING: smr_alloc - ", m)
        message("NOTE: smr_alloc - rows produced = ", nrow(alloc),
                " (matches N); max deviation = ", format(dev_max))
    }
    res
}

## ---------------------------------------------------------------------------
## Derive n independent L'Ecuyer-CMRG sub-streams from one base seed.
##
## This is the simple, efficient way to call smr_alloc many times without
## ever re-using the same random stream: the streams are non-overlapping
## (each about 2^127 draws long), so the k-th call is statistically
## independent of every other even when they all share the same base seed.
## ---------------------------------------------------------------------------
#' Independent random streams for high-frequency calling
#'
#' Derives \code{n} non-overlapping L'Ecuyer-CMRG sub-streams from one base
#' seed. Each returned value is the raw RNG state of a stream and can be passed
#' directly to \code{\link{smr_alloc}} via its \code{stream} argument. Because
#' the streams do not overlap, a loop that hands one stream per call never
#' re-uses the same random numbers, which is what keeps high-frequency calling
#' safe and reproducible.
#'
#' @param seed Integer. Base random seed.
#' @param n Integer. Number of independent streams to derive.
#'
#' @return A list of \code{n} elements, each the numeric state of an
#'   independent L'Ecuyer-CMRG stream.
#'
#' @examples
#' ss <- smr_streams(seed = 1, n = 3)
#' a <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[1]])
#' b <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[2]])
#' ## the two sequences are different even though they share the base seed
#' !identical(a$alloc, b$alloc)
#'
#' @seealso \code{\link{smr_alloc}}, \code{\link{smr_simulate}}.
#'
#' @export
smr_streams <- function(seed, n) {
    seed <- suppressWarnings(as.numeric(seed))
    if (length(seed) != 1L || !is.finite(seed)) {
        stop("`seed` must be a single number (received ",
             paste(format(seed), collapse = " "), ").", call. = FALSE)
    }
    n <- suppressWarnings(as.numeric(n))
    if (length(n) != 1L || !is.finite(n) || n < 1) {
        stop("`n` must be a positive number (received ",
             paste(format(n), collapse = " "), ").", call. = FALSE)
    }
    .smr_streams(as.integer(round(seed)), as.integer(round(n)))
}

## ---------------------------------------------------------------------------
## Regenerate a sequence from its stored specification.
## ---------------------------------------------------------------------------

#' Reproduce an allocation sequence
#'
#' Re-runs the generator with the specification and the seed stored in an
#' \code{"sandwichrand"} object. The result is identical to the original object,
#' which is the property an auditor needs: the seed plus the package version
#' fully determine the list.
#'
#' @param x An object of class \code{"sandwichrand"}.
#'
#' @return A new object of class \code{"sandwichrand"}, identical to \code{x}.
#'
#' @examples
#' x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 2026)
#' identical(smr_reproduce(x)$alloc, x$alloc)
#'
#' @export
smr_reproduce <- function(x) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    s <- x$spec
    ## `arms_input` is the arm specification exactly as the user supplied it;
    ## `arms` holds only the arm labels and would lose the allocation ratio.
    y <- smr_alloc(N = s$N, arms = s$arms_input, seed = s$seed,
                   block.sizes = s$block.sizes,
                   imbalance = s$imbalance, max.cand = s$max.cand,
                   max.iter = s$max.iter,
                   keep.blocks = isTRUE(x$design$keep.blocks),
                   verbose = FALSE)
    y$call <- x$call
    y
}

## ---------------------------------------------------------------------------
## One row per block (the statistician's copy of the block structure)
## ---------------------------------------------------------------------------
.smr_block_table <- function(a) {
    if (!all(c("block", "block_size") %in% names(a))) return(NULL)
    b <- unique(a$block)
    first <- vapply(b, function(z) min(a$id[a$block == z]), integer(1L))
    last <- vapply(b, function(z) max(a$id[a$block == z]), integer(1L))
    size <- vapply(b, function(z) a$block_size[a$block == z][1L], integer(1L))
    arms <- sort(unique(a$arm))
    tab <- data.frame(block = as.integer(b),
                      size = as.integer(size),
                      first_id = as.integer(first),
                      last_id = as.integer(last),
                      stringsAsFactors = FALSE)
    cnt <- matrix(0L, nrow = length(b), ncol = length(arms),
                  dimnames = list(NULL, paste0("n_", arms)))
    for (i in seq_along(b)) {
        sub <- a$arm[a$block == b[i]]
        cnt[i, ] <- .smr_tabulate(sub, arms)
    }
    cbind(tab, as.data.frame(cnt, check.names = FALSE))
}
