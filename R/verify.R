## ---------------------------------------------------------------------------
## Verification, balance diagnostics and assessment.
## ---------------------------------------------------------------------------

#' Verify the structural invariants of an allocation sequence
#'
#' Re-checks every invariant that the SMR scheme promises: the length and the
#' numbering of the sequence, the final balance, and -- when the block
#' structure is available -- the consistency of the blocks and the
#' complementarity of the first and the last block. The function works on any
#' \code{"sandwichrand"} object, so a list received from a third party can be
#' re-checked independently of how it was produced.
#'
#' @param x An object of class \code{"sandwichrand"}.
#' @param stop.on.fail Logical. Stop with an error when at least one check
#'   fails, instead of returning the table (default \code{FALSE}).
#'
#' @return A data frame of class \code{"smr_check"} with one row per check and
#'   the columns \code{check}, \code{status} (\code{"PASS"} or \code{"FAIL"})
#'   and \code{detail}. Invisibly the same value when \code{stop.on.fail =
#'   TRUE} and every check passes.
#'
#' @examples
#' x <- smr_alloc(N = 48, arms = c("TestTrt", "Placebo"), seed = 84529052)
#' v <- smr_verify(x)
#' all(v$status == "PASS")
#'
#' @export
smr_verify <- function(x, stop.on.fail = FALSE) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    if (!is.logical(stop.on.fail) || length(stop.on.fail) != 1L ||
        is.na(stop.on.fail)) {
        stop("`stop.on.fail` must be TRUE or FALSE.", call. = FALSE)
    }

    a <- x$alloc
    d <- x$design
    N <- d$N
    res <- list()
    add <- function(id, ok, detail) {
        res[[length(res) + 1L]] <<- data.frame(
            check = id,
            status = if (isTRUE(ok)) "PASS" else "FAIL",
            detail = detail,
            stringsAsFactors = FALSE)
    }

    add("rows", nrow(a) == N,
        sprintf("rows = %d, expected %d", nrow(a), N))
    add("arms_complete",
        !any(is.na(a$arm)) && all(nzchar(as.character(a$arm))) &&
            all(a$arm %in% d$arms),
        sprintf("unknown or missing arm values = %d",
                sum(is.na(a$arm) | !nzchar(as.character(a$arm)) |
                        !(a$arm %in% d$arms))))
    add("ids", !any(is.na(a$id)) && length(unique(a$id)) == N &&
            min(a$id) == 1L && max(a$id) == N,
        sprintf("distinct ids = %d, range [%s, %s], expected %d",
                length(unique(a$id)), format(min(a$id)), format(max(a$id)), N))
    add("balance", max(abs(x$balance$deviation)) <= 1 + 1e-6,
        sprintf("max |deviation| = %s",
                format(max(abs(x$balance$deviation)))))
    add("seed_or_stream",
        isTRUE(d$stream_used) ||
            (length(d$seed) == 1L && is.finite(d$seed) && !is.na(d$seed)),
        sprintf("seed = %s, stream_used = %s",
                format(d$seed), isTRUE(d$stream_used)))

    has_blocks <- !is.null(x$blocks)

    if (has_blocks && d$path == 1L) {
        add("single_block", nrow(x$blocks) == 1L,
            sprintf("number of blocks = %d, expected 1 (path 1)",
                    nrow(x$blocks)))
    } else if (has_blocks) {
        bt <- x$blocks
        add("blocks_sum", sum(bt$size) == N,
            sprintf("sum of block sizes = %d, expected %d", sum(bt$size), N))
        add("blocks_contiguous",
            all(bt$first_id[-1] == bt$last_id[-nrow(bt)] + 1L) &&
                bt$first_id[1L] == 1L && bt$last_id[nrow(bt)] == N,
            sprintf("blocks cover ids %d..%d contiguously",
                    bt$first_id[1L], bt$last_id[nrow(bt)]))
        b1 <- bt$block[1L]
        bK <- bt$block[nrow(bt)]
        add("n1_equal", bt$size[1L] == bt$size[nrow(bt)],
            sprintf("head block size = %d, tail block size = %d",
                    bt$size[1L], bt$size[nrow(bt)]))
        add("n1_min", bt$size[1L] >= 2L * d$ratio_sum,
            sprintf("head block size = %d, minimum = %d",
                    bt$size[1L], 2L * d$ratio_sum))

        ## head + tail complementarity (path 2 or 3)
        cols <- paste0("n_", d$arms)
        if (all(cols %in% names(bt))) {
            cnt1 <- unlist(bt[1L, cols])
            cntK <- unlist(bt[nrow(bt), cols])
            combined <- cnt1 + cntK
            n1 <- bt$size[1L]
            target2 <- 2L * n1 * d$mult / d$ratio_sum
            dev <- combined - target2
            add("complementarity",
                all(abs(dev) < 0.5 + 1e-9),
                sprintf("combined head+tail deviation from target = [%s]",
                        paste(sprintf("%.3f", dev), collapse = ", ")))
        } else {
            add("complementarity", FALSE, "arm-count columns missing")
        }
    } else {
        add("blocks_present", FALSE, "block structure not stored")
    }

    df <- do.call(rbind, res)
    rownames(df) <- NULL
    class(df) <- c("smr_check", "data.frame")
    if (isTRUE(stop.on.fail) && any(df$status == "FAIL")) {
        stop("smr_verify: the following checks FAILED: ",
             paste(df$check[df$status == "FAIL"], collapse = ", "),
             call. = FALSE)
    }
    invisible(df)
}

#' @export
print.smr_check <- function(x, ...) {
    cat("SMR verification report\n")
    for (i in seq_len(nrow(x))) {
        cat(sprintf("  [%s] %s : %s\n", x$status[i], x$check[i], x$detail[i]))
    }
    invisible(x)
}

#' Balance summary of an allocation sequence
#'
#' Returns the arm counts, the expected counts under the target allocation and
#' the deviation of each arm from its target.
#'
#' @param x An object of class \code{"sandwichrand"}.
#'
#' @return A data frame with columns \code{arm}, \code{n}, \code{expected} and
#'   \code{deviation}.
#'
#' @examples
#' x <- smr_alloc(N = 150, arms = c("A", "A", "B", "B", "C"), seed = 3)
#' smr_balance(x)
#'
#' @export
smr_balance <- function(x) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    x$balance
}

#' Block structure of an allocation sequence
#'
#' Returns one row per block with its size, the first and last allocation ids
#' and the arm counts. This is the statistician's copy of the block structure
#' (kept even when \code{keep.blocks = FALSE} was used at generation time).
#'
#' @param x An object of class \code{"sandwichrand"}.
#'
#' @return A data frame, or \code{NULL} if no block structure is stored.
#'
#' @examples
#' x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
#' smr_blocks(x)
#'
#' @export
smr_blocks <- function(x) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    x$blocks
}

#' Running imbalance of an allocation sequence
#'
#' For each allocation position, computes the spread of the cumulative
#' deviations of the arms from their targets. For a 1:1 allocation this is the
#' running difference between the two arms; it is the quantity the SMR scheme
#' is designed to keep small and that starts (by construction) at the head
#' block's deliberate imbalance.
#'
#' @param x An object of class \code{"sandwichrand"}.
#'
#' @return A numeric vector of length \code{N}, the running imbalance.
#'
#' @examples
#' x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
#' tail(smr_imbalance(x))
#'
#' @export
smr_imbalance <- function(x) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    a <- x$alloc
    d <- x$design
    N <- d$N
    mult <- d$mult
    rs <- d$ratio_sum
    ng <- length(mult)
    mi <- match(a$arm, d$arms)
    ## one-hot matrix of the allocations, then cumulative counts per arm:
    ## apply(..., 2, cumsum) yields an N x ng matrix whose i-th row holds the
    ## cumulative counts of each arm up to (and including) position i.
    M <- matrix(0L, nrow = N, ncol = ng)
    M[cbind(seq_len(N), mi)] <- 1L
    cum <- apply(M, 2L, cumsum)
    out <- numeric(N)
    for (i in seq_len(N)) {
        out[i] <- .smr_spread(cum[i, ], i, mult, rs)
    }
    out
}

#' The spread of deviations from the target allocation
#'
#' \code{spread} is the difference between the largest and the smallest
#' deviation of the arm counts from the target allocation
#' \code{count_g - n * mult_g / ratio_sum}. For a 1:1 allocation it is exactly
#' the difference of the two arm sizes (the imbalance used in the SMR paper);
#' for unequal ratios it is the natural generalisation (\code{spread = 0}
#' means exactly proportional).
#'
#' @param cnt Numeric vector of arm counts, in the declared arm order.
#' @param n Numeric. Number of slots the counts are expressed for.
#' @param mult Integer vector. Allocation-ratio units of each arm.
#' @param ratio_sum Numeric. Sum of \code{mult}.
#'
#' @return A single numeric value, the spread.
#'
#' @examples
#' ## 1:1, 10 slots, counts 7 / 3 -> spread 4
#' smr_spread(c(7, 3), 10, c(1, 1), 2)
#' ## 2:2:1, 10 slots, counts 5 / 4 / 1 -> spread 2
#' smr_spread(c(5, 4, 1), 10, c(2, 2, 1), 5)
#'
#' @export
smr_spread <- function(cnt, n, mult, ratio_sum) {
    cnt <- as.numeric(cnt)
    dev <- cnt - n * mult / ratio_sum
    max(dev) - min(dev)
}

#' Assess an SMR design by simulation
#'
#' Runs \code{reps} independent replications of \code{\link{smr_alloc}} with
#' independent random streams and summarises the two key performance metrics of
#' a randomization procedure: the final and the maximum running imbalance
#' (allocation balance) and the proportion of correct guesses, PCG (allocation
#' unpredictability). This is the routine to use when a design (sample size,
#' arms, imbalance window) needs to be assessed before a trial is started.
#'
#' @param N Integer. Total sample size.
#' @param arms Treatment-arm specification, see \code{\link{smr_alloc}}.
#' @param reps Integer. Number of replications.
#' @param seed Integer. Base seed; the \code{reps} independent streams are
#'   derived from it, so the assessment is fully reproducible.
#' @param ... Further arguments passed to \code{\link{smr_alloc}}.
#'
#' @return An object of class \code{"smr_simulation"}, a list with
#'   \code{reps}, \code{final_imbalance}, \code{max_imbalance}, \code{pcg}
#'   (the proportion of correct guesses in each replication), the mean running
#'   curves \code{imbalance_curve} and \code{pcg_curve} (one value per
#'   allocation step), and \code{n1}.
#'
#' @examples
#' s <- smr_simulate(N = 48, arms = c("A", "B"), reps = 200, seed = 11)
#' s
#'
#' @export
smr_simulate <- function(N, arms, reps = 100L, seed = NULL, ...) {
    if (is.null(seed)) seed <- as.integer(as.double(Sys.time()) %% 100000000)
    reps <- as.integer(reps)
    streams <- .smr_streams(as.integer(round(seed)), reps)
    final_imb <- numeric(reps)
    max_imb <- numeric(reps)
    pcg <- numeric(reps)
    n1_used <- integer(reps)
    imb_sum <- NULL
    pcg_sum <- NULL
    n_ok <- 0L
    for (i in seq_len(reps)) {
        x <- try(smr_alloc(N = N, arms = arms, stream = streams[[i]], ...),
                 silent = TRUE)
        if (inherits(x, "try-error")) {
            final_imb[i] <- NA_real_
            max_imb[i] <- NA_real_
            pcg[i] <- NA_real_
            n1_used[i] <- NA_integer_
        } else {
            imb <- smr_imbalance(x)
            cg <- smr_correct_guesses(x)
            final_imb[i] <- imb[length(imb)]
            max_imb[i] <- max(imb)
            pcg[i] <- mean(cg)
            n1_used[i] <- x$design$n1
            if (is.null(imb_sum)) {
                imb_sum <- imb
                pcg_sum <- cg
            } else {
                imb_sum <- imb_sum + imb
                pcg_sum <- pcg_sum + cg
            }
            n_ok <- n_ok + 1L
        }
    }
    structure(list(reps = reps,
                   final_imbalance = final_imb,
                   max_imbalance = max_imb,
                   pcg = pcg,
                   imbalance_curve = if (n_ok > 0L) imb_sum / n_ok else NA_real_,
                   pcg_curve = if (n_ok > 0L) pcg_sum / n_ok else NA_real_,
                   n1 = n1_used,
                   seed = seed,
                   call = match.call()),
              class = "smr_simulation")
}

#' @export
print.smr_simulation <- function(x, ...) {
    cat("SMR design assessment by simulation\n")
    cat(sprintf("  replications        = %d\n", x$reps))
    f <- x$final_imbalance
    m <- x$max_imbalance
    if (all(is.na(f))) {
        cat("  every replication failed\n")
        return(invisible(x))
    }
    cat(sprintf("  final imbalance     : mean %s, median %s, max %s, P(|imb| <= 1) %s\n",
                format(mean(f, na.rm = TRUE)),
                format(stats::median(f, na.rm = TRUE)),
                format(max(f, na.rm = TRUE)),
                format(mean(f <= 1, na.rm = TRUE))))
    cat(sprintf("  maximum imbalance   : mean %s, median %s, 95%% quantile %s\n",
                format(mean(m, na.rm = TRUE)),
                format(stats::median(m, na.rm = TRUE)),
                format(quantile(m, 0.95, na.rm = TRUE))))
    if (!is.null(x$pcg) && !all(is.na(x$pcg))) {
        cat(sprintf("  correct-guess PCG   : mean %s, median %s, min %s, max %s\n",
                    format(mean(x$pcg, na.rm = TRUE)),
                    format(stats::median(x$pcg, na.rm = TRUE)),
                    format(min(x$pcg, na.rm = TRUE)),
                    format(max(x$pcg, na.rm = TRUE))))
    }
    if (!all(is.na(x$n1))) {
        cat(sprintf("  first block size n1 : %s\n",
                    paste(range(x$n1, na.rm = TRUE), collapse = " .. ")))
    }
    invisible(x)
}
