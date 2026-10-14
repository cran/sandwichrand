## ---------------------------------------------------------------------------
## The three generation paths of the sandwich mixed randomization scheme.
##
##   Path 1 : N < 4 * ratio_sum
##            Complete randomization with a pre-determined N, one single block.
##            The schedule is short enough that no block structure is needed.
##
##   Path 2 : Module I  -> head block of size n1, deliberately imbalanced
##            Module II -> middle section of permuted blocks with varying sizes
##            head + tail are two halves of ONE balanced pool of 2*n1 slots,
##            so the tail is the exact complement of the head.
##
##   Path 3 : N == 2 * n1 -> head + tail only, no middle section.
##
## Each module draws from its own L'Ecuyer-CMRG sub-stream, so the modules are
## independent and the whole sequence is a deterministic function of the stream
## (or, equivalently, of the seed).
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Candidate values of n1.
##
##   Multiples of ratio_sum in [2 * ratio_sum, N / 2], at most max.cand of
##   them, preferring the values that leave a middle section of 0 or of at
##   least the smallest block size -- so that the sandwich never needs a
##   degenerate middle block unless that is unavoidable.
## ---------------------------------------------------------------------------
.smr_n1_candidates <- function(spec) {
    kmax <- floor(spec$N / (2 * spec$ratio_sum))
    if (kmax < 2L) return(integer(0))
    cand <- (2L:kmax) * spec$ratio_sum
    if (length(cand) > spec$max.cand) cand <- cand[seq_len(spec$max.cand)]
    mid <- spec$N - 2L * cand
    keep <- cand[mid == 0L | mid >= spec$vmin]
    if (length(keep)) cand <- keep
    as.integer(cand)
}

## ---------------------------------------------------------------------------
## Module I : the two complementary head / tail blocks.
##
##   One balanced pool of 2*n1 slots is shuffled; the first n1 slots become the
##   head block and the remaining n1 slots ARE the complement of the head, so
##   the head and the tail together are exactly proportional.  The shuffle is
##   repeated until the head block is imbalanced by an amount inside
##   [ndiffmin, ndiffmax] AND is not exactly proportional; the closest
##   candidate is kept as a fallback so that a sequence is always produced.
## ---------------------------------------------------------------------------
.smr_module1 <- function(spec, n1, state) {
    n2 <- 2L * n1
    rs <- spec$ratio_sum
    mult <- spec$mult
    ng <- spec$n_arms
    ## The balanced pool reproduces the target allocation: cycling through the
    ## token map tmap (one token per allocation-ratio unit) exactly r times
    ## builds 2*n1 slots of which group g appears 2*(n1/r)*mult_g times.
    pool <- rep(spec$tmap, length.out = n2)

    .with_rng(state, {
        best_pen <- Inf
        best_pool <- pool
        best_spr <- NA_real_
        ok <- FALSE
        spr <- NA_real_
        it <- 0L
        for (it in seq_len(spec$max.iter)) {
            pool <- .smr_shuffle(pool)
            cnt <- tabulate(pool[seq_len(n1)], nbins = ng)
            dev <- cnt - n1 * mult / rs
            spr <- max(dev) - min(dev)
            ## The head block must NOT be exactly proportional, otherwise the
            ## scheme degenerates into a plain permuted-block design.
            viol <- as.integer(any(cnt * rs != mult * n1))

            if (spr >= spec$ndiffmin && spr <= spec$ndiffmax && viol == 1L) {
                ok <- TRUE
                break
            }
            pen <- 0
            if (spr < spec$ndiffmin) {
                pen <- spec$ndiffmin - spr
            } else if (spr > spec$ndiffmax) {
                pen <- spr - spec$ndiffmax
            }
            if (viol == 0L) pen <- pen + 1000
            if (pen < best_pen) {
                best_pen <- pen
                best_spr <- spr
                best_pool <- pool
            }
        }
        if (!ok) pool <- best_pool
        list(head = as.integer(pool[seq_len(n1)]),
             tail = as.integer(pool[(n1 + 1L):n2]),
             spread = if (ok) spr else best_spr,
             ok = ok,
             iter = it)
    })
}

## ---------------------------------------------------------------------------
## Path 1 : complete randomization with a pre-determined N, single block.
##
##   The pool holds ceil(N / ratio_sum) complete copies of the allocation
##   pattern; only the first N slots are released. Rejected draws are retried
##   until the spread of the whole (short) block is at most ndiffmin.
## ---------------------------------------------------------------------------
.smr_path1 <- function(spec, state) {
    N <- spec$N
    rs <- spec$ratio_sum
    ng <- spec$n_arms
    npool <- as.integer(ceiling(N / rs) * rs)
    pool <- rep(spec$tmap, length.out = npool)

    .with_rng(state, {
        best <- Inf
        best_pool <- pool
        ok <- FALSE
        spr <- NA_real_
        it <- 0L
        for (it in seq_len(spec$max.iter)) {
            pool <- .smr_shuffle(pool)
            cnt <- tabulate(pool[seq_len(N)], nbins = ng)
            spr <- .smr_spread(cnt, N, spec$mult, rs)
            if (spr <= spec$ndiffmin) {
                ok <- TRUE
                break
            }
            if (spr < best) {
                best <- spr
                best_pool <- pool
            }
        }
        if (!ok) {
            pool <- best_pool
            spr <- best
        }
        list(slot = pool, spread = spr, ok = ok, iter = it)
    })
}

## ---------------------------------------------------------------------------
## Module II : permuted blocks with varying block sizes (middle section).
##
##   At every step only those candidate sizes are drawn from that leave a
##   residue of 0 or of at least the smallest block size, so a "crumb" block
##   can never be created.  When no candidate qualifies the middle section is
##   closed with a single block (larger than max(block.sizes), never smaller
##   than min(block.sizes)).  The content of each block is apportioned by the
##   largest-remainder rule and then shuffled.
## ---------------------------------------------------------------------------
.smr_module2 <- function(spec, n_mid, state) {
    vs <- spec$block.sizes
    vmin <- spec$vmin
    vmax <- spec$vmax
    rs <- spec$ratio_sum
    mult <- spec$mult

    .with_rng(state, {
        blk <- integer(n_mid)
        bs <- integer(n_mid)
        sq <- integer(n_mid)
        gi <- integer(n_mid)
        rem <- n_mid
        pos <- 1L
        nblk <- 0L
        while (rem > 0L) {
            nblk <- nblk + 1L
            if (rem <= vmax) {
                b <- rem                      # finish with a single block
            } else {
                fea <- vs[(rem - vs) == 0L | (rem - vs) >= vmin]
                if (length(fea)) {
                    k <- 1L + floor(runif(1L) * length(fea))
                    if (k > length(fea)) k <- length(fea)
                    b <- fea[k]
                } else {
                    b <- rem
                }
            }
            b <- as.integer(b)
            idx <- seq.int(pos, length.out = b)
            cnt <- .smr_apportion(b, mult, rs)
            pat <- .smr_shuffle(rep.int(seq_along(mult), cnt))
            blk[idx] <- nblk
            bs[idx] <- b
            sq[idx] <- seq_len(b)
            gi[idx] <- pat
            pos <- pos + b
            rem <- rem - b
        }
        list(blk = blk, seq = sq, bs = bs, gidx = gi, nblk = nblk)
    })
}

## ---------------------------------------------------------------------------
## Orchestrator: assemble one sequence from a specification and the four
## independent module sub-streams.
##
## This is a genuine function (not a brace expression) so that `return()` works
## normally; it returns the allocation data frame, the raw design fields and
## the notes/warnings accumulated during generation. The caller runs it under
## the base RNG state and merges the notes/warnings into its own.
## ---------------------------------------------------------------------------
.smr_generate <- function(spec, st_sel, st_m1, st_m2, st_p1) {
    notes <- character(0)
    warns <- character(0)

    ## =======================================================================
    ## Path 1 : small sample, complete randomization in one block
    ## =======================================================================
    if (spec$N < 4L * spec$ratio_sum) {
        p1 <- .smr_path1(spec, st_p1)
        if (!p1$ok) {
            warns <- c(warns,
                       sprintf(paste("Path 1: no sequence with imbalance",
                                     "<= %s within %d attempts, so the",
                                     "closest candidate (imbalance = %s)",
                                     "was accepted."),
                               format(spec$ndiffmin), spec$max.iter,
                               format(p1$spread)))
        }
        notes <- c(notes,
                   sprintf(paste("path = 1 (N < 4*ratio_sum, complete",
                                 "randomization with fixed N; imbalance = %s)"),
                           format(p1$spread)))

        slot <- p1$slot[seq_len(spec$N)]
        alloc <- data.frame(id = seq_len(spec$N),
                            block = rep(1L, spec$N),
                            block_size = rep(spec$N, spec$N),
                            arm = spec$arms[slot],
                            stringsAsFactors = FALSE)
        design <- list(N = spec$N, path = 1L, n1 = 0L, n_mid = 0L,
                       n_mid_blocks = 0L, imbalance = p1$spread,
                       iter = p1$iter)
        return(list(alloc = alloc, design = design,
                    notes = notes, warns = warns))
    }

    ## =======================================================================
    ## Module I : choose n1 at random, then build the two complementary blocks
    ## =======================================================================
    cand <- .smr_n1_candidates(spec)
    if (!length(cand)) {
        stop("no admissible first-block size for N = ", spec$N,
             " and ratio_sum = ", spec$ratio_sum, ".", call. = FALSE)
    }
    n1 <- .with_rng(st_sel, {
        k <- 1L + floor(runif(1L) * length(cand))
        if (k > length(cand)) k <- length(cand)
        cand[k]
    })
    n1 <- as.integer(n1)
    n_mid <- as.integer(spec$N - 2L * n1)

    if (n_mid > 0L && n_mid < spec$vmin) {
        notes <- c(notes,
                   sprintf(paste("the middle section (%d slots) is smaller",
                                 "than the smallest block size (%d), so it",
                                 "is generated as a single permuted block."),
                            n_mid, spec$vmin))
    }

    m1 <- .smr_module1(spec, n1, st_m1)
    if (!m1$ok) {
        warns <- c(warns,
                   sprintf(paste("Module I: no first block with imbalance in",
                                 "[%s,%s] within %d attempts, so the closest",
                                 "candidate (imbalance = %s) was used.",
                                 "Consider a larger n1 or a wider",
                                 "imbalance window."),
                           format(spec$ndiffmin), format(spec$ndiffmax),
                           spec$max.iter, format(m1$spread)))
    }

    ## =======================================================================
    ## Module II : permuted blocks with varying sizes (middle section)
    ## =======================================================================
    if (n_mid > 0L) {
        mid <- .smr_module2(spec, n_mid, st_m2)
        nmidblk <- mid$nblk
        block <- c(rep(1L, n1), mid$blk + 1L, rep(nmidblk + 2L, n1))
        bsize <- c(rep(n1, n1), mid$bs, rep(n1, n1))
        gidx <- c(m1$head, mid$gidx, m1$tail)
        path <- 2L
    } else {
        nmidblk <- 0L
        block <- c(rep(1L, n1), rep(2L, n1))
        bsize <- c(rep(n1, n1), rep(n1, n1))
        gidx <- c(m1$head, m1$tail)
        path <- 3L
    }

    notes <- c(notes,
               if (n_mid == 0L)
                   "path = 3 (N = 2*n1, complementary blocks only)"
               else
                   "path = 2 (n1 + permuted blocks with varying sizes + n1)")
    notes <- c(notes,
               sprintf("first block n1 = %d (imbalance = %s, attempts = %d)",
                       n1, format(m1$spread), m1$iter),
               sprintf("middle blocks = %d (%d slots, candidate sizes %s)",
                       nmidblk, n_mid,
                       paste(spec$block.sizes, collapse = " ")))

    alloc <- data.frame(id = seq_len(spec$N),
                        block = as.integer(block),
                        block_size = as.integer(bsize),
                        arm = spec$arms[gidx],
                        stringsAsFactors = FALSE)
    design <- list(N = spec$N, path = path, n1 = n1, n_mid = n_mid,
                   n_mid_blocks = nmidblk, imbalance = m1$spread,
                   iter = m1$iter)
    list(alloc = alloc, design = design, notes = notes, warns = warns)
}
