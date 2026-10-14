## ---------------------------------------------------------------------------
## Argument checking and normalisation.
##
## Every error message names the offending argument and the admissible range,
## so that a mis-specified design fails loudly instead of silently producing a
## wrong randomization list. The result of .smr_spec() is the single source of
## truth used by the three generation paths.
## ---------------------------------------------------------------------------

## Parse the arm specification.
##
##   c("TestTrt", "Placebo")        -> 1:1,   arms TestTrt / Placebo
##   c("A", "A", "B", "B", "C")     -> 2:2:1, arms A / B / C
##   "A A B B C"                    -> 2:2:1  (single blank-separated string)
##   c(TestTrt = 1, Placebo = 1)    -> 1:1    (named allocation-ratio vector)
##   c(TestTrt = 2, Placebo = 2, C = 1) -> 2:2:1
##
## One token per allocation-ratio unit; a repeated name raises the
## multiplicity of that arm. Arm identity is case-insensitive, as in the SAS
## macro (which compares with UPCASE).
.smr_parse_arms <- function(arms) {
    if (is.null(arms) || length(arms) == 0L) {
        stop("`arms` is required, e.g. c(\"TestTrt\", \"Placebo\") for 1:1 ",
             "or c(\"A\", \"A\", \"B\", \"B\", \"C\") for 2:2:1.",
             call. = FALSE)
    }

    if (is.character(arms)) {
        if (length(arms) == 1L && grepl("\\s", arms)) {
            tok <- unlist(strsplit(trimws(arms), "[[:space:]]+"))
        } else {
            tok <- as.character(arms)
        }
        tok <- trimws(tok)
        tok <- tok[nzchar(tok)]
        if (!length(tok)) {
            stop("`arms` must contain at least one non-empty arm name.",
                 call. = FALSE)
        }
        key <- toupper(tok)
        ukey <- unique(key)
        tmap <- match(key, ukey)
        mult <- tabulate(tmap, nbins = length(ukey))
        labels <- tok[match(ukey, key)]
        ratio_sum <- length(tok)
    } else {
        v <- suppressWarnings(as.numeric(arms))
        if (length(v) < 1L || any(!is.finite(v)) || any(v <= 0) ||
            !.smr_is_whole_set(v)) {
            stop("`arms` must be a character vector of arm names or a vector ",
                 "of positive whole allocation-ratio units (got ",
                 paste(format(arms), collapse = " "), ").", call. = FALSE)
        }
        mult <- as.integer(round(v))
        nm <- names(arms)
        if (is.null(nm)) {
            nm <- if (length(mult) <= 26L) LETTERS[seq_along(mult)]
                  else paste0("Arm", seq_along(mult))
        }
        nm <- trimws(as.character(nm))
        if (any(!nzchar(nm)) || anyDuplicated(toupper(nm))) {
            stop("the names of `arms` must be non-empty and unique (got ",
                 paste(nm, collapse = ", "), ").", call. = FALSE)
        }
        labels <- nm
        tmap <- rep.int(seq_along(mult), mult)
        ratio_sum <- sum(mult)
    }

    if (ratio_sum < 2L) {
        stop("`arms` must contain at least 2 allocation-ratio units (received ",
             paste(as.character(arms), collapse = " "), ").", call. = FALSE)
    }
    if (length(labels) < 2L) {
        stop("`arms` must define at least 2 DISTINCT arms (received ",
             paste(as.character(arms), collapse = " "), ").", call. = FALSE)
    }
    list(tmap = as.integer(tmap),
         mult = as.integer(mult),
         labels = labels,
         ratio_sum = as.integer(ratio_sum),
         n_arms = length(labels))
}

## Validate and normalise the design specification.
.smr_spec <- function(N, arms, seed = NULL, stream = NULL,
                      block.sizes = NULL, imbalance = c(1, 4),
                      max.cand = 10L, max.iter = 1000L) {

    notes <- character(0)
    warnings <- character(0)

    ## -- N : total sample size --------------------------------------------
    if (missing(N) || is.null(N)) {
        stop("`N` is required (total sample size).", call. = FALSE)
    }
    n_total <- suppressWarnings(as.numeric(N))
    if (length(n_total) != 1L || !is.finite(n_total)) {
        stop("`N` must be a single number (received ",
             paste(format(N), collapse = " "), ").", call. = FALSE)
    }
    if (n_total <= 0) {
        stop("`N` must be > 0 (received ", format(N), ").", call. = FALSE)
    }
    if (!.smr_is_whole(n_total)) {
        stop("`N` must be a whole number (received ", format(N), ").",
             call. = FALSE)
    }
    n_total <- as.integer(round(n_total))

    ## -- arms : allocation ratio ------------------------------------------
    g <- .smr_parse_arms(arms)

    ## -- imbalance window --------------------------------------------------
    if (is.null(imbalance) || length(imbalance) == 0L) imbalance <- c(1, 4)
    if (length(imbalance) == 1L) imbalance <- c(0, imbalance)
    imb <- suppressWarnings(as.numeric(imbalance))
    if (length(imb) != 2L || any(!is.finite(imb))) {
        stop("`imbalance` must be two numbers c(min, max) (received ",
             paste(format(imbalance), collapse = " "), ").", call. = FALSE)
    }
    ndiffmin <- imb[1L]
    ndiffmax <- imb[2L]
    if (ndiffmin < 0 || ndiffmin > ndiffmax) {
        stop("require 0 <= imbalance[1] <= imbalance[2] (got ",
             paste(format(imbalance), collapse = ", "), ").", call. = FALSE)
    }

    ## -- max.cand / max.iter ----------------------------------------------
    if (is.null(max.cand) || length(max.cand) != 1L ||
        !is.finite(as.numeric(max.cand)) || as.numeric(max.cand) < 1) {
        max.cand <- 10L
    }
    if (is.null(max.iter) || length(max.iter) != 1L ||
        !is.finite(as.numeric(max.iter)) || as.numeric(max.iter) < 1) {
        max.iter <- 1000L
    }
    max.cand <- as.integer(max.cand)
    max.iter <- as.integer(max.iter)

    ## -- seed / stream ------------------------------------------------------
    ## At most one of `seed` / `stream` should be supplied. A `stream` is the
    ## raw RNG state of a L'Ecuyer-CMRG sub-stream (e.g. from smr_streams());
    ## it is the most efficient way to give one call an independent stream.
    if (!is.null(stream) && !is.null(seed)) {
        stop("supply at most one of `seed` and `stream`.", call. = FALSE)
    }
    if (!is.null(stream) &&
        (is.null(stream) || !is.numeric(stream) ||
            length(stream) == 0L)) {
        stop("`stream` must be a numeric vector (an L'Ecuyer-CMRG state).",
             call. = FALSE)
    }
    if (is.null(stream)) {
        if (is.null(seed) || (is.character(seed) && !nzchar(seed))) {
            seed <- as.integer(as.double(Sys.time()) %% 100000000)
            seed_note <- "The seed was generated automatically from the system clock."
        } else {
            s <- suppressWarnings(as.numeric(seed))
            if (length(s) != 1L || !is.finite(s)) {
                stop("`seed` must be a single number (received ",
                     paste(format(seed), collapse = " "), ").", call. = FALSE)
            }
            seed <- as.integer(round(s))
            seed_note <- "The seed was specified by the user."
        }
    } else {
        seed <- NA_integer_
        seed_note <- "An independent stream was supplied by the caller."
    }

    ## -- middle block sizes -------------------------------------------------
    if (is.null(block.sizes) || length(block.sizes) == 0L ||
        (is.character(block.sizes) && !any(nzchar(block.sizes)))) {
        block.sizes <- c(2L, 3L, 4L) * g$ratio_sum
    } else {
        if (is.character(block.sizes)) {
            block.sizes <- suppressWarnings(as.numeric(
                unlist(strsplit(trimws(paste(block.sizes, collapse = " ")),
                                "[[:space:]]+"))))
        }
        block.sizes <- suppressWarnings(as.numeric(block.sizes))
        if (length(block.sizes) == 0L || any(!is.finite(block.sizes))) {
            stop("`block.sizes` must be a vector of positive whole numbers ",
                 "(got ", paste(format(block.sizes), collapse = " "), ").",
                 call. = FALSE)
        }
        if (!.smr_is_whole_set(block.sizes) || any(block.sizes < 1)) {
            stop("`block.sizes` must be positive whole numbers (got ",
                 paste(format(block.sizes), collapse = " "), ").",
                 call. = FALSE)
        }
        block.sizes <- as.integer(round(block.sizes))
    }
    small <- block.sizes[block.sizes < g$ratio_sum]
    if (length(small)) {
        notes <- c(notes,
                   sprintf(paste("block size %s is below ratio_sum (%d), so such",
                                 "a block cannot contain every treatment arm."),
                           paste(small, collapse = ", "), g$ratio_sum))
    }
    odd <- block.sizes[block.sizes %% g$ratio_sum != 0L]
    if (length(odd)) {
        notes <- c(notes,
                   sprintf(paste("block size %s is not a multiple of ratio_sum",
                                 "(%d), so such a block cannot be exactly",
                                 "proportional; the residual imbalance is kept",
                                 "as small as possible but adds up over the",
                                 "middle blocks."),
                         paste(odd, collapse = ", "), g$ratio_sum))
    }

    ## -- n1 : the head-block size is always sampled at random --------------
    ## (no user override; see .smr_generate() in core.R for the sampling.)

    if (n_total %% g$ratio_sum != 0L) {
        warnings <- c(warnings,
                      sprintf(paste("N (%d) is not a multiple of ratio_sum (%d),",
                                    "so exact proportional balance is",
                                    "impossible; the deviation stays <= 1 per",
                                    "arm."), n_total, g$ratio_sum))
    }

    list(N = n_total,
         arms = g$labels,
         arms_input = arms,
         tmap = g$tmap,
         mult = g$mult,
         ratio_sum = g$ratio_sum,
         n_arms = g$n_arms,
         imbalance = c(ndiffmin, ndiffmax),
         ndiffmin = ndiffmin,
         ndiffmax = ndiffmax,
         max.cand = max.cand,
         max.iter = max.iter,
         block.sizes = block.sizes,
         vmin = min(block.sizes),
         vmax = max(block.sizes),
         seed = seed,
         stream = stream,
         seed_note = seed_note,
         notes = notes,
         warnings = warnings)
}
