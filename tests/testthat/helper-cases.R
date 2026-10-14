## ---------------------------------------------------------------------------
## Shared helpers and the design matrix used by several test files
## ---------------------------------------------------------------------------

## Number of allocation-ratio units (sum of the ratio) implied by `arms`.
ratio_sum_of <- function(arms) {
    if (is.character(arms)) {
        if (length(arms) == 1L && grepl("\\s", arms)) {
            length(unlist(strsplit(trimws(arms), "[[:space:]]+")))
        } else {
            length(arms)
        }
    } else {
        sum(as.integer(arms))
    }
}

## Structural invariants that must hold for EVERY generated sequence.
expect_smr_ok <- function(x, N, arms, keep.blocks = TRUE) {
    rs <- ratio_sum_of(arms)
    d <- x$design

    expect_true(inherits(x, "sandwichrand"), label = "class")
    expect_identical(d$N, N, label = "N recorded")

    a <- x$alloc
    expect_equal(nrow(a), N, label = "rows")
    expect_true(all(!is.na(a$arm)), label = "no missing arm")
    expect_identical(as.integer(a$id), seq_len(N), label = "ids 1..N")
    expect_true(all(a$arm %in% d$arms), label = "arms declared")

    ## balance
    expect_lte(max(abs(x$balance$deviation)), 1 + 1e-6, label = "deviation")
    expect_equal(sum(x$balance$n), N, label = "counts add up to N")
    expect_equal(x$balance$expected, N * d$mult / rs, label = "expected")

    ## columns
    if (isTRUE(keep.blocks)) {
        expect_named(a, c("id", "block", "block_size", "arm"),
                     label = "columns with blocks")
    } else {
        expect_named(a, c("id", "arm"), label = "columns released")
    }

    ## block structure (always internally available)
    bt <- x$blocks
    expect_true(!is.null(bt), label = "blocks stored")
    expect_equal(sum(bt$size), N, label = "block sizes add up to N")
    expect_true(all(bt$first_id[-1] == bt$last_id[-nrow(bt)] + 1L),
                label = "blocks contiguous")
    expect_equal(bt$first_id[1L], 1L, label = "first id")
    expect_equal(bt$last_id[nrow(bt)], N, label = "last id")

    if (d$path == 1L) {
        expect_equal(nrow(bt), 1L, label = "path 1: one block")
    } else {
        expect_equal(bt$size[1L], bt$size[nrow(bt)], label = "n1 equal")
        expect_gte(bt$size[1L], 2L * rs, label = "n1 >= 2*ratio_sum")
        ## head + tail complementarity
        cols <- paste0("n_", d$arms)
        cnt <- unlist(bt[1L, cols]) + unlist(bt[nrow(bt), cols])
        expect_equal(max(abs(cnt * rs - d$mult * (bt$size[1L] +
                                                      bt$size[nrow(bt)]))),
                     0, label = "head+tail complementary")
    }

    ## running imbalance
    imb <- smr_imbalance(x)
    expect_length(imb, N)
    expect_true(all(imb >= 0), label = "imbalance non-negative")
    expect_lte(imb[N], 1 + 1e-6, label = "final imbalance <= 1")

    ## verification
    v <- smr_verify(x)
    expect_true(all(v$status == "PASS"),
                label = paste0("verify: ",
                               paste(v$check[v$status == "FAIL"],
                                     collapse = ", ")))
    invisible(TRUE)
}

## The design matrix: two 1:1 arms, three 2:2:1 arms, 3:1, three-arm 1:1:1
## and five-arm 1:1:1:1:1, over a range of sample sizes.
design_matrix <- function() {
    out <- list()
    for (N in c(6L, 8L, 9L, 11L, 12L, 13L, 24L, 36L, 37L, 48L, 49L, 100L)) {
        out[[length(out) + 1L]] <- list(N = N, arms = c("TestTrt", "Placebo"))
    }
    for (N in c(10L, 20L, 25L, 60L, 150L, 151L)) {
        out[[length(out) + 1L]] <- list(N = N, arms = c("A", "A", "B", "B", "C"))
    }
    for (N in c(16L, 40L, 41L)) {
        out[[length(out) + 1L]] <- list(N = N, arms = c("A", "A", "A", "B"))
    }
    for (N in c(12L, 30L, 99L)) {
        out[[length(out) + 1L]] <- list(N = N, arms = c("X", "Y", "Z"))
    }
    for (N in c(20L, 50L)) {
        out[[length(out) + 1L]] <- list(N = N, arms = c("P", "Q", "R", "S", "T"))
    }
    out
}
