## ---------------------------------------------------------------------------
## Structural invariants of the generated sequences
## ---------------------------------------------------------------------------

test_that("the design matrix produces valid sequences", {
    for (cs in design_matrix()) {
        for (sd in c(1L, 2026L)) {
            x <- smr_alloc(N = cs$N, arms = cs$arms, seed = sd)
            expect_smr_ok(x, cs$N, cs$arms)
        }
    }
})

test_that("keep.blocks = FALSE releases a list without block information", {
    x <- smr_alloc(N = 48, arms = c("TestTrt", "Placebo"), seed = 7L,
                   keep.blocks = FALSE)
    expect_named(x$alloc, c("id", "arm"))
    expect_true(!is.null(x$blocks))          # statistician's copy is retained
    expect_false(x$design$keep.blocks)
})

test_that("the three paths are reached", {
    ## path 1 : N < 4 * ratio_sum
    expect_equal(smr_alloc(N = 6, arms = c("A", "B"), seed = 1L)$design$path, 1L)
    expect_equal(smr_alloc(N = 7, arms = c("A", "B"), seed = 1L)$design$path, 1L)
    expect_equal(smr_alloc(N = 10, arms = c("A", "A", "B", "B", "C"),
                           seed = 1L)$design$path, 1L)
    ## path 2 : sandwich with a middle section
    expect_equal(smr_alloc(N = 48, arms = c("A", "B"), seed = 1L)$design$path, 2L)
    ## path 3 : N == 2*n1 -> head + tail only (for 1:1 with N = 8, the only
    ## admissible n1 is 4, so the middle section is empty)
    z <- smr_alloc(N = 8, arms = c("A", "B"), seed = 1L)
    expect_equal(z$design$path, 3L)
    expect_equal(z$design$n_mid, 0L)
    expect_equal(nrow(z$blocks), 2L)
})

test_that("n1 is sampled at random (never pre-specified)", {
    ## repeated calls with different seeds should not always pick the same n1
    n1s <- vapply(1:50, function(s) {
        smr_alloc(N = 48, arms = c("A", "B"), seed = s)$design$n1
    }, integer(1L))
    expect_true(length(unique(n1s)) > 1L)
    expect_true(all(n1s %% 2L == 0L))
    expect_true(all(n1s >= 4L) && all(n1s <= 24L))
    ## there is no n1 argument in the public interface
    expect_false("n1" %in% names(formals(smr_alloc)))
})

test_that("custom block sizes are used for the middle section", {
    x <- smr_alloc(N = 96, arms = c("A", "B"),
                   block.sizes = c(4L, 6L, 8L, 10L), seed = 5L)
    bt <- x$blocks
    mid <- bt[bt$block > 1L & bt$block < nrow(bt), , drop = FALSE]
    expect_true(all(mid$size %in% c(4L, 6L, 8L, 10L)) ||
                    (nrow(mid) == 1L && mid$size[1L] == x$design$n_mid))
    expect_equal(x$design$block.sizes, c(4L, 6L, 8L, 10L))
})

test_that("unequal allocation ratios are respected", {
    x <- smr_alloc(N = 150, arms = c("A", "A", "B", "B", "C"), seed = 3L)
    expect_equal(x$design$ratio, "2:2:1")
    expect_equal(x$design$ratio_sum, 5L)
    expect_equal(x$design$n_arms, 3L)
    expect_equal(x$balance$n, c(60L, 60L, 30L))
    expect_equal(x$balance$expected, c(60, 60, 30))
})

test_that("named ratio vectors are accepted", {
    x <- smr_alloc(N = 60, arms = c(TestTrt = 2, Placebo = 2, Ctrl = 1),
                   seed = 11L)
    expect_equal(x$design$arms, c("TestTrt", "Placebo", "Ctrl"))
    expect_equal(x$design$ratio, "2:2:1")
    expect_equal(x$balance$n, c(24L, 24L, 12L))
})

test_that("a blank-separated arm string is accepted", {
    x <- smr_alloc(N = 50, arms = "A A B B C", seed = 11L)
    expect_equal(x$design$arms, c("A", "B", "C"))
    expect_equal(x$design$ratio_sum, 5L)
})

test_that("arm identity is case-insensitive", {
    x <- smr_alloc(N = 60, arms = c("a", "A", "b"), seed = 1L)
    expect_equal(x$design$n_arms, 2L)
    expect_equal(x$design$arms, c("a", "b"))
    expect_equal(x$design$ratio, "2:1")
    expect_error(smr_alloc(N = 48, arms = c("a", "A")), "at least 2 DISTINCT")
})

test_that("N that is not a multiple of ratio_sum is flagged", {
    x <- smr_alloc(N = 49, arms = c("A", "B"), seed = 1L)
    expect_true(any(grepl("not a multiple of ratio_sum", x$warnings)))
    expect_lte(max(abs(x$balance$deviation)), 1)
})

test_that("large N works", {
    for (N in c(500L, 1000L, 2000L)) {
        x <- smr_alloc(N = N, arms = c("A", "B"), seed = 99L)
        expect_smr_ok(x, N, c("A", "B"))
    }
})

test_that("verbose emits messages and returns the same object", {
    expect_message(smr_alloc(N = 48, arms = c("A", "B"), seed = 1L,
                             verbose = TRUE))
    expect_silent(smr_alloc(N = 48, arms = c("A", "B"), seed = 1L,
                            verbose = FALSE))
})

test_that("the first block is imbalanced, not proportional", {
    for (sd in 1:20) {
        x <- smr_alloc(N = 48, arms = c("A", "B"), seed = sd)
        bt <- x$blocks
        cols <- paste0("n_", x$design$arms)
        cnt <- unlist(bt[1L, cols])
        n1 <- bt$size[1L]
        expect_true(any(cnt * x$design$ratio_sum != x$design$mult * n1),
                    label = paste("seed", sd))
        expect_gte(x$design$imbalance, 1)
        expect_lte(x$design$imbalance, 4)
    }
})
