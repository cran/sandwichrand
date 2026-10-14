## ---------------------------------------------------------------------------
## Argument validation: mis-specified designs must fail loudly
## ---------------------------------------------------------------------------

test_that("missing or invalid N is rejected", {
    expect_error(smr_alloc(arms = c("A", "B")), "N.*required")
    expect_error(smr_alloc(N = 0, arms = c("A", "B")), "> 0")
    expect_error(smr_alloc(N = -5, arms = c("A", "B")), "> 0")
    expect_error(smr_alloc(N = 10.5, arms = c("A", "B")), "whole number")
    expect_error(smr_alloc(N = "x", arms = c("A", "B")), "single number")
})

test_that("missing or invalid arms is rejected", {
    expect_error(smr_alloc(N = 48), "arms.*required")
    expect_error(smr_alloc(N = 48, arms = c("A")), "at least 2")
    expect_error(smr_alloc(N = 48, arms = c("a", "A")), "at least 2 DISTINCT")
    expect_error(smr_alloc(N = 48, arms = c(1, -1)), "positive whole")
})

test_that("invalid imbalance is rejected", {
    expect_error(smr_alloc(N = 48, arms = c("A", "B"), imbalance = c(4, 1)),
                 "imbalance\\[1\\] <= imbalance\\[2\\]")
    expect_error(smr_alloc(N = 48, arms = c("A", "B"), imbalance = c(-1, 2)),
                 "imbalance\\[1\\] <= imbalance\\[2\\]")
})

test_that("invalid block.sizes is rejected", {
    expect_error(smr_alloc(N = 48, arms = c("A", "B"), block.sizes = c(0, 4)),
                 "positive whole")
    expect_error(smr_alloc(N = 48, arms = c("A", "B"), block.sizes = c(-2)),
                 "positive whole")
})

test_that("invalid seed is rejected", {
    expect_error(smr_alloc(N = 48, arms = c("A", "B"), seed = "not-a-number"),
                 "single number")
})

test_that("a custom imbalance window is honoured", {
    x <- smr_alloc(N = 96, arms = c("A", "B"), imbalance = c(2, 3), seed = 5)
    expect_true(x$design$imbalance >= 2)
    expect_true(x$design$imbalance <= 3)
})

test_that("block sizes below ratio_sum produce a note", {
    x <- smr_alloc(N = 96, arms = c("A", "B"), block.sizes = c(2, 3, 4),
                   seed = 1)
    expect_true(any(grepl("below ratio_sum", x$notes)) ||
                    any(grepl("not a multiple", x$notes)))
})
