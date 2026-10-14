## ---------------------------------------------------------------------------
## Reproducibility, seeds and independent streams
## ---------------------------------------------------------------------------

test_that("the same seed reproduces the same sequence", {
    a <- smr_alloc(N = 48, arms = c("A", "B"), seed = 2026)
    b <- smr_alloc(N = 48, arms = c("A", "B"), seed = 2026)
    expect_identical(a$alloc, b$alloc)
    expect_identical(a$blocks, b$blocks)
    expect_equal(a$design$n1, b$design$n1)
})

test_that("different seeds give different sequences", {
    a <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
    b <- smr_alloc(N = 48, arms = c("A", "B"), seed = 2)
    expect_false(identical(a$alloc, b$alloc))
})

test_that("smr_reproduce regenerates an identical object", {
    x <- smr_alloc(N = 150, arms = c("A", "A", "B", "B", "C"), seed = 123)
    y <- smr_reproduce(x)
    expect_identical(x$alloc, y$alloc)
    expect_identical(x$blocks, y$blocks)
})

test_that("smr_streams gives independent, reproducible streams", {
    ss <- smr_streams(seed = 42, n = 5)
    expect_length(ss, 5L)
    ## same seed -> same streams
    expect_identical(ss, smr_streams(seed = 42, n = 5))
    ## a given stream is reproducible
    a1 <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[3]])
    a2 <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[3]])
    expect_identical(a1$alloc, a2$alloc)
    ## different streams give different sequences
    b1 <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[1]])
    b2 <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[2]])
    expect_false(identical(b1$alloc, b2$alloc))
})

test_that("high-frequency calling never re-uses the same stream", {
    ## the key property for a simulation loop: 200 calls, all distinct
    ss <- smr_streams(seed = 7, n = 200)
    allocs <- lapply(ss, function(s)
        smr_alloc(N = 48, arms = c("A", "B"), stream = s)$alloc)
    ids <- vapply(allocs, function(a) paste(a$arm, collapse = ""), character(1))
    expect_equal(length(unique(ids)), 200L)
})

test_that("calls never disturb the session RNG state", {
    set.seed(123)
    before <- .Random.seed
    invisible(smr_alloc(N = 48, arms = c("A", "B"), seed = 999))
    after <- .Random.seed
    expect_identical(before, after)
})

test_that("a NULL seed draws a seed from the system clock", {
    x <- smr_alloc(N = 48, arms = c("A", "B"))
    expect_true(is.finite(x$design$seed))
    expect_match(x$design$seed_note, "system clock")
})

test_that("seed and stream are mutually exclusive", {
    expect_error(smr_alloc(N = 48, arms = c("A", "B"), seed = 1,
                           stream = smr_streams(1, 1)[[1]]),
                 "at most one of `seed` and `stream`")
})
