## ---------------------------------------------------------------------------
## print / summary / as.data.frame / plot methods
## ---------------------------------------------------------------------------

test_that("print and format return the object invisibly", {
    x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
    expect_output(print(x), "SMR allocation sequence")
    expect_silent(f <- format(x))
    expect_true(is.character(f))
    expect_true(any(grepl("N", f)))
    expect_invisible(print(x))
})

test_that("print shows warnings when present", {
    x <- smr_alloc(N = 49, arms = c("A", "B"), seed = 1)
    expect_output(print(x), "Warnings")
})

test_that("summary returns a summary.sandwichrand object", {
    x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
    s <- summary(x)
    expect_s3_class(s, "summary.sandwichrand")
    expect_output(print(s), "Verification")
})

test_that("as.data.frame extracts the allocation", {
    x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
    d <- as.data.frame(x)
    expect_equal(nrow(d), 48)
    expect_true("arm" %in% names(d))
})

test_that("plot methods run without error", {
    x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
    f <- tempfile(fileext = ".pdf")
    pdf(f)
    expect_silent(plot(x, type = "imbalance"))
    expect_silent(plot(x, type = "cumulative"))
    dev.off()
    expect_true(file.exists(f))
    file.remove(f)
})

test_that("methods reject non-sandwichrand input", {
    expect_error(print.sandwichrand(1), "sandwichrand")
    expect_error(summary.sandwichrand(1), "sandwichrand")
    expect_error(plot.sandwichrand(1), "sandwichrand")
    expect_error(as.data.frame.sandwichrand(1))
})
