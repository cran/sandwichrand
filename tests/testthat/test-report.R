## ---------------------------------------------------------------------------
## RTF and CSV release
## ---------------------------------------------------------------------------

test_that("smr_write_csv writes a released list without block columns", {
    x <- smr_alloc(N = 12, arms = c("Drug", "Placebo"), seed = 777,
                   keep.blocks = FALSE)
    f <- tempfile(fileext = ".csv")
    invisible(smr_write_csv(x, f))
    d <- read.csv(f, stringsAsFactors = FALSE)
    expect_named(d, c("id", "arm"))
    expect_equal(nrow(d), 12)
    file.remove(f)
})

test_that("smr_write_csv writes block columns when keep.blocks = TRUE", {
    x <- smr_alloc(N = 12, arms = c("Drug", "Placebo"), seed = 777,
                   keep.blocks = TRUE)
    f <- tempfile(fileext = ".csv")
    invisible(smr_write_csv(x, f, keep.blocks = TRUE))
    d <- read.csv(f, stringsAsFactors = FALSE)
    expect_true(all(c("block", "block_size") %in% names(d)))
    file.remove(f)
})

test_that("smr_rtf writes a readable RTF document", {
    x <- smr_alloc(N = 12, arms = c("Drug", "Placebo"), seed = 777)
    f <- tempfile(fileext = ".rtf")
    invisible(smr_rtf(x, f))
    txt <- readLines(f, warn = FALSE)
    expect_true(any(grepl("rtf1", txt[1])))
    expect_true(any(grepl("Drug", txt)))
    file.remove(f)
})

test_that("smr_rtf honours keep.blocks", {
    x <- smr_alloc(N = 12, arms = c("Drug", "Placebo"), seed = 777)
    f1 <- tempfile(fileext = ".rtf")
    f2 <- tempfile(fileext = ".rtf")
    smr_rtf(x, f1, keep.blocks = TRUE)
    smr_rtf(x, f2, keep.blocks = FALSE)
    t1 <- paste(readLines(f1, warn = FALSE), collapse = "\n")
    t2 <- paste(readLines(f2, warn = FALSE), collapse = "\n")
    expect_true(grepl("Block size", t1))
    expect_false(grepl("Block size", t2))
    file.remove(f1, f2)
})
