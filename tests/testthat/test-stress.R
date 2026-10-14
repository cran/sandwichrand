## ---------------------------------------------------------------------------
## Stress tests: many seeds, many arms, many designs, always valid
## ---------------------------------------------------------------------------

test_that("every seed in a sweep produces a valid sequence", {
    for (sd in 1:40) {
        x <- smr_alloc(N = 48, arms = c("A", "B"), seed = sd)
        expect_smr_ok(x, 48L, c("A", "B"))
    }
})

test_that("many designs and ratios are valid", {
    cases <- list(
        list(N = 48L,  arms = c("A", "B")),
        list(N = 150L, arms = c("A", "A", "B", "B", "C")),
        list(N = 40L,  arms = c("A", "A", "A", "B")),
        list(N = 99L,  arms = c("X", "Y", "Z")),
        list(N = 50L,  arms = c("P", "Q", "R", "S", "T")),
        list(N = 20L,  arms = c("P", "Q", "R", "S", "T"))
    )
    for (cs in cases) {
        for (sd in c(1L, 99L)) {
            x <- smr_alloc(N = cs$N, arms = cs$arms, seed = sd)
            expect_smr_ok(x, cs$N, cs$arms)
        }
    }
})

test_that("simulation summarises a design", {
    s <- smr_simulate(N = 48, arms = c("A", "B"), reps = 30, seed = 11)
    expect_s3_class(s, "smr_simulation")
    expect_equal(s$reps, 30L)
    expect_true(mean(s$final_imbalance, na.rm = TRUE) <= 1 + 1e-6)
    expect_output(print(s), "replications")
})

test_that("N = 4*ratio_sum - 1 goes to path 1 and stays valid", {
    ## 1:1 -> 4*2 = 8; N = 7 < 8 -> path 1
    x <- smr_alloc(N = 7, arms = c("A", "B"), seed = 3)
    expect_equal(x$design$path, 1L)
    expect_smr_ok(x, 7L, c("A", "B"))
})

test_that("N = 4*ratio_sum goes to path 2/3 and stays valid", {
    ## 1:1 -> N = 8 -> path 2 (n1=4 -> middle 0 -> path 3 actually)
    x <- smr_alloc(N = 8, arms = c("A", "B"), seed = 3)
    expect_true(x$design$path %in% c(2L, 3L))
    expect_smr_ok(x, 8L, c("A", "B"))
})
