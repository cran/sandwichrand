## ---------------------------------------------------------------------------
## Evaluation of allocation balance and unpredictability
## ---------------------------------------------------------------------------

test_that("correct guesses match the 1:1 Blackwell-Hodges definition", {
    fake <- structure(
        list(alloc = data.frame(id = 1:4, arm = c("A", "B", "A", "B"),
                                stringsAsFactors = FALSE),
             design = list(N = 4L, arms = c("A", "B"), mult = c(1L, 1L),
                           ratio_sum = 2L, n_arms = 2L)),
        class = "sandwichrand")
    ## i=1 tie -> 0.5; i=2 B has fewer -> guess B (correct);
    ## i=3 tie -> 0.5; i=4 B has fewer -> guess B (correct)
    expect_equal(smr_correct_guesses(fake), c(0.5, 1, 0.5, 1))
    expect_equal(smr_pcg(fake), 0.75)

    fake2 <- structure(
        list(alloc = data.frame(id = 1:4, arm = c("A", "A", "B", "B"),
                                stringsAsFactors = FALSE),
             design = list(N = 4L, arms = c("A", "B"), mult = c(1L, 1L),
                           ratio_sum = 2L, n_arms = 2L)),
        class = "sandwichrand")
    expect_equal(smr_correct_guesses(fake2), c(0.5, 0, 1, 1))
    expect_equal(smr_pcg(fake2), 0.625)
})

test_that("multi-arm tie handling follows the convergence strategy", {
    fake3 <- structure(
        list(alloc = data.frame(id = 1:3, arm = c("A", "B", "C"),
                                stringsAsFactors = FALSE),
             design = list(N = 3L, arms = c("A", "B", "C"),
                           mult = c(1L, 1L, 1L),
                           ratio_sum = 3L, n_arms = 3L)),
        class = "sandwichrand")
    ## i=1 tie among 3 -> 1/3; i=2 B/C tied for fewest -> 1/2; i=3 C -> 1
    expect_equal(smr_correct_guesses(fake3), c(1 / 3, 1 / 2, 1))
})

test_that("correct guesses on real sequences stay in [0,1] and start at 1/n_arms", {
    x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 1)
    cg <- smr_correct_guesses(x)
    expect_length(cg, 48)
    expect_true(all(cg >= 0 & cg <= 1))
    expect_equal(cg[1L], 0.5)
    expect_equal(smr_pcg(x), mean(cg))

    y <- smr_alloc(N = 150, arms = c("A", "A", "B", "B", "C"), seed = 3)
    cgy <- smr_correct_guesses(y)
    expect_length(cgy, 150)
    expect_true(all(cgy >= 0 & cgy <= 1))
    expect_equal(cgy[1L], 1 / 3)
})

test_that("smr_simulate reports PCG and the running curves", {
    s <- smr_simulate(N = 48, arms = c("A", "B"), reps = 20, seed = 11)
    expect_true(all(c("pcg", "imbalance_curve", "pcg_curve") %in% names(s)))
    expect_length(s$pcg, 20)
    expect_length(s$imbalance_curve, 48)
    expect_length(s$pcg_curve, 48)
    expect_true(all(s$pcg >= 0 & s$pcg <= 1))
    expect_true(all(s$pcg_curve >= 0 & s$pcg_curve <= 1))
    ## both are averages over the same set of correct-guess indicators
    expect_equal(mean(s$pcg), mean(s$pcg_curve), tolerance = 1e-12)
})
