## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

test_that(".smr_shuffle is a permutation and consumes the stream predictably", {
    set.seed(1L, kind = "L'Ecuyer-CMRG")
    x <- 1:20
    y1 <- sandwichrand:::.smr_shuffle(x)
    set.seed(1L, kind = "L'Ecuyer-CMRG")
    y2 <- sandwichrand:::.smr_shuffle(x)
    expect_identical(y1, y2)
    expect_setequal(y1, x)
    expect_equal(length(y1), 20L)
    for (n in c(2L, 3L, 5L, 17L, 50L)) {
        expect_setequal(sandwichrand:::.smr_shuffle(seq_len(n)), seq_len(n))
    }
    expect_identical(sandwichrand:::.smr_shuffle(integer(0)), integer(0))
    expect_identical(sandwichrand:::.smr_shuffle(5L), 5L)
})

test_that(".smr_apportion distributes exactly b slots by largest remainder", {
    for (b in 1:12) {
        a <- sandwichrand:::.smr_apportion(b, c(1L, 1L), 2L)
        expect_equal(sum(a), b)
        expect_equal(max(abs(a - b / 2)), if (b %% 2L == 0L) 0 else 0.5)
    }
    for (b in c(5L, 10L, 12L, 13L, 20L, 23L)) {
        a <- sandwichrand:::.smr_apportion(b, c(2L, 2L, 1L), 5L)
        expect_equal(sum(a), b)
        expect_equal(a[1L], a[2L])
        expect_lte(max(abs(a - b * c(2, 2, 1) / 5)), 1)
    }
    for (b in 1:30) {
        a <- sandwichrand:::.smr_apportion(b, c(3L, 1L), 4L)
        expect_equal(sum(a), b)
        expect_lte(max(abs(a - b * c(3, 1) / 4)), 1)
    }
})

test_that(".smr_spread is the spread of the deviations from the target", {
    expect_equal(sandwichrand:::.smr_spread(c(5L, 5L), 10L, c(1L, 1L), 2L), 0)
    expect_equal(sandwichrand:::.smr_spread(c(7L, 3L), 10L, c(1L, 1L), 2L), 4)
    expect_equal(sandwichrand:::.smr_spread(c(4L, 4L, 2L), 10L, c(2L, 2L, 1L), 5L), 0)
    expect_equal(sandwichrand:::.smr_spread(c(5L, 4L, 1L), 10L, c(2L, 2L, 1L), 5L), 2)
})

test_that(".smr_parse_arms understands every accepted form", {
    p <- sandwichrand:::.smr_parse_arms
    expect_equal(p(c("A", "B"))$mult, c(1L, 1L))
    expect_equal(p(c("A", "A", "B", "B", "C"))$mult, c(2L, 2L, 1L))
    expect_equal(p(c("A", "A", "B", "B", "C"))$labels, c("A", "B", "C"))
    expect_equal(p("A A B B C")$ratio_sum, 5L)
    expect_equal(p(c(Test = 2, Placebo = 2, C = 1))$mult, c(2L, 2L, 1L))
    expect_equal(p(c(Test = 2, Placebo = 2, C = 1))$labels,
                 c("Test", "Placebo", "C"))
    expect_equal(p(c(2, 2, 1))$labels, c("A", "B", "C"))
    expect_error(p(c("a", "A")), "at least 2 DISTINCT")
    expect_error(p(c(1, -1)), "positive whole")
    expect_error(p(c(A = 1, A = 1)), "unique")
})

test_that(".smr_n1_candidates only returns admissible multiples of ratio_sum", {
    spec <- sandwichrand:::.smr_spec(N = 48L, arms = c("A", "B"))
    cand <- sandwichrand:::.smr_n1_candidates(spec)
    expect_true(all(cand %% 2L == 0L))
    expect_true(all(cand >= 4L))
    expect_true(all(cand <= 24L))
    expect_true(all(spec$N - 2L * cand == 0L |
                        spec$N - 2L * cand >= spec$vmin))
})

test_that(".smr_streams gives n distinct sub-streams from one seed", {
    s <- sandwichrand:::.smr_streams(123L, 4L)
    expect_length(s, 4L)
    for (i in 1:4) expect_true(length(s[[i]]) > 1L)
    expect_true(!identical(s[[1L]], s[[2L]]))
    expect_true(!identical(s[[2L]], s[[3L]]))
    expect_true(!identical(s[[3L]], s[[4L]]))
    expect_identical(s, sandwichrand:::.smr_streams(123L, 4L))
})

test_that(".smr_tabulate counts in the declared arm order", {
    expect_equal(sandwichrand:::.smr_tabulate(c("A", "B", "A"), c("A", "B")),
                 c(2L, 1L))
    expect_equal(sandwichrand:::.smr_tabulate(c("B", "B"), c("A", "B")), c(0L, 2L))
})

test_that(".smr_rtf_escape escapes the RTF special characters", {
    e <- sandwichrand:::.smr_rtf_escape
    expect_equal(e("a\\b"), "a\\\\b")
    expect_equal(e("{x}"), "\\{x\\}")
    expect_equal(e("plain"), "plain")
    expect_equal(length(e(character(0))), 0L)
})
