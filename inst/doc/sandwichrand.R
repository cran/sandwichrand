## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(sandwichrand)

## -----------------------------------------------------------------------------
library(sandwichrand)

## -----------------------------------------------------------------------------
x <- smr_alloc(N = 48, arms = c("TestTrt", "Placebo"), seed = 84529052)
print(x, n = 8)

## -----------------------------------------------------------------------------
x$balance
x$blocks

## -----------------------------------------------------------------------------
y <- smr_alloc(N = 150, arms = c("A", "A", "B", "B", "C"), seed = 123456)
y$balance

## -----------------------------------------------------------------------------
z <- smr_alloc(N = 60, arms = c(TestTrt = 2, Placebo = 2, Control = 1), seed = 7)
z$design[c("arms", "ratio", "ratio_sum")]

## ----eval = FALSE-------------------------------------------------------------
# x <- smr_alloc(N = 48, arms = c("TestTrt", "Placebo"), seed = 84529052,
#                keep.blocks = FALSE)
# names(x$alloc)             # "id" "arm" only -- no block information
# smr_write_csv(x, "randomization.csv")
# smr_rtf(x, "randomization.rtf")

## -----------------------------------------------------------------------------
centers <- data.frame(
  center = c("Site 01", "Site 02", "Site 03"),
  n      = c(48, 36, 24)            # planned enrolment per centre
)
streams <- smr_streams(seed = 20260301, n = nrow(centers))

## One independent SMR run per centre, then one merged list
rand_list <- do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
  x <- smr_alloc(N = centers$n[i], arms = c("TestTrt", "Placebo"),
                 stream = streams[[i]])
  cbind(center = centers$center[i], x$alloc)
}))
rand_list$subject <- seq_len(nrow(rand_list))   # global sequential id

head(rand_list[, c("center", "subject", "id", "arm")], 10)

## -----------------------------------------------------------------------------
## Balance is exact within every centre (and therefore overall)
with(rand_list, table(center, arm))

## ----eval = FALSE-------------------------------------------------------------
# release <- rand_list[, c("center", "subject", "id", "arm")]
# write.csv(release, "randomization_multicenter.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
v <- smr_verify(x)
v
all(v$status == "PASS")

## ----fig.width = 6, fig.height = 4--------------------------------------------
plot(x, type = "imbalance")

## -----------------------------------------------------------------------------
s <- smr_simulate(N = 48, arms = c("A", "B"), reps = 1000, seed = 11)
s

## -----------------------------------------------------------------------------
cg <- smr_correct_guesses(x)
smr_pcg(x)            # proportion of correct guesses (0.5 = complete randomization)

## -----------------------------------------------------------------------------
x1 <- smr_alloc(N = 48, arms = c("A", "B"), seed = 2026)
x2 <- smr_alloc(N = 48, arms = c("A", "B"), seed = 2026)
identical(x1$alloc, x2$alloc)      # TRUE
identical(smr_reproduce(x1)$alloc, x1$alloc)   # TRUE

## -----------------------------------------------------------------------------
ss <- smr_streams(seed = 42, n = 3)
a <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[1]])
b <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[2]])
!identical(a$alloc, b$alloc)       # independent streams -> different sequences

