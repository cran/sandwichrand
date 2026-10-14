# sandwichrand

**Allocation sequences by Sandwich Mixed Randomization (SMR) for randomized
controlled trials.**

`sandwichrand` generates treatment allocation sequences (randomization lists) with
the Sandwich Mixed Randomization scheme described in Wang et al. (2026),
*BMC Medical Research Methodology* **26**(1), 152,
[doi:10.1186/s12874-026-02877-2](https://doi.org/10.1186/s12874-026-02877-2).

## The sandwich

An SMR sequence is a sandwich:

```
|<--- n1 --->|<-------- N - 2*n1 -------->|<--- n1 --->|
+------------+----------------------------+------------+
| block 1    | blocks 2 .. K              | block K+1  |
| imbalanced | permuted blocks of         | complement |
|            | varying sizes              | of block 1 |
+------------+----------------------------+------------+
```

* The **head block** (size `n1`) is deliberately imbalanced: it is drawn by
  repeated complete randomization until its imbalance (the *spread* of the
  deviations from the target allocation) lies in `[1, 4]`. `n1` is sampled at
  random from the admissible multiples of `ratio_sum`, so there is no need to
  pre-specify it.
* The **middle** is filled with permuted blocks of varying sizes (the PBVB
  procedure).
* The **tail block** (size `n1`) is the exact complement of the head block, so
  overall balance is guaranteed while the next allocation stays unpredictable.

## Installation

```r
# from the source tarball
install.packages("sandwichrand_1.0.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

```r
library(sandwichrand)

## Two-arm 1:1 trial, N = 48 -- the illustration of the SMR paper
x <- smr_alloc(N = 48, arms = c("TestTrt", "Placebo"), seed = 84529052)
x
head(x$alloc)
x$balance

## Three-arm 2:2:1 trial, N = 150, released without the block structure
y <- smr_alloc(N = 150, arms = c("A", "A", "B", "B", "C"), seed = 123456,
               keep.blocks = FALSE)
smr_balance(y)

## Independent streams for high-frequency calling (simulation loops)
ss <- smr_streams(seed = 42, n = 3)
a <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[1]])
b <- smr_alloc(N = 48, arms = c("A", "B"), stream = ss[[2]])

## Verify, document and release
smr_verify(x)
smr_write_csv(x, file = "randomization.csv")
smr_rtf(x, file = "randomization.rtf")
```

## Evaluating the design

A randomization procedure is assessed on two performance metrics (Wang et al.
2026):

* **Allocation balance** — the absolute imbalance `E|D(i)|`,
  `D(i) = n_A(i) - n_B(i)`:

  ```r
  smr_balance(x)            # final arm counts and deviations
  smr_imbalance(x)          # running |D(i)| at every allocation step
  ```

* **Allocation unpredictability** — the Blackwell–Hodges proportion of correct
  guesses (PCG). Under the convergence strategy the next assignment is
  predicted to go to the arm with the fewer prior allocations (`0.5` on ties):

  ```r
  smr_correct_guesses(x)    # per-step correct-guess indicators
  smr_pcg(x)                # proportion of correct guesses (0.5 = CRD)
  ```

Both are summarised across replications by `smr_simulate()`:

```r
smr_simulate(N = 48, arms = c("A", "B"), reps = 1000, seed = 11)
# final imbalance (balance) and correct-guess PCG (unpredictability)
```

## Design highlights

* **Random head-block size.** `n1` is always sampled at random, so the block
  structure is not revealed in advance and there is no `n1` argument to specify.
* **Spread-based imbalance.** The imbalance is the *spread*
  `max_g(count_g - n*mult_g/ratio_sum) - min_g(...)`, which equals the
  between-group count difference for 1:1 and is zero for an exactly
  proportional block.
* **Independent streams.** Each call draws from an independent
  L'Ecuyer-CMRG sub-stream; `smr_streams()` derives a family of
  non-overlapping streams up front so high-frequency loops never reuse the same
  random numbers.
* **Reproducible.** Every sequence is a deterministic function of a single seed
  (or of a supplied stream); `smr_reproduce()` regenerates it.
* **ICH E9 / GCP aware.** `keep.blocks = FALSE` releases a list without any
  block information, while the statistician's copy stays in `$blocks`.

## Documentation

The package vignette (`vignette("sandwichrand")`) follows the structure
recommended for a *Journal of Statistical Software* contribution:
motivation → methods → implementation → usage → validation → reproducibility →
computational details → summary.

## Clinical use and validation disclaimer

`sandwichrand` is research software provided **"AS IS"**, without warranty of
any kind, express or implied. It has **not been validated** as a computerised
system under Good Clinical Practice (GCP), Good Clinical/Manufacturing
Practice (GxP), 21 CFR Part 11, EU GMP Annex 11, or any other regulatory
framework.

Before using `sandwichrand` to generate allocation sequences for a
registration trial, or for any trial whose results may support a regulatory
submission, the sponsor is solely responsible for:

* independently validating the software and the generated sequences under its
  own quality system;
* maintaining allocation concealment (ICH E9), including keeping the seed,
  the block structure and any stream information confidential; and
* ensuring compliance with all applicable regulations.

The authors accept no liability for any consequences arising from the use of
this software.

## License

Apache License (>= 2.0). The SMR scheme was designed by Bingshun Wang.
