# sandwichrand 1.0.0

## New features

* `smr_alloc()` generates allocation sequences by Sandwich Mixed Randomization
  (SMR), with the three generation paths of the methodology paper:
  * Path 1: `N < 4 * ratio_sum` — complete randomization in one block;
  * Path 2: `n1 + permuted blocks with varying sizes + n1` — the full sandwich;
  * Path 3: `N == 2 * n1` — head and tail blocks only.

* The head-block size `n1` is always sampled at random (no need to
  pre-specify it, and no fixed `n1` override is accepted).

* The imbalance is the **spread** of the deviations from the target allocation
  (`spread = max_g(dev_g) - min_g(dev_g)`), which equals the between-group count
  difference for 1:1 and is zero for an exactly proportional block.

* `smr_streams()` derives non-overlapping L'Ecuyer-CMRG sub-streams from one
  base seed, so high-frequency loops never re-use the same random numbers.

* `smr_reproduce()` regenerates a sequence from its stored specification.

* `smr_verify()`, `smr_balance()`, `smr_blocks()`, `smr_imbalance()`,
  `smr_spread()` and `smr_simulate()` provide verification and diagnostics.

* `smr_correct_guesses()` and `smr_pcg()` evaluate allocation
  unpredictability with the Blackwell-Hodges proportion of correct guesses;
  `smr_simulate()` now reports this alongside the imbalance metrics.

* `smr_rtf()` and `smr_write_csv()` release the list, honouring `keep.blocks`
  for ICH E9 / GCP compliance.

## Notes

* This package is an optimization and refactor of the earlier `sandwichrand`
  package and of the SAS macro `SanMR2`. The changes focus on (i) random
  head-block size, (ii) the spread-based imbalance measure, (iii) independent
  per-call random streams, and (iv) a cleaner and faster core.
