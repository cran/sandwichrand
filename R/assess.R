## ---------------------------------------------------------------------------
## Evaluation of allocation unpredictability (Blackwell-Hodges model).
##
## Allocation balance is evaluated with smr_balance(), smr_imbalance() and
## smr_simulate(). Allocation unpredictability is evaluated here with the
## "convergence" guessing strategy: at every step the next assignment is
## predicted to go to the arm that is most under-allocated relative to the
## target ratio. The proportion of correct guesses (PCG) is the standard
## Blackwell-Hodges measure of allocation predictability [Blackwell & Hodges
## 1957; Wang et al. 2026].
## ---------------------------------------------------------------------------

#' Per-step correct guesses under the convergence strategy
#'
#' Evaluates allocation unpredictability of a generated sequence with the
#' Blackwell-Hodges \dQuote{convergence} guessing strategy. At each allocation
#' step the next assignment is predicted to go to the treatment arm that is
#' most under-allocated relative to the target allocation ratio; ties are
#' broken uniformly at random. \code{CG(i)} equals 1 when the guess is right,
#' 0 when it is wrong, and \code{1/k} when \code{k} arms are tied for the
#' largest deficit. For a two-arm 1:1 trial this is exactly the definition of
#' Wang et al. (2026): guess the arm with the fewer prior allocations, with
#' \code{CG(i) = 0.5} when the two counts are equal.
#'
#' @param x An object of class \code{"sandwichrand"}.
#'
#' @return A numeric vector of length \code{N} with the per-step correct-guess
#'   indicators \code{CG(1), ..., CG(N)}. The first element is always
#'   \code{1/n_arms} (no prior information), which is \code{0.5} for 1:1.
#'
#' @references
#' Blackwell D, Hodges JL (1957). Design for the control of selection bias.
#' \emph{Annals of Mathematical Statistics}, \bold{28}(2), 449--460.
#'
#' Wang B, Wang X, Ni C, Xie J, Wu S, Wu Y, Shen Y, Dong Z, Lv W, He Y (2026).
#' \dQuote{Reducing selection bias risk to enhance RCT validity: sandwich mixed randomization outperforms permuted block design.}
#' \emph{BMC Medical Research Methodology}, \bold{26}(1), 152.
#' \doi{10.1186/s12874-026-02877-2}.
#'
#' @seealso \code{\link{smr_pcg}}, \code{\link{smr_imbalance}},
#'   \code{\link{smr_simulate}}.
#'
#' @examples
#' x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 84529052)
#' smr_correct_guesses(x)
#'
#' @export
smr_correct_guesses <- function(x) {
    if (!inherits(x, "sandwichrand")) {
        stop("`x` must be an object of class \"sandwichrand\".", call. = FALSE)
    }
    a <- x$alloc
    d <- x$design
    N <- d$N
    mult <- d$mult
    rs <- d$ratio_sum
    ng <- d$n_arms
    mi <- match(a$arm, d$arms)

    cg <- numeric(N)
    cnt <- integer(ng)                 # cumulative counts BEFORE step i
    for (i in seq_len(N)) {
        ## Scaled deficit = (i-1)*mult_g - cnt_g*ratio_sum. Positive when arm g
        ## is under-allocated relative to the target after i-1 allocations. It
        ## is an integer, so ties (== max) are exact.
        sd <- (i - 1L) * mult - cnt * rs
        guess_set <- which(sd == max(sd))
        k <- length(guess_set)
        cg[i] <- if (mi[i] %in% guess_set) 1 / k else 0
        cnt[mi[i]] <- cnt[mi[i]] + 1L
    }
    cg
}

#' Proportion of correct guesses (allocation predictability)
#'
#' The mean of \code{\link{smr_correct_guesses}}, i.e. the expected proportion
#' of correct guesses (PCG). A value close to \code{0.5} means the sequence is
#' about as unpredictable as complete randomization; values well above
#' \code{0.5} mean the next allocation can be anticipated from the running
#' imbalance, which is the predictability that the SMR scheme is designed to
#' reduce.
#'
#' @param x An object of class \code{"sandwichrand"}.
#'
#' @return A single number in \code{[0, 1]}, the proportion of correct guesses.
#'
#' @examples
#' x <- smr_alloc(N = 48, arms = c("A", "B"), seed = 84529052)
#' smr_pcg(x)
#'
#' ## Expected PCG over many replications
#' s <- smr_simulate(N = 48, arms = c("A", "B"), reps = 200, seed = 11)
#' s
#'
#' @seealso \code{\link{smr_correct_guesses}}, \code{\link{smr_simulate}}.
#'
#' @export
smr_pcg <- function(x) {
    mean(smr_correct_guesses(x))
}
