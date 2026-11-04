library(BCalm)

# Tests for mpra_treat()
#
# mpra_treat() is a variant of limma::treat() whose null hypothesis is derived from negative controls:
#   1. all coefficients are centred on the mean of the negative controls,
#   2. R = `percentile` quantile and L = `1 - percentile` quantile of the centred controls,
#   3. the null is "true centred effect in [L, R]". For skewed controls |L| != R.
#
# With C the centred coefficient, se = stdev.unscaled * sqrt(s2.post), T ~ t(df.total) and mid = (L + R) / 2:
#   C >= mid: p = P(T >= (C - R) / se) + P(T >= (C - L) / se)   (Up side, plus the far tail of the other side)
#   C <  mid: p = P(T <= (C - L) / se) + P(T <= (C - R) / se)   (Down side, plus the far tail of the other side)
# i.e. the first term tests the boundary on the coefficient's own side, the second term is the probability
# that the opposite boundary would give an equally extreme result. t is (C - R) / se for C > R,
# (C - L) / se for C < L (negative) and 0 inside [L, R].
#
# Background: a previous version used abs(C) with the signed thresholds for every coefficient. That is only
# correct for C > 0. Negative coefficients were tested against R instead of |L| (and vice versa), which is
# invisible for symmetric controls (|L| == R) but wrong for skewed ones. These tests guard against that.

## Helpers -----------------------------------------------------------------------------------------------------

# Builds a minimal MArrayLM as produced by fit_elements(); `label` marks the negative controls.
make_fit <- function(coef, label, se = 0.1, df = 20, sigma = rep(1, length(coef))) {
    n <- length(coef)
    new("MArrayLM", list(
        coefficients = matrix(coef, ncol = 1),
        stdev.unscaled = matrix(se, n, 1),
        sigma = sigma,
        df.residual = rep(df, n),
        label = label
    ))
}

# Reference implementation, written independently of mpra_treat(): signed distances to the two null
# boundaries, no abs(), explicit upper/lower tails. `centred` is the coefficient minus the control mean.
reference_p <- function(centred, se, df, R, L) {
    up <- centred >= (R + L) / 2
    ifelse(
        up,
        pt((centred - R) / se, df, lower.tail = FALSE) + pt((centred - L) / se, df, lower.tail = FALSE),
        pt((centred - L) / se, df, lower.tail = TRUE) + pt((centred - R) / se, df, lower.tail = TRUE)
    )
}

reference_t <- function(centred, se, R, L) {
    ifelse(centred > R, (centred - R) / se, ifelse(centred < L, (centred - L) / se, 0))
}

# Negative-control distributions. Only the first one is symmetric.
control_generators <- list(
    symmetric = function(n) rnorm(n, 0, 0.2),
    right_skewed = function(n) rgamma(n, 2, scale = 0.15),
    left_skewed = function(n) -rgamma(n, 2, scale = 0.15),
    heavy_tailed = function(n) 0.1 * rt(n, 3),
    bimodal = function(n) c(rnorm(n / 2, -0.3, 0.05), rnorm(n / 2, 0.5, 0.05))
)
percentiles <- c(0.5, 0.85, 0.95, 0.975)

# One simulated data set: controls first, then tested elements spread over and beyond the null region.
simulate_fit <- function(generator, n_ctrl = 400, n_test = 300, shift = 0) {
    ctrl <- generator(n_ctrl)
    coef <- c(ctrl, runif(n_test, -1.5, 2)) + shift
    n <- length(coef)
    label <- c(rep("neg", n_ctrl), rep("test", n_test))
    # varying sigma so that variance moderation (squeezeVar) is actually exercised
    make_fit(coef, label, sigma = sqrt(rchisq(n, 20) / 20))
}

## 1. Agreement with the independent reference implementation ---------------------------------------------------
# Covers symmetric and skewed controls and a range of percentiles (including the extreme 0.5, where for skewed
# controls both thresholds can even be on the same side of 0, which is why the side is chosen by the midpoint).

for (scenario in names(control_generators)) {
    for (pct in percentiles) {
        test_that(sprintf("mpra_treat matches the reference p-values and t (%s controls, percentile %s)", scenario, pct), {
            set.seed(1)
            fit <- simulate_fit(control_generators[[scenario]])
            res <- mpra_treat(fit, percentile = pct, neg_label = "neg")

            is_neg <- fit$label == "neg"
            centred <- fit$coefficients[, 1] - mean(fit$coefficients[is_neg, 1])
            R <- unname(quantile(centred[is_neg], pct))
            L <- unname(quantile(centred[is_neg], 1 - pct))
            se <- fit$stdev.unscaled[, 1] * sqrt(res$s2.post)

            expect_equal(unname(res$treat.lfc_right), R)
            expect_equal(unname(res$treat.lfc_left), L)
            expect_equal(unname(res$p.value[, 1]), reference_p(centred, se, res$df.total, R, L))
            expect_equal(unname(res$t[, 1]), reference_t(centred, se, R, L))
        })
    }
}

## 2. Exact hand-computed example ------------------------------------------------------------------------------
# 11 controls with mean exactly 0, sorted x1..x11. For percentile 0.9 the type-7 quantile index is 1 + 10 * 0.9 = 10
# and for 0.1 it is 2, so R = x10 = 0.6 and L = x2 = -0.2 exactly (skewed: |L| = 0.2 << R = 0.6). mid = 0.2.
# With the old code, a coefficient at -0.2 was compared with 0.6 and looked hopelessly non-significant.

test_that("mpra_treat reproduces hand-computed values for skewed controls", {
    ctrl <- c(-1.0, -0.2, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.6, 0.6)
    expect_equal(sum(ctrl), 0)
    tested <- c(1.0, 0.6, 0.2, 0, -0.2, -0.5, -1.0)
    fit <- make_fit(c(ctrl, tested), c(rep("neg", 11), rep("test", 7)))
    res <- mpra_treat(fit, percentile = 0.9, neg_label = "neg")

    expect_equal(unname(res$treat.lfc_right), 0.6)
    expect_equal(unname(res$treat.lfc_left), -0.2)

    se <- 0.1 * sqrt(res$s2.post[1]) # identical sigma for all rows, so one scalar
    df <- res$df.total[1]
    upper <- function(x) pt(x / se, df, lower.tail = FALSE)
    rows <- 12:18

    # Each line: own-side distance beyond its boundary + distance to the far boundary (as in the formula header).
    expected_p <- c(
        upper(1.0 - 0.6) + upper(1.0 + 0.2), # C = 1.0:  Up,   beyond R by 0.4; far boundary L at distance 1.2
        upper(0.6 - 0.6) + upper(0.6 + 0.2), # C = 0.6:  exactly on R
        upper(0.2 - 0.6) + upper(0.2 + 0.2), # C = 0.2:  exactly the midpoint, p = 1
        (1 - upper(0.2 - 0)) + upper(0.6 - 0), # C = 0: below the midpoint, lower tails: P(T<=(C-L)/se) + P(T<=(C-R)/se)
        upper(-0.2 + 0.2) + upper(0.2 + 0.6), # C = -0.2: exactly on L
        upper(0.5 - 0.2) + upper(0.5 + 0.6), # C = -0.5: Down, beyond L by 0.3
        upper(1.0 - 0.2) + upper(1.0 + 0.6) # C = -1.0: Down, beyond L by 0.8
    )
    expect_equal(unname(res$p.value[rows, 1]), expected_p)

    # Spot checks that do not depend on s2.post: boundary elements have p ~ 0.5, the midpoint has p = 1.
    expect_equal(unname(res$p.value[13, 1]), 0.5, tolerance = 1e-3)
    expect_equal(unname(res$p.value[16, 1]), 0.5, tolerance = 1e-3)
    expect_equal(unname(res$p.value[14, 1]), 1)

    expected_t <- c(0.4 / se, 0, 0, 0, 0, -0.3 / se, -0.8 / se)
    expect_equal(unname(res$t[rows, 1]), expected_t)
})

## 3. Mirror symmetry ------------------------------------------------------------------------------------------
# Flipping the sign of all coefficients must exactly swap Up and Down: same p-values, opposite t, swapped
# thresholds. This cannot hold for an implementation that treats the two sides differently by mistake,
# and it was violated by the old code for every skewed control distribution.

for (scenario in names(control_generators)) {
    test_that(sprintf("mpra_treat is mirror-symmetric (%s controls)", scenario), {
        set.seed(2)
        fit <- simulate_fit(control_generators[[scenario]])
        mirrored <- fit
        mirrored$coefficients <- -fit$coefficients
        for (pct in percentiles) {
            res <- mpra_treat(fit, percentile = pct, neg_label = "neg")
            res_m <- mpra_treat(mirrored, percentile = pct, neg_label = "neg")
            expect_equal(res$p.value, res_m$p.value)
            expect_equal(res$t, -res_m$t)
            expect_equal(unname(res$treat.lfc_right), -unname(res_m$treat.lfc_left))
            expect_equal(unname(res$treat.lfc_left), -unname(res_m$treat.lfc_right))
        }
    })
}

## 4. Equivalence with limma::treat for symmetric thresholds ----------------------------------------------------
# If the centred controls are exactly symmetric (R = -L), mpra_treat must reduce to limma::treat(lfc = R). This
# anchors the sign convention (down = negative t) and the variance moderation to the original implementation.

test_that("mpra_treat equals limma::treat when the controls are exactly symmetric", {
    set.seed(3)
    x <- rgamma(200, 2, scale = 0.15)
    ctrl <- c(x, -x) # mean 0, quantiles symmetric by construction
    fit <- make_fit(c(ctrl, runif(300, -1.5, 1.5)), c(rep("neg", 400), rep("test", 300)),
        sigma = sqrt(rchisq(700, 20) / 20)
    )
    for (pct in c(0.85, 0.95, 0.975)) {
        res <- mpra_treat(fit, percentile = pct, neg_label = "neg")
        expect_equal(unname(res$treat.lfc_right), -unname(res$treat.lfc_left))
        lim <- limma::treat(fit, lfc = unname(res$treat.lfc_right))
        expect_equal(res$p.value, lim$p.value)
        expect_equal(res$t, lim$t)
    }
})

test_that("mpra_treat uses the same variance moderation as limma::treat, also for skewed controls", {
    set.seed(4)
    fit <- simulate_fit(control_generators$right_skewed)
    res <- mpra_treat(fit, percentile = 0.95, neg_label = "neg")
    lim <- limma::treat(fit, lfc = 0.1)
    expect_equal(res$s2.post, lim$s2.post)
    expect_equal(res$df.total, lim$df.total)
    expect_equal(res$df.prior, lim$df.prior)
})

## 5. Structural properties of the p-values and t statistics ---------------------------------------------------

for (scenario in names(control_generators)) {
    test_that(sprintf("p-values and t statistics are internally consistent (%s controls)", scenario), {
        set.seed(5)
        fit <- simulate_fit(control_generators[[scenario]])
        for (pct in percentiles) {
            res <- mpra_treat(fit, percentile = pct, neg_label = "neg")
            R <- unname(res$treat.lfc_right)
            L <- unname(res$treat.lfc_left)
            centred <- fit$coefficients[, 1] - mean(fit$coefficients[fit$label == "neg", 1])
            p <- res$p.value[, 1]
            tt <- res$t[, 1]
            eps <- 1e-12

            # valid probabilities
            expect_true(all(is.finite(p)))
            expect_true(all(p >= 0 & p <= 1 + eps))
            # thresholds are ordered
            expect_true(L <= R)
            # sign convention of t: positive iff above R, negative iff below L, zero inside the null
            expect_true(all(tt[centred > R] > 0))
            expect_true(all(tt[centred < L] < 0))
            expect_true(all(tt[centred >= L & centred <= R] == 0))
            # inside the null interval an element can never look significant
            inside <- centred >= L & centred <= R
            expect_true(all(p[inside] >= 0.5 - eps))
            # elements at / beyond a boundary: p < 0.5 iff strictly outside the null (up to the tiny far-tail term)
            expect_true(all(p[centred > R + 1e-6 & tt > 3] < 0.01))
            expect_true(all(p[centred < L - 1e-6 & tt < -3] < 0.01))
        }
    })
}

test_that("p-value decreases monotonically moving away from the null on each side", {
    set.seed(6)
    ctrl <- rgamma(400, 2, scale = 0.15)
    grid <- seq(-2, 2, length.out = 401)
    fit <- make_fit(c(ctrl, grid + mean(ctrl)), c(rep("neg", 400), rep("test", 401)))
    res <- mpra_treat(fit, percentile = 0.95, neg_label = "neg")
    p <- res$p.value[401:801, 1] # ordered by increasing centred coefficient (= grid)
    R <- unname(res$treat.lfc_right)
    L <- unname(res$treat.lfc_left)
    mid <- (R + L) / 2
    # p falls to the left of the midpoint and rises towards it from the left, mirrored on the right
    expect_true(all(diff(p[grid <= L]) >= -1e-12)) # grid increases towards L: p increases (decreases going left)
    expect_true(all(diff(p[grid >= R]) <= 1e-12)) # moving right, away from R
    left_part <- p[grid <= mid]
    right_part <- p[grid >= mid]
    expect_true(all(diff(left_part) >= -1e-12)) # increases towards the midpoint
    expect_true(all(diff(right_part) <= 1e-12)) # decreases away from the midpoint
})

test_that("a higher percentile never makes any element look more significant (|t| is non-increasing)", {
    set.seed(7)
    for (scenario in names(control_generators)) {
        fit <- simulate_fit(control_generators[[scenario]])
        t_low <- mpra_treat(fit, percentile = 0.85, neg_label = "neg")$t
        t_high <- mpra_treat(fit, percentile = 0.975, neg_label = "neg")$t
        expect_true(all(abs(t_high) <= abs(t_low) + 1e-12))
    }
})

## 6. Invariances ----------------------------------------------------------------------------------------------

test_that("results do not depend on the row order", {
    set.seed(8)
    fit <- simulate_fit(control_generators$right_skewed)
    n <- nrow(fit$coefficients)
    perm <- sample(n)
    permuted <- fit
    permuted$coefficients <- fit$coefficients[perm, , drop = FALSE]
    permuted$stdev.unscaled <- fit$stdev.unscaled[perm, , drop = FALSE]
    permuted$sigma <- fit$sigma[perm]
    permuted$df.residual <- fit$df.residual[perm]
    permuted$label <- fit$label[perm]
    res <- mpra_treat(fit, 0.95, "neg")
    res_p <- mpra_treat(permuted, 0.95, "neg")
    expect_equal(res_p$p.value, res$p.value[perm, , drop = FALSE])
    expect_equal(res_p$t, res$t[perm, , drop = FALSE])
    expect_equal(res_p$treat.lfc_right, res$treat.lfc_right)
})

test_that("results do not depend on the overall offset of the coefficients", {
    # centring on the control mean removes any global shift (e.g. a reporter-wide effect)
    set.seed(9)
    fit <- simulate_fit(control_generators$left_skewed)
    shifted <- fit
    shifted$coefficients <- fit$coefficients + 3.7
    res <- mpra_treat(fit, 0.95, "neg")
    res_s <- mpra_treat(shifted, 0.95, "neg")
    expect_equal(res_s$p.value, res$p.value)
    expect_equal(res_s$t, res$t)
    expect_equal(res_s$treat.lfc_right, res$treat.lfc_right)
})

test_that("results are scale-equivariant (changing the units of logFC and its se changes nothing)", {
    set.seed(10)
    fit <- simulate_fit(control_generators$right_skewed)
    scaled <- fit
    scaled$coefficients <- 2.5 * fit$coefficients
    scaled$stdev.unscaled <- 2.5 * fit$stdev.unscaled
    res <- mpra_treat(fit, 0.95, "neg")
    res_s <- mpra_treat(scaled, 0.95, "neg")
    expect_equal(res_s$p.value, res$p.value)
    expect_equal(res_s$t, res$t)
    expect_equal(res_s$treat.lfc_right, 2.5 * res$treat.lfc_right)
    expect_equal(res_s$treat.lfc_left, 2.5 * res$treat.lfc_left)
})

test_that("only the negative controls define the thresholds", {
    set.seed(11)
    fit <- simulate_fit(control_generators$right_skewed)
    other <- fit
    test_rows <- which(fit$label == "test")
    other$coefficients[test_rows, 1] <- other$coefficients[test_rows, 1] + 10 * rnorm(length(test_rows))
    res <- mpra_treat(fit, 0.95, "neg")
    res_o <- mpra_treat(other, 0.95, "neg")
    expect_equal(res_o$treat.lfc_right, res$treat.lfc_right)
    expect_equal(res_o$treat.lfc_left, res$treat.lfc_left)
    # the controls' own p-values are unchanged too (the variance moderation uses all rows, so we use the thresholds)
})

## 7. Statistical calibration ----------------------------------------------------------------------------------
# For an element whose TRUE effect lies exactly on a boundary of the null, the test must reject at p < 0.05 in
# about 5% of repetitions, on BOTH sides, for skewed controls too. Noise has sd = stdev.unscaled (sigma = 1).
# The old code gave about 0.1% (right-skewed) and about 58% (left-skewed) on the left boundary.
# With 5000 repetitions the binomial sd is 0.003, so a tolerance of 0.02 is more than 6 sd.

simulate_boundary_rate <- function(controls, side, pct = 0.95, n_rep = 5000, se = 0.1, df = 20) {
    mu <- mean(controls)
    centred <- controls - mu
    boundary <- if (side == "right") quantile(centred, pct) else quantile(centred, 1 - pct)
    true_effect <- mu + unname(boundary)
    coef <- c(controls, true_effect + rnorm(n_rep, 0, se))
    n <- length(coef)
    fit <- make_fit(coef, c(rep("neg", length(controls)), rep("test", n_rep)),
        se = se, df = df, sigma = sqrt(rchisq(n, df) / df)
    )
    res <- mpra_treat(fit, percentile = pct, neg_label = "neg")
    mean(res$p.value[fit$label == "test", 1] < 0.05)
}

for (scenario in c("symmetric", "right_skewed", "left_skewed", "heavy_tailed")) {
    for (side in c("right", "left")) {
        test_that(sprintf("type I error is ~5%% on the %s boundary (%s controls)", side, scenario), {
            set.seed(12)
            ctrl <- control_generators[[scenario]](1000)
            rate <- simulate_boundary_rate(ctrl, side)
            expect_gt(rate, 0.03)
            expect_lt(rate, 0.07)
        })
    }
}

test_that("elements well inside the null are (almost) never called, elements well outside (almost) always", {
    set.seed(13)
    ctrl <- rgamma(1000, 2, scale = 0.15)
    mu <- mean(ctrl)
    centred <- ctrl - mu
    R <- unname(quantile(centred, 0.95))
    L <- unname(quantile(centred, 0.05))
    n_rep <- 2000
    coef <- c(
        ctrl,
        mu + (L + R) / 2 + rnorm(n_rep, 0, 0.02), # centre of the null
        mu + R + 1 + rnorm(n_rep, 0, 0.02), # far above
        mu + L - 1 + rnorm(n_rep, 0, 0.02) # far below
    )
    label <- c(rep("neg", 1000), rep(c("centre", "up", "down"), each = n_rep))
    fit <- make_fit(coef, label, se = 0.02, df = 20, sigma = sqrt(rchisq(length(coef), 20) / 20))
    res <- mpra_treat(fit, 0.95, "neg")
    p <- res$p.value[, 1]
    expect_lt(mean(p[label == "centre"] < 0.05), 0.001)
    expect_gt(mean(p[label == "up"] < 0.05), 0.999)
    expect_gt(mean(p[label == "down"] < 0.05), 0.999)
})

## 8. Input handling -------------------------------------------------------------------------------------------

test_that("mpra_treat validates its input", {
    set.seed(14)
    fit <- simulate_fit(control_generators$symmetric)
    expect_error(mpra_treat(list(), neg_label = "neg"), "MArrayLM")
    no_label <- fit
    no_label$label <- NULL
    expect_error(mpra_treat(no_label, neg_label = "neg"), "label")
    expect_error(mpra_treat(fit, neg_label = "not_a_label"))
    no_coef <- fit
    no_coef$coefficients <- NULL
    expect_error(mpra_treat(no_coef, neg_label = "neg"), "coefficients")
})

test_that("mpra_treat keeps dimensions and names of the coefficients", {
    set.seed(15)
    fit <- simulate_fit(control_generators$symmetric)
    rownames(fit$coefficients) <- paste0("e", seq_len(nrow(fit$coefficients)))
    rownames(fit$stdev.unscaled) <- rownames(fit$coefficients)
    res <- mpra_treat(fit, 0.95, "neg")
    expect_equal(dim(res$p.value), dim(fit$coefficients))
    expect_equal(dim(res$t), dim(fit$coefficients))
    expect_equal(rownames(res$p.value), rownames(fit$coefficients))
    expect_equal(rownames(res$t), rownames(fit$coefficients))
})
