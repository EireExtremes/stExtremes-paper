##==============================================================================
## Verify the bGEV return-level computation used in 02_compare_barrier.R.
##
## Source this AFTER 02_compare_barrier.R: it reuses fits, model_pars,
## idat, dc_max, last_year, p_a and p_b. It demonstrates the parametrisation
## bug and confirms the fix. Three diagnostics per model:
##   - bug_z2_vs_loc : max |2-yr level - location| using return_level_bgev()
##     fed the NEW-parametrisation inputs (q_alpha, s_beta) as if they were
##     (mu, sigma). This is the bug; it is ~0.10 here, NOT 0.
##   - fix_z2_vs_loc : same, using return_level_bgev2(). Must be ~0, because
##     with q.location = 0.5 the 2-yr level equals the median location.
##   - fix_vs_ref_50 : max |return_level_bgev2 - closed-form reference| at
##     50 yr. Must be ~0.
##==============================================================================

## Closed-form upper-tail reference (equals new_to_old() + qgev(); the q.mix
## probabilities do not enter the upper tail). alpha = q.location, beta =
## q.spread; s_beta is the central (1 - beta) interval width.
rl_bgev_ref <- function(period, location, spread, tail,
    alpha = 0.5, beta = 0.5) {
    p <- 1 - 1 / period
    if(abs(tail) < 1e-8) {
        num <- log(-log(alpha)) - log(-log(p))
        den <- log(-log(beta / 2)) - log(-log(1 - beta / 2))
    } else {
        av <- function(q) (-log(q))^(-tail)
        num <- av(p) - av(alpha)
        den <- av(1 - beta / 2) - av(beta / 2)
    }
    location + spread * num / den
}

rl_check <- imap(fits, \(fit, nm) {
    pr <- filter(model_pars, model == nm)
    loc <- fit$summary.fitted.values$mean[idat][dc_max$year == last_year]
    tibble(loc = loc,
        z2_old = map_dbl(loc, \(l)
            return_level_bgev(2, l, pr$spread, pr$tail,
                p_a = p_a, p_b = p_b)),
        z2_new = map_dbl(loc, \(l)
            return_level_bgev2(2, l, pr$spread, pr$tail,
                p_a = p_a, p_b = p_b)),
        z50_new = map_dbl(loc, \(l)
            return_level_bgev2(50, l, pr$spread, pr$tail,
                p_a = p_a, p_b = p_b)),
        z50_ref = map_dbl(loc, \(l)
            rl_bgev_ref(50, l, pr$spread, pr$tail))) |>
        summarise(
            bug_z2_vs_loc = max(abs(z2_old - loc)),
            fix_z2_vs_loc = max(abs(z2_new - loc)),
            fix_vs_ref_50 = max(abs(z50_new - z50_ref))) |>
        mutate(model = nm, .before = 1)
}) |>
    list_rbind()
print(rl_check)
