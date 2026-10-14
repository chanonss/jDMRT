# Core computation (no JASP or jamovi dependency, so it can be unit-tested in plain R).
#
# Report formatting -------------------------------------------------------------
# Unicode superscript for group letters (there is no superscript q; uppercase is left as is).
.duncanSuper <- function(x) {
  glyph <- c(a = "\u1d43", b = "\u1d47", c = "\u1d9c", d = "\u1d48", e = "\u1d49", f = "\u1da0",
             g = "\u1d4d", h = "\u02b0", i = "\u2071", j = "\u02b2", k = "\u1d4f", l = "\u02e1",
             m = "\u1d50", n = "\u207f", o = "\u1d52", p = "\u1d56", r = "\u02b3", s = "\u02e2",
             t = "\u1d57", u = "\u1d58", v = "\u1d5b", w = "\u02b7", x = "\u02e3", y = "\u02b8",
             z = "\u1dbb")
  vapply(as.character(x), function(s) {
    ch <- strsplit(s, "", fixed = TRUE)[[1]]
    paste0(ifelse(ch %in% names(glyph), glyph[ch], ch), collapse = "")
  }, "", USE.NAMES = FALSE)
}

# "mean \u00b1 SD" + superscript letter. digits = "auto" (2 decimals, 3 when |value| < 1) or an integer.
.duncanReport <- function(mean, sd, n, letter = "", digits = "auto") {
  fmt <- function(v) {
    dg <- if (identical(digits, "auto")) ifelse(abs(v) < 1, 3L, 2L) else as.integer(digits)
    mapply(function(x, k) formatC(x, format = "f", digits = k, big.mark = ","), v, dg)
  }
  one <- ifelse(n < 2 | !is.finite(sd), paste0(fmt(mean), " (n = ", n, ")"),
                paste0(fmt(mean), " \u00b1 ", fmt(sd)))
  paste0(one, .duncanSuper(ifelse(is.na(letter), "", letter)))
}

# Duncan-Kramer letters --------------------------------------------------------
#   est    estimated means; vd  matrix of Var(mean_i - mean_j); dfe  error df.
#   Returns the descending order of the means and the letter string for each (in that order).
.duncanLetters <- function(est, vd, dfe, alpha) {
  k <- length(est)
  ord <- order(est, decreasing = TRUE)
  est_s <- est[ord]
  vd_s  <- vd[ord, ord, drop = FALSE]

  # Pair (a, b) with a < b in sorted order spans r = b - a + 1 means.
  # Significant when |diff| > q((1-alpha)^(r-1), r, dfe) * sqrt(Var(diff) / 2)
  nonsig <- matrix(TRUE, k, k)
  if (k > 1) for (a in seq_len(k - 1)) for (b in (a + 1):k) {
    r <- b - a + 1
    crit <- stats::qtukey((1 - alpha)^(r - 1), r, dfe) * sqrt(vd_s[a, b] / 2)
    nonsig[a, b] <- nonsig[b, a] <- (est_s[a] - est_s[b]) <= crit
  }

  # Maximal runs of consecutive means with no significant pair -> letters
  last <- integer(k)
  for (i in seq_len(k)) {
    j <- i
    while (j < k && all(nonsig[i:(j + 1), i:(j + 1)])) j <- j + 1
    last[i] <- j
  }
  mat <- matrix(FALSE, k, 0)
  prev_end <- 0
  for (i in seq_len(k)) if (last[i] > prev_end) {
    col <- rep(FALSE, k); col[i:last[i]] <- TRUE
    mat <- cbind(mat, col); prev_end <- last[i]
  }
  lab <- vapply(seq_len(k), function(i) paste(c(letters, LETTERS)[which(mat[i, ])], collapse = ""), "")
  list(ord = ord, letter = lab)
}

# Least-squares means of the groups in `key` (one entry per row of `grid`, which lists every
# design cell the means are averaged over), with the variance matrix of their differences.
.duncanLsm <- function(fit, grid, key, mse) {
  ok  <- !is.na(stats::coef(fit))
  X   <- stats::model.matrix(stats::delete.response(stats::terms(fit)), grid,
                             contrasts.arg = fit$contrasts, xlev = fit$xlevels)[, ok, drop = FALSE]
  lev <- levels(key)
  L   <- t(vapply(lev, function(l) colMeans(X[key == l, , drop = FALSE]), numeric(sum(ok))))
  V   <- summary(fit)$cov.unscaled[ok, ok, drop = FALSE] * mse
  list(est = as.numeric(L %*% stats::coef(fit)[ok]),
       vd  = outer(seq_along(lev), seq_along(lev), Vectorize(function(i, j) {
         cvec <- L[i, ] - L[j, ]
         as.numeric(t(cvec) %*% V %*% cvec)
       })))
}

# ANOVA (any factorial model, optional additive block) + Duncan post hoc -----------
#
#   y         numeric response
#   factors   named list of factor vectors (names = the labels shown to the user)
#   block     optional block factor, entered additively (RCBD)
#   terms     model terms, a list of character vectors of factor names, e.g.
#             list("A", "B", c("A", "B")); NULL = full factorial of `factors`
#   posthoc   terms to run Duncan's test on (subset of `terms`)
#   describe  terms to give N / raw mean / SD for (should include `posthoc`; default = posthoc)
#   adjusted  TRUE : least-squares means and Duncan-Kramer ranges (correct for unbalanced data)
#             FALSE: raw means and the harmonic mean of the group sizes (agricolae method)
#   Both are identical for balanced designs.
#
# Sums of squares are Type III (sum-to-zero contrasts), equal to the sequential values when the
# design is balanced. Every comparison uses the error term of the full model.

.duncanTermLabel <- function(term, fn) paste(fn[sort(match(term, fn))], collapse = " \u00d7 ")

# all main effects and interactions of the factors `fn`, lowest order first
.duncanFullFactorial <- function(fn) {
  out <- list()
  for (m in seq_along(fn))
    out <- c(out, if (m == length(fn)) list(fn) else utils::combn(fn, m, simplify = FALSE))
  out
}

.duncanAnova <- function(y, factors, block = NULL, terms = NULL, posthoc = list(),
                         describe = posthoc, alpha = 0.05, adjusted = TRUE) {
  fn <- names(factors)
  k  <- length(fn)
  d  <- data.frame(y = y)
  for (i in seq_len(k)) d[[paste0("f", i)]] <- factors[[i]]
  if (!is.null(block)) d$blk <- block
  d <- d[stats::complete.cases(d), , drop = FALSE]
  for (col in setdiff(names(d), "y")) d[[col]] <- droplevels(as.factor(d[[col]]))

  if (is.null(terms) || !length(terms)) {
    terms <- .duncanFullFactorial(fn)
  }
  # sort the factors inside each term by their position, drop duplicate terms
  terms <- lapply(terms, function(t) fn[sort(match(t, fn))])
  terms <- terms[!duplicated(vapply(terms, paste, "", collapse = "\r"))]
  terms <- terms[order(lengths(terms))]
  internal <- function(nm) paste0("f", match(nm, fn))
  tint  <- lapply(terms, internal)
  tlab  <- vapply(tint, paste, "", collapse = ":")
  mvars <- unique(unlist(tint))
  if (!length(mvars)) stop("Choose at least one model term.", call. = FALSE)
  for (v in c(mvars, if (!is.null(block)) "blk"))
    if (nlevels(d[[v]]) < 2)
      stop("'", if (v == "blk") "block" else fn[match(v, paste0("f", seq_len(k)))],
           "' needs at least 2 levels after removing missing values.", call. = FALSE)

  disp <- function(v) fn[match(v, paste0("f", seq_len(k)))]
  for (t in tint) if (length(t) > 1) {
    cnt <- table(d[t])
    if (any(cnt == 0)) {
      e <- which(cnt == 0, arr.ind = TRUE)[1, ]
      stop("No observations in the cell ", paste(vapply(seq_along(t), function(j)
             dimnames(cnt)[[j]][e[j]], ""), collapse = " \u00d7 "), " of the term ",
           paste(disp(t), collapse = " \u00d7 "), "; the interaction cannot be estimated.", call. = FALSE)
    }
  }

  contr <- stats::setNames(rep(list(stats::contr.sum), length(mvars) + !is.null(block)),
                           c(mvars, if (!is.null(block)) "blk"))
  rhs <- c(if (!is.null(block)) "blk", tlab)
  fit <- stats::lm(stats::reformulate(rhs, response = "y"), d, contrasts = contr)
  dfe <- fit$df.residual
  if (dfe < 1)
    stop("Not enough data to estimate the error term; the model needs replication ",
         "(for an RCBD, remove the block interactions from the model terms).", call. = FALSE)
  rss <- sum(stats::residuals(fit)^2)
  mse <- rss / dfe

  # Type III table
  dr  <- suppressWarnings(stats::drop1(fit, scope = rhs, test = "F"))
  ss  <- dr[rhs, "Sum of Sq"]; df <- dr[rhs, "Df"]
  ms  <- ss / df; f <- ms / mse
  p   <- stats::pf(f, df, dfe, lower.tail = FALSE)
  labs <- c(if (!is.null(block)) "(block)", vapply(terms, paste, "", collapse = " \u00d7 "))
  anova_tab <- data.frame(term = c(labs, "Residuals"), df = c(df, dfe), ss = c(ss, rss),
                          ms = c(ms, mse), f = c(f, NA_real_), p = c(p, NA_real_),
                          stringsAsFactors = FALSE)
  term_p <- stats::setNames(p, rhs)

  balanced <- {
    cnt <- table(d[c(mvars, if (!is.null(block)) "blk")])
    length(unique(as.vector(cnt))) == 1
  }

  # Group tables -----------------------------------------------------------------
  sep <- " \u00d7 "
  # every design cell the least-squares means are averaged over
  grid <- NULL
  gridOf <- function() {
    if (is.null(grid)) {
      vs <- c(mvars, if (!is.null(block)) "blk")
      grid <<- expand.grid(stats::setNames(lapply(vs, function(v) levels(d[[v]])), vs),
                           stringsAsFactors = TRUE)
    }
    grid
  }
  one <- function(term, doPost) {
    fac <- internal(term)
    key <- interaction(d[fac], sep = sep, lex.order = TRUE)
    lev <- levels(key)
    n   <- as.numeric(table(key)[lev])
    raw <- as.numeric(tapply(d$y, key, mean)[lev])
    sdv <- as.numeric(tapply(d$y, key, stats::sd)[lev])
    lvl <- lapply(fac, function(v) levels(d[[v]]))
    ai  <- arrayInd(seq_along(lev), rev(lengths(lvl)))
    ai  <- ai[, rev(seq_along(fac)), drop = FALSE]
    cols <- lapply(seq_along(fac), function(j) lvl[[j]][ai[, j]])
    names(cols) <- paste0("f", seq_along(fac))
    out <- data.frame(cols, level = lev, n = as.integer(n), mean = raw, rawMean = raw, sd = sdv,
                      pos = seq_along(lev), letter = "", stringsAsFactors = FALSE)
    if (doPost) {
      if (adjusted) {
        g <- gridOf()
        lsm <- .duncanLsm(fit, g, interaction(g[fac], sep = sep, lex.order = TRUE), mse)
        est <- lsm$est; vd <- lsm$vd
      } else {
        est <- raw
        vd  <- matrix(2 * mse / (1 / mean(1 / n)), length(lev), length(lev))
      }
      let <- .duncanLetters(est, vd, dfe, alpha)
      out$mean <- est
      out <- out[let$ord, , drop = FALSE]
      out$letter <- let$letter
      rownames(out) <- NULL
    }
    out
  }
  lab_of <- function(t) .duncanTermLabel(t, fn)
  plabs <- vapply(posthoc, lab_of, "")
  dlabs <- vapply(describe, lab_of, "")
  groups <- list()
  for (i in seq_along(describe)) {
    t <- describe[[i]]
    if (!(dlabs[i] %in% vapply(terms, paste, "", collapse = " \u00d7 ")))
      stop("The term ", dlabs[i], " is not in the model.", call. = FALSE)
    if (!is.null(groups[[dlabs[i]]])) next
    groups[[dlabs[i]]] <- list(factors = fn[sort(match(t, fn))],
                               p = unname(term_p[paste(internal(fn[sort(match(t, fn))]), collapse = ":")]),
                               posthoc = dlabs[i] %in% plabs)
    groups[[dlabs[i]]]$groups <- one(t, dlabs[i] %in% plabs)
  }

  list(anova = anova_tab, dfe = dfe, mse = mse, alpha = alpha, adjusted = adjusted,
       unbalanced = !balanced, terms = vapply(terms, paste, "", collapse = " \u00d7 "),
       term_p = stats::setNames(p, vapply(terms, paste, "", collapse = " \u00d7 ")),
       groups = groups)
}
