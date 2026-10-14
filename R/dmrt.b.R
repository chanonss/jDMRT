# jamovi wrapper around .duncanAnova() (see R/duncanCompute.R)
dmrtClass <- R6::R6Class(
  "dmrtClass",
  inherit = dmrtBase,
  private = list(

    # factors of a term in the order of the Fixed Factors list
    .allFactors = function() {
      opts <- self$options
      c(opts$factors, if (!is.null(opts$blk) && nzchar(opts$blk)) opts$blk)
    },
    .sortTerm = function(term) {
      fn <- private$.allFactors()
      fn[sort(match(term, fn))]
    },

    .init = function() {
      opts <- self$options
      if (!isTRUE(opts$postHocDuncan) || length(opts$factors) == 0) return()
      for (term in opts$postHoc) {
        if (!all(term %in% private$.allFactors())) next
        ft <- private$.sortTerm(term)
        tb <- self$results$postHocTables$get(key = term)
        tb$setTitle(paste0("Duncan's Multiple Range Test: ", paste(ft, collapse = " \u00d7 ")))
        for (i in seq_along(ft))
          tb$addColumn(name = paste0("f", i), title = ft[i], type = "text", index = i)
      }
    },

    .run = function() {
      opts <- self$options
      if (is.null(opts$dep) || length(opts$factors) == 0) return()
      hasBlock <- !is.null(opts$blk) && nzchar(opts$blk)
      if (hasBlock && opts$blk %in% opts$factors)
        jmvcore::reject("The block variable must differ from the fixed factors.", code = "")
      # the block is an ordinary factor of the model; it only enters additively unless the
      # user adds interaction terms with it
      fn <- private$.allFactors()

      d <- self$data
      facs <- lapply(fn, function(v) as.factor(d[[v]]))
      names(facs) <- fn
      inModelVars <- function(t) length(t) > 0 && all(t %in% fn)
      terms <- Filter(inModelVars, opts$modelTerms)
      if (!length(terms)) terms <- .duncanFullFactorial(opts$factors)
      if (hasBlock && !any(vapply(terms, function(t) identical(t, opts$blk), TRUE)))
        terms <- c(terms, list(opts$blk))
      blkTerm <- function(t) hasBlock && opts$blk %in% t
      lab <- function(t) .duncanTermLabel(t, fn)
      termLabs <- vapply(terms, lab, "")
      inModel <- function(t) inModelVars(t) && lab(t) %in% termLabs
      keyOf <- function(lst) if (length(lst)) vapply(lst, lab, "") else character()
      uniq  <- function(lst) lst[!duplicated(keyOf(lst))]

      ph <- if (opts$postHocDuncan) Filter(inModel, opts$postHoc) else list()
      repTerms <- if (opts$reportTable)
        (if (length(ph)) ph else Filter(function(t) !blkTerm(t), terms)) else list()
      plotTerm <- NULL
      if (opts$plot) {
        two <- Filter(function(t) length(t) == 2 && !blkTerm(t), c(ph, terms))
        if (length(two)) plotTerm <- two[[1]]
      }
      describe <- uniq(c(ph, repTerms, if (!is.null(plotTerm)) list(plotTerm)))

      res <- tryCatch(
        .duncanAnova(jmvcore::toNumeric(d[[opts$dep]]), facs,
                     NULL, terms = terms, posthoc = uniq(ph), describe = describe,
                     alpha = opts$alpha, adjusted = identical(opts$meanType, "adjusted")),
        error = function(e) jmvcore::reject(conditionMessage(e), code = ""))

      # ANOVA table ----------------------------------------------------------------
      tab <- self$results$anova
      a <- res$anova
      if (tab$rowCount == 0) for (i in seq_len(nrow(a))) {
        vals <- list(term = if (hasBlock && a$term[i] == opts$blk) paste0(opts$blk, " (block)") else a$term[i],
                     df = a$df[i], ss = a$ss[i], ms = a$ms[i])
        if (is.finite(a$f[i])) { vals$f <- a$f[i]; vals$p <- a$p[i] }
        tab$addRow(rowKey = i, values = vals)
      }
      if (res$unbalanced)
        tab$setNote("unb", "Unbalanced design: each term is adjusted for all the others (Type III).")

      # letters for a term: "ns" when its own F-test is not significant (if requested)
      lettersOf <- function(g) {
        gated <- opts$requireSignificantF && is.finite(g$p) && g$p > opts$alpha
        out <- g$groups
        if (g$posthoc && gated) out$letter <- "ns"
        list(groups = out, gated = gated)
      }
      tl <- lapply(terms, function(t) fn[sort(match(t, fn))])

      # Post hoc tables ------------------------------------------------------------
      if (opts$postHocDuncan) {
        info <- sprintf("alpha = %s, error df = %d, MSE = %s (full model).",
                        format(res$alpha), res$dfe, format(signif(res$mse, 5)))
        unb <- if (res$unbalanced) {
          if (res$adjusted) "Unbalanced design: least-squares means; Duncan-Kramer pairwise ranges."
          else "Unbalanced design: raw means; harmonic mean of group sizes (agricolae method)."
        }
        for (t in ph) {
          g  <- res$groups[[lab(t)]]
          lt <- lettersOf(g)
          tb <- self$results$postHocTables$get(key = t)
          tb$getColumn("mean")$setTitle(if (res$adjusted) "Adjusted mean" else "Mean")
          gr <- lt$groups
          nf <- length(g$factors)
          for (i in seq_len(nrow(gr))) {
            vals <- list(n = gr$n[i], sd = gr$sd[i], mean = gr$mean[i], letter = gr$letter[i])
            for (j in seq_len(nf)) vals[[paste0("f", j)]] <- gr[[paste0("f", j)]][i]
            tb$addRow(rowKey = i, values = vals)
          }
          if (lt$gated)
            tb$setNote("gate", "The F-test for this term is not significant; group letters are suppressed (ns).")
          else
            tb$setNote("share", "Means sharing a letter are not significantly different (Duncan's test).")
          # a significant higher-order term containing this one makes the comparison questionable
          S <- g$factors
          hi <- termLabs[vapply(tl, function(u) length(u) > length(S) && all(S %in% u), TRUE)]
          hi <- hi[is.finite(res$term_p[hi]) & res$term_p[hi] <= opts$alpha]
          if (length(hi))
            tb$setNote("int", paste0("Significant interaction (", paste(hi, collapse = ", "),
                                     "): these marginal comparisons may be misleading."))
          if (!is.null(unb)) tb$setNote("unb", unb)
          tb$setNote("info", info)
        }
      }

      # Report table ---------------------------------------------------------------
      if (opts$reportTable) {
        dg <- if (identical(opts$decimals, "auto")) "auto" else as.integer(opts$decimals)
        rt <- self$results$report
        i <- 0
        for (t in uniq(repTerms)) {
          g  <- res$groups[[lab(t)]]
          gr <- lettersOf(g)$groups
          gr <- gr[order(gr$pos), , drop = FALSE]
          txt <- .duncanReport(gr$rawMean, gr$sd, gr$n, gr$letter, dg)
          for (r in seq_len(nrow(gr))) {
            i <- i + 1
            rt$addRow(rowKey = i, values = list(part = lab(t), level = gr$level[r], n = gr$n[r],
                                                value = txt[r]))
          }
        }
        rt$setNote("raw", "Raw (unadjusted) group means \u00b1 standard deviation.")
        if (length(ph))
          rt$setNote("share", "Within a term, means sharing a letter are not significantly different (Duncan's test).")
        if (length(ph) && res$unbalanced && res$adjusted)
          rt$setNote("unb", "Unbalanced design: letters are based on least-squares means, the table shows raw means.")
      }

      # Plot state -----------------------------------------------------------------
      if (!is.null(plotTerm)) {
        g  <- res$groups[[lab(plotTerm)]]
        gr <- lettersOf(g)$groups
        self$results$plot$setState(list(
          cells = data.frame(a = gr$f1, b = gr$f2, mean = gr$rawMean,
                             se = gr$sd / sqrt(gr$n), letter = gr$letter, stringsAsFactors = FALSE),
          facA = g$factors[1], facB = g$factors[2], dep = opts$dep,
          levA = levels(droplevels(facs[[g$factors[1]]])),
          levB = levels(droplevels(facs[[g$factors[2]]]))))
      } else if (opts$plot) {
        self$results$plot$setState(NULL)
      }
    },

    # Means of the cells of a two-factor term; one line per level of B; bars = +/- SE of the raw
    # cell mean; letters are drawn when the term was compared.
    .plot = function(image, ...) {
      st <- image$state
      if (is.null(st)) return(FALSE)
      cl <- st$cells
      la <- st$levA; lb <- st$levB
      col <- stats::setNames(rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                                   "#56B4E9", "#F0E442", "#000000"), length.out = length(lb)), lb)
      pch <- stats::setNames(rep(c(16, 17, 15, 18, 1, 2, 0, 5), length.out = length(lb)), lb)
      se <- ifelse(is.finite(cl$se), cl$se, 0)
      yr <- range(c(cl$mean - se, cl$mean + se))
      pad <- 0.12 * diff(yr); if (!is.finite(pad) || pad == 0) pad <- 1
      graphics::plot(NA, xlim = c(0.6, length(la) + 0.4), ylim = c(yr[1] - pad, yr[2] + pad),
                     xaxt = "n", xlab = st$facA, ylab = paste("Mean of", st$dep), las = 1)
      graphics::axis(1, at = seq_along(la), labels = la)
      for (l in lb) {
        s <- cl[cl$b == l, , drop = FALSE]
        x <- match(s$a, la); o <- order(x); x <- x[o]; s <- s[o, , drop = FALSE]
        sev <- ifelse(is.finite(s$se), s$se, 0)
        graphics::lines(x, s$mean, col = col[[l]], lwd = 2)
        ok <- sev > 0
        if (any(ok))
          graphics::arrows(x[ok], s$mean[ok] - sev[ok], x[ok], s$mean[ok] + sev[ok],
                           angle = 90, code = 3, length = 0.04, col = col[[l]])
        graphics::points(x, s$mean, pch = pch[[l]], col = col[[l]], cex = 1.3)
        # letter above the bar for the highest line at each x, below it for the others
        top <- vapply(seq_len(nrow(s)), function(i) s$mean[i] >= max(cl$mean[cl$a == s$a[i]]), TRUE)
        lt <- !(s$letter %in% c("", "ns"))
        if (any(lt))
          graphics::text(x[lt], ifelse(top, s$mean + sev, s$mean - sev)[lt], labels = s$letter[lt],
                         pos = ifelse(top, 3, 1)[lt], offset = 0.5, cex = 0.85, col = col[[l]], font = 2)
      }
      graphics::legend("topleft", legend = lb, col = col[lb], pch = pch[lb], lty = 1, lwd = 2,
                       title = st$facB, bty = "n", cex = 0.85)
      TRUE
    })
)
