# Input-robustness checks that go through the generated jamovi analysis function dmrt().
# Needs the module installed in jamovi and jamovi's bundled R 4.5 (it loads jmvcore from jamovi):
#   R_HOME=<jamovi>/Frameworks/R.framework/Versions/4.5-arm64/Resources \
#   DYLD_FALLBACK_LIBRARY_PATH=$R_HOME/lib:/usr/lib  $R_HOME/bin/exec/R --vanilla --slave -f tests/robustness_test.R
.libPaths(c(path.expand("~/Library/Application Support/jamovi/modules/jDMRT/R"),
            "/Applications/jamovi.app/Contents/Resources/modules/base/R", .libPaths()))
suppressMessages({library(jmvcore); library(jDMRT)})

set.seed(7)
base <- data.frame(trt = gl(4, 6, labels = paste0("T", 1:4)), blk = factor(rep(paste0("B", 1:6), 4)),
                   y = rnorm(24, rep(c(5, 6, 8, 8.2), each = 6)))
fails <- 0
# expect = NULL  -> must run; otherwise a regular expression the (friendly) error message must match
check <- function(name, d, dep = "y", factors = "trt", blk = NULL, expect = NULL, ...) {
  args <- list(data = d, dep = dep, factors = factors, postHoc = as.list(factors), postHocDuncan = TRUE, ...)
  if (!is.null(blk)) args$blk <- blk
  msg <- tryCatch({ do.call(dmrt, args); NULL }, error = function(e) conditionMessage(e))
  ok <- if (is.null(expect)) is.null(msg) else (!is.null(msg) && grepl(expect, msg))
  cat(sprintf("%s  %-40s %s\n", if (ok) "OK  " else "FAIL", name, if (is.null(msg)) "" else msg))
  if (!ok) fails <<- fails + 1
}
mod <- function(f) { d <- base; f(d) }

check("baseline RCBD", base, blk = "blk")
check("NA in response",  mod(function(d) { d$y[c(3, 10)] <- NA; d }), blk = "blk")
check("NA in factor",    mod(function(d) { d$trt[c(2, 9)] <- NA; d }), blk = "blk")
check("NA in block",     mod(function(d) { d$blk[5] <- NA; d }), blk = "blk")
check("empty dataset",   base[0, ], expect = "Fewer than 2 complete observations")
check("all-NA response", mod(function(d) { d$y <- NA_real_; d }), expect = "Fewer than 2 complete observations")
check("constant response", mod(function(d) { d$y <- 5; d }), expect = "constant")
check("one-level factor", mod(function(d) { d$trt <- factor("T1"); d }), expect = "at least 2 levels")
check("one obs per group", base[!duplicated(base$trt), ], expect = "error term")
check("response as text", mod(function(d) { d$y <- as.character(d$y); d }), expect = "numeric")
check("Inf in response", mod(function(d) { d$y[1] <- Inf; d }), blk = "blk", expect = "Inf")
check("heavily unbalanced RCBD", base[base$trt != "T3" | base$blk %in% c("B1", "B2"), ], blk = "blk")
d <- base; names(d) <- c("สูตร", "ผู้บริโภค", "สี")
check("Thai variable names", d, dep = names(d)[3], factors = names(d)[1], blk = names(d)[2])
d <- base; names(d) <- c("treat (%)", "block no.", "y-value")
check("spaces/symbols in names", d, dep = "y-value", factors = "treat (%)", blk = "block no.")
d <- base; levels(d$trt) <- c("a b", "c/d", "e,f", "ก ข")
check("special characters in levels", d, blk = "blk")
check("integer-coded factor", mod(function(d) { d$trt <- as.integer(d$trt); d }), blk = "blk")
check("unused factor level", mod(function(d) { d$trt <- factor(d$trt, levels = c(levels(d$trt), "T5")); d }), blk = "blk")
check("block same as factor", base, blk = "trt", expect = "must differ")
check("alpha = 0.01", base, blk = "blk", alpha = 0.01)

# proper two-factor + block design (4 x 2, 3 blocks)
f <- expand.grid(A = paste0("a", 1:4), B = c("b1", "b2"), blk = paste0("R", 1:3)); f$y <- rnorm(nrow(f), 10)
check("two factors + block", f, factors = c("A", "B"), blk = "blk")
f2 <- f; f2$B <- f2$A    # second factor identical to the first: aliased terms
check("aliased factors refused", f2, factors = c("A", "B"), blk = "blk", expect = "confounded|empty|cell|error term")

cat(if (fails == 0) "All robustness checks passed\n" else sprintf("%d robustness check(s) FAILED\n", fails))
if (fails) quit(status = 1)
