# Run from the module directory: Rscript tests/test_duncan.R   (needs agricolae, car, emmeans)
source("R/duncanCompute.R")
set.seed(42)

cp <- function(x) paste(utf8ToInt(enc2utf8(x)), collapse = " ")
same_text <- function(x, y) identical(cp(x), cp(y))
agri <- function(y, key, dfe, mse) {
  r <- agricolae::duncan.test(y, key, dfe, mse, group = TRUE, console = FALSE)$groups
  setNames(as.character(r$groups), rownames(r))
}
letters_ok <- function(tab, y, key, res, what) {
  ref  <- agri(y, key, res$dfe, res$mse)
  mine <- setNames(tab$letter, tab$level)
  stopifnot(identical(unname(mine[names(ref)]), unname(ref)))
  cat("OK  ", what, "\n")
}
run <- function(y, f, block = NULL, terms = NULL, post, adj) {
  .duncanAnova(y, f, block, terms, posthoc = post, adjusted = adj)
}

# 1. one-way balanced / unbalanced: letters equal agricolae::duncan.test (aov y ~ g)
g <- gl(5, 6, labels = LETTERS[1:5]); y <- rnorm(30, rep(c(10, 12, 12.5, 15, 15.4), each = 6))
for (adj in c(TRUE, FALSE)) {
  r <- run(y, list(G = g), post = list("G"), adj = adj)
  letters_ok(r$groups$G$groups, y, g, r, paste("one-way balanced, adjusted =", adj))
}
g2 <- factor(sample(letters[1:4], 40, TRUE)); y2 <- rnorm(40, as.integer(g2) * 1.2)
for (adj in c(TRUE, FALSE)) {
  r <- run(y2, list(G = g2), post = list("G"), adj = adj)
  # agricolae's group tables for unbalanced one-way use the harmonic mean = raw mode
  if (!adj) letters_ok(r$groups$G$groups, y2, g2, r, "one-way unbalanced, raw")
  fit <- lm(y2 ~ g2)
  stopifnot(isTRUE(all.equal(r$dfe, fit$df.residual)),
            isTRUE(all.equal(r$mse, sum(resid(fit)^2) / fit$df.residual)))
}
cat("OK   one-way unbalanced error terms\n")

# 2. RCBD, additive model
trt <- gl(5, 1, 30, labels = LETTERS[1:5])[rep(1:5, 6)]; blk <- gl(6, 5, 30)
yb <- rnorm(30, rep(c(10, 12, 12.5, 15, 15.4), 6)) + rep(rnorm(6, 0, 2), each = 5)
for (adj in c(TRUE, FALSE)) {
  r <- run(yb, list(T = trt), blk, post = list("T"), adj = adj)
  ref <- agricolae::duncan.test(aov(yb ~ blk + trt), "trt", group = TRUE, console = FALSE)$groups
  mine <- setNames(r$groups$T$groups$letter, r$groups$T$groups$level)
  stopifnot(identical(unname(mine[rownames(ref)]), as.character(ref$groups)))
  cat("OK   RCBD balanced vs agricolae, adjusted =", adj, "\n")
}
keep <- setdiff(seq_along(yb), c(3, 9, 14, 22, 27))
yu <- yb[keep]; tu <- trt[keep]; bu <- blk[keep]
r <- run(yu, list(T = tu), bu, post = list("T"), adj = TRUE)
fit <- lm(yu ~ bu + tu, contrasts = list(bu = contr.sum, tu = contr.sum))
em <- as.data.frame(emmeans::emmeans(fit, "tu"))
got <- setNames(r$groups$T$groups$mean, r$groups$T$groups$level)[as.character(em$tu)]
stopifnot(isTRUE(all.equal(unname(got), em$emmean)), r$unbalanced)
rr <- run(yu, list(T = tu), bu, post = list("T"), adj = FALSE)
ref <- agricolae::duncan.test(aov(yu ~ bu + tu), "tu", group = TRUE, console = FALSE)$groups
stopifnot(identical(unname(setNames(rr$groups$T$groups$letter, rr$groups$T$groups$level)[rownames(ref)]),
                    as.character(ref$groups)))
cat("OK   RCBD unbalanced: adjusted = emmeans, raw = agricolae\n")

# 3. two-way factorial in an RCBD
A <- gl(3, 1, 72, labels = c("a1", "a2", "a3"))[rep(1:3, 24)]
B <- gl(4, 1, 72, labels = c("b1", "b2", "b3", "b4"))[rep(rep(1:4, each = 3), 6)]
K <- gl(6, 12, 72)
y <- rnorm(72, 10 + as.integer(A) * 1.5 + as.integer(B) * 0.8 + (A == "a3" & B == "b4") * 3) +
  rep(rnorm(6, 0, 1.5), each = 12)
for (adj in c(TRUE, FALSE)) {
  r <- run(y, list(A = A, B = B), K, post = list("A", "B", c("A", "B")), adj = adj)
  letters_ok(r$groups$A$groups, y, A, r, paste("two-way RCBD, A, adjusted =", adj))
  letters_ok(r$groups$B$groups, y, B, r, paste("two-way RCBD, B, adjusted =", adj))
  cg <- r$groups[["A \u00d7 B"]]$groups
  cg$level <- paste(cg$f1, cg$f2, sep = "|")   # ASCII keys: agricolae may re-encode the labels
  letters_ok(cg, y, interaction(A, B, sep = "|", lex.order = TRUE), r,
             paste("two-way RCBD, cells, adjusted =", adj))
}
old <- options(contrasts = c("contr.sum", "contr.poly")); sq <- summary(aov(y ~ K + A * B))[[1]]; options(old)
r <- run(y, list(A = A, B = B), K, post = list(), adj = TRUE)
stopifnot(isTRUE(all.equal(r$anova$ss, sq[["Sum Sq"]], check.attributes = FALSE)))
cat("OK   Type III = sequential SS when balanced\n")

# 4. two-way unbalanced: Type III vs car, adjusted means vs emmeans
kp <- setdiff(seq_along(y), c(2, 9, 17, 30, 31, 44, 58, 66))
yu <- y[kp]; Au <- A[kp]; Bu <- B[kp]; Ku <- K[kp]
r <- run(yu, list(A = Au, B = Bu), Ku, post = list("A", c("A", "B")), adj = TRUE)
stopifnot(r$unbalanced)
fit <- lm(yu ~ Ku + Au * Bu, contrasts = list(Ku = contr.sum, Au = contr.sum, Bu = contr.sum))
ct <- as.data.frame(car::Anova(fit, type = 3))
stopifnot(isTRUE(all.equal(r$anova$ss[r$anova$term == "A"], ct["Au", "Sum Sq"])),
          isTRUE(all.equal(r$anova$ss[r$anova$term == "B"], ct["Bu", "Sum Sq"])),
          isTRUE(all.equal(r$anova$ss[r$anova$term == "A \u00d7 B"], ct["Au:Bu", "Sum Sq"])))
cat("OK   Type III SS match car::Anova on unbalanced data\n")
emc <- as.data.frame(emmeans::emmeans(fit, ~ Au * Bu))
cg <- r$groups[["A \u00d7 B"]]$groups
got <- setNames(cg$mean, paste(cg$f1, cg$f2, sep = "|"))[paste(emc$Au, emc$Bu, sep = "|")]
stopifnot(isTRUE(all.equal(unname(got), emc$emmean)))
ema <- as.data.frame(emmeans::emmeans(fit, "Au"))
ga <- r$groups$A$groups
stopifnot(isTRUE(all.equal(unname(setNames(ga$mean, ga$level)[as.character(ema$Au)]), ema$emmean)))
cat("OK   adjusted means match emmeans (A, A x B)\n")

# 5. three factors, user-chosen terms, describe-only terms, errors
C <- factor(sample(c("c1", "c2"), 72, TRUE))
r3 <- .duncanAnova(y, list(A = A, B = B, C = C), terms = list("A", "B", "C", c("A", "B")),
                   posthoc = list("A"), describe = list("A", "C"))
stopifnot(identical(names(r3$groups), c("A", "C")), all(r3$groups$C$groups$letter == ""),
          r3$groups$A$posthoc, !r3$groups$C$posthoc, nrow(r3$anova) == 5)
bad <- !(A == "a1" & B == "b1")
msg <- tryCatch(.duncanAnova(y[bad], list(A = A[bad], B = B[bad])), error = function(e) conditionMessage(e))
stopifnot(grepl("No observations", msg))
msg <- tryCatch(.duncanAnova(y, list(A = A, B = B), K, terms = list("A", "B", c("A", "B"), "K")),
                error = function(e) conditionMessage(e))
cat("OK   three factors / describe-only terms / empty cell refused\n")

# 6. report formatting (built from code points, so the check does not depend on the locale)
u <- function(...) intToUtf8(c(...))
stopifnot(same_text(.duncanSuper("abc"), u(0x1d43, 0x1d47, 0x1d9c)),
          same_text(.duncanSuper("ns"), u(0x207f, 0x02e2)),
          same_text(.duncanReport(66.754, 1.5234, 4, "a"), paste0("66.75 ", u(0xb1), " 1.52", u(0x1d43))),
          same_text(.duncanReport(0.5, 0.25, 4, ""), paste0("0.500 ", u(0xb1), " 0.250")),
          same_text(.duncanReport(0.5, 0.25, 4, "", 1), paste0("0.5 ", u(0xb1), " 0.2")),
          same_text(.duncanReport(3, NA, 1, "b", 1), paste0("3.0 (n = 1)", u(0x1d47))))
cat("OK   report formatting\nAll checks passed\n")
