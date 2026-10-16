# jDMRT — ANOVA with Duncan's multiple range test for jamovi

Adds **DMRT → ANOVA with Duncan Post Hoc** to jamovi: ANOVA (one-way, RCBD, factorial; balanced or
unbalanced) followed by Duncan's multiple range test (DMRT) group letters, as used in agricultural
and food/sensory research.

## How to use
1. **Dependent Variable**, **Fixed Factors** (one or more), optional **Block** (additive, for RCBD).
2. **Model** (collapsed): terms default to the full factorial; edit them as in jamovi's ANOVA.
   For an RCBD put the block in the Block box so it is not crossed with the factors.
3. **Post Hoc Tests**: move terms to the box and tick **Duncan's multiple range test**.
   Nothing post hoc runs until you do.
4. **Report**: *Report table* gives `Mean ± S.D.` with superscript letters (raw means);
   optional interaction plot for the first two-factor term.

Means sharing a letter are not significantly different at the chosen α.
If the term's ANOVA F-test is not significant, letters are suppressed (`ns`) — this can be switched off.

## Example data
Two simulated datasets ship with the module (jamovi: Open → Data Library → *jDMRT*):
- **Rice RCBD** — 5 varieties × 6 blocks with two plots missing (unbalanced). Dependent *Yield*,
  fixed factor *Variety*, block *Block*, tick Duncan for *Variety*.
- **Fertilizer × Variety factorial** — 3 × 2 in 4 blocks, with a significant interaction.
  Dependent *Yield*, fixed factors *Fertilizer* and *Variety*, block *Block*.

## Methods
- Type III sums of squares for the ANOVA table; letters use the error term of the full model.
- *Means (unbalanced designs)*:
  - **Least-squares (recommended)**: block-adjusted means with Duncan-Kramer pairwise ranges
    (Kramer, 1956); correct when plots are missing.
  - **Raw means (agricolae)**: raw means and the harmonic mean of replications, identical to
    `agricolae::duncan.test`.
  For balanced designs both give the same result.
- Compared against `agricolae::duncan.test`, `car` Type III and `emmeans` in `tests/test_duncan.R`.

## Limitations
- One dependent variable per run.
- Duncan's test controls the family-wise error rate weakly; for strict control use Tukey HSD.
- Duncan-Kramer is an approximation for unequal replication.

## References
Duncan, D. B. (1955). Multiple range and multiple F tests. *Biometrics, 11*(1), 1-42.
Kramer, C. Y. (1956). Extension of multiple range tests to group means with unequal numbers of
replications. *Biometrics, 12*(3), 307-310.

## Author and licence
Chanon Sarasuk, Faculty of Home Economics Technology, Rajamangala University of Technology Krungthep.
GPL (>= 2), see `LICENSE`.

## For developers
### Build (needs R, the `jmvtools` package and jamovi installed)
    install.packages("jmvtools", repos = c("https://repo.jamovi.org", "https://cran.r-project.org"))
    jmvtools::check()                 # locate jamovi
    jmvtools::install("jDMRT")      # builds and installs into jamovi

### Shared code
`R/duncanCompute.R` is a copy of `../../shared/duncanCompute.R`. Edit it there and run `../../scripts/sync.sh`.

### Tests
`tests/test_duncan.R` (plain R, needs agricolae, car, emmeans) checks the statistics;
`tests/robustness_test.R` runs the installed module through missing data, odd variable names,
empty/constant data etc. (see its header for how to run it with jamovi's R).

### Test (statistics)
    Rscript tests/test_duncan.R
