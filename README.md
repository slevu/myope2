---
title: "myope2"
output: 
  html_document:
    keep_md: yes
---



Code and simulated data for the article **Age and sex-specific risks of myocarditis and pericarditis following Covid-19 messenger RNA Vaccines** by Le Vu et al., available at https://www.nature.com/articles/s41467-022-31401-5


```r
source("functions.R")
```

# Population
- We simulate a population of size `npop` with a base prevalence of the outcome `p0` and exposed to two vaccines, V1 and V2, according to the ranking of the first two doses, D1 or D2.
- A vector `pe` of four probabilities defines exposures levels to the combination of vaccines and doses. 
- Note that in our study, exposure refers to being within a window of 21 days from the last dose of vaccine when included in the case-control sample.
- Relative risks are parameterized like so :
    - `risk_v1` risk of outcome associated with exposure to V1 relative to being unexposed.
    - `risk_v2_v1` risk of outcome associated with exposure to V2 relative to v1.
    - `risk_d2` risk of outcome associated with exposure to the second dose relative to the first, applied for both vaccines.


```r
set.seed(1234)
n <-  32e6
pe <- setNames(c(0.08, 0.08, 0.01, 0.01), c("V1_D1", "V1_D2", "V2_D1", "V2_D2"))
pop0 <- make_population(npop = n, 
                        p0 = 1e-4,
                        pe = pe,
                        risk_v1 = 1.5,
                        risk_d2 = 2,
                        risk_v2_v1 = 2)
```

- Description of the population

```r
with(pop0, table(expo, case))
```

```
##        case
## expo           0        1
##   UNEXP 26237228     2640
##   V1_D1  2558275      392
##   V1_D2  2559883      801
##   V2_D1   320638       86
##   V2_D2   319891      166
```

# Case-control study
- We now sample `ncase` cases and `ra * ncase` (matched) controls.


```r
df_cc <- make_cc_study(pop0, ncase = 1500, ra = 10)
```

- Description of the study sample

```r
(t0 <- with(df_cc, table(expo, case)) )
```

```
##        case
## expo        0     1
##   UNEXP 12305   965
##   V1_D1  1197   139
##   V1_D2  1208   301
##   V2_D1   150    32
##   V2_D2   140    63
```

## Measure of association
- We estimate the odds ratios by conditional logistic regression

```r
require(survival)
```

```
## Le chargement a nécessité le package : survival
```

```r
mod <- clogit(case ~ expo + strata(set), data = df_cc)
r0 <- cbind( OR = round(exp(coef(mod)), 3),
               round(exp(confint(mod)), 3) )
rownames(r0) <- gsub("expo", "", rownames(r0))
res <- merge(unclass(t0), r0, by = 0, all.x = TRUE)
colnames(res) <- c("expo", "ctl", "case", "or", "lo", "up")
res
```

```
##    expo   ctl case    or    lo    up
## 1 UNEXP 12305  965    NA    NA    NA
## 2 V1_D1  1197  139 1.473 1.221 1.777
## 3 V1_D2  1208  301 3.178 2.754 3.667
## 4 V2_D1   150   32 2.718 1.844 4.006
## 5 V2_D2   140   63 5.645 4.157 7.666
```

## Vaccine-associated cases
- Given a vector `pd` of proportion of doses received, we then estimate
    - the population attributable fraction `paf`,
    - the excess of vaccine-associated cases `ec` per `nper` doses, and
    - the number of doses needed for the occurrence of one vaccine-associated case `nnh`.
- Here, we assume that the associations are causal, that the relative risks can reasonably be estimated by the odds ratios,  and that each vaccinee has experienced the (short) window of exposure, so that exposure is given by the number of respective doses received by the population.
- We derive the confidence intervals according to Greenland, 1987 [1].


```r
pd <- setNames( c(0.65, 0.59, 0.09, 0.08), names(pe) )
res <- cbind(res, nd = c(NA, pd * n) )
nper <- 1e5 # denominator
```



```r
u_case <- res[1, "case"] ## unexposed cases
ar_ic0 <- with(res[-1,], ar_delta_greenland(RR = or, RR_lo = lo,
                                              RR_up = up, n1 = case,
                                              n0 = u_case, type = 2))

paf_ic0 <- matrix(ar_ic0$ci_logitbased_greenland, ncol = 3, byrow = FALSE)
ec_ic0 <- matrix(
    with(res[-1,], ar_ic0$ci_logitbased_greenland * (case + u_case) / nd * nper),
    ncol = 3, byrow = FALSE)
nnh_ic0 <- (nper / ec_ic0)
nnh_ic0 <- nnh_ic0[, c(1,3,2)] ## reorder
a <- matrix(cbind(paf_ic0, ec_ic0, nnh_ic0), ncol = 9, byrow = FALSE,
               dimnames = list(res[-1, 1], c("paf", "paf_lo", "paf_up",
                                       "ec", "ec_lo", "ec_up",
                                       "nnh", "nnh_lo", "nnh_up")) )
signif(a, 4)
```

```
##           paf  paf_lo  paf_up     ec  ec_lo  ec_up    nnh nnh_lo nnh_up
## V1_D1 0.04043 0.02329 0.06929 0.2146 0.1236 0.3678 466000 271900 808900
## V1_D2 0.16290 0.13840 0.19080 1.0930 0.9284 1.2800  91520  78160 107700
## V2_D1 0.02029 0.01163 0.03517 0.7023 0.4025 1.2170 142400  82140 248400
## V2_D2 0.05043 0.03741 0.06765 2.0250 1.5020 2.7170  49380  36810  66560
```

# Reference
1. Greenland S. Variance estimators for attributable fraction estimates consistent in both large strata and sparse data. Statistics in Medicine. 1987;6(6):701–8. 
