
##---- make_population ----
make_population <- function(
  npop = 1e3,
  # pop exposure (within window) to vaccine and dose
  pe = setNames(c(0.08, 0.08, 0.01, 0.01),
                c("V1_D1", "V1_D2", "v2_D1", "V2_D2")),
  # base prevalence (arbitrarily relative to pop size)
  p0 = 1e-2,
  # risk of vaccine 1
  risk_v1 = 1.5,
  # added risk of dose 2 for both vaccines
  risk_d2 = 1,
  # RR of v2 relative to v1
  risk_v2_v1 = 1,
  outcome = "case", unexposed = "UNEXP"
){
  expos <- c(unexposed, names(pe))
  ## exposure
  expo <- as.factor( sample(expos , npop,
                            prob = c(1 - sum(pe), pe),
                            replace = TRUE) )
  pop <- data.frame(expo)
  ## risks
  ps <- setNames(c(p0,
                   p0 * risk_v1,
                   p0 * risk_v1 * risk_d2,
                   p0 * risk_v1 * risk_v2_v1,
                   p0 * risk_v1 * risk_v2_v1 * risk_d2), expos)
  ## case or not
  pop[, outcome] <- vector("integer", npop)
  for (e in expos){
    pop[pop$expo == e, outcome] <- rbinom(sum(pop$expo == e), 1, ps[e])
  }
  pop
}


####---- make_cc_study ----
make_cc_study <- function(
  pop,
  ncase,
  ra = 1 # ctl / case ratio
){
  sample_case <- pop[pop$case == 1,][sample(c(1:sum(pop$case == 1)),
                                            ncase, replace = FALSE),]
  sample_ctl <- pop[pop$case == 0,][sample(c(1:sum(pop$case == 0)),
                                           ra * ncase, replace = FALSE),]
  df <- data.frame(rbind(sample_case, sample_ctl),
                   set = c(1:ncase, rep(1:ncase, each = ra)))
}


##---- ar_delta_greenland ----
#' @param RR relative risk
#' @param RR_lo lower limit of RR
#' @param RR_up upper limit of RR
#' @param n1 number exposed cases
#' @param n0 number unexposed cases
#' @references
#' Greenland S. Variance estimators for attributable fraction
#' estimates consistent in both large strata and sparse data.
#' Statistics in Medicine 1987;6(6):701–8.
#'
#' Steenland K, Armstrong B. An Overview of Methods for Calculating
#' the Burden of Disease Due to Specific Risk Factors: Epidemiology.
#' 2006 Sep;17(5):512–9.
#'
ar_delta_greenland <- function(RR, RR_lo, RR_up, n1, n0, type = "Steenland"){
  n <- n1 + n0
  var_log_RR <- ( (log(RR_up) - log(RR))/qnorm(0.975) )^2
  # exp(log(RR) + qnorm(0.975)*sqrt(var_log_RR))
  AF <- n1/n * (RR-1)/RR ## AF_p in Greenland
  ## Greenland
  var_log_AF_p <- var_log_RR / (RR - 1)^2 +
    n0 / (n1 * n) +
    2/(n1 *  (RR - 1))
  ## logit -based
  x1 <- log(AF / (1 - AF)) + qnorm(0.025) * sqrt(var_log_AF_p) / (1 - AF)
  x2 <- log(AF / (1 - AF)) + qnorm(0.975) * sqrt(var_log_AF_p) / (1 - AF)
  ci_logitbased_greenland <- 1/(1 + exp(-c(x1, x2)))
  ## direct
  ci_direct_greenland <- c(
    AF + qnorm(0.025) * AF * sqrt(var_log_AF_p), # var(log(x)) = (1/x)^2 var(x)
    AF + qnorm(0.975) * AF * sqrt(var_log_AF_p))
  ## Steenland, p. 515 Formula (7)
  var_log_1minusAF <- (AF^2 / (1 - AF)^2)*(
    (var_log_RR / (RR - 1)^2) +
      2 / (n1 * (RR - 1)) +
      n0 / (n1 * n) )
  LL <- log(1 - AF) + qnorm(0.025)*sqrt(var_log_1minusAF)
  UL <- log(1 - AF) + qnorm(0.975)*sqrt(var_log_1minusAF)
  LLAF <- 1 - exp(UL)
  ULAF <- 1 - exp(LL)
  ci_steenland <- c(AF, LLAF, ULAF)
  if ( type == "Steenland") {
    ci_steenland
  } else {
    list(ci_steenland = ci_steenland,
         ci_logitbased_greenland = c(AF, ci_logitbased_greenland),
         ci_direct_greenland = c(AF, ci_direct_greenland),
         AF = AF,
         var_log_AF_p = var_log_AF_p)
  }
}

