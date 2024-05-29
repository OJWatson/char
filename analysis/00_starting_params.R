library(magenta)
library(tidyverse)

EIRs <- c(0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,
          15,20,30,50,100,200)
fts <- seq(0.1,0.9,0.1)

equilibrium_init_create <- function(par) {

  ## EIR
  EIRY_eq <- par$EIR  # initial annual EIR
  EIRd_eq <- EIR_eq <- EIRY_eq/365

  # FOI
  FOI_eq <- EIR_eq * par$b

  # FOI to T and D
  aT <- FOI_eq * par$phi * par$ft/par$rT_S
  aD <- FOI_eq * par$phi * (1 - par$ft)/par$rD

  Z_eq <- rep(NA, 3)
  Z_eq[1] <- 1/(1 + aT + aD)
  Z_eq[2] <- aT * Z_eq[1]
  Z_eq[3] <- aD * Z_eq[1]

  Y_eq <- Z_eq[1]
  T_eq <- Z_eq[2]
  D_eq <- Z_eq[3]

  betaS <- FOI_eq
  betaA <- FOI_eq * par$phi + par$rA

  A_eq <- (FOI_eq * (1 - par$phi) * Y_eq + par$rD * D_eq)/(betaA + FOI_eq * (1 - par$phi))
  S_eq <- Y_eq - A_eq
  FOIv_eq <- par$a * (par$cT*T_eq + par$cD*D_eq + par$cA*A_eq)

  # mosquito states
  Iv_eq <- FOIv_eq * exp(-par$mu * par$n)/(FOIv_eq + par$mu)
  Sv_eq <- par$mu * Iv_eq/(FOIv_eq * exp(-par$mu * par$n))
  Ev_eq <- 1 - Sv_eq - Iv_eq

  # mosquito density needed to give this EIR
  mv0 <- EIRd_eq/(Iv_eq * par$a)

  ## collate init
  list(
    EIR = par$EIR, ft = par$ft,
    S = S_eq, D = D_eq, A = A_eq, T = T_eq,
    phi = par$phi, b = par$b,
    m = mv0, Sv = Sv_eq, Ev = Ev_eq, Iv = Iv_eq, a = par$a,
    cA = par$cA, cD = par$cD, cT = par$cT,
    n = par$n,
    mu = par$mu,
    rD = par$rD,
    rA = 1/(1/par$rA +1/par$rU),
    rT_S = par$rT_S,
    rT_R = par$rT_R
  ) %>%
    as.data.frame()

}

phi_eir_rel <- function(eir, ft){

  mpl <- ICDMM::model_param_list_create(rho=0, rU = Inf, rP = Inf, sigma2 = 0)
  eq <- ICDMM::equilibrium_init_create(
    age_vector=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
    EIR=eir,ft=ft,
    model_param_list = mpl, het_brackets=2,
    country = NULL,
    admin_unit = NULL)

  phi <- weighted.mean(
    apply(eq$phi_eq, 1, weighted.mean, eq$het_wt),
    eq$den
  )

  c_A <- weighted.mean(
    apply(eq$cA_eq, 1, weighted.mean, eq$het_wt),
    rowMeans(eq$init_A)
  )


  b <- weighted.mean(rowMeans(eq$b0 * ((1 - eq$b1)/(1 + (eq$init_IB/eq$IB0)^eq$kB) + eq$b1)), eq$den)

  S <- sum(eq$init_S) + sum(eq$init_P)
  D <- sum(eq$init_D)
  A <- sum(eq$init_A + eq$init_U)
  T <- sum(eq$init_T)

  lambda_v_scale <- ((eq$av0 * (c_A*A + eq$cD*D + eq$cT*T))/eq$FOIv_eq)

  par <- list(
    EIR = eir, ft = ft,
    S = S, D = D, A = A, T = T, phi = phi, b = b,
    m = eq$mv0, Sv = eq$init_Sv, Ev = eq$init_Ev, Iv = eq$init_Iv, a = eq$av0,
    cA = c_A, cD = eq$cD, cT = eq$cT,
    n = eq$delayMos,
    mu = eq$mu0,
    rD = eq$rD,
    rU = eq$rU,
    rA = eq$rA,
    rT_S = eq$rT,
    rT_R = eq$rT,
    lambda_v_scale = lambda_v_scale
  ) %>%
    as.data.frame()

  equilibrium_init_create(par)

}

pars <- expand.grid("eir" = EIRs, "ft" = fts)
pars$n <- seq_along(pars$eir)

starting_params <- lapply(split(pars, pars$n),
                          function(x){
                            phi_eir_rel(x$eir, x$ft)
                          }) %>% do.call(rbind, .)

saveRDS(starting_params, "analysis/data/derived/starting_params.rds")
write.csv(starting_params, "analysis/data/derived/starting_params.csv")

starting_params <- readRDS("analysis/data/derived/starting_params.rds")

# Example plots and relationships
starting_params %>% ggplot(aes(EIR, D+A+T, color = as.factor(ft), group = as.factor(ft))) + geom_line() + theme_bw() + ylim(c(0,1)) + scale_color_discrete("ft")
starting_params %>% ggplot(aes(EIR, D+A+T, color = as.factor(ft), group = as.factor(ft))) + geom_line() + theme_bw() + ylim(c(0,1)) + scale_color_discrete("ft") + scale_x_log10()
