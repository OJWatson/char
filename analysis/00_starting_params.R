library(magenta)
library(tidyverse)

EIRs <- c(0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,
          15,20,30,50,100,200)
fts <- seq(0.1,0.9,0.1)

phi_eir_rel <- function(eir, ft){

  age_vector <- magenta:::age_brackets(100, 40, TRUE)
  eq <- magenta::equilibrium_init_create(age_vector = age_vector,
                                         het_brackets = 5,
                                         ft = ft,
                                         EIR = eir,
                                         model_param_list = magenta::model_param_list_create())

  phi <- weighted.mean(
    apply(eq$phi_eq, 1, weighted.mean, eq$het_wt),
    eq$den
  )

  b <- weighted.mean(
    apply(eq$init_IB[,,1], 1, weighted.mean, eq$het_wt),
    eq$den
  )

  b <- eq$b0 * ((1 - eq$b1)/(1 + (b/eq$IB0)^eq$kB) + eq$b1)

  S <- sum(eq$init_S)
  D <- sum(eq$init_D)
  A <- sum(eq$init_A + eq$init_U)
  T <- sum(eq$init_T)

  list(
    EIR = eir, ft = ft,
    S = S, D = D, A = A, T = T, phi = phi, b = b, n = eq$delayGam,
    m = eq$mv0, Sv = eq$init_Sv, Ev = eq$init_Ev, Iv = eq$init_Iv, mu = eq$mu0,
    n = eq$delayGam,
    mu = eq$mu0,
    rD = eq$rD,
    rA = eq$rA,
    rT_S = eq$rT,
    rT_R = eq$rT*2
  ) %>%
    as.data.frame()

}

pars <- expand.grid("eir" = EIRs, "ft" = fts)
pars$n <- seq_along(pars$eir)

starting_params <- lapply(split(pars, pars$n),
                          function(x){
                            phi_eir_rel(x$eir, x$ft)
                          }) %>% do.call(rbind, .)

saveRDS(starting_params, "analysis/data/derived/starting_params.rds")
write.csv(starting_params, "analysis/data/derived/starting_params.csv")
