#' Run the SEEIR model
#'
#' @param pars Parameter List
#' @param time_period Duration
#'
#' @return Simulation output
#' @export
#'

run_simple_SEEIR_model <- function(pars,
                                   time_period = 60000) {


  # Running the Model
  mod <- model$new(user = form_pars(pars, res_start = 0.5, ton = 500,toff=1500), unused_user_action = "ignore")
  t <- seq(from = 0, to = time_period, 1)
  results <- mod$run(t) %>% as.data.frame()
  plot(results$res)

  # Summarise inputs
  out <- list(output = results, parameters = parameters, model = mod)
  out <- structure(out, class = "squire_simulation")
  return(out)
}


form_pars <- function(pars, res_start = 0,ton=2000,toff=4000){

  list(a= pars$a,
       AR0= pars$A*res_start,
       As0= pars$A*(1-res_start),
       b= pars$b,
       c_A= pars$cA,
       c_D= pars$cD,
       c_T= pars$cT,
       DR0= pars$D*res_start,
       Ds0= pars$D*(1-res_start),
       e= pars$mu,
       Ev_s0= pars$Ev*(1-res_start),
       Ev_r0= pars$Ev*(res_start),
       fT= pars$ft,
       Iv_s0= pars$Iv*(1-res_start),
       Iv_r0= pars$Iv*(res_start),
       lambda_v_scale= pars$lambda_v_scale,
       m= pars$m,
       mu= pars$mu,
       n= pars$n,
       Phi= pars$phi,
       rA= pars$rA,
       rD= pars$rD,
       res_start=res_start,
       rTR_true= pars$rT_R,
       rTs= pars$rT_S,
       S0= pars$S,
       Sv0= pars$Sv,
       ton= ton,
       toff= toff,
       TR0= pars$T*res_start,
       Ts0= pars$T*(1-res_start))


}
