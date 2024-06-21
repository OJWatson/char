## Human Equations

deriv(S) <- (-S * Lambda_s * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) -
               S * Lambda_R * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) +
               Ts * rTs + As * rA + AR * rA + TR * rTR)

deriv(Ds) <- (S * Lambda_s * Phi * (1 - fT) +
                Lambda_s * AR * Phi * (1 - fT) +
                Lambda_s * As * Phi * (1 - fT) -
                Ds * rD)

deriv(As) <- (S * Lambda_s * (1 - Phi) +
                Lambda_s * AR * (1 - Phi) +
                Ds * rD -
                Lambda_s * As * Phi * (1 - fT) -
                Lambda_s * As * Phi * fT -
                Lambda_R * As * Phi * (1 - fT) -
                Lambda_R * As * Phi * fT -
                Lambda_R * As * (1 - Phi) -
                As * rA) - invading_AR

deriv(Ts) <- (S * Lambda_s * Phi * fT +
                Lambda_s * AR * Phi * fT +
                Lambda_s * As * Phi * fT -
                Ts * rTs)

deriv(DR) <- (S * Lambda_R * Phi * (1 - fT) +
                Lambda_R * As * Phi * (1 - fT) +
                Lambda_R * AR * Phi * (1 - fT) -
                DR * rD)

deriv(AR) <- invading_AR + (S * Lambda_R * (1 - Phi) +
                              Lambda_R * As * (1 - Phi) +
                              DR * rD -
                              Lambda_R * AR * Phi * (1 - fT) -
                              Lambda_R * AR * Phi * fT -
                              Lambda_s * AR * Phi * (1 - fT) -
                              Lambda_s * AR * Phi * fT -
                              Lambda_s * AR * (1 - Phi) -
                              AR * rA)

deriv(TR) <- (S * Lambda_R * Phi * fT +
                Lambda_R * AR * Phi * fT +
                Lambda_R * As * Phi * fT -
                TR * rTR)

## Mosquito Equations
deriv(Sv) <- e - (Lambda_v_s + Lambda_v_r) * Sv - mu * Sv

delayed_Lambda_v_s_Sv <- delay(Lambda_v_s * Sv * exp(-mu * n), n)

deriv(Ev_s) <- Lambda_v_s * Sv - delayed_Lambda_v_s_Sv - mu * Ev_s
deriv(Iv_s) <- delayed_Lambda_v_s_Sv - mu * Iv_s

delayed_Lambda_v_r_Sv <- delay(Lambda_v_r * Sv * exp(-mu * n), n)

deriv(Ev_r) <- Lambda_v_r * Sv - delayed_Lambda_v_r_Sv - mu * Ev_r
deriv(Iv_r) <- delayed_Lambda_v_r_Sv - mu * Iv_r

Lambda_v_s <- a * (c_A * As + c_D * Ds + c_T * Ts)
Lambda_v_r <- a * (c_A * AR + c_D * DR + c_T * TR)


## Outputs
output(prevalence) <- As + Ds + Ts + AR + DR + TR
output(prevalence_res) <- (AR + DR + TR) / (As + Ds + Ts + AR + DR + TR)
output(population) <- S + Ds + As + Ts + DR + AR + TR
output(population_v) <- Sv + Ev_s + Iv_s + Ev_r + Iv_r

# Resistance introduction
invading_AR <- if(t < res_time || t > (res_time+0.05)) 0 else As*init_res
output(invading_AR_out) <- invading_AR

## Initial conditions
initial(S) <- S0
initial(Ds) <- Ds0
initial(As) <- As0
initial(Ts) <- Ts0
initial(DR) <- DR0
initial(AR) <- AR0
initial(TR) <- TR0
initial(Sv) <- Sv0
initial(Ev_s) <- Ev_s0
initial(Iv_s) <- Iv_s0
initial(Ev_r) <- Ev_r0
initial(Iv_r) <- Iv_r0


## User-defined parameters
S0 <- user()
Ds0 <- user()
As0 <- user()
Ts0 <- user()
DR0 <- user()
AR0 <- user()
TR0 <- user()
m <- user()
a <- user()
b <- user()
Lambda_s <- m * a * b * Iv_s
Lambda_R <- m * a * b * Iv_r
Phi <- user()
fT <- user()
rD <- user()
rA <- user()
rTs <- user()
rTR_true <- user()
rTR <- if (t > ton && t < toff) rTR_true else rTs
Sv0 <- user()
Ev_s0 <- user()
Iv_s0 <- user()
Ev_r0 <- user()
Iv_r0 <- user()
e <- user()
mu <- user()
n <- user()
c_A <- user()
c_D <- user()
c_T <- user()
ton <- user()
toff <- user()
res_time <- user()
init_res <- user()
