## Human Equations
deriv(S) <- -S * Lambda_s * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) -
  S * Lambda_R * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) +
  Ts * rTs + As * rA + AR * rA + TR * rTR

deriv(Ds) <- S * Lambda_s * Phi * (1 - fT) +
  Lambda_s * As * Phi * (1 - fT) -
  Ds * rD

deriv(As) <- S * Lambda_s * (1 - Phi) +
  Ds * rD -
  Lambda_s * As * Phi * (1 - fT) -
  Lambda_s * As * Phi * fT -
  As * rA

deriv(Ts) <- S * Lambda_s * Phi * fT +
  Lambda_s * As * Phi * fT -
  Ts * rTs

deriv(DR) <- S * Lambda_R * Phi * (1 - fT) +
  Lambda_R * AR * Phi * (1 - fT) -
  DR * rD

deriv(AR) <- S * Lambda_R * (1 - Phi) +
  DR * rD -
  Lambda_R * AR * Phi * (1 - fT) -
  Lambda_R * AR * Phi * fT -
  AR * rA

deriv(TR) <- S * Lambda_R * Phi * fT +
  Lambda_R * AR * Phi * fT -
  TR * rTR

## Mosquito Equations
deriv(Sv) <- e - (Lambda_v_s + Lambda_v_r) * Sv - mu * Sv

#delayed_Lambda_v_s_Sv <- delay(Lambda_v_s * Sv * exp(-mu * n), n)
delayed_Lambda_v_s_Sv <- Lambda_v_s * Sv * exp(-mu * n)

deriv(Ev_s) <- Lambda_v_s * Sv - delayed_Lambda_v_s_Sv - mu * Ev_s
deriv(Iv_s) <- delayed_Lambda_v_s_Sv - mu * Iv_s

#delayed_Lambda_v_r_Sv <- delay(Lambda_v_r * Sv * exp(-mu * n), n)
delayed_Lambda_v_r_Sv <- Lambda_v_r * Sv * exp(-mu * n)

deriv(Ev_r) <- Lambda_v_r * Sv - delayed_Lambda_v_r_Sv - mu * Ev_r
deriv(Iv_r) <- delayed_Lambda_v_r_Sv - mu * Iv_r

output(mos) <- Sv + Ev_s + Iv_s + Ev_r + Iv_r
output(pop) <- S + As + Ds + Ts + AR + DR + TR
output(res) <- (AR + DR + TR)/(As + Ds + Ts + AR + DR + TR)

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
S0 <- user(0.8380032)

Ds0 <-  user()
As0 <-  user()
Ts0 <-  user()
DR0 <-  user()
AR0 <-  user()
TR0 <-  user()
#Lambda_s <- user(0.001214977)
m <- user(1.271483)
b <- user(0.5876259)
Lambda_s <- m*a*b*Iv_s
Lambda_R <- m*a*b*Iv_r
Phi <- user(0.6998258)
fT <- user(0.1)
rD <- user(0.2)
rA <- user(0.00512821)

rTs <- user(0.2)
rTR_true <- user(0.2)
# ton is a time where we are happy the simulation is at equilibrium
ton <- user(10000)
toff <- user(20000)
rTR <- if(t > ton && t < toff) rTR_true else rTs


Sv0 <- user(0.961586)
Ev_s0 <- user()
Iv_s0 <- user()
Ev_r0 <- user()
Iv_r0 <- user()
e <- user(0.132)
mu <- user(0.132)
n <- user(10)
a <- user(0.3066667)
c_A <- user(0.05047807)
c_D <- user(0.0676909)
c_T <- user(0.02179647)
Lambda_v_s <- a * (c_A*As + c_D*Ds + c_T*Ts)
Lambda_v_r <- a * (c_A*AR + c_D*DR + c_T*TR)
output(Lambda_v_s) <- Lambda_v_s
output(Lambda_v_r) <- Lambda_v_r
output(Lambda_s) <- Lambda_s
output(Lambda_R) <- Lambda_R
output(delayed_Lambda_v_s_Sv) <- delayed_Lambda_v_s_Sv
output(delayed_Lambda_v_r_Sv) <- delayed_Lambda_v_r_Sv

