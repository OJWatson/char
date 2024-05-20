f1 <- 0.1
t <- seq(0,40,1)
s_values <- c(-0.1, 0, 0.05, 0.1, 0.2)

library(tidyverse)
dat <- expand.grid(t = t, s = s_values)
dat <- dat %>% mutate(f = (exp(t * s) * (f1 / (1 - f1))) / (1 + exp(t * s) * (f1 / (1 - f1))))

dat %>% ggplot(aes(t, f, color = as.factor(s), group = s)) + geom_line(lwd = 2) + ylab("Resistance Prevalence") + xlab("Years") + ylim(c(0,1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  scale_color_viridis_d("Selection Coefficient") +
  theme_bw()
