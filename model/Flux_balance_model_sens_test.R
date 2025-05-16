rm(list = ls())
pacman::p_load(tidyverse, ggpubr)
source('model/Flux_balance_model.R')
vars = ctrl()
results = snail(vars)
theme = theme(axis.text.x = element_text(margin = margin(t = 0.1, unit = "cm")),
              axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm")),
              axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())

#### tau ----
vars$theta = seq(0, 1, .2)
for (i in 1:100) {
  vars$F_evap = vars$M_bw * i * 1e-5
  results = snail(vars)
  results$theta = vars$theta
  results$F_evap_over_M_bw = i * 1e-5
  if(i == 1){
    sims = results
  }else{
    sims = rbind(sims, results)
  }
}
p1 = ggplot(sims) +
  geom_vline(xintercept = (2.7e-6/2e-2)*1e4, linetype = "dashed", linewidth = .5) +
  geom_path(aes(x = F_evap_over_M_bw*1e4, y = tau18 / 60, color = theta, group = theta)) +
  scale_color_viridis_c(option = "turbo") +
  labs(x = expression("F"[evap]*"/M"[bw]*" (10"^"-4"*")"), 
       y = expression(tau*" (min)"),
       color = expression(italic(theta))) +
  theme_bw() + theme +
  scale_y_continuous(limits = c(0, 1e2))

p2 = ggplot(sims) +
  geom_vline(xintercept = (2.7e-6/2e-2)*1e4, linetype = "dashed", linewidth = .5) +
  geom_path(aes(x = F_evap_over_M_bw*1e4, y = d18_l, color = theta, group = theta)) +
  scale_color_viridis_c(option = "turbo") +
  labs(x = expression("F"[evap]*"/M"[bw]*" (10"^"-4"*")"), 
       y = expression(delta^"18"*"O"[bw]*" (\u2030)"),
       color = expression(italic(theta))) +
  theme_bw() + theme 
plot_theta = ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", labels = c("a", "b"),
          common.legend = TRUE, legend = "right")
plot_theta

#### RH ----
vars = ctrl()
vars$RH = seq(.75, 1, .05)
for (i in 1:100) {
  vars$F_evap = vars$M_bw * i * 1e-5
  results = snail(vars)
  results$RH = vars$RH
  results$F_evap_over_M_bw = i * 1e-5
  if(i == 1){
    sims = results
  }else{
    sims = rbind(sims, results)
  }
}

p1 = ggplot(sims) +
  geom_vline(xintercept = (2.7e-6/2e-2)*1e4, linetype = "dashed", linewidth = .5) +
  geom_path(aes(x = F_evap_over_M_bw*1e4, y = tau18 / 60, color = RH * 1e2, group = RH)) +
  scale_color_viridis_c(option = "turbo") +
  labs(x = expression("F"[evap]*"/M"[bw]), 
       y = expression(tau*" (min)"),
       color = "RH (%)") +
  theme_bw() + theme +
  scale_y_continuous(limits = c(0, 1e2))

p2 = ggplot(sims) +
  geom_vline(xintercept = (2.7e-6/2e-2)*1e4, linetype = "dashed", linewidth = .5) +
  geom_path(aes(x = F_evap_over_M_bw*1e4, y = d18_l, color = RH, group = RH)) +
  scale_color_viridis_c(option = "turbo") +
  labs(x = expression("F"[evap]*"/M"[bw]), 
       y = expression(delta^"18"*"O"[bw]*" (\u2030)"),
       color = "RH (%)") +
  theme_bw() + theme
plot_RH = ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", labels = c("c", "d"),
                       common.legend = TRUE, legend = "right")
plot_RH
plot_tau_d18 = ggarrange(plot_theta, plot_RH, nrow = 2, ncol = 1, align = "hv")
plot_tau_d18
ggsave(plot_tau_d18, file = "figures/steady_state_time.jpg",
       width = 6, height = 5.4,)

#### sensitivity analyses on body water isotopes ----
# temperature
vars = ctrl() 
vars$temp = seq(10, 30, 5)
sims_temp = snail(vars)
sims_temp$temp = vars$temp
p1 = ggplot(sims_temp) +
  geom_abline(slope = 8, intercept = 10) +
  annotate("point", x = -5, y = -30, size = 4, shape = 22, fill = "white") +
  annotate("text", x = -9, y = 40, label = "a", fontface = "bold", size = 5) +
  geom_point(aes(x = d18_ss, y = d2_ss,
                 fill = temp),
             shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako") +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression(delta*"D (\u2030)"), 
       fill = expression(paste("T (", degree, "C)"))) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-10, 5)) +
  scale_y_continuous(limits = c(-100, 50))
p2 = ggplot(sims_temp) +
  geom_abline(slope = 0, intercept = 10) +
  annotate("point", x = -5, y = 10, size = 4, shape = 22, fill = "white") +
  annotate("text", x = -9, y = -45, label = "d", fontface = "bold", size = 5) +
  geom_point(aes(x = d18_ss, y = d2_ss - 8*d18_ss,
                 fill = temp),
             shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako") +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression("d-excess (\u2030)"), 
       fill = expression(paste("T (", degree, "C)"))) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-10, 5)) +
  scale_y_continuous(limits = c(-50, 20))
plot_temp = ggarrange(p1, p2, nrow = 2, ncol = 1, align = "hv", common.legend = TRUE, legend = "top")
plot_temp

# humidity
vars = ctrl()
vars$RH = seq(.7, 1, .05)
sims_RH = snail(vars)
sims_RH$RH = vars$RH
p1 = ggplot(sims_RH) +
  geom_abline(slope = 8, intercept = 10) +
  annotate("point", x = -5, y = -30, size = 4, shape = 22, fill = "white") +
  annotate("text", x = -9, y = 40, label = "b", fontface = "bold", size = 5) +
  geom_point(aes(x = d18_ss, y = d2_ss,
                 fill = RH),
             shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako",
                       direction = -1) +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression(delta*"D (\u2030)"), 
       fill = "RH (%)") +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-10, 5)) +
  scale_y_continuous(limits = c(-100, 50))
p2 = ggplot(sims_RH) +
  geom_abline(slope = 0, intercept = 10) +
  annotate("point", x = -5, y = 10, size = 4, shape = 22, fill = "white") +
  annotate("text", x = -9, y = -45, label = "e", fontface = "bold", size = 5) +
  geom_point(aes(x = d18_ss, y = d2_ss - 8*d18_ss,
                 fill = RH),
             shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako") +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression("d-excess (\u2030)"), 
       fill = "RH (%)") +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-10, 5)) +
  scale_y_continuous(limits = c(-50, 20))
plot_RH = ggarrange(p1, p2, nrow = 2, ncol = 1, align = "hv", common.legend = TRUE, legend = "top")
plot_RH

### theta (F_out / F_in for body water) ---- 
vars = ctrl()
vars$theta = seq(0, 1, .1)
sims_theta = snail(vars)
sims_theta$theta = vars$theta
p1 = ggplot(sims_theta) +
  geom_abline(slope = 8, intercept = 10) +
  annotate("point", x = -5, y = -30, size = 4, shape = 22, fill = "white") +
  annotate("text", x = -9, y = 40, label = "c", fontface = "bold", size = 5) +
  geom_point(aes(x = d18_ss, y = d2_ss,
                 fill = theta),
             shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako",
                       direction = -1) +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression(delta*"D (\u2030)"), 
       fill = expression(theta)) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-10, 5)) +
  scale_y_continuous(limits = c(-100, 50))
p2 = ggplot(sims_theta) +
  geom_abline(slope = 0, intercept = 10) +
  annotate("point", x = -5, y = 10, size = 4, shape = 22, fill = "white") +
  annotate("text", x = -9, y = -45, label = "f", fontface = "bold", size = 5) +
  geom_point(aes(x = d18_ss, y = d2_ss - 8*d18_ss,
                 fill = theta),
             shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako") +
  labs(x = expression(delta^"18"*"O (\u2030)"), 
       y = expression("d-excess (\u2030)"), 
       fill = expression(theta)) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-10, 5)) +
  scale_y_continuous(limits = c(-50, 20))
plot_theta = ggarrange(p1, p2, nrow = 2, ncol = 1, align = "hv", common.legend = TRUE, legend = "top")
plot_theta
ggarrange(plot_temp, plot_RH, plot_theta, nrow = 1, ncol = 3, align = "hv")
ggsave("figures/sens_test.jpg", width = 8, height = 6)
