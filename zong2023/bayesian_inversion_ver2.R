rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl, R2jags)
set.seed(42)
theme = theme(axis.text.x = element_text(margin = margin(t = 0.1, unit = "cm")),
              axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm")),
              axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 11),
              legend.key.width = unit(.5, "cm"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())

#### load and groom data ----
meteorological_data = read_csv("output/meteorological_data_zong2023.csv")
sims = read_csv("output/forward_model_results.csv") # forward modeling results
monitor_data = read_csv("output/monitor_data_zong2023.csv") |>
  drop_na(d18_sw)
# avoid 100% RH
for (i in 1:nrow(monitor_data)) {
  if(monitor_data$RH_max[i] == 100){
    monitor_data$RH_max[i] = 99.9
  }
  if(monitor_data$RH_min[i] >= 99.9){
    monitor_data$RH_min[i] = monitor_data$RH_max[i] - .1
  }
  if(monitor_data$RH_night[i] >= 99.9){
    monitor_data$RH_night[i] = monitor_data$RH_night[i] - .1
  }
}

#### BAYESIAN INVERSION ----
d18_bw = data.frame(monitor_data$d18_bw, "d18_bw_sd" = .2)
d2_bw = data.frame(monitor_data$dD_bw, "d2_bw_sd" = .2)
temp = data.frame("temp_min" = monitor_data$T_min,
                  "temp_max" = monitor_data$T_max)
RH = data.frame("RH_min" = monitor_data$RH_min / 100,
                "RH_max" = monitor_data$RH_max / 100)
data_obs = list(
  "length" = nrow(d18_bw),
  "d18_bw_obs" = d18_bw,
  "d2_bw_obs" = d2_bw,
  # "d18_w_obs" = monitor_data$d18_sw,
  # "d2_w_obs" = monitor_data$dD_sw,
  "temp_obs" = temp,
  "RH_obs" = RH
)
parms = c("RH", "temp", "theta",
          "d18_in", "d2_in", 
          "d18_bw", "d2_bw")
system.time({post.clp = jags.parallel(data_obs, NULL, parms, "model/flux_balance_model_bayesian_inversion_ver2.R",
                                      n.iter = 1e5, n.chains = 3, n.burnin = 5e4)})

view(post.clp$BUGSoutput$summary)
# traceplot(post.clp, varname = "theta")
post_data = data.frame("date" = monitor_data$date)
for (i in 1:length(parms)) {
  post_data$Rhat = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "Rhat"]
  post_data$n.eff = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "n.eff"]
  post_data[, paste0("post_", parms[i])] = post.clp$BUGSoutput$mean[parms[i]]
  post_data[, paste0("post_", parms[i], "_sd")] = post.clp$BUGSoutput$sd[parms[i]]
  for (p in 1:nrow(post_data)) {
    if (post_data$Rhat[p] <= 1.05 & post_data$n.eff[p] >= 200){
      post_data[p, paste0(parms[i], "_eff")] = "positive"
    } else {
      post_data[p, paste0(parms[i], "_eff")] = "negative"
    }
  }
}
data = left_join(monitor_data, post_data, by = "date")
write.csv(data, file = "output/posteriors_scenario2.csv")

#### PLOT ----
data = read_csv("output/posteriors_scenario2.csv")
data$DRH = data$post_RH - data$RH_mean/1e2
p1 = ggplot(data) +
  geom_abline(slope = 0, intercept = 0) +
  # geom_errorbar(aes(x = RH_mean, y = DRH * 1e2,
  #                   ymin = (DRH - post_RH_sd)*1e2, ymax = (DRH + post_RH_sd)*1e2),
  #               linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(x = RH_mean, y = DRH * 1e2, fill = T_mean, alpha = RH_eff), 
             shape = 21, size = 2) +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.5)) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "top") +
  labs(x = expression("RH"[mean]*" (%)"),
       y = expression(Delta*"RH"[post-mean]*" (%)"), 
       fill = expression(paste("T"[mean]*" (", degree, "C)")))

data$DT = data$post_temp - data$T_mean
p2 = ggplot(data) +
  geom_abline(slope = 0, intercept = 0) +
  # geom_errorbar(aes(x = T_mean, y = post_temp, xmin = T_min, xmax = T_max), 
  #               linewidth = .2, width = 0, color = "grey80") +
  geom_errorbar(aes(x = RH_mean, y = DT, 
                    ymin = DT - post_temp_sd, ymax = DT + post_temp_sd), 
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(x = RH_mean, y = DT, fill = T_mean, alpha = temp_eff), 
             shape = 21, size = 2) +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.5)) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "top") +
  labs(x = expression("RH"[mean]*" (%)"),
       y = expression(paste(Delta*"T (", degree, "C)")),
       fill = expression(paste("T"[mean]*" (", degree, "C)")))

p3 = ggplot(data) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(data = sims, aes(x = d18_bw, y = d18_ss), shape = 22, size = 2, color = "grey80") +
  # geom_errorbar(aes(x = d18_bw, y = post_d18_bw, 
  #                   ymin = post_d18_bw - post_d18_bw_sd, ymax = post_d18_bw + post_d18_bw_sd), 
  #               linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(x = d18_bw, y = post_d18_bw, fill = RH_mean, alpha = d18_bw_eff), shape = 21, size = 2) +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.5)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "top") +
  labs(x = expression(paste(delta^"18"*"O"[bw])),
       y = expression(paste(delta^"18"*"O"[bw_post])),
       fill = expression("RH"[mean]*" (%)"))

p4 = ggplot(data) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(data = sims, aes(x = d2_bw, y = d2_ss), shape = 22, size = 2, color = "grey80") +
  geom_point(aes(x = dD_bw, y = post_d2_bw, fill = RH_mean, alpha = d2_bw_eff), shape = 21, size = 2) +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.5)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "top") +
  labs(x = expression(paste(delta*"D"[bw])),
       y = expression(paste(delta*"D"[bw_post])),
       fill = expression("RH"[mean]*" (%)"))

p5 = ggplot(data) +
  geom_errorbar(aes(x = RH_mean, y = post_theta,
                    ymin = post_theta - post_theta_sd, ymax = post_theta + post_theta_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(x = RH_mean, y = post_theta, fill = T_mean, alpha = theta_eff),
             shape = 21, size = 2) +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.5)) +
  scale_fill_distiller(palette = "RdBu") +
  guides(alpha = "none") +
  theme_bw() + theme +
  theme(legend.position = "top") +
  labs(fill = expression(paste("T"[mean]*" (", degree, "C)")),
       y = expression(italic(theta)*" (F"[out]*"/ F"["in"]*")"),
       x = expression("RH"[mean]*" (%)"))

ggarrange(p1, p2, p5, p3, p4, nrow = 2, ncol = 3, align = "hv", labels = c("a", "b", "c", "d", "e"))
ggarrange(p1, p2, p5, p3, p4, p6, nrow = 2, ncol = 3, align = "hv", labels = c("a", "b", "c", "d", "e", "f"))
ggsave("figures/scenario1_rainfall.jpg", width = 7.9, height = 6.6, dpi = 500)

plot_vapor = ggplot(data) +
  geom_errorbar(aes(x = RH_mean, y = post_dew_influx_ratio,
                    ymin = post_dew_influx_ratio - post_dew_influx_ratio_sd, ymax = post_dew_influx_ratio + post_dew_influx_ratio_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(x = RH_mean, y = post_dew_influx_ratio, fill = T_mean),
             shape = 21, size = 2) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(fill = expression(paste("T"[mean]*" (", degree, "C)")),
       y = expression("F"[dew]*"/F"["in"]),
       x = expression("RH"[mean]*" (%)"))
ggarrange(plot_theta, plot_vapor, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE, legend = "top")

summary(lm(data = meteorological_data, dD_p ~ d18_p))
summary(lm(data = monitor_data, dD_sw ~ d18_sw))
p6 = ggplot(data) +
  geom_abline(slope = 7.8, intercept = 12) +
  geom_abline(slope = 7.3, intercept = .9, linetype = "dashed") +
  # geom_point(data = meteorological_data,
  #            aes(x = d18_p, y = dD_p), color = "grey80", shape = 21) +
  geom_point(aes(x = d18_sw, y = dD_sw, color = RH_mean), shape = 21, size = 2) +
  geom_point(aes(x = post_d18_in, y = post_d2_in, fill = RH_mean), shape = 22, size = 2) +
  # scale_alpha_manual(values = c("positive" = 1, "negative" = 0.5)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  theme_bw() + theme +
  theme(legend.position = "top") +
  guides(color = "none") +
  labs(x = expression(delta^"18"*"O (\u2030)"),
       y = expression(delta*"D (\u2030)"),
       fill = expression("RH"[mean]*" (%)"))
p6

#### comparing din posteriors with dp ----
p1 = ggplot(data, aes(x = RH_mean, y = d18_p - post_d18_in)) +
  geom_abline(slope = 0, intercept = 0) +
  geom_errorbar(aes(ymin = d18_p - post_d18_in - post_d18_in_sd,
                    ymax = d18_p - post_d18_in + post_d18_in_sd),
                linewidth = .2, color = "grey80") +
  geom_point(aes(fill = T_mean),
             shape = 21, size = 2.5) +
  annotate("text", x = 100, y = 18, hjust = 1,
           label = "a", fontface = "bold", size = 5) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_bw() + theme +
  labs(x = expression("RH"[mean]*" (%)"),
       y = expression(Delta^"18"*"O"[p-"in"]*" (\u2030)"),
       fill = expression(paste("T"[mean]*" (", degree, "C)")))
  
p2 = ggplot(data, aes(x = RH_mean, y = d18_sw - post_d18_in)) +
  geom_abline(slope = 0, intercept = 0) +
  geom_errorbar(aes(ymin = d18_sw - post_d18_in - post_d18_in_sd,
                    ymax = d18_sw - post_d18_in + post_d18_in_sd),
                linewidth = .2, color = "grey80") +
  geom_point(aes(fill = T_mean),
             shape = 21, size = 2.5) +
  annotate("text", x = 100, y = 20, hjust = 1,
           label = "b", fontface = "bold", size = 5) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_bw() + theme +
  labs(x = expression("RH"[mean]*" (%)"),
       y = expression(Delta^"18"*"O"[sw-"in"]*" (\u2030)"),
       fill = expression(paste("T"[mean]*" (", degree, "C)")))
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE, legend = "right")
ggsave("figures/D18_difference_scenario2.jpg", width = 7, height = 3.5, dpi = 500)







