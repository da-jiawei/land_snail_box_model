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
  "d18_w_obs" = monitor_data$d18_sw,
  "d2_w_obs" = monitor_data$dD_sw,
  "temp_obs" = temp,
  "RH_obs" = RH
)
parms = c("RH", "temp", "theta",
          "d18_bw", "d2_bw")
system.time({post.clp = jags.parallel(data_obs, NULL, parms, "model/flux_balance_model_bayesian_inversion_ver1.R",
                                      n.iter = 5e5, n.chains = 5, n.burnin = 1e5)})
save(post.clp, file = "output/ver1.rda")

view(post.clp$BUGSoutput$summary)
# traceplot(post.clp, varname = "theta")
post_data = data.frame("date" = monitor_data$date)
for (i in 1:length(parms)) {
  post_data$Rhat = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "Rhat"]
  post_data$n.eff = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "n.eff"]
  post_data[, paste0("post_", parms[i])] = post.clp$BUGSoutput$mean[parms[i]]
  post_data[, paste0("post_", parms[i], "_sd")] = post.clp$BUGSoutput$sd[parms[i]]
  for (p in 1:nrow(post_data)) {
    if (post_data$Rhat[p] < 1.05 & post_data$n.eff[p] >= 1000){
      post_data[p, paste0(parms[i], "_eff")] = "positive"
    } else {
      post_data[p, paste0(parms[i], "_eff")] = "negative"
    }
  }
}
data = left_join(monitor_data, post_data, by = "date")
write.csv(data, file = "output/posteriors_scenario1.csv")

#### PLOT ----
data = read_csv("output/posteriors_scenario1.csv")
data$DRH = data$post_RH - data$RH_mean/1e2
p1 = ggplot(data, aes(x = RH_mean, y = DRH * 1e2)) +
  geom_abline(slope = 0, intercept = 0) +
  geom_errorbar(aes(ymin = (DRH - post_RH_sd)*1e2, 
                    ymax = (DRH + post_RH_sd)*1e2,
                    alpha = RH_eff),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = T_mean, alpha = RH_eff), 
             shape = 21, size = 2) +
  annotate("text", x = 95, y = 35, label = "a", size = 5, fontface = "bold") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.2)) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "right") +
  labs(x = expression("RH"[mean]*" (%)"),
       y = expression(Delta*"RH"[post-mean]*" (%)"), 
       fill = expression(paste("T"[mean]*" (", degree, "C)")))

data$DT = data$post_temp - data$T_mean
p2 = ggplot(data, aes(x = RH_mean, y = DT)) +
  geom_abline(slope = 0, intercept = 0) +
  geom_errorbar(aes(ymin = DT - post_temp_sd, 
                    ymax = DT + post_temp_sd,
                    alpha = temp_eff), 
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = T_mean, alpha = temp_eff), 
             shape = 21, size = 2) +
  annotate("text", x = 95, y = -7, label = "b", size = 5, fontface = "bold") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.2)) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "right") +
  labs(x = expression("RH"[mean]*" (%)"),
       y = expression(paste(Delta*"T (", degree, "C)")),
       fill = expression(paste("T"[mean]*" (", degree, "C)")))

p3 = ggplot(data, aes(x = RH_mean, y = post_theta)) +
  geom_errorbar(aes(ymin = post_theta - post_theta_sd, 
                    ymax = post_theta + post_theta_sd,
                    alpha = theta_eff),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = T_mean, alpha = theta_eff),
             shape = 21, size = 2) +
  annotate("text", x = 35, y = .04, label = "c", size = 5, fontface = "bold") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.2)) +
  scale_fill_distiller(palette = "RdBu") +
  guides(alpha = "none") +
  theme_bw() + theme +
  theme(legend.position = "right") +
  labs(fill = expression(paste("T"[mean]*" (", degree, "C)")),
       y = expression(italic(theta)*" (F"[out]*"/ F"["in"]*")"),
       x = expression("RH"[mean]*" (%)"))
upper = ggarrange(p1, p2, p3, nrow = 1, ncol = 3, align = "hv", 
          common.legend = TRUE, legend = "right")
upper

p4 = ggplot(data, aes(x = d18_bw, y = post_d18_bw)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(data = sims, aes(x = d18_bw, y = d18_ss), shape = 22, size = 2, color = "grey80") +
  geom_errorbar(aes(ymin = post_d18_bw - post_d18_bw_sd, 
                    ymax = post_d18_bw + post_d18_bw_sd,
                    alpha = d18_bw_eff),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = RH_mean, alpha = d18_bw_eff), 
             shape = 21, size = 2) +
  annotate("text", x = -15, y = 14, label = "d", size = 5, fontface = "bold") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.2)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "right") +
  labs(x = expression(paste("measured "*delta^"18"*"O"[bw]*" (\u2030)")),
       y = expression(paste("modeled "*delta^"18"*"O"[bw]*" (\u2030)")),
       fill = expression("RH"[mean]*" (%)"))

p5 = ggplot(data, aes(x = dD_bw, y = post_d2_bw)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(data = sims, aes(x = d2_bw, y = d2_ss), shape = 22, size = 2, color = "grey80") +
  geom_errorbar(aes(ymin = post_d2_bw - post_d2_bw_sd, 
                    ymax = post_d2_bw + post_d2_bw_sd,
                    alpha = d2_bw_eff),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = RH_mean, alpha = d2_bw_eff), shape = 21, size = 2) +
  annotate("text", x = -110, y = 65, label = "e", size = 5, fontface = "bold") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 0.2)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_bw() + theme +
  guides(alpha = "none") +
  theme(legend.position = "right") +
  labs(x = expression(paste("measured "*delta*"D"[bw]*" (\u2030)")),
       y = expression(paste("modeled "*delta*"D"[bw]*" (\u2030)")),
       fill = expression("RH"[mean]*" (%)"))
lower = ggarrange(p4, p5, nrow = 1, ncol = 2,
                  common.legend = TRUE, legend = "right")
lower
ggarrange(upper, lower, nrow = 2, ncol = 1, common.legend = TRUE)
ggsave("figures/scenario1.jpg", width = 7.7, height = 4.7, dpi = 500)






