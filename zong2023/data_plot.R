rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
theme = theme(axis.text.x = element_text(margin = margin(t = 0.1, unit = "cm")),
              axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm")),
              axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())

#### LOAD AND GROOM DATA ----
raw_data = read_xlsx("data/continuous_monitoring/Zong-2023-QSR.xlsx", sheet = 1, skip = 1)
meteorological_data = raw_data[, 1:4] |>
  drop_na() 
names(meteorological_data) = c("date", "P_mm", "d18_p", "dD_p")
meteorological_data = meteorological_data |>
  mutate(date = as.Date(date)) |>
  group_by(date) |>
  summarise(P_mm = mean(P_mm),
            d18_p = mean(d18_p),
            dD_p = mean(dD_p))

monitor_data = raw_data[, 5:ncol(raw_data)] |>
  drop_na()
names(monitor_data) = c("date", "T_mean", "RH_mean", "T_max", "RH_min", "T_min",
                        "RH_max", "T_day", "RH_day", "T_night", "RH_night",
                        "d18_sw", "dD_sw", "d18_bw", "dD_bw",
                        "bwc", "swc", "snail_activity")
monitor_data = monitor_data |>
  mutate(date = as.Date(date)) |>
  mutate(across(-1, as.numeric)) |>
  mutate(d_sw_excess = dD_sw - 8 * d18_sw,
         d_bw_excess = dD_bw - 8 * d18_bw)

# find the days after rainfall events
for (i in 1:nrow(monitor_data)) {
  target_date = monitor_data$date[i]
  meteo_filtered = meteorological_data[meteorological_data$date <= target_date, ]
  closest_date = meteo_filtered[which.max(meteo_filtered$date), ]
  monitor_data$drought_days[i] = monitor_data$date[i] - closest_date$date
  monitor_data$d18_p[i] = closest_date$d18_p
  monitor_data$d2_p[i] = closest_date$dD_p
}
write_csv(meteorological_data, file = "output/meteorological_data_zong2023.csv")
write_csv(monitor_data, file = "output/monitor_data_zong2023.csv")

#### d18O-dD plot ----
m1 = lm(dD_p ~ d18_p, data = meteorological_data)
summary(m1)
m2 = lm(dD_sw ~ d18_sw, data = monitor_data)
summary(m2)
m3 = lm(dD_bw ~ d18_bw, data = monitor_data)
summary(m3)
p1 = ggplot() +
  geom_abline(slope = m1$coefficients[2],
              intercept = m1$coefficients[1],
              color = "grey50") +
  geom_abline(slope = m2$coefficients[2],
              intercept = m2$coefficients[1],
              color = "#8c510a") +
  geom_abline(slope = m3$coefficients[2],
              intercept = m3$coefficients[1],
              color = "#01665e") +
  geom_point(data = meteorological_data, 
             aes(x = d18_p, y = dD_p),
             size = 2, color = "grey80", shape = 21) +
  geom_point(data = monitor_data,
             aes(x = d18_sw, y = dD_sw), na.rm = TRUE,
             size = 2, shape = 21, color = "#d8b365") +
  geom_point(data = monitor_data,
             aes(x = d18_bw, y = dD_bw),
             size = 2, shape = 22, color = "#5ab4ac") +
  annotate("text", x = -18, y = 10, hjust = 0,
           label = expression(delta*"D"[p]*" = 7.8 × "*delta^"18"*"O"[p]*" + 12.0"),
           color = "grey50", size = 4) +
  annotate("text", x = -12, y = -120, hjust = 0,
           label = expression(delta*"D"[sw]*" = 7.3 × "*delta^"18"*"O"[sw]*" + 9.0"),
           color = "#8c510a", size = 4) +
  annotate("text", x = -12, y = -140, hjust = 0,
           label = expression(delta*"D"[bw]*" = 7.4 × "*delta^"18"*"O"[bw]*" + 2.8"),
           color = "#01665e", size = 4) +
  annotate("text", x = -17, y = 30, hjust = 1,
           label = "a", size = 5, fontface = "bold") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O (\u2030, VSMOW)"),
       y = expression(delta*"D (\u2030, VSMOW)"))
p1
p2 = ggplot(monitor_data) +
  geom_point(aes(x = drought_days, y = d_sw_excess), shape = 21, size = 2, color = "#d8b365") +
  geom_point(aes(x = drought_days, y = d_bw_excess), shape = 21, size = 2, color = "#5ab4ac") +
  annotate("text", x = 17, y = 18, hjust = 1,
           label = "b", size = 5, fontface = "bold") +
  theme_bw() + theme +
  labs(x = "days after rain",
       y = "d-excess (\u2030)")
p2
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv")
ggsave("figures/zong2023_data.jpg", width = 7.3, height = 3.8)
