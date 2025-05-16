rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
source("model/Flux_balance_model.R")
monitor_data = read_csv("output/monitor_data_zong2023.csv")
meteorological_data = read_csv("output/meteorological_data_zong2023.csv")
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
cat("\014")

#### ISOTOPE STEADY STATE MODEL ----
vars = ctrl()
for (i in 1:nrow(monitor_data)) {
  vars$RH = min(monitor_data$RH_max[i] / 100, .99)
  vars$air_RH = vars$RH
  vars$temp = monitor_data$T_min[i]
  vars$air_temp = vars$temp
  vars$d18_mw = monitor_data$d18_p[i]
  vars$d2_mw = monitor_data$d2_p[i]
  results = snail(vars)
  results$date = monitor_data$date[i]
  results$d18_sw = monitor_data$d18_sw[i]
  results$d18_bw = monitor_data$d18_bw[i]
  results$d2_bw = monitor_data$dD_bw[i]
  results$d2_sw = monitor_data$dD_sw[i]
  if(i == 1){
    sims = results
  }else{
    sims = rbind(sims, results)
  }
}
sims = sims |>
  mutate(D18 = d18_ss - d18_sw) |>
  drop_na(d18_sw)
write_csv(sims, file = "output/forward_model_results.csv")

# plot model-data comparison
pal = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
jpeg("figures/time_series.jpg", width = 4.4, height = 5.2, units = "in", res = 600)
par(mar = c(4,4,1,4))
plot(0, 0, xlim = c(range(sims$date)), ylim = c(0, 3.5), axes = FALSE,
     xlab = "", ylab = "")
meteorological_data = meteorological_data |>
  filter(date >= min(sims$date) & date <= max(sims$date))
# d18O
yext = range(sims$d18_ss, sims$d18_bw, sims$d18_sw, meteorological_data$d18_p)
tix = seq(floor(min(yext)-1), ceiling(max(yext)+2), by = 5)
d18_ss.rs = cbind(sims$date, 
                  2.5 + (sims$d18_ss - min(tix)) / diff(range(tix)))
d18_bw.rs = cbind(sims$date, 
                  2.5 + (sims$d18_bw - min(tix)) / diff(range(tix)))
d18_sw.rs = cbind(sims$date, 
                  2.5 + (sims$d18_sw - min(tix)) / diff(range(tix)))
d18_p.rs = cbind(meteorological_data$date, 
                 2.5 + (meteorological_data$d18_p - min(tix)) / diff(range(tix)))
lines(d18_ss.rs[, 1], d18_ss.rs[, 2], col = pal[1])
lines(d18_bw.rs[, 1], d18_bw.rs[, 2], col = pal[2])
lines(d18_sw.rs[, 1], d18_sw.rs[, 2], col = pal[3])
points(d18_p.rs[, 1], d18_p.rs[, 2], col = "grey")
axis(2, 2.5 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste(delta^"18"*"O (\u2030)")), 2, line = 2.5, at = 3)

# dD
yext = range(sims$d2_ss, sims$d2_bw, sims$d2_sw, meteorological_data$dD_p)
tix = seq(floor(min(yext)-8), ceiling(max(yext)), by = 50)
d2_ss.rs = cbind(sims$date, 
                 1.5 + (sims$d2_ss - min(tix)) / diff(range(tix)))
d2_bw.rs = cbind(sims$date, 
                 1.5 + (sims$d2_bw - min(tix)) / diff(range(tix)))
d2_sw.rs = cbind(sims$date, 
                 1.5 + (sims$d2_sw - min(tix)) / diff(range(tix)))
d2_p.rs = cbind(meteorological_data$date, 
                1.5 + (meteorological_data$dD_p - min(tix)) / diff(range(tix)))
lines(d2_ss.rs[, 1], d2_ss.rs[, 2], col = pal[1])
lines(d2_bw.rs[, 1], d2_bw.rs[, 2], col = pal[2])
lines(d2_sw.rs[, 1], d2_sw.rs[, 2], col = pal[3])
points(d2_p.rs[, 1], d2_p.rs[, 2], col = "grey")
axis(4, 1.5 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste(delta*"D (\u2030)")), 4, line = 2.5, at = 2)

# rainfall
yext = range(meteorological_data$P_mm)
tix = seq(floor(min(yext)), ceiling(max(yext)), by = 5)
P.rs = cbind(meteorological_data$date,
             0 + (meteorological_data$P_mm - min(tix)) / diff(range(tix)))
for (i in seq_len(nrow(P.rs))) {
  rect(xleft = P.rs[i, 1] - 0.5,
       xright = P.rs[i, 1] + 0.5,
       ybottom = 0,
       ytop = P.rs[i, 2],
       col = "grey90",
       border = NA)
}
axis(2, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext("P (mm)", 2, line = 2.5, at = 0.5)

yext = range(sims$D18)
tix = seq(floor(min(yext)-1), ceiling(max(yext)+3), by = 5)
D18.rs = cbind(sims$date,
               0 + (sims$D18 - min(tix)) / diff(range(tix)))
lines(D18.rs[, 1], D18.rs[, 2], col = pal[2], lwd = 1)
axis(4, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste(Delta^"18"*"O (\u2030)")), 4, line = 2.5, at = 0.5)

tix = seq(from = min(sims$date), to = max(sims$date), by = "1 month")
axis(1, at = tix, labels = format(tix, "%Y-%m"))

legend(x = min(sims$date), y = 1.9, 
       legend = c(expression(delta^"18"*"O"[bw_modeled]),
                  expression(delta^"18"*"O"[bw]),
                  expression(delta^"18"*"O"[sw])),
       col = pal, lty = 1, lwd = 2, cex = 1, bty = "n", bg = NA)
dev.off()

