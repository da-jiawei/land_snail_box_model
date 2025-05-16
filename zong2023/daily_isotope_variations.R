rm(list = ls())
pacman::p_load(tidyverse, ggpubr, readxl)
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
monitor_data = read_csv("output/monitor_data_zong2023.csv") |>
  drop_na(d18_sw) |>
  select(date, drought_days, d18_bw, dD_bw, d18_sw, dD_sw, d18_p, d2_p)
# function to pick target dates
rainfall_isotope_effect = function(start_date, end_date) {
  body_water = monitor_data |>
    filter(date >= start_date & date <= end_date) |>
    select(date, d18_bw, dD_bw) |>
    mutate(type = "body water")
  names(body_water) = c("date", "d18", "d2", "type") 
  soil_water = monitor_data |>
    filter(date >= start_date & date <= end_date) |>
    select(date, d18_sw, dD_sw) |>
    mutate(type = "soil water")
  names(soil_water) = c("date", "d18", "d2", "type") 
  rainfall = monitor_data |>
    filter(date >= start_date & date <= end_date & drought_days == 0) |>
    select(date, d18_p, d2_p) |>
    mutate(type = "rainfall")
  names(rainfall) = c("date", "d18", "d2", "type") 
  sample = rbind(body_water, soil_water, rainfall)
  sample$date2 = format(sample$date, "%m-%d")
  p1 = ggplot(sample) +
    geom_path(aes(x = d18, y = d2, group = type, linetype = type), color = "gray80") +
    # geom_text(aes(x = d18, y = d2, label = date2),
    #           nudge_y = 1, size = 3) +
    geom_point(aes(x = d18, y = d2, fill = as.numeric(date), shape = type),
               size = 3) +
    scale_fill_viridis_c(option = "plasma",
                         breaks = as.numeric(pretty(sample$date, n = 5)),
                         labels = format(pretty(sample$date, n = 5), "%m-%d")) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    theme_bw() + theme +
    guides(linetype = "none",
           shape = "none") +
    labs(x = expression(delta^"18"*"O (\u2030)"),
         y = expression(delta*"D (\u2030)"),
         shape = "", fill = "")
  return(p1)
}

p1 = rainfall_isotope_effect("2021-03-30", "2021-04-07")
p2 = rainfall_isotope_effect("2021-04-19", "2021-05-02")
p3 = rainfall_isotope_effect("2021-05-13", "2021-05-21")
# rainfall_isotope_effect("2021-05-29", "2021-06-12")
# rainfall_isotope_effect("2021-06-12", "2021-06-23")
rainfall_isotope_effect("2021-07-17", "2021-07-22")
rainfall_isotope_effect("2021-08-01", "2021-08-06")
rainfall_isotope_effect("2021-08-19", "2021-08-29")
p4 = rainfall_isotope_effect("2021-07-05", "2021-07-17")
ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")
