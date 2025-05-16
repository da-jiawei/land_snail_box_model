ctrl = function(){
  vars = list(
    t = 60 * 60,
    M_bw = 2e-2, # water molar mass per wet body weight - unit: mol/g
    M_b = 1.6e-7, # HCO3- molar mass per wet body weight - unit:mol/g
    Po = 15800, # partial pressure of CO2 in the lungs at the interface between the snail body and the diffusive boundary layer
    Pa = 400, # partial pressure of CO2 in the open atmosphere
    F_evap = 2.7e-6, # evaporation rate per wet body weight - unit: mol/g/s (Dainton, 1954)
    # F_evap might depend on the size of the surface area of the soft body 
    F_rs = 9.7e-10, # CO2 respiration rate per wet body weight - unit: mol/g/s (Barnhart, 1986)
    theta = .4, # F_out/F_in for oxygen (i.e, water)
    phi = .4, # F_out/F_in for carbon (input is food, output is HCO3- in body fluid)
    temp = 15, # air temperature during snail activity
    RH = .8, # air humidity during snail activity
    air_temp = 15, # mean daily air temp
    air_RH = .6, # mean daily air humidity
    d18_mw = -5,
    d13_in = -15,
    d13_a = -8,
    dew_influx_ratio = 0
  )
}

snail = function(vars) {
  list2env(vars, environment())
  
  #### CONSTANTS
  R18_smow = 0.0020052
  R17_smow = 0.0003799
  R2_smow = 0.00015576
  R18_vpdb = 0.0020672
  R17_vpdb = 0.0003860
  R13_vpdb = 0.011237
  theta_eq = 0.529 
  theta_diff = 0.5185
  alpha18_diff = 1/0.9723 # kinetic fractionation factor - H216O / H218O (Merlivat 1978)
  alpha17_diff = alpha18_diff ^ theta_diff
  alpha2_diff = 1/0.9755
  alpha13_eq_a_b = 1.0027 # equilibrium isotopic fractionation factor between aragonite and HCO3- at 10-25 degrees 
  alpha13_diff = 1.0044 # kinetic fractionation factor between 12CO2 / 13CO2 - Cerling and Quade (1993)
  temp.kelvin = temp + 273.15
  air_temp.kelvin = air_temp + 273.15 # evaporation mainly occur during the day
  
  #### ATMOSPHERIC OXYGEN
  # d18_o2 = 23.881 # atmospheric O2 - Barkan and Luz (2011)
  # Dp17_o2 = -506 # Barkan and Luz (2011)
  # R18_o2 = (d18_o2 / 1e3 + 1) * R18_smow
  # dp18_o2 = log(R18_o2 / R18_smow) * 1e3
  # dp17_o2 = Dp17_o2 / 1e3 + 0.528 * dp18_o2
  # R17_o2 = exp(dp17_o2 / 1e3) * R17_smow
  # alpha18_o2 = .993
  # lambda_o2 = .5179
  # alpha17_o2 = alpha18_o2 ^ lambda_o2
  # R18_o2_lung = R18_o2 * alpha18_o2
  # R17_o2_lung = R17_o2 * alpha17_o2
  
  
  #### OXYGEN ISOTOPE ----
  # equilibrium fractionation factor between liquid and vapor (Majoube 1971)
  alpha18_eq_l_v = exp(-2.0667e-3 - 0.4156/(temp.kelvin) + 1.137e3/((temp.kelvin)^2))
  alpha17_eq_l_v = alpha18_eq_l_v ^ theta_eq 
  alpha2_eq_l_v = exp(52.612e-3 - 76.248/(temp.kelvin) + (24.844e3)/(temp.kelvin)^2)
  alpha18_eq_air_l_v = exp(-2.0667e-3 - 0.4156/(air_temp.kelvin) + 1.137e3/((air_temp.kelvin)^2))
  alpha17_eq_air_l_v = alpha18_eq_air_l_v ^ theta_eq 
  alpha2_eq_air_l_v = exp(52.612e-3 - 76.248/(air_temp.kelvin) + (24.844e3)/(air_temp.kelvin)^2)
  
  # isotopic values of meteroic water
  R18_mw = (d18_mw / 1e3 + 1) * R18_smow 
  dp18_mw = 1e3 * log(d18_mw / 1e3 + 1)
  dp17_mw = 0.5268 * dp18_mw + 0.015 # Global meteoric water - Aron (2021)
  R17_mw = exp(dp17_mw / 1e3) * R17_smow
  d2_mw = 8 * d18_mw + 10
  R2_mw = (d2_mw / 1e3 + 1) * R2_smow
  
  # isotopic values of atmospheric water vapor
  # assume to be mainly from evaporation
  # R18_a = R18_mw / (alpha18_eq_air_l_v * (alpha18_diff * (1 - air_RH) + air_RH)) # assume local atmospheric wate vapor in isotopic equilibrium with meteoric water
  # R17_a = R17_mw / (alpha17_eq_air_l_v * (alpha17_diff * (1 - air_RH) + air_RH))
  # R2_a = R2_mw / (alpha2_eq_air_l_v * (alpha2_diff * (1 - air_RH) + air_RH))
  # assume to be equilibrium with meteoric water
  R18_a = R18_mw / alpha18_eq_air_l_v
  R17_a = R17_mw / alpha17_eq_air_l_v
  R2_a = R2_mw / alpha2_eq_air_l_v

  # isotopic values of the dew 
  R18_d = R18_a * alpha18_eq_l_v
  R17_d = R17_a * alpha17_eq_l_v
  R2_d = R2_a * alpha2_eq_l_v
  
  R18_in = R18_mw * (1 - dew_influx_ratio) + R18_d * dew_influx_ratio
  R17_in = R17_mw * (1 - dew_influx_ratio) + R17_d * dew_influx_ratio
  R2_in = R2_mw * (1 - dew_influx_ratio) + R2_d * dew_influx_ratio
  # F_o2 = F_rs * 2 # atmospheric O2 as HCO3-
  # F_mw = F_in - F_o2
  # R18_in = (R18_o2_lung * F_o2 + R18_mw * F_mw) / F_in
  # R17_in = (R17_o2_lung * F_o2 + R17_mw * F_mw) / F_in
  
  # Influx assuming steady state
  ratio = F_evap / M_bw
  A18 = ratio * (alpha18_diff * alpha18_eq_l_v * (1 - RH) * R18_in + alpha18_eq_l_v * R18_a * RH * (1 - theta)) / (alpha18_diff * alpha18_eq_l_v * (1 - RH) * (1 - theta))
  A17 = ratio * (alpha17_diff * alpha17_eq_l_v * (1 - RH) * R17_in + alpha17_eq_l_v * R17_a * RH * (1 - theta)) / (alpha17_diff * alpha17_eq_l_v * (1 - RH) * (1 - theta))
  A2 = ratio * (alpha2_diff * alpha2_eq_l_v * (1 - RH) * R2_in + alpha2_eq_l_v * R2_a * RH * (1 - theta)) / (alpha2_diff * alpha2_eq_l_v * (1 - RH) * (1 - theta))
  
  B18 = ratio * ((1 - theta) + alpha18_diff * alpha18_eq_l_v * (1 - RH) * theta) / (alpha18_diff * alpha18_eq_l_v * (1 - RH) * (1 - theta)) 
  B17 = ratio * ((1- theta) + alpha17_diff * alpha17_eq_l_v * (1 - RH) * theta) / (alpha17_diff * alpha17_eq_l_v * (1 - RH) * (1 - theta)) 
  B2 = ratio * ((1- theta) + alpha2_diff * alpha2_eq_l_v * (1 - RH) * theta) / (alpha2_diff * alpha2_eq_l_v * (1 - RH) * (1 - theta))
  
  tau18 = 1 / B18 # time constant for 63% approach to steady state
  tau17 = 1 / B17
  tau2 = 1 / B2
  
  R18_ss = A18 / B18
  d18_ss = (R18_ss / R18_smow - 1) * 1e3 # steady state d18O of body fluid
  R17_ss = A17 / B17
  d17_ss = (R17_ss / R17_smow - 1) * 1e3
  Dp17_ss = (1e3 * log(R17_ss / R17_smow) - .528 * 1e3 * log(R18_ss / R18_smow)) * 1e3 # per meg
  R2_ss = A2 / B2
  d2_ss = (R2_ss / R2_smow - 1) * 1e3
  
  R18_l = (R18_in - R18_ss) * exp(-t / tau18) + R18_ss
  R17_l = (R17_in - R17_ss) * exp(-t / tau17) + R17_ss
  R2_l = (R2_in - R2_ss) * exp(-t / tau2) + R2_ss
  d18_l = (R18_l / R18_smow - 1) * 1e3
  d17_l = (R17_l / R17_smow - 1) * 1e3
  Dp17_l = (1e3 * log(R17_l / R17_smow) - .528 * 1e3 * log(R18_l / R18_smow)) * 1e3 # per meg
  d2_l = (R2_l / R2_smow - 1) * 1e3
  d18_c = d18_l + (19.7 - temp) / 4.34 # Grossman and Ku (1986)
  
  # alpha18_eq_a_b = exp((18.672e3 / (temp + 273.15) - 63.44) / 1e3) # equilibrium fractionation factor between aragonite shell and water - Grossman (1986)
  # alpha18_eq_a_b = exp((16.74e3 / (temp + 273.15) - 26.39) / 1e3) # White (1999)
  # alpha18_eq_a_b = exp((18.45e3 / (temp + 273.15) - 32.54) / 1e3) # Bohm (2000)
  # lamda_eq = 0.5252 # sharp experimental data
  # alpha17_eq_a_b = alpha18_eq_a_b ^ lamda_eq
  
  # R18_c = R18_l * alpha18_eq_a_b
  # R17_c = R17_l * alpha17_eq_a_b
  # d18_c = (R18_c / R18_vpdb - 1) * 1e3
  # dp18_c = 1e3 * log(R18_c / R18_vpdb)
  # dp17_c = 1e3 * log(R17_c / R17_vpdb)
  # Dp17_c = (dp17_c - 0.528 * dp18_c) * 1e3
  
  #### CARBON ISOTOPE  ----
  alpha13_eq_b_g = exp((9.552e3 / temp.kelvin - 24.1) / 1e3) # equilibrium isotopic fractionation factor between HCO3- and CO2 - Mook (1974)
  eta = Pa / Po
  
  d13_a = -8
  R13_a = (d13_a / 1e3 + 1) * R13_vpdb
  R13_in = (d13_in / 1e3 + 1) * R13_vpdb
  G13 = ((alpha13_diff * alpha13_eq_b_g * (1 - eta) * R13_in + alpha13_eq_b_g * eta * R13_a * (1 - phi)) / (alpha13_diff * alpha13_eq_b_g * (1 - eta) * (1 - phi))) * (F_rs / M_b)       
  H13 = (((1 - phi) + alpha13_diff * alpha13_eq_b_g * (1 - eta) * phi) / (alpha13_diff * alpha13_eq_b_g * (1 - eta) * (1 - phi))) * (F_rs / M_b)

  tau13 = 1 / H13
  d13_ss = ((G13 / H13) / R13_vpdb - 1) * 1e3
  d13_b = (d13_in - d13_ss) * exp(-t * tau13) + d13_ss
  d13_c = (1e3 + d13_b) * alpha13_eq_a_b - 1e3
  
  d18_a = (R18_a / R18_smow - 1) * 1e3
  d2_a = (R2_a / R2_smow - 1) * 1e3
  
  results = data.frame(tau18 = tau18, tau13 = tau13,
                       d18_l = d18_l, Dp17_l = Dp17_l, d2_l = d2_l,
                       d18_ss = d18_ss, Dp17_ss = Dp17_ss, d2_ss = d2_ss,
                       # dp18_c = dp18_c, Dp17_c = Dp17_c,
                       d18_a = d18_a, d2_a = d2_a,
                       d18_c = d18_c, d13_c = d13_c
                       )
  return(results)
}








