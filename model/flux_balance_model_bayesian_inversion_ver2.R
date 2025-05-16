model{
  for (i in 1:length) {
    d18_bw_obs[i, 1] ~ dnorm(d18_bw[i], d18_bw_pre[i])
    d18_bw_pre[i] = 1 / d18_bw_obs[i, 2] ^ 2
    d2_bw_obs[i, 1] ~ dnorm(d2_bw[i], d2_bw_pre[i])
    d2_bw_pre[i] = 1 / d2_bw_obs[i, 2] ^ 2
    temp[i] ~ dunif(temp_obs[i, 1], temp_obs[i, 2]) # temperature during snail activity
    RH[i] ~ dunif(RH_obs[i, 1], RH_obs[i, 2])
    d18_w[i] ~ dunif(-20, 5)
    d2_w[i] ~ dunif(-150, 50)
    theta[i] ~ dunif(0, 1)

    #### OXYGEN ISOTOPE ----
    # temperature during the snail activity
    temp.kelvin[i] = temp[i] + 273.15
    alpha18_eq_l_v[i] = exp(-2.0667e-3 - 0.4156/(temp.kelvin[i]) + 1.137e3/((temp.kelvin[i])^2))
    alpha17_eq_l_v[i] = alpha18_eq_l_v[i] ^ theta_eq 
    alpha2_eq_l_v[i] = exp(52.612e-3 - 76.248/(temp.kelvin[i]) + (24.844e3)/(temp.kelvin[i])^2)
    
    # meteoric water
    R18_w[i] = (d18_w[i] / 1e3 + 1) * R18_smow
    dp18_w[i] = 1e3 * log(R18_w[i] / R18_smow)
    dp17_w[i] = 0.5268 * dp18_w[i] + 0.015 # Global meteoric water - Aron (2021)
    R17_w[i] = exp(dp17_w[i] / 1e3) * R17_smow
    R2_w[i] = (d2_w[i] / 1e3 + 1) * R2_smow
    
    # assume local water vapor is predominantly evaporated - different RH and temp
    # R18_a[i] = R18_w[i] / (alpha18_eq_l_v[i] * (alpha18_diff * (1 - RH[i]) + RH[i]))
    # R17_a[i] = R17_w[i] / (alpha17_eq_l_v[i] * (alpha17_diff * (1 - RH[i]) + RH[i]))
    # R2_a[i] = R2_w[i] / (alpha2_eq_l_v[i] * (alpha2_diff * (1 - RH[i]) + RH[i]))
    # assume local water vapor is in isotope equilibrium with meteoric water
    R18_a[i] = R18_w[i] / alpha18_eq_l_v[i]
    R17_a[i] = R17_w[i] / alpha17_eq_l_v[i]
    R2_a[i] = R2_w[i] / alpha2_eq_l_v[i]

    # influx
    R18_in[i] = R18_w[i] 
    R17_in[i] = R17_w[i] 
    R2_in[i] = R2_w[i] 
    
    d18_in[i] = (R18_in[i] / R18_smow - 1) * 1e3
    d2_in[i] = (R2_in[i] / R2_smow - 1) * 1e3
    
    # assuming steady state
    A18[i] = (alpha18_diff * alpha18_eq_l_v[i] * (1 - RH[i]) * R18_in[i] + alpha18_eq_l_v[i] * R18_a[i] * RH[i] * (1 - theta[i])) / (alpha18_diff * alpha18_eq_l_v[i] * (1 - RH[i]) * (1 - theta[i])) * (F_evap / M_bw)
    B18[i] = ((1 - theta[i]) + alpha18_diff * alpha18_eq_l_v[i] * (1 - RH[i]) * theta[i]) / (alpha18_diff * alpha18_eq_l_v[i] * (1 - RH[i]) * (1 - theta[i])) * (F_evap / M_bw)
    A17[i] = (alpha17_diff * alpha17_eq_l_v[i] * (1 - RH[i]) * R17_in[i] + alpha17_eq_l_v[i] * R17_a[i] * RH[i] * (1 - theta[i])) / (alpha17_diff * alpha17_eq_l_v[i] * (1 - RH[i]) * (1 - theta[i])) * (F_evap / M_bw)
    B17[i] = ((1- theta[i]) + alpha17_diff * alpha17_eq_l_v[i] * (1 - RH[i]) * theta[i]) / (alpha17_diff * alpha17_eq_l_v[i] * (1 - RH[i]) * (1 - theta[i])) * (F_evap / M_bw)
    A2[i] = (alpha2_diff * alpha2_eq_l_v[i] * (1 - RH[i]) * R2_in[i] + alpha2_eq_l_v[i] * R2_a[i] * RH[i] * (1 - theta[i])) / (alpha2_diff * alpha2_eq_l_v[i] * (1 - RH[i]) * (1 - theta[i])) * (F_evap / M_bw)
    B2[i] = ((1- theta[i]) + alpha2_diff * alpha2_eq_l_v[i] * (1 - RH[i]) * theta[i]) / (alpha2_diff * alpha2_eq_l_v[i] * (1 - RH[i]) * (1 - theta[i])) * (F_evap / M_bw)
    
    tau18[i] = 1 / B18[i] # time constant for 63% approach to steady state
    tau17[i] = 1 / B17[i]
    tau2[i] = 1 / B2[i]
    
    R18_ss[i] = A18[i] / B18[i]
    R17_ss[i] = A17[i] / B17[i]
    R2_ss[i] = A2[i] / B2[i]
    
    # R18_l[i] = (R18_in[i] - R18_ss[i]) * exp(-time[i] / tau18[i]) + R18_ss[i]
    # R17_l[i] = (R17_in[i] - R17_ss[i]) * exp(-time[i] / tau17[i]) + R17_ss[i]
    # R2_l[i] = (R2_in[i] - R2_ss[i]) * exp(-time[i] / tau2[i]) + R2_ss[i]
    
    d18_bw[i] = (R18_ss[i] / R18_smow - 1) * 1e3 # steady state d18O of body fluid
    d17_bw[i] = (R17_ss[i] / R17_smow - 1) * 1e3
    Dp17_bw[i] = (1e3 * log(R17_ss[i] / R17_smow) - .528 * 1e3 * log(R18_ss[i] / R18_smow)) * 1e3 # per meg
    d2_bw[i] = (R2_ss[i] / R2_smow - 1) * 1e3
  }
  
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
  alpha17_diff = alpha18_diff^theta_diff
  alpha2_diff = 1/0.9755
  alpha13_eq_a_b = 1.0027 # equilibrium isotopic fractionation factor between aragonite and HCO3- at 10-25 degrees 
  alpha13_diff = 1.0044 # kinetic fractionation factor between 12CO2 / 13CO2 - Cerling and Quade (1993)
  M_bw = 2e-2 # water molar mass per wet body weight - unit: mol/g
  F_evap = 2.7e-6 # evaporation rate per wet body weight - unit: mol/g/s (Dainton, 1954)
}