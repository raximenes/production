###############################################################################
# PSA Analysis: Reading Parameters from Excel, Constructing Base-Case Params,
# Generating PSA Samples, and Combining Them into a List of Simulations
#
# Files used:
#   - data/IMD Data.xls (main parameters)
#   - data/confidential/vaccine_costs.xlsx (vaccine costs)
###############################################################################

#----------------------------------------------------------------------------
# 1) HELPER FUNCTIONS
#----------------------------------------------------------------------------

# calculate_all_ve: Computes a linear decay from the base effectiveness (starting at cycle 1)
# to 0 over max_year cycles.
calculate_all_ve <- function(data, max_year = 10, n_cycles = 89) {
  years <- 1:n_cycles
  out_list <- lapply(data$vaccine, function(vac) {
    base_ve <- data[data$vaccine == vac, "effectiveness"]
    ve_vec <- sapply(years, function(yr) {
      if (yr >= 1 && yr < (max_year + 1)) {
        x0 <- 1
        base_ve * (1 - (yr - x0) / max_year)
      } else {
        0
      }
    })
    data.frame(Cycle = years, Vaccine = vac, Effectiveness = ve_vec)
  })
  do.call(rbind, out_list)
}

# calculate_all_ve_shifted: For a delayed single vaccination (vaccination_time > 0),
# effectiveness remains 0 until cycle (shift+1), then decays linearly over max_year cycles.
calculate_all_ve_shifted <- function(data, max_year = 10, n_cycles = 89, shift = 0) {
  years <- 1:n_cycles
  out_list <- lapply(data$vaccine, function(vac) {
    base_ve <- data[data$vaccine == vac, "effectiveness"]
    ve_vec <- sapply(years, function(yr) {
      if (yr <= shift) {
        0
      } else if (yr > shift && yr < (shift + max_year + 1)) {
        x0 <- shift + 1
        base_ve * (1 - (yr - x0) / max_year)
      } else {
        0
      }
    })
    data.frame(Cycle = years, Vaccine = vac, Effectiveness = ve_vec)
  })
  do.call(rbind, out_list)
}

# calculate_boosted_ve: For multiple vaccination events (e.g. a booster),
# computes, for each cycle, the maximum effectiveness contribution from any event.
calculate_boosted_ve <- function(base_ve, vac_times, doses, max_year, n_cycles, coverage) {
  ve_vector <- numeric(n_cycles)
  for (t in 1:n_cycles) {
    contributions <- sapply(seq_along(vac_times), function(i) {
      vt <- vac_times[i]
      if (t >= (vt + 1) && t < (vt + max_year + 1)) {
        # Note: The dose multiplier is applied in cost calculations;
        # here we assume the base VE is already defined per dose.
        base_ve * (1 - (t - (vt + 1)) / max_year)
      } else {
        0
      }
    })
    ve_vector[t] <- max(contributions)
  }
  return(ve_vector * coverage)
}

#----------------------------------------------------------------------------
# 1-A) Define scenarios
#----------------------------------------------------------------------------
define_scenarios <- function() {
  scenarios <- list(
    "Scenario A" = list(
      vaccination_time = 12,    # Vaccination at cycle 12
      doses = 1,                # Single dose
      n_cycles = 89
    ),
    "Scenario B" = list(
      vaccination_time = 0,   # Vaccination at cycle 0
      doses = 2,                # Double dose
      n_cycles = 100
    ),
    "Scenario C" = list(
      vaccination_time = c(0,12),  # Vaccination at cycles 0 and 12
      doses = c(2,1),               # Two doses initially; one booster at cycle 12
      n_cycles = 100
    )
  )
  return(scenarios)
}

#----------------------------------------------------------------------------
# 2) READING DATA AND BUILDING BASE-CASE PARAMETERS
#----------------------------------------------------------------------------

# Define file paths.
file_path       <- "data/IMD Data complete.xls"
file_path_costs <- "data/confidential/vaccine_costs.xlsx"

# get_basecase_params: Reads Excel files and merges with scenario-specific settings.
get_basecase_params <- function(scenario_name = "Scenario A") {
  scenarios <- define_scenarios()
  if (!scenario_name %in% names(scenarios))
    stop("Invalid scenario name. Please choose one of: ", paste(names(scenarios), collapse = ", "))
  scenario_params <- scenarios[[scenario_name]]
  
  # sheet_names <- excel_sheets(file_path)
  # data_list <- lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet))
  # names(data_list) <- sheet_names
  
  sheet_names <- excel_sheets(file_path)
  data_list <- lapply(sheet_names, function(sheet) {
    if (scenario_name == "Scenario A" && sheet %in% c("infection", "mortality", "cost_IMD")) {
      # Read only the first row to get the header
      header <- suppressMessages(read_excel(file_path, sheet = sheet, n_max = 1, col_names = TRUE))
      # Read the data starting from row 13 (i.e., skipping rows 2 to 12), without headers
      df <- suppressMessages(read_excel(file_path, sheet = sheet, skip = 12, col_names = FALSE))
      # Assign the header names to the data frame
      names(df) <- names(header)
      return(df)
    } else {
      read_excel(file_path, sheet = sheet)
    }
  })
  names(data_list) <- sheet_names
  
  vaccine_costs <- read_excel(file_path_costs)
  
  n_cycles <- scenario_params$n_cycles # with(data_list$model_settings, Value[Name == "n_cycles"])
  cycle_length <- with(data_list$model_settings, Value[Name == "cycle_len"])
  d_c <- with(data_list$model_settings, Value[Name == "d_c"])
  d_e <- with(data_list$model_settings, Value[Name == "d_e"])
  coverage <- with(data_list$model_settings, Value[Name == "coverage"])
  
  c_admin <- with(data_list$costs, Value[Name == "c_admin"])
  c_MenABCWY <- with(vaccine_costs, Value[Name == "c_MenABCWY"])
  c_MenACWY <- with(vaccine_costs, Value[Name == "c_MenACWY"])
  c_MenC <- with(vaccine_costs, Value[Name == "c_MenC"])
  
  p_bg_mort <- rate_to_prob(with(data_list$mortality, Background_Mortality), t = cycle_length)
  p_B <- rate_to_prob(with(data_list$infection, SerogroupB_Infection), t = cycle_length)
  p_C <- rate_to_prob(with(data_list$infection, SerogroupC_Infection), t = cycle_length)
  p_W <- rate_to_prob(with(data_list$infection, SerogroupW_Infection), t = cycle_length)
  p_Y <- rate_to_prob(with(data_list$infection, SerogroupY_Infection), t = cycle_length)
  
  p_B_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupB_Dead), t = cycle_length)
  p_C_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupC_Dead), t = cycle_length)
  p_W_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupW_Dead), t = cycle_length)
  p_Y_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupY_Dead), t = cycle_length)
  
  p_IMD_Scarring <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Scarring"])
  p_IMD_Single_Amput <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Single_Amput"])
  p_IMD_Multiple_Amput <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Multiple_Amput"])
  p_IMD_Neuro_Disability <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Neuro_Disability"])
  p_IMD_Hearing_Loss <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Hearing_Loss"])
  p_IMD_Renal_Failure <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Renal_Failure"])
  p_IMD_Seizure <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Seizure"])
  p_IMD_Paralysis <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Paralysis"])
  
  ve_data <- data.frame(
    vaccine = c("MenABCWY_for_SeroACWY", "MenABCWY_for_SeroB", "MenACWY", "MenC"),
    effectiveness = with(data_list$vac_effect,
                         Value[Name %in% c("MenABCWY_for_SeroACWY", 
                                           "MenABCWY_for_SeroB", 
                                           "MenACWY", 
                                           "MenC")])
  )
  
  # Calculate VE according to the scenario.
  if (length(scenario_params$vaccination_time) == 1) {
    vt <- scenario_params$vaccination_time
    if (vt == 0) {
      ve_all <- calculate_all_ve(ve_data, max_year = 10, n_cycles = n_cycles)
    } else {
      ve_all <- calculate_all_ve_shifted(ve_data, max_year = 10, n_cycles = n_cycles, shift = vt)
    }
    ve_MenABCWY_forACWY <- ve_all[ve_all$Vaccine == "MenABCWY_for_SeroACWY", "Effectiveness"] * coverage
    ve_MenABCWY_forB    <- ve_all[ve_all$Vaccine == "MenABCWY_for_SeroB",    "Effectiveness"] * coverage 
    ve_MenACWY          <- ve_all[ve_all$Vaccine == "MenACWY",               "Effectiveness"] * coverage
    ve_MenC             <- ve_all[ve_all$Vaccine == "MenC",                  "Effectiveness"] * coverage 
  } else {
    ve_MenABCWY_forACWY <- calculate_boosted_ve(
      base_ve = ve_data$effectiveness[ve_data$vaccine == "MenABCWY_for_SeroACWY"],
      vac_times = scenario_params$vaccination_time,
      doses = scenario_params$doses,
      max_year = 10,
      n_cycles = n_cycles,
      coverage = coverage
    )
    ve_MenABCWY_forB <- calculate_boosted_ve(
      base_ve = ve_data$effectiveness[ve_data$vaccine == "MenABCWY_for_SeroB"],
      vac_times = scenario_params$vaccination_time,
      doses = scenario_params$doses,
      max_year = 10,
      n_cycles = n_cycles,
      coverage = coverage
    )
    ve_MenACWY <- calculate_boosted_ve(
      base_ve = ve_data$effectiveness[ve_data$vaccine == "MenACWY"],
      vac_times = scenario_params$vaccination_time,
      doses = scenario_params$doses,
      max_year = 10,
      n_cycles = n_cycles,
      coverage = coverage
    )
    ve_MenC <- calculate_boosted_ve(
      base_ve = ve_data$effectiveness[ve_data$vaccine == "MenC"],
      vac_times = scenario_params$vaccination_time,
      doses = scenario_params$doses,
      max_year = 10,
      n_cycles = n_cycles,
      coverage = coverage
    )
  }
  
  # Extract health state costs.
  if ("cost_IMD" %in% names(data_list)) {
    c_IMD_infection <- data_list$cost_IMD$cost
  } else {
    c_IMD_infection <- with(data_list$costs, Value[Name == "c_IMD_infection"])
  }
  c_Scarring <- with(data_list$costs, Value[Name == "c_Scarring"])
  c_Single_Amput <- with(data_list$costs, Value[Name == "c_Single_Amput"])
  c_Multiple_Amput <- with(data_list$costs, Value[Name == "c_Multiple_Amput"])
  c_Neuro_Disab <- with(data_list$costs, Value[Name == "c_Neuro_Disab"])
  c_Hearing_Loss <- with(data_list$costs, Value[Name == "c_Hearing_Loss"])
  c_Renal_Failure <- with(data_list$costs, Value[Name == "c_Renal_Failure"])
  c_Seizure <- with(data_list$costs, Value[Name == "c_Seizure"])
  c_Paralysis <- with(data_list$costs, Value[Name == "c_Paralysis"])
  c_Dead <- with(data_list$costs, Value[Name == "c_Dead"])
  c_Healthy <- with(data_list$costs, Value[Name == "c_Healthy"])
  
  # Extract health state utilities.
  u_Healthy <- with(data_list$utilities, Value[Name == "u_Healthy"])
  u_IMD <- with(data_list$utilities, Value[Name == "u_IMD"])
  u_Scarring <- with(data_list$utilities, Value[Name == "u_Scarring"])
  u_Single_Amput <- with(data_list$utilities, Value[Name == "u_Single_Amput"])
  u_Multiple_Amput <- with(data_list$utilities, Value[Name == "u_Multiple_Amput"])
  u_Neuro_Disability <- with(data_list$utilities, Value[Name == "u_Neuro_Disability"])
  u_Hearing_Loss <- with(data_list$utilities, Value[Name == "u_Hearing_Loss"])
  u_Renal_Failure <- with(data_list$utilities, Value[Name == "u_Renal_Failure"])
  u_Seizure <- with(data_list$utilities, Value[Name == "u_Seizure"])
  u_Paralysis <- with(data_list$utilities, Value[Name == "u_Paralysis"])
  u_Dead <- with(data_list$utilities, Value[Name == "u_Dead"])
  
  params_base <- list(
    n_cycles = n_cycles,
    cycle_length = cycle_length,
    d_c = d_c,
    d_e = d_e,
    coverage = coverage,
    c_MenABCWY = c_MenABCWY,
    c_MenACWY = c_MenACWY,
    c_MenC = c_MenC,
    c_admin = c_admin,
    p_bg_mort = p_bg_mort,
    p_B = p_B,
    p_C = p_C,
    p_W = p_W,
    p_Y = p_Y,
    p_B_DeadIMD = p_B_DeadIMD,
    p_C_DeadIMD = p_C_DeadIMD,
    p_W_DeadIMD = p_W_DeadIMD,
    p_Y_DeadIMD = p_Y_DeadIMD,
    p_IMD_Scarring = p_IMD_Scarring,
    p_IMD_Single_Amput = p_IMD_Single_Amput,
    p_IMD_Multiple_Amput = p_IMD_Multiple_Amput,
    p_IMD_Neuro_Disability = p_IMD_Neuro_Disability,
    p_IMD_Hearing_Loss = p_IMD_Hearing_Loss,
    p_IMD_Renal_Failure = p_IMD_Renal_Failure,
    p_IMD_Seizure = p_IMD_Seizure,
    p_IMD_Paralysis = p_IMD_Paralysis,
    ve_MenABCWY_forACWY = ve_MenABCWY_forACWY,
    ve_MenABCWY_forB = ve_MenABCWY_forB,
    ve_MenACWY = ve_MenACWY,
    ve_MenC = ve_MenC,
    c_IMD_infection = c_IMD_infection,
    c_Scarring = c_Scarring,
    c_Single_Amput = c_Single_Amput,
    c_Multiple_Amput = c_Multiple_Amput,
    c_Neuro_Disab = c_Neuro_Disab,
    c_Hearing_Loss = c_Hearing_Loss,
    c_Renal_Failure = c_Renal_Failure,
    c_Seizure = c_Seizure,
    c_Paralysis = c_Paralysis,
    c_Dead = c_Dead,
    c_Healthy = c_Healthy,
    u_Healthy = u_Healthy,
    u_IMD = u_IMD,
    u_Scarring = u_Scarring,
    u_Single_Amput = u_Single_Amput,
    u_Multiple_Amput = u_Multiple_Amput,
    u_Neuro_Disability = u_Neuro_Disability,
    u_Hearing_Loss = u_Hearing_Loss,
    u_Renal_Failure = u_Renal_Failure,
    u_Seizure = u_Seizure,
    u_Paralysis = u_Paralysis,
    u_Dead = u_Dead
  )
  
  params <- c(params_base, scenario_params)
  return(params)
}

#----------------------------------------------------------------------------
# 3) BUILDING THE STANDARD DEVIATION (SD) VALUES LIST USING THE BASE-CASE STRUCTURE
#----------------------------------------------------------------------------
get_basecase_sd_values <- function() {
  sheet_names <- excel_sheets(file_path)
  data_list <- lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet))
  names(data_list) <- sheet_names
  vaccine_costs <- read_excel(file_path_costs)
  
  extract_sd <- function(df) {
    if ("sd" %in% names(df)) as.list(setNames(as.numeric(df$sd), df$Name)) else list()
  }
  
  sd_costs         <- extract_sd(data_list$costs)
  sd_vaccine_costs <- extract_sd(vaccine_costs)
  sd_utilities     <- extract_sd(data_list$utilities)
  sd_sequelae      <- extract_sd(data_list$sequelae_probp_IMD)
  sd_vac_effect    <- extract_sd(data_list$vac_effect)
  
  sd_values <- c(sd_costs, sd_vaccine_costs, sd_utilities, sd_sequelae, sd_vac_effect)
  
  if ("cost_IMD" %in% names(data_list))
    sd_values[["c_IMD_infection"]] <- as.numeric(data_list$cost_IMD$sd)
  
  if (!is.null(data_list$mortality)) {
    mortality_sd_cols <- grep("^sd_", names(data_list$mortality), value = TRUE)
    for (col in mortality_sd_cols) {
      param_name <- sub("^sd_", "p_", col)
      sd_values[[param_name]] <- as.numeric(data_list$mortality[[col]])
    }
  }
  
  if (!is.null(data_list$infection)) {
    infection_sd_cols <- grep("^sd_", names(data_list$infection), value = TRUE)
    for (col in infection_sd_cols) {
      param_name <- sub("^sd_", "p_", col)
      sd_values[[param_name]] <- as.numeric(data_list$infection[[col]])
    }
  }
  
  cycle_length <- with(data_list$model_settings, Value[Name == "cycle_len"])
  if ("sd_Background_Mortality" %in% names(data_list$mortality)) {
    sd_values$p_bg_mort <- rate_to_prob(as.numeric(data_list$mortality$sd_Background_Mortality), t = cycle_length)
  }
  
  sd_values <- sd_values[!sapply(sd_values, function(x) any(is.na(x)))]
  
  rename_map <- c(
    "p_SerogroupB_Dead"      = "p_B_DeadIMD",
    "p_SerogroupC_Dead"      = "p_C_DeadIMD",
    "p_SerogroupW_Dead"      = "p_W_DeadIMD",
    "p_SerogroupY_Dead"      = "p_Y_DeadIMD",
    "p_Background_Mortality" = "p_bg_mort",
    "p_SerogroupB_Infection" = "p_B",
    "p_SerogroupC_Infection" = "p_C",
    "p_SerogroupW_Infection" = "p_W",
    "p_SerogroupY_Infection" = "p_Y",
    "MenABCWY_for_SeroACWY"  = "ve_MenABCWY_forACWY",
    "MenABCWY_for_SeroB"     = "ve_MenABCWY_forB",
    "MenACWY"                = "ve_MenACWY",
    "MenC"                   = "ve_MenC"
  )
  
  names(sd_values) <- sapply(names(sd_values), function(name) {
    if (name %in% names(rename_map)) rename_map[[name]] else name
  })
  
  coverage <- with(data_list$model_settings, Value[Name == "coverage"])
  vaccine_names <- c("ve_MenABCWY_forACWY", "ve_MenABCWY_forB", "ve_MenACWY", "ve_MenC")
  ve_data <- data.frame(
    vaccine = vaccine_names,
    effectiveness = sapply(vaccine_names, function(vac) {
      if (!is.null(sd_values[[vac]])) sd_values[[vac]] else NA
    })
  )
  # Compute deterministic VE curves.
  ve_all <- calculate_all_ve(ve_data, max_year = 10, n_cycles = 89)
  # For each VE parameter, if no SD was provided, do not include it in sampling.
  for (vac in vaccine_names) {
    if (is.null(sd_values[[vac]])) {
      sd_values[[vac]] <- NULL
    }
    # Otherwise, leave the provided SD value unchanged.
  }
  
  return(sd_values)
}

sd_values <- get_basecase_sd_values()

#----------------------------------------------------------------------------
# 4) PSA SAMPLE GENERATION
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# 1) Define scenarios
#----------------------------------------------------------------------------
# define_scenarios <- function() {
#   scenarios <- list(
#     "Scenario A" = list(
#       vaccination_time = 0,    # Vacinação no ciclo 0
#       doses = 1                # Single dose
#     ),
#     "Scenario B" = list(
#       vaccination_time = 12,   # Vacinação no ciclo 12
#       doses = 1                # Single dose
#     ),
#     "Scenario C" = list(
#       vaccination_time = c(0, 12),  # Vacinação nos ciclos 0 e 12
#       doses = c(2, 1)               # Duas doses inicialmente; uma dose de reforço no ciclo 12
#     )
#   )
#   return(scenarios)
# }

#----------------------------------------------------------------------------
# 2) Generate PSA samples for a specific scenario
#----------------------------------------------------------------------------
# generate_psa_samples <- function(nsamp, scenario_name = "Scenario A") {
#   # 1) Retrieve the base-case parameters for the specified scenario
#   basecase_params <- get_basecase_params(scenario_name)
#   
#   # 2) Filter only those parameters for which SD values are available
#   filtered_list_for_psa <- basecase_params[names(basecase_params) %in% names(sd_values)]
#   
#   # 3) Separate zero and non-zero parameters
#   params_for_psa <- unlist(filtered_list_for_psa)
#   non_zero_params <- params_for_psa[sapply(params_for_psa, function(x) !all(x == 0))]
#   zero_params <- params_for_psa[sapply(params_for_psa, function(x) all(x == 0))]
#   
#   # 4) Flatten the filtered list into a named vector
#   base_params <- non_zero_params
#   my_params <- names(base_params)
#   
#   # 5) Define distribution types based on parameter names
#   my_dists <- sapply(my_params, function(param) {
#     if (grepl("^c_", param)) {
#       "gamma"
#     } else if (grepl("^(u_|p_|ve_)", param)) {
#       "beta"
#     } else {
#       NA
#     }
#   })
#   
#   # 6) Define parameterization types for each distribution
#   my_parameterization_types <- sapply(my_dists, function(dist) {
#     if (is.na(dist)) {
#       NA
#     } else if (dist == "beta") {
#       "mean, sd"
#     } else if (dist == "gamma") {
#       "mean, sd"
#     } else {
#       NA
#     }
#   })
#   
#   # 7) Create a list pairing each parameter's mean and SD
#   unl_sd_val <- unlist(sd_values)
#   mean_sd_list <- list()
#   for (param in names(base_params)) {
#     if (param %in% names(unl_sd_val)) {
#       mean_sd_list[[param]] <- c(mean = base_params[[param]], sd = unl_sd_val[[param]])
#     }
#   }
#   
#   # 8) Generate distribution parameters for use with gen_psa_samp
#   generate_my_dists_params <- function(mean_sd_list, my_parameterization_types) {
#     dists_params <- list()
#     for (param in names(mean_sd_list)) {
#       if (param %in% names(my_parameterization_types)) {
#         if (my_parameterization_types[[param]] == "mean, sd") {
#           # For Beta distribution, use the mean and SD directly
#           dists_params[[param]] <- mean_sd_list[[param]]
#         } else {
#           NA
#         }
#       }
#     }
#     return(dists_params)
#   }
#   
#   my_dists_params <- generate_my_dists_params(mean_sd_list, my_parameterization_types)
#   
#   # 9) Remove parameters that have NA for distribution types
#   valid_params_idx <- !is.na(my_parameterization_types)
#   valid_params <- my_params[valid_params_idx]
#   valid_dists <- my_dists[valid_params_idx]
#   valid_param_types <- my_parameterization_types[valid_params_idx]
#   
#   # 10) Generate PSA samples using dampack::gen_psa_samp
#   df_psa <- gen_psa_samp(
#     params = valid_params,
#     dists = valid_dists,
#     parameterization_types = valid_param_types,
#     dists_params = my_dists_params,
#     nsamp = nsamp
#   )
#   
#   # 11) Add zero parameters back to the data frame
#   if (length(zero_params) > 0) {
#     for (param in names(zero_params)) {
#       df_psa[[param]] <- rep(0, nrow(df_psa))
#     }
#   }
#   
#   # 12) Convert the PSA data frame into a list of parameter sets
#   scalar_params <- names(filtered_list_for_psa)[sapply(filtered_list_for_psa, length) == 1]
#   vector_params <- names(filtered_list_for_psa)[sapply(filtered_list_for_psa, length) > 1]
#   
#   vector_param_cols <- map(vector_params, function(param) {
#     pattern <- paste0("^", param, "\\d+$")  # e.g., "p_bg_mort" followed by digits
#     cols <- grep(pattern, names(df_psa), value = TRUE)
#     cols_sorted <- cols[order(as.numeric(str_extract(cols, "\\d+$")))]
#     return(cols_sorted)
#   })
#   names(vector_param_cols) <- vector_params
#   
#   samples_list <- df_psa %>%
#     select(-nsamp) %>%
#     split(1:nrow(.)) %>%
#     map(function(row) {
#       sample_list <- list()
#       
#       # Scalar parameters: use the drawn value if available; otherwise, revert to base-case
#       for (param in scalar_params) {
#         sample_list[[param]] <- if (!is.null(row[[param]])) row[[param]] else filtered_list_for_psa[[param]]
#       }
#       
#       # Vector parameters: combine the corresponding columns into a numeric vector
#       for (param in vector_params) {
#         cols <- vector_param_cols[[param]]
#         if (length(cols) > 0) {
#           sample_list[[param]] <- as.numeric(row[cols])
#         } else {
#           sample_list[[param]] <- filtered_list_for_psa[[param]]
#         }
#       }
#       
#       # Add any missing parameters (with their base-case values) to each sample
#       missing_params <- setdiff(names(basecase_params), names(sample_list))
#       for (param in missing_params) {
#         sample_list[[param]] <- basecase_params[[param]]
#       }
#       
#       return(sample_list)
#     })
#   
#   return(samples_list)
# }

# ----------------------------------------------------------------------------
# Function: generate_psa_samples
#
# This function generates probabilistic sensitivity analysis (PSA) samples.
# It uses the base-case parameters for a given scenario (controlled by the 
# "scenario_name" argument) and the corresponding standard deviations (SDs)
# to sample from defined distributions.
#
# For cost parameters (starting with "c_") a Gamma distribution is used,
# while probabilities, utilities, and vaccine effectiveness (starting with "p_", "u_", or "ve_")
# use a Beta distribution.
#
# The function returns a list where each element is a parameter set corresponding
# to one PSA sample.
#
# Arguments:
#   nsamp         - The number of PSA samples to generate.
#   scenario_name - (Optional) The scenario to use ("Scenario A", "Scenario B", etc.).
#                   Defaults to "Scenario A".
# ----------------------------------------------------------------------------
generate_psa_samples <- function(nsamp, scenario_name = "Scenario A") {
  
  # 1) Retrieve base-case parameters for the selected scenario.
  basecase_params <- get_basecase_params(scenario_name)
  
  # 2) Filter parameters that have an associated SD value (from the global sd_values).
  #    This step excludes scenario-specific parameters like vaccination_time and doses.
  filtered_params <- basecase_params[names(basecase_params) %in% names(sd_values)]
  
  # 3) Flatten the filtered list into a named vector.
  params_vector <- unlist(filtered_params)
  
  # Separate non-zero and zero parameters.
  non_zero_params <- params_vector[sapply(params_vector, function(x) !all(x == 0))]
  zero_params     <- params_vector[sapply(params_vector, function(x) all(x == 0))]
  
  # 4) Define the distribution type for each parameter based on its name.
  #    - Cost parameters (starting with "c_") use the Gamma distribution.
  #    - Probabilities, utilities, and vaccine effectiveness (starting with "p_", "u_", or "ve_")
  #      use the Beta distribution.
  parameter_names <- names(non_zero_params)
  dist_types <- sapply(parameter_names, function(param) {
    if (grepl("^c_", param)) {
      "gamma"
    } else if (grepl("^(p_|u_|ve_)", param)) {
      "beta"
    } else {
      NA
    }
  })
  
  # 5) Define parameterization types.
  #     Here we assume both Beta and Gamma distributions are parameterized by "mean, sd".
  param_types <- sapply(dist_types, function(dist) {
    if (is.na(dist)) {
      NA
    } else {
      "mean, sd"
    }
  })
  
  # 6) Create a list of mean and SD pairs for each parameter.
  #    Only include those parameters that are available in the global sd_values.
  sd_vector <- unlist(sd_values)
  mean_sd_list <- list()
  for (param in parameter_names) {
    # Remove trailing digits to get the base name.
    base_name <- sub("\\d+$", "", param)
    if (base_name %in% names(sd_vector)) {
      mean_sd_list[[param]] <- c(mean = non_zero_params[[param]], sd = sd_vector[[base_name]])
    }
  }
  
  # 7) Filter valid parameters: only those for which we have both a mean and an SD.
  valid_params <- names(mean_sd_list)
  valid_dist_types <- dist_types[valid_params]
  valid_param_types <- param_types[valid_params]
  
  # 8) Generate distribution parameters for each valid parameter.
  #     For both Beta and Gamma we assume a "mean, sd" parameterization.
  generate_dist_params <- function(mean_sd_list, param_types) {
    dist_params <- list()
    for (param in names(mean_sd_list)) {
      if (!is.na(param_types[[param]]) && param_types[[param]] == "mean, sd") {
        dist_params[[param]] <- mean_sd_list[[param]]
      }
    }
    return(dist_params)
  }
  
  dist_params <- generate_dist_params(mean_sd_list, param_types)
  
  # 9) Generate the PSA samples using the sampling function (e.g., gen_psa_samp from dampack).
  df_psa <- gen_psa_samp(
    params = valid_params,
    dists = valid_dist_types,
    parameterization_types = valid_param_types,
    dists_params = dist_params,
    nsamp = nsamp
  )
  
  # 10) For parameters that had a base-case value of zero, add them back to the data frame.
  if (length(zero_params) > 0) {
    for (param in names(zero_params)) {
      df_psa[[param]] <- rep(0, nrow(df_psa))
    }
  }
  
  # 11) Separate scalar and vector parameters based on their length in the filtered parameters.
  scalar_params <- names(filtered_params)[sapply(filtered_params, length) == 1]
  vector_params <- names(filtered_params)[sapply(filtered_params, length) > 1]
  
  # 12) Identify the columns in the PSA data frame corresponding to each vector parameter.
  #     For instance, vector parameters might appear as p_bg_mort1, p_bg_mort2, etc.
  vector_param_cols <- lapply(vector_params, function(param) {
    pattern <- paste0("^", param, "\\d+$")
    cols <- grep(pattern, names(df_psa), value = TRUE)
    # Sort columns by their numeric suffix to preserve order.
    cols[order(as.numeric(gsub(".*?(\\d+)$", "\\1", cols)))]
  })
  names(vector_param_cols) <- vector_params
  
  # 13) Convert the PSA data frame into a list of parameter sets (one per PSA draw).
  samples_list <- split(df_psa, seq(nrow(df_psa)))
  samples_list <- lapply(samples_list, function(row_df) {
    sample_list <- list()
    
    # For scalar parameters, use the drawn value (or revert to the base-case if missing).
    for (param in scalar_params) {
      sample_list[[param]] <- if (!is.null(row_df[[param]])) row_df[[param]] else filtered_params[[param]]
    }
    
    # For vector parameters, combine the corresponding columns into a numeric vector.
    for (param in vector_params) {
      cols <- vector_param_cols[[param]]
      if (length(cols) > 0) {
        sample_list[[param]] <- as.numeric(row_df[cols])
      } else {
        sample_list[[param]] <- filtered_params[[param]]
      }
    }
    return(sample_list)
  })
  
  # 14) Add any missing global parameters from the base-case (if they were not sampled).
  existing_psa_params <- unique(unlist(lapply(samples_list, names)))
  all_params <- names(basecase_params)
  missing_params <- setdiff(all_params, existing_psa_params)
  if (length(missing_params) > 0) {
    samples_list <- lapply(samples_list, function(sample) {
      for (param in missing_params) {
        sample[[param]] <- basecase_params[[param]]
      }
      return(sample)
    })
  }
  
  # 15) Return the final list of PSA sample parameter sets.
  return(samples_list)
}


###############################################################################
# EXAMPLE USAGE:
#
# 1) Generate 3 PSA samples for a specified scenario (e.g., "Scenario A", "Scenario B", or "Scenario C"):
#      samples <- generate_psa_samples(nsamp = 10, scenario_name = "Scenario A")
#
# 2) Inspect one sample set:
#      str(samples[[1]])
#
# 3) Check for any NA values:
#      na_parameters_per_sim <- lapply(samples, function(sim) {
#         na_flags <- sapply(sim, function(param) any(is.na(param)))
#         names(na_flags[na_flags])
#      })
#      na_parameters_per_sim
###############################################################################
