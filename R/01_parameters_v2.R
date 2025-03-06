###############################################################################
# PSA Analysis: Reading Parameters from Excel, Constructing Base-Case Params,
# Generating PSA Samples, and Combining Them into a List of Simulations
#
# Files used:
#   - data/IMD Data complete.xls (main parameters)
#   - data/confidential/vaccine_costs.xlsx (vaccine costs)
###############################################################################

#----------------------------------------------------------------------------
# 1) HELPER FUNCTIONS FOR CALCULATING VACCINE EFFECTIVENESS (VE)
#----------------------------------------------------------------------------

# calculate_all_ve: Computes a linear decay from the base effectiveness 
# (starting at cycle 1) to 0 over max_year cycles.
calculate_all_ve <- function(data, max_year = 10, n_cycles = n_cycles) {
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

# calculate_all_ve_shifted: For a delayed vaccination (vaccination_time > 0),
# effectiveness remains 0 until cycle (shift+1), then decays linearly over max_year cycles.
calculate_all_ve_shifted <- function(data, max_year = 10, n_cycles = n_cycles, shift = 0) {
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

# calculate_boosted_ve: For multiple vaccination events (e.g., a booster),
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
# 1-A) DEFINE SCENARIOS (INCLUDING NEW FIELDS FOR INITIAL AND BOOSTER VACCINES)
#----------------------------------------------------------------------------
define_scenarios <- function() {
  scenarios <- list(
    "Scenario A" = list(
      vaccination_time = 12,    # Vaccination at cycle 12
      doses = 1,                # Single dose
      n_cycles = 89,
      allowed_vaccines = c("MenABCWY_for_SeroACWY", "MenACWY")
    ),
    "Scenario B" = list(
      vaccination_time = 0,     # Vaccination at cycle 0
      doses = 2,                # Double dose
      n_cycles = 100,
      allowed_vaccines = c("MenACWY", "MenC")
    ),
    "Scenario C" = list(
      vaccination_time = c(0, 12),  # Vaccination at cycles 0 and 12
      doses = c(2, 1),               # Two doses initially; one booster at cycle 12
      n_cycles = 100,
      allowed_vaccines_init = c("MenACWY", "MenC"),       # Vaccines for initial vaccination
      allowed_vaccines_boost = c("MenABCWY_for_SeroACWY")  # Booster vaccine
    ),
    "Scenario D" = list(
      vaccination_time = c(0, 12),  # Vaccination at cycles 0 and 12
      doses = c(2, 1),               # Two doses initially; one booster at cycle 12
      n_cycles = 100,
      allowed_vaccines_init = c("MenACWY", "MenC"),       # Vaccines for initial vaccination
      allowed_vaccines_boost = c("MenACWY")               # Booster vaccine alternative
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
  
  # Read all sheets from the Excel file.
  sheet_names <- excel_sheets(file_path)
  data_list <- lapply(sheet_names, function(sheet) {
    # For specific sheets in Scenario A, skip the first 12 rows (if needed)
    if (scenario_name == "Scenario A" && sheet %in% c("infection", "mortality", "cost_IMD")) {
      header <- suppressMessages(read_excel(file_path, sheet = sheet, n_max = 1, col_names = TRUE))
      df <- suppressMessages(read_excel(file_path, sheet = sheet, skip = 12, col_names = FALSE))
      names(df) <- names(header)
      return(df)
    } else {
      read_excel(file_path, sheet = sheet)
    }
  })
  names(data_list) <- sheet_names
  
  # Read vaccine cost data.
  vaccine_costs <- read_excel(file_path_costs)
  
  # Extract general parameters.
  n_cycles    <- scenario_params$n_cycles
  cycle_length<- with(data_list$model_settings, Value[Name == "cycle_len"])
  d_c         <- with(data_list$model_settings, Value[Name == "d_c"])
  d_e         <- with(data_list$model_settings, Value[Name == "d_e"])
  coverage    <- with(data_list$model_settings, Value[Name == "coverage"])
  
  c_admin     <- with(data_list$costs, Value[Name == "c_admin"])
  c_MenABCWY  <- with(vaccine_costs, Value[Name == "c_MenABCWY"])
  c_MenACWY   <- with(vaccine_costs, Value[Name == "c_MenACWY"])
  c_MenC      <- with(vaccine_costs, Value[Name == "c_MenC"])
  
  p_bg_mort   <- rate_to_prob(with(data_list$mortality, Background_Mortality), t = cycle_length)
  p_B         <- rate_to_prob(with(data_list$infection, SerogroupB_Infection), t = cycle_length)
  p_C         <- rate_to_prob(with(data_list$infection, SerogroupC_Infection), t = cycle_length)
  p_W         <- rate_to_prob(with(data_list$infection, SerogroupW_Infection), t = cycle_length)
  p_Y         <- rate_to_prob(with(data_list$infection, SerogroupY_Infection), t = cycle_length)
  
  p_B_DeadIMD<- rate_to_prob(with(data_list$mortality, SerogroupB_Dead), t = cycle_length)
  p_C_DeadIMD<- rate_to_prob(with(data_list$mortality, SerogroupC_Dead), t = cycle_length)
  p_W_DeadIMD<- rate_to_prob(with(data_list$mortality, SerogroupW_Dead), t = cycle_length)
  p_Y_DeadIMD<- rate_to_prob(with(data_list$mortality, SerogroupY_Dead), t = cycle_length)
  
  p_IMD_Scarring       <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Scarring"])
  p_IMD_Single_Amput   <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Single_Amput"])
  p_IMD_Multiple_Amput <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Multiple_Amput"])
  p_IMD_Neuro_Disability<- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Neuro_Disability"])
  p_IMD_Hearing_Loss   <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Hearing_Loss"])
  p_IMD_Renal_Failure  <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Renal_Failure"])
  p_IMD_Seizure        <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Seizure"])
  p_IMD_Paralysis      <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Paralysis"])
  
  # Read vaccine effectiveness data from the "vac_effect" sheet.
  ve_data <- data.frame(
    vaccine = c("MenABCWY_for_SeroACWY", "MenABCWY_for_SeroB", "MenACWY", "MenC"),
    effectiveness = with(data_list$vac_effect,
                         Value[Name %in% c("MenABCWY_for_SeroACWY", 
                                           "MenABCWY_for_SeroB", 
                                           "MenACWY", 
                                           "MenC")])
  )
  
  # Calculate VE based on the scenario.
  # If the scenario contains both allowed_vaccines_init and allowed_vaccines_boost,
  # we are in a two-dose (initial + booster) setting.
  if (!is.null(scenario_params$allowed_vaccines_init) &&
      !is.null(scenario_params$allowed_vaccines_boost)) {
    
    # --- Initial Vaccination Calculation ---
    vt_init <- scenario_params$vaccination_time[1]  # typically 0
    allowed_init <- scenario_params$allowed_vaccines_init
    ve_data_init <- ve_data[ve_data$vaccine %in% allowed_init, ]
    
    if (vt_init == 0) {
      ve_init_df <- calculate_all_ve(ve_data_init, max_year = 10, n_cycles = n_cycles)
    } else {
      ve_init_df <- calculate_all_ve_shifted(ve_data_init, max_year = 10, n_cycles = n_cycles, shift = vt_init)
    }
    # Here we assume that for the initial vaccination the effectiveness is based on the maximum
    # of the allowed vaccines at each cycle.
    ve_init <- aggregate(Effectiveness ~ Cycle, data = ve_init_df, max)$Effectiveness
    
    # --- Booster Vaccination Calculation ---
    vt_boost <- scenario_params$vaccination_time[2]  # e.g., cycle 12
    allowed_boost <- scenario_params$allowed_vaccines_boost
    ve_data_boost <- ve_data[ve_data$vaccine %in% allowed_boost, ]
    
    # For booster, we use a shifted function (assuming booster is delayed until vt_boost)
    ve_boost_df <- calculate_all_ve_shifted(ve_data_boost, max_year = 10, n_cycles = n_cycles, shift = vt_boost)
    ve_boost <- aggregate(Effectiveness ~ Cycle, data = ve_boost_df, max)$Effectiveness
    
    # --- Combine Initial and Booster VE Curves ---
    # Here we take the maximum VE at each cycle (this can be adjusted as needed).
    ve_combined <- pmax(ve_init, ve_boost)
    
    # Map the combined VE to the appropriate vaccine parameter.
    # For example, if the booster vaccine is MenABCWY_for_SeroACWY (Scenario C), assign:
    if (scenario_name == "Scenario C") {
      ve_MenABCWY_forACWY <- ve_combined * coverage
      # For the other vaccine, assign the initial effectiveness curve for reference.
      ve_MenACWY <- ve_init * coverage
      # (Adjust naming as needed by the simulation model)
    } else if (scenario_name == "Scenario D") {
      # In Scenario D, the booster is MenACWY.
      ve_MenACWY <- ve_combined * coverage
      # Similarly, you may decide how to treat the other allowed vaccine.
    }
    
  } else {
    # Single vaccination event scenario.
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
  }
  
  # --- Extract health state costs ---
  if ("cost_IMD" %in% names(data_list)) {
    c_IMD_infection <- data_list$cost_IMD$cost
  } else {
    c_IMD_infection <- with(data_list$costs, Value[Name == "c_IMD_infection"])
  }
  c_Scarring       <- with(data_list$costs, Value[Name == "c_Scarring"])
  c_Single_Amput   <- with(data_list$costs, Value[Name == "c_Single_Amput"])
  c_Multiple_Amput <- with(data_list$costs, Value[Name == "c_Multiple_Amput"])
  c_Neuro_Disab    <- with(data_list$costs, Value[Name == "c_Neuro_Disab"])
  c_Hearing_Loss   <- with(data_list$costs, Value[Name == "c_Hearing_Loss"])
  c_Renal_Failure  <- with(data_list$costs, Value[Name == "c_Renal_Failure"])
  c_Seizure        <- with(data_list$costs, Value[Name == "c_Seizure"])
  c_Paralysis      <- with(data_list$costs, Value[Name == "c_Paralysis"])
  c_Dead           <- with(data_list$costs, Value[Name == "c_Dead"])
  c_Healthy        <- with(data_list$costs, Value[Name == "c_Healthy"])
  
  # --- Extract health state utilities ---
  u_Healthy        <- with(data_list$utilities, Value[Name == "u_Healthy"])
  u_IMD            <- with(data_list$utilities, Value[Name == "u_IMD"])
  u_Scarring       <- with(data_list$utilities, Value[Name == "u_Scarring"])
  u_Single_Amput   <- with(data_list$utilities, Value[Name == "u_Single_Amput"])
  u_Multiple_Amput <- with(data_list$utilities, Value[Name == "u_Multiple_Amput"])
  u_Neuro_Disability<- with(data_list$utilities, Value[Name == "u_Neuro_Disability"])
  u_Hearing_Loss   <- with(data_list$utilities, Value[Name == "u_Hearing_Loss"])
  u_Renal_Failure  <- with(data_list$utilities, Value[Name == "u_Renal_Failure"])
  u_Seizure        <- with(data_list$utilities, Value[Name == "u_Seizure"])
  u_Paralysis      <- with(data_list$utilities, Value[Name == "u_Paralysis"])
  u_Dead           <- with(data_list$utilities, Value[Name == "u_Dead"])
  
  # --- Build the base parameter list ---
  params_base <- list(
    n_cycles         = n_cycles,
    cycle_length     = cycle_length,
    d_c              = d_c,
    d_e              = d_e,
    coverage         = coverage,
    c_MenABCWY       = c_MenABCWY,
    c_MenACWY        = c_MenACWY,
    c_MenC           = c_MenC,
    c_admin          = c_admin,
    p_bg_mort        = p_bg_mort,
    p_B              = p_B,
    p_C              = p_C,
    p_W              = p_W,
    p_Y              = p_Y,
    p_B_DeadIMD      = p_B_DeadIMD,
    p_C_DeadIMD      = p_C_DeadIMD,
    p_W_DeadIMD      = p_W_DeadIMD,
    p_Y_DeadIMD      = p_Y_DeadIMD,
    p_IMD_Scarring   = p_IMD_Scarring,
    p_IMD_Single_Amput = p_IMD_Single_Amput,
    p_IMD_Multiple_Amput = p_IMD_Multiple_Amput,
    p_IMD_Neuro_Disability = p_IMD_Neuro_Disability,
    p_IMD_Hearing_Loss = p_IMD_Hearing_Loss,
    p_IMD_Renal_Failure = p_IMD_Renal_Failure,
    p_IMD_Seizure    = p_IMD_Seizure,
    p_IMD_Paralysis  = p_IMD_Paralysis,
    c_IMD_infection  = c_IMD_infection,
    c_Scarring       = c_Scarring,
    c_Single_Amput   = c_Single_Amput,
    c_Multiple_Amput = c_Multiple_Amput,
    c_Neuro_Disab    = c_Neuro_Disab,
    c_Hearing_Loss   = c_Hearing_Loss,
    c_Renal_Failure  = c_Renal_Failure,
    c_Seizure        = c_Seizure,
    c_Paralysis      = c_Paralysis,
    c_Dead           = c_Dead,
    c_Healthy        = c_Healthy,
    u_Healthy        = u_Healthy,
    u_IMD            = u_IMD,
    u_Scarring       = u_Scarring,
    u_Single_Amput   = u_Single_Amput,
    u_Multiple_Amput = u_Multiple_Amput,
    u_Neuro_Disability = u_Neuro_Disability,
    u_Hearing_Loss   = u_Hearing_Loss,
    u_Renal_Failure  = u_Renal_Failure,
    u_Seizure        = u_Seizure,
    u_Paralysis      = u_Paralysis,
    u_Dead           = u_Dead
  )
  
  # --- Combine base parameters with scenario-specific settings ---
  # Also add the VE parameters calculated above.
  if (!is.null(scenario_params$allowed_vaccines_init) &&
      !is.null(scenario_params$allowed_vaccines_boost)) {
    if (scenario_name == "Scenario C") {
      params_base$ve_MenABCWY_forACWY <- ve_MenABCWY_forACWY
      params_base$ve_MenACWY <- ve_MenACWY
    } else if (scenario_name == "Scenario D") {
      params_base$ve_MenACWY <- ve_MenACWY
      # If needed, add additional VE parameters for scenario D.
    }
  } else {
    params_base$ve_MenABCWY_forACWY <- ve_MenABCWY_forACWY
    params_base$ve_MenABCWY_forB    <- ve_MenABCWY_forB
    params_base$ve_MenACWY          <- ve_MenACWY
    params_base$ve_MenC             <- ve_MenC
  }
  
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
  for (vac in vaccine_names) {
    if (is.null(sd_values[[vac]])) {
      sd_values[[vac]] <- NULL
    }
  }
  
  return(sd_values)
}

sd_values <- get_basecase_sd_values()

#----------------------------------------------------------------------------
# 4) PSA SAMPLE GENERATION
#----------------------------------------------------------------------------
generate_psa_samples <- function(nsamp, scenario_name = "Scenario A") {
  # 1) Retrieve the base-case parameters for the selected scenario.
  basecase_params <- get_basecase_params(scenario_name)
  
  # 2) Filter parameters that have an associated SD value.
  filtered_params <- basecase_params[names(basecase_params) %in% names(sd_values)]
  
  # 3) Flatten the filtered list into a named vector.
  params_vector <- unlist(filtered_params)
  
  # Separate non-zero and zero parameters.
  non_zero_params <- params_vector[sapply(params_vector, function(x) !all(x == 0))]
  zero_params     <- params_vector[sapply(params_vector, function(x) all(x == 0))]
  
  # 4) Define the distribution type for each parameter based on its name.
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
  
  # 5) Define parameterization types ("mean, sd") for all valid parameters.
  param_types <- sapply(dist_types, function(dist) {
    if (is.na(dist)) NA else "mean, sd"
  })
  
  # 6) Create a list pairing each parameter's mean and SD.
  sd_vector <- unlist(sd_values)
  mean_sd_list <- list()
  for (param in parameter_names) {
    base_name <- sub("\\d+$", "", param)
    if (base_name %in% names(sd_vector)) {
      mean_sd_list[[param]] <- c(mean = non_zero_params[[param]], sd = sd_vector[[base_name]])
    }
  }
  
  # 7) Filter valid parameters (those with both mean and SD).
  valid_params <- names(mean_sd_list)
  valid_dist_types <- dist_types[valid_params]
  valid_param_types <- param_types[valid_params]
  
  # 8) Generate distribution parameters for each valid parameter.
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
  
  # 9) Generate the PSA samples using a sampling function (e.g., gen_psa_samp from dampack).
  df_psa <- gen_psa_samp(
    params = valid_params,
    dists = valid_dist_types,
    parameterization_types = valid_param_types,
    dists_params = dist_params,
    nsamp = nsamp
  )
  
  # 10) For parameters that had a base-case value of zero, add them back.
  if (length(zero_params) > 0) {
    for (param in names(zero_params)) {
      df_psa[[param]] <- rep(0, nrow(df_psa))
    }
  }
  
  # 11) Separate scalar and vector parameters.
  scalar_params <- names(filtered_params)[sapply(filtered_params, length) == 1]
  vector_params <- names(filtered_params)[sapply(filtered_params, length) > 1]
  
  # 12) Identify the columns corresponding to each vector parameter.
  vector_param_cols <- lapply(vector_params, function(param) {
    pattern <- paste0("^", param, "\\d+$")
    cols <- grep(pattern, names(df_psa), value = TRUE)
    cols[order(as.numeric(gsub(".*?(\\d+)$", "\\1", cols)))]
  })
  names(vector_param_cols) <- vector_params
  
  # 13) Convert the PSA data frame into a list of parameter sets.
  samples_list <- split(df_psa, seq(nrow(df_psa)))
  samples_list <- lapply(samples_list, function(row_df) {
    sample_list <- list()
    for (param in scalar_params) {
      sample_list[[param]] <- if (!is.null(row_df[[param]])) row_df[[param]] else filtered_params[[param]]
    }
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
  
  # 14) Add any missing global parameters from the base-case.
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
# 1) Generate 10 PSA samples for a specified scenario (e.g., "Scenario A", "Scenario B", "Scenario C", or "Scenario D"):
#      samples <- generate_psa_samples(nsamp = 10, scenario_name = "Scenario C")
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
