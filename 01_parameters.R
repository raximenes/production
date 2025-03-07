# 01_parameters.R
###############################################################################
# PSA Analysis: Reading Parameters from Excel, Constructing Base-Case Params,
# Generating PSA Samples, and Combining Them into a List of Simulations
#
# This script includes:
#   1) Helper functions (e.g., for calculating Gamma parameters)
#   2) A function to read Excel sheets and build a base‐case parameter list
#      (get_basecase_params)
#   3) A function to calculate vaccine effectiveness decay over time
#      (calculate_all_ve)
#   4) A function to read Excel sheets and build a list of standard deviation 
#      (SD) values (get_basecase_sd_values) using the same structure as 
#      get_basecase_params
#   5) A function to generate PSA samples (generate_psa_samples)
#
# Files used:
#   - data/IMD Data.xls (main parameters)
#   - data/confidential/vaccine_costs.xlsx (vaccine costs)
###############################################################################

#----------------------------------------------------------------------------
# 1) HELPER FUNCTIONS
#----------------------------------------------------------------------------

# Function to calculate the linear decay of vaccine effectiveness over time.
# For each vaccine, effectiveness starts at a base value (at cycle 1) and decays
# linearly to 0 by the specified max_year.
calculate_all_ve <- function(data, max_year = 4, n_cycles = 89) {
  years <- 1:n_cycles
  out_list <- lapply(data$vaccine, function(vac) {
    # Extract the base effectiveness for this vaccine.
    base_ve <- data[data$vaccine == vac, "effectiveness"]
    ve_vec <- sapply(years, function(yr) {
      if (yr >= 1 && yr <= max_year) {
        # Linear interpolation between base effectiveness (at year 1)
        # and 0 (at max_year).
        x0 <- 1
        x1 <- max_year
        y0 <- base_ve
        y1 <- 0
        y0 + (y1 - y0) / (x1 - x0) * (yr - x0)
      } else {
        0
      }
    })
    data.frame(Cycle = years, Vaccine = vac, Effectiveness = ve_vec)
  })
  do.call(rbind, out_list)
}

#----------------------------------------------------------------------------
# 2) READING DATA AND BUILDING BASE-CASE PARAMETERS
#----------------------------------------------------------------------------

# Define the file paths for the Excel data files.
file_path       <- "data/IMD Data.xls"
file_path_costs <- "data/confidential/vaccine_costs.xlsx"

# Function to read the Excel sheets and construct the base-case parameter list.
get_basecase_params <- function() {
  # 1) Read all sheets from the main Excel file.
  sheet_names <- excel_sheets(file_path)
  data_list   <- lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet))
  names(data_list) <- sheet_names
  
  # 2) Read the separate vaccine costs file.
  vaccine_costs <- read_excel(file_path_costs)
  
  # 3) Extract model settings from the 'model_settings' sheet.
  n_cycles     <- with(data_list$model_settings, Value[Name == "n_cycles"])
  cycle_length <- with(data_list$model_settings, Value[Name == "cycle_len"])
  d_c          <- with(data_list$model_settings, Value[Name == "d_c"])  # Discount rate for costs
  d_e          <- with(data_list$model_settings, Value[Name == "d_e"])  # Discount rate for effectiveness
  coverage     <- with(data_list$model_settings, Value[Name == "coverage"])
  
  # 4) Extract the administrative cost from the 'costs' sheet.
  c_admin <- with(data_list$costs, Value[Name == "c_admin"])
  
  # 5) Extract vaccine costs from the vaccine_costs file.
  c_MenABCWY <- with(vaccine_costs, Value[Name == "c_MenABCWY"])
  c_MenACWY  <- with(vaccine_costs, Value[Name == "c_MenACWY"])
  c_MenC     <- with(vaccine_costs, Value[Name == "c_MenC"])
  
  # 6) Convert rates to probabilities for mortality and infection.
  p_bg_mort <- rate_to_prob(with(data_list$mortality, Background_Mortality), t = cycle_length)
  p_B       <- rate_to_prob(with(data_list$infection, SerogroupB_Infection), t = cycle_length)
  p_C       <- rate_to_prob(with(data_list$infection, SerogroupC_Infection), t = cycle_length)
  p_W       <- rate_to_prob(with(data_list$infection, SerogroupW_Infection), t = cycle_length)
  p_Y       <- rate_to_prob(with(data_list$infection, SerogroupY_Infection), t = cycle_length)
  
  # 7) Convert rates to probabilities for death from infection.
  p_B_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupB_Dead), t = cycle_length)
  p_C_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupC_Dead), t = cycle_length)
  p_W_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupW_Dead), t = cycle_length)
  p_Y_DeadIMD <- rate_to_prob(with(data_list$mortality, SerogroupY_Dead), t = cycle_length)
  
  # 8) Extract sequela probabilities from the 'sequelae_probp_IMD' sheet.
  p_IMD_Scarring         <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Scarring"])
  p_IMD_Single_Amput     <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Single_Amput"])
  p_IMD_Multiple_Amput   <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Multiple_Amput"])
  p_IMD_Neuro_Disability <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Neuro_Disability"])
  p_IMD_Hearing_Loss     <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Hearing_Loss"])
  p_IMD_Renal_Failure    <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Renal_Failure"])
  p_IMD_Seizure          <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Seizure"])
  p_IMD_Paralysis        <- with(data_list$sequelae_probp_IMD, Value[Name == "p_IMD_Paralysis"])
  
  # 9) Extract vaccine effectiveness from the 'vac_effect' sheet and compute its linear decay.
  ve_data <- data.frame(
    vaccine       = c("MenABCWY_for_SeroACWY", "MenABCWY_for_SeroB", "MenACWY", "MenC"),
    effectiveness = with(data_list$vac_effect,
                         Value[Name %in% c("MenABCWY_for_SeroACWY", 
                                           "MenABCWY_for_SeroB", 
                                           "MenACWY", 
                                           "MenC")])
  )
  # Compute the time-varying vaccine effectiveness using linear decay.
  ve_all <- calculate_all_ve(ve_data, max_year = 10, n_cycles = n_cycles)
  
  # Multiply each vaccine’s effectiveness by the coverage rate.
  ve_MenABCWY_forACWY <- ve_all[ve_all$Vaccine == "MenABCWY_for_SeroACWY", "Effectiveness"] 
  ve_MenABCWY_forB    <- ve_all[ve_all$Vaccine == "MenABCWY_for_SeroB",    "Effectiveness"]
  ve_MenACWY          <- ve_all[ve_all$Vaccine == "MenACWY",               "Effectiveness"] 
  ve_MenC             <- ve_all[ve_all$Vaccine == "MenC",                  "Effectiveness"] 
  
  # 10) Extract health state costs.
  # For c_IMD_infection, use the cost_IMD sheet if available (age-specific cost),
  # otherwise use the 'costs' sheet.
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
  
  # 11) Extract health state utilities from the 'utilities' sheet.
  u_Healthy          <- with(data_list$utilities, Value[Name == "u_Healthy"])
  u_IMD              <- with(data_list$utilities, Value[Name == "u_IMD"])
  u_Scarring         <- with(data_list$utilities, Value[Name == "u_Scarring"])
  u_Single_Amput     <- with(data_list$utilities, Value[Name == "u_Single_Amput"])
  u_Multiple_Amput   <- with(data_list$utilities, Value[Name == "u_Multiple_Amput"])
  u_Neuro_Disability <- with(data_list$utilities, Value[Name == "u_Neuro_Disability"])
  u_Hearing_Loss     <- with(data_list$utilities, Value[Name == "u_Hearing_Loss"])
  u_Renal_Failure    <- with(data_list$utilities, Value[Name == "u_Renal_Failure"])
  u_Seizure          <- with(data_list$utilities, Value[Name == "u_Seizure"])
  u_Paralysis        <- with(data_list$utilities, Value[Name == "u_Paralysis"])
  u_Dead             <- with(data_list$utilities, Value[Name == "u_Dead"])
  
  # 12) Combine all the extracted parameters into a single list and return it.
  params <- list(
    n_cycles     = n_cycles,
    cycle_length = cycle_length,
    d_c = d_c,
    d_e = d_e,
    coverage = coverage,
    
    # Vaccine costs
    c_MenABCWY = c_MenABCWY,
    c_MenACWY  = c_MenACWY,
    c_MenC     = c_MenC,
    c_admin    = c_admin,
    
    # Infection and mortality probabilities
    p_bg_mort = p_bg_mort,
    p_B       = p_B,
    p_C       = p_C,
    p_W       = p_W,
    p_Y       = p_Y,
    
    # Mortality from infection probabilities
    p_B_DeadIMD = p_B_DeadIMD,
    p_C_DeadIMD = p_C_DeadIMD,
    p_W_DeadIMD = p_W_DeadIMD,
    p_Y_DeadIMD = p_Y_DeadIMD,
    
    # Sequela probabilities
    p_IMD_Scarring         = p_IMD_Scarring,
    p_IMD_Single_Amput     = p_IMD_Single_Amput,
    p_IMD_Multiple_Amput   = p_IMD_Multiple_Amput,
    p_IMD_Neuro_Disability = p_IMD_Neuro_Disability,
    p_IMD_Hearing_Loss     = p_IMD_Hearing_Loss,
    p_IMD_Renal_Failure    = p_IMD_Renal_Failure,
    p_IMD_Seizure          = p_IMD_Seizure,
    p_IMD_Paralysis        = p_IMD_Paralysis,
    
    # Vaccine effectiveness (time-varying)
    ve_MenABCWY_forACWY = ve_MenABCWY_forACWY,
    ve_MenABCWY_forB    = ve_MenABCWY_forB,
    ve_MenACWY          = ve_MenACWY,
    ve_MenC             = ve_MenC,
    
    # Health state costs
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
    
    # Health state utilities
    u_Healthy          = u_Healthy,
    u_IMD              = u_IMD,
    u_Scarring         = u_Scarring,
    u_Single_Amput     = u_Single_Amput,
    u_Multiple_Amput   = u_Multiple_Amput,
    u_Neuro_Disability = u_Neuro_Disability,
    u_Hearing_Loss     = u_Hearing_Loss,
    u_Renal_Failure    = u_Renal_Failure,
    u_Seizure          = u_Seizure,
    u_Paralysis        = u_Paralysis,
    u_Dead             = u_Dead
  )
  
  return(params)
}

#----------------------------------------------------------------------------
# 3) BUILDING THE STANDARD DEVIATION (SD) VALUES LIST USING THE BASE-CASE STRUCTURE
#----------------------------------------------------------------------------

# This function reads the Excel sheets in the same way as get_basecase_params and
# extracts the SD values for the parameters. The SD values are stored in a list so
# that vector-valued parameters (e.g., age-specific costs) are preserved.
get_basecase_sd_values <- function() {
  # 1) Read all sheets from the main Excel file.
  sheet_names <- excel_sheets(file_path)
  data_list <- lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet))
  names(data_list) <- sheet_names
  
  # 2) Read the separate vaccine costs file.
  vaccine_costs <- read_excel(file_path_costs)
  
  # 3) Define a helper function to extract SD values from a sheet.
  #    Returns a list mapping parameter names (from the 'Name' column) to their SD.
  extract_sd <- function(df) {
    if ("sd" %in% names(df)) {
      return(as.list(setNames(as.numeric(df$sd), df$Name)))
    } else {
      return(list())
    }
  }
  
  # 4) Extract SD values from sheets that have one row per parameter.
  sd_costs         <- extract_sd(data_list$costs)
  sd_vaccine_costs <- extract_sd(vaccine_costs)
  sd_utilities     <- extract_sd(data_list$utilities)
  sd_sequelae      <- extract_sd(data_list$sequelae_probp_IMD)
  sd_vac_effect    <- extract_sd(data_list$vac_effect)
  
  # 5) Combine the above SD values into one list.
  sd_values <- c(sd_costs, sd_vaccine_costs, sd_utilities, sd_sequelae, sd_vac_effect)
  
  # 5.5) Override the SD for c_IMD_infection using the cost_IMD sheet if available.
  if ("cost_IMD" %in% names(data_list)) {
    # Extract the entire 'sd' column (a vector) and assign it to parameter "c_IMD_infection"
    sd_values[["c_IMD_infection"]] <- as.numeric(data_list$cost_IMD$sd)
  }
  
  # 6) From the 'mortality' sheet, extract SD values from columns starting with "sd_".
  if (!is.null(data_list$mortality)) {
    mortality_sd_cols <- grep("^sd_", names(data_list$mortality), value = TRUE)
    for (col in mortality_sd_cols) {
      # Rename the parameter by replacing "sd_" with "p_"
      param_name <- sub("^sd_", "p_", col)
      # Store the entire column (which may be a vector) in the list.
      sd_values[[param_name]] <- as.numeric(data_list$mortality[[col]])
    }
  }
  
  # 7) From the 'infection' sheet, extract SD values from columns starting with "sd_".
  if (!is.null(data_list$infection)) {
    infection_sd_cols <- grep("^sd_", names(data_list$infection), value = TRUE)
    for (col in infection_sd_cols) {
      param_name <- sub("^sd_", "p_", col)
      sd_values[[param_name]] <- as.numeric(data_list$infection[[col]])
    }
  }
  
  # 8) For background mortality, convert its SD from a rate to a probability.
  cycle_length <- with(data_list$model_settings, Value[Name == "cycle_len"])
  if ("sd_Background_Mortality" %in% names(data_list$mortality)) {
    sd_values$p_bg_mort <- rate_to_prob(as.numeric(data_list$mortality$sd_Background_Mortality), t = cycle_length)
  }
  
  # 9) Remove any entries that are NA.
  sd_values <- sd_values[!sapply(sd_values, function(x) any(is.na(x)))]
  
  # 10) Rename parameters to match those in the base-case list using a mapping.
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
    if (name %in% names(rename_map)) {
      rename_map[[name]]
    } else {
      name
    }
  })
  
  # 11) For vaccine effectiveness parameters, calculate their time-varying SD values.
  #      Retrieve the coverage rate from model settings.
  coverage <- with(data_list$model_settings, Value[Name == "coverage"])
  vaccine_names <- c("ve_MenABCWY_forACWY", "ve_MenABCWY_forB", "ve_MenACWY", "ve_MenC")
  
  # Create a data frame with vaccine names and their base SD values (if available).
  ve_data <- data.frame(
    vaccine = vaccine_names,
    effectiveness = sapply(vaccine_names, function(vac) {
      if (!is.null(sd_values[[vac]])) { 
        sd_values[[vac]] 
      } else { 
        NA 
      }
    })
  )
  
  # Compute the time-varying vaccine effectiveness SD values using the same linear decay.
  ve_all <- calculate_all_ve(ve_data, max_year = 10, n_cycles = 89)
  for (vac in vaccine_names) {
    sd_values[[vac]] <- ve_all[ve_all$Vaccine == vac, "Effectiveness"]
  }
  
  return(sd_values)
}

# Precompute the global SD values using the new function.
sd_values <- get_basecase_sd_values()

#----------------------------------------------------------------------------
# 4) PSA SAMPLE GENERATION
#----------------------------------------------------------------------------

# This function generates PSA samples. It creates a list of parameter sets
# (one per PSA draw) using the base-case parameters and the SD values for sampling.
generate_psa_samples <- function(nsamp) {
  
  # 1) Retrieve the base-case parameters.
  basecase_params <- get_basecase_params()
  
  # 2) Filter only those parameters for which SD values are available.
  filtered_list_for_psa <- basecase_params[names(basecase_params) %in% names(sd_values)]
  
  params_for_psa <- unlist(filtered_list_for_psa)
  
  non_zero_params <- params_for_psa[sapply(params_for_psa, function(x) !all(x == 0))]
  zero_params     <- params_for_psa[sapply(params_for_psa, function(x) all(x == 0))]
  
  
  # 3) Flatten the filtered list into a named vector.
  base_params <- non_zero_params
  my_params   <- names(base_params)
  
  
  # 4) Define distribution types based on parameter names:
  #    - Costs (starting with "c_") use the Gamma distribution.
  #    - Probabilities (starting with "p_"), utilities (starting with "u_"),
  #      and vaccine effectiveness (starting with "ve_") use the Beta distribution.
  my_dists <- sapply(my_params, function(param) {
    if (grepl("^c_", param)) {
      "gamma"
    } else if (grepl("^(u_|p_|ve_)", param)) {
      "beta"
    } else {
      NA
    }
  })
  
  # 5) Define parameterization types for each distribution:
  #    For Beta: "mean, sd" and for Gamma: "shape, scale"
  # Using mean, sd for all parameters since gen_pas_sample supports it.
  my_parameterization_types <- sapply(my_dists, function(dist) {
    if (is.na(dist)) {
      NA
    } else if (dist == "beta") {
      "mean, sd"
    } else if (dist == "gamma") {
      "mean, sd"
    } else {
      NA
    }
  })
  
  # # 6)
  
  # 7) Create a list pairing each parameter's mean and SD.
  
  
  unl_sd_val <- unlist(sd_values)
  mean_sd_list <- list()
  for (param in names(base_params)) {
    if (param %in% names(unl_sd_val)) {
      mean_sd_list[[param]] <- c(mean = base_params[[param]], sd = unl_sd_val[[param]])
    }
  }
  
  
  # 8) Generate distribution parameters for use with gen_psa_samp.
  #    (For Beta: use (mean, sd) directly; for Gamma: convert (mean, sd) to (shape, scale)).
  generate_my_dists_params <- function(mean_sd_list, my_parameterization_types) {
    dists_params <- list()
    for (param in names(mean_sd_list)) {
      
      if (param %in% names(my_parameterization_types)) {
        if (my_parameterization_types[[param]] == "mean, sd") {
          # For Beta distribution, use the mean and SD directly.
          dists_params[[param]] <- mean_sd_list[[param]]
        } else {
          NA
        }
        
      }
    }
    return(dists_params)
  }
  
  my_dists_params <- generate_my_dists_params(mean_sd_list, my_parameterization_types)
  
  # Remove parameters that have NA for distribution types. # Make sure this is only a test and exclude if not necessary
  valid_params_idx <- !is.na(my_parameterization_types)
  valid_params     <- my_params[valid_params_idx]
  valid_dists      <- my_dists[valid_params_idx]
  valid_param_types<- my_parameterization_types[valid_params_idx]
  
  # 9) Generate PSA samples using dampack::gen_psa_samp.
  # Note: Ensure that the function gen_psa_samp is available (e.g., from the dampack package).
  df_psa <- gen_psa_samp(
    params = valid_params,
    dists = valid_dists,
    parameterization_types = valid_param_types,
    dists_params = my_dists_params, # need to fix the values!
    nsamp = nsamp
  )
  # including zero parameters back to the data frame
  if (length(zero_params) > 0) {
    for (param in names(zero_params)) {
      df_psa[[param]] <- rep(0, nrow(df_psa))
    }
  }
  
  # 10) Convert the PSA data frame into a list of parameter sets.
  #     Distinguish between scalar parameters (length = 1) and vector parameters (length > 1).
  scalar_params <- names(filtered_list_for_psa)[sapply(filtered_list_for_psa, length) == 1]
  vector_params <- names(filtered_list_for_psa)[sapply(filtered_list_for_psa, length) > 1]
  
  # 11) Identify columns corresponding to each vector parameter.
  #     For example, p_bg_mort1, p_bg_mort2, etc.
  vector_param_cols <- map(vector_params, function(param) {
    pattern <- paste0("^", param, "\\d+$")  # e.g., "p_bg_mort" followed by digits.
    cols <- grep(pattern, names(df_psa), value = TRUE)
    cols_sorted <- cols[order(as.numeric(str_extract(cols, "\\d+$")))]
    return(cols_sorted)
  })
  names(vector_param_cols) <- vector_params
  
  
  # 12) Build a list of samples (one element per PSA draw).
  samples_list <- df_psa %>%
    select(-nsamp) %>%
    split(1:nrow(.)) %>%
    map(function(row) {
      sample_list <- list()
      
      # Scalar parameters: use the drawn value if available; otherwise, revert to base-case.
      for (param in scalar_params) {
        sample_list[[param]] <- if (!is.null(row[[param]])) row[[param]] else filtered_list_for_psa[[param]]
      }
      
      # Vector parameters: combine the corresponding columns into a numeric vector.
      for (param in vector_params) {
        cols <- vector_param_cols[[param]]
        if (length(cols) > 0) {
          sample_list[[param]] <- as.numeric(row[cols])
        } else {
          sample_list[[param]] <- filtered_list_for_psa[[param]]
        }
      }
      return(sample_list)
    })
  
  
  
  
  # 13) Identify any global parameters present in the base-case that are missing from the PSA samples.
  existing_params_psa <- unique(unlist(lapply(samples_list, names)))
  all_params <- names(basecase_params)
  missing_params <- setdiff(all_params, existing_params_psa)
  
  # 14) Add any missing parameters (with their base-case values) to each sample.
  samples_list <- lapply(samples_list, function(sample) {
    for (param in missing_params) {
      sample[[param]] <- basecase_params[[param]]
    }
    return(sample)
  })
  
  # 15) Return the final list of PSA samples.
  return(samples_list)
}

###############################################################################
# EXAMPLE USAGE:
#
# 1) Generate 3 PSA samples:
#      samples <- generate_psa_samples(nsamp = 3)
#
# 2) Inspect one of the sample sets:
#      str(samples[[1]])
#
# 3) Check if any parameters ended up as NA:
#      na_parameters_per_sim <- lapply(samples, function(sim) {
#        na_flags <- sapply(sim, function(param) any(is.na(param)))
#        names(na_flags[na_flags])
#      })
#      na_parameters_per_sim
###############################################################################
