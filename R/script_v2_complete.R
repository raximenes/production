# --------------------------------------------------------
# Deterministic Analysis for All Scenarios of the IMD Economic Model
# --------------------------------------------------------

# 1) Load necessary scripts (ensure that file paths are correct)
source("00_global_setup.R")
source("R/01_parameters_v2.R")
source("R/02_transition_array_v2.R")
source("R/03_markov_model_v2.R")
source("R/04_analysis_v2.R")  # Make sure the updated run_scenario() is used

# 2) Define the scenarios to run
scenarios_to_run <- c("Scenario A", "Scenario B", "Scenario C")

# 3) Loop through each scenario and run the analysis
results <- lapply(scenarios_to_run, function(scenario_name) {
  cat("Running scenario:", scenario_name, "\n")
  
  # Load base-case parameters for the current scenario
  params_det <- get_basecase_params(scenario_name)
  params_det$scenario_name <- scenario_name
  # Run the deterministic simulation to evaluate cost and QALYs for each strategy
  df_det <- simulate_strategies(params_det, wtp = 50000)
  
  # Run the complete scenario analysis (including ICER calculation)
  res_det <- run_scenario(params_det, strategies = c("MenABCWY", "MenACWY", "MenC"), wtp = 50000)
  
  # Format the cost-effectiveness analysis (CEA) table
  table_cea <- format_table_cea(res_det$cea)
  
  # Generate the Markov trace for the "MenABCWY" strategy
  transition_matrix <- build_transition_array(params_det, strategy = "MenABCWY")
  markov_results <- run_markov_model(transition_matrix, params_det)
  
  # Calculate overall survival for "MenABCWY"
  # (Overall survival is computed as 1 minus the sum of the death state probabilities)
  survival_prob <- 1 - rowSums(markov_results[, c("Dead_IMD", "Background_mortality")])
  
  # Return all results for the current scenario
  list(
    scenario = scenario_name,
    deterministic_results = df_det,
    analysis = res_det,
    cea_table = table_cea,
    markov_trace = markov_results,
    survival_probability = survival_prob
    
  )
})

# 4) Print results for all scenarios
print(results)
# Loop through the results list and print each scenario's results with its name
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Scenario:", res$scenario, "\n\n")
  
  cat("Analysis - Raw Results:\n")
  print(res$analysis$raw)
  cat("\n")
  
  cat("Analysis - CEA Results:\n")
  print(res$analysis$cea)
  cat("\n")
  
  cat("CEA Table:\n")
  print(res$cea_table)
  cat("\n--------------------------------------------------------\n")
}))

#results[[1]]$analysis$cea
#results[[2]]$analysis$cea
#results[[3]]$analysis$cea

# 5) Plot the cost-effectiveness frontier for each scenario with a title including the scenario name
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting CE frontier for scenario:", res$scenario, "\n")
  frontier_plot <- plot(res$analysis$cea, label = "all", txtsize = 12) +
    theme(legend.position = c(0.8, 0.3)) +
    ggtitle(paste("Cost-effectiveness Frontier -", res$scenario))
  print(frontier_plot)
}))

# Loop through the results list and plot the Markov trace for each scenario
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting Markov trace for scenario:", res$scenario, "\n")
  
  # Extract the stored markov trace from the results
  markov_results <- res$markov_trace
  
  # Define n_cycles based on the markov trace dimensions and assign it globally
  n_cycles <- nrow(markov_results) - 1
  assign("n_cycles", n_cycles, envir = .GlobalEnv)
  
  # Define v_names_states using the desired column ordering and assign globally
  v_names_states <- colnames(markov_results)[c(1,14,15,2,3,4,5,6,7,8,9,10,11,12,13)]
  assign("v_names_states", v_names_states, envir = .GlobalEnv)
  
  # Generate the Markov trace plot using plot_trace() and add a title
  trace_plot <- plot_trace(markov_results) +
    ggtitle(paste("Markov Trace for MenABCWY -", res$scenario))
  
  print(trace_plot)
}))




# Loop through the results list and plot overall survival for each scenario
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting overall survival for scenario:", res$scenario, "\n")
  
  # Plot overall survival using the stored survival_probability from each scenario
  plot(res$survival_probability, type = 'l', ylim = c(0, 1),
       xlab = "Cycle (Year)", ylab = "Survival Probability",
       main = paste("Overall Survival - MenABCWY -", res$scenario))
  grid()
  
  # Print the approximate life expectancy (sum of survival probabilities)
  cat("Approximate life expectancy for", res$scenario, ":", sum(res$survival_probability), "\n")
  cat("--------------------------------------------------------\n")
}))


#############################################
# Probabilistic Sensitivity Analysis (PSA) - Extended Analyses for Each Scenario
#############################################
n_samp <- 100
scenarios_to_run <- c("Scenario A", "Scenario B", "Scenario C")
set.seed(123)

# ----- Step A: Loop over each scenario to run the PSA -------------------

psa_results <- lapply(scenarios_to_run, function(scenario_name) {
  cat("Running PSA for scenario:", scenario_name, "\n")
  
  # Load base-case parameters for this scenario and tag with scenario name
  params_det <- get_basecase_params(scenario_name)
  params_det$scenario_name <- scenario_name
  # Make sure scenario_name is available globally for any internal functions that require it
  assign("scenario_name", scenario_name, envir = .GlobalEnv)
  
  # Generate PSA samples (ensure your generate_psa_samples() uses the scenario as needed)
  samples_list <- generate_psa_samples(n_samp)
  
  # Run the PSA using the scenario-specific parameters
  df_psa_results <- run_psa_custom(
    psa_samp = samples_list,
    params_basecase = params_det,
    strategies = c("MenABCWY", "MenACWY", "MenC"),
    wtp = 50000,
    parallel = TRUE,
    progress = FALSE
  )
  
  # Summarize PSA outcomes by strategy (mean cost and QALYs)
  df_psa_summary <- df_psa_results %>%
    group_by(Strategy) %>%
    summarise(
      Cost = mean(Cost),
      QALYs = mean(QALYs)
    ) %>%
    arrange(Cost)
  
  # Compute a cost-effectiveness analysis (CEA) based on the summary (if desired)
  df_psa_cea <- calculate_icers(
    cost = df_psa_summary$Cost,
    effect = df_psa_summary$QALYs,
    strategies = df_psa_summary$Strategy
  )
  
  # Tag the summaries with the scenario name
  df_psa_summary$Scenario <- scenario_name
  df_psa_cea$Scenario <- scenario_name
  
  list(
    scenario = scenario_name,
    psa_results = df_psa_results,   # Raw PSA simulation results (with a "sim" column)
    psa_samples = samples_list,       # The list of parameter samples used
    psa_summary = df_psa_summary,
    psa_cea = df_psa_cea
  )
})

# Print a summary of PSA CEA results for each scenario
invisible(lapply(psa_results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("PSA CEA Results for scenario:", res$scenario, "\n")
  print(res$psa_cea)
  cat("--------------------------------------------------------\n")
}))


# 10) Print PSA CEA results for each scenario
invisible(lapply(psa_results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("PSA Results for scenario:", res$scenario, "\n")
  cat("PSA CEA Results:\n")
  print(res$psa_cea)
  cat("--------------------------------------------------------\n")
}))

# 11) Plot the cost-effectiveness frontier for each scenario (PSA) with a title including the scenario name
invisible(lapply(psa_results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting PSA CE frontier for scenario:", res$scenario, "\n")
  frontier_plot <- plot(res$psa_cea, label = "all", txtsize = 12) +
    theme(legend.position = c(0.8, 0.3)) +
    ggtitle(paste("PSA Cost-effectiveness Frontier -", res$scenario))
  print(frontier_plot)
}))


# Loop over each scenario's PSA results and generate the ICEP plot
invisible(lapply(psa_results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting ICEP for scenario:", res$scenario, "\n")
  
  # Extract the raw PSA results for the current scenario
  df_psa_results <- res$psa_results
  
  # Define the willingness-to-pay (WTP) threshold
  ce_threshold <- 50000  # $50,000 per QALY
  
  # Choose the reference strategy (the least expensive) based on average cost
  reference_strategy <- df_psa_results %>%
    group_by(Strategy) %>%
    summarise(Cost = mean(Cost)) %>%
    arrange(Cost) %>%
    slice(1)
  
  # For each simulation, calculate incremental costs and QALYs relative to the reference
  df_psa_incremental <- df_psa_results %>%
    group_by(sim) %>%
    mutate(
      Incremental_Cost = Cost - Cost[Strategy == reference_strategy$Strategy],
      Incremental_QALYs = QALYs - QALYs[Strategy == reference_strategy$Strategy]
    ) %>%
    ungroup()
  
  # Create the ICEP plot with all PSA points
  icep_plot <- ggplot(df_psa_incremental, aes(x = Incremental_QALYs, y = Incremental_Cost, color = Strategy)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = ce_threshold, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste("Incremental Cost-Effectiveness Plane (ICEP) - PSA -", res$scenario),
      x = "Incremental QALYs",
      y = "Incremental Cost",
      color = "Strategy"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Calculate strategy means to highlight on the plot
  df_psa_means <- df_psa_incremental %>%
    group_by(Strategy) %>%
    summarise(
      Mean_Incremental_Cost = mean(Incremental_Cost),
      Mean_Incremental_QALYs = mean(Incremental_QALYs)
    )
  
  # Add highlighted mean points and text labels to the plot
  icep_plot <- icep_plot +
    geom_point(data = df_psa_means, aes(x = Mean_Incremental_QALYs, y = Mean_Incremental_Cost, color = Strategy),
               size = 4, shape = 18) +
    geom_text(data = df_psa_means, aes(x = Mean_Incremental_QALYs, y = Mean_Incremental_Cost, label = Strategy),
              vjust = -1, hjust = 1, size = 4, color = "black") +
    annotate(
      "text",
      label = paste("Cost-Effectiveness Threshold = $", format(ce_threshold, big.mark = ","), "/QALY"),
      color = "red",
      size = 4,
      hjust = 0.5,
      vjust = -0.5,
      angle = 9
    )
  
  print(icep_plot)
}))

# -------------------------------
# PSA Analysis Plots for Each Scenario
# -------------------------------

invisible(lapply(psa_results, function(res) {
  # Extract the raw PSA simulation results for the current scenario
  df_psa_results <- res$psa_results
  
  # ====== Plot 2: Cost Distribution by Strategy ======
  cost_density_plot <- ggplot(df_psa_results, aes(x = Cost, fill = Strategy)) +
    geom_density(alpha = 0.5) +  # Density plot to show cost variability
    labs(
      title = paste("Cost Distribution Across Strategies -", res$scenario),
      x = "Cost ($)",
      y = "Density",
      fill = "Strategy"
    ) +
    theme_minimal()
  print(cost_density_plot)
  
  # ====== Plot 3: QALY Distribution by Strategy ======
  qaly_density_plot <- ggplot(df_psa_results, aes(x = QALYs, fill = Strategy)) +
    geom_density(alpha = 0.5) +  # Density plot to show QALY variability
    labs(
      title = paste("QALY Distribution Across Strategies -", res$scenario),
      x = "QALYs",
      y = "Density",
      fill = "Strategy"
    ) +
    theme_minimal()
  print(qaly_density_plot)
  
  # ====== Cost-Effectiveness Acceptability Curve (CEAC) ======
  # Define the WTP threshold
  ce_threshold <- 50000  # Example: $50,000 per QALY
  
  # Define a range of cost-effectiveness thresholds to evaluate
  threshold_values <- seq(0, 200000, by = 1000)
  
  # Initialize an empty data frame to store the CEAC results for this scenario
  ceac_data <- data.frame()
  
  # Loop over each threshold value
  for (threshold in threshold_values) {
    # Calculate the Net Monetary Benefit (NMB) for each simulation
    df_nmb <- df_psa_results %>%
      mutate(NMB = QALYs * threshold - Cost)
    
    # Identify the best (cost-effective) strategy for each simulation
    df_cost_effective <- df_nmb %>%
      group_by(sim) %>%
      summarise(BestStrategy = Strategy[which.max(NMB)]) %>%
      ungroup()
    
    # Calculate the probability of each strategy being cost-effective
    df_prob <- df_cost_effective %>%
      group_by(BestStrategy) %>%
      summarise(Probability = n() / nrow(df_cost_effective)) %>%
      ungroup() %>%
      # Ensure all strategies are included even if they have zero probability
      complete(BestStrategy = unique(df_psa_results$Strategy), fill = list(Probability = 0)) %>%
      mutate(CostEffectivenessThreshold = threshold)
    
    # Append to the CEAC data frame
    ceac_data <- bind_rows(ceac_data, df_prob)
  }
  
  # Rename column for clarity
  ceac_data <- ceac_data %>% rename(Strategy = BestStrategy)
  
  # Generate the CEAC plot
  ceac_plot <- ggplot(ceac_data, 
                      aes(x = CostEffectivenessThreshold, y = Probability, color = Strategy)) +
    geom_line(size = 1.2) +
    labs(
      title = paste("Cost-Effectiveness Acceptability Curve (CEAC) -", res$scenario),
      x = "Cost-Effectiveness Threshold per QALY",
      y = "Probability of Being Cost-Effective"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  print(ceac_plot)
}))

################################################################################
# ----- Step B: Extended PSA Analyses for Each Scenario -------------------
################################################################################
invisible(lapply(psa_results, function(res) {
  cat("========================================================\n")
  cat("Performing extended PSA analyses for scenario:", res$scenario, "\n")
  
  # (1) Convert raw PSA results to wide format
  df_psa_wide <- res$psa_results %>%
    select(sim, Strategy, Cost, QALYs) %>%
    pivot_wider(
      id_cols = sim,
      names_from = Strategy,
      values_from = c(Cost, QALYs),
      names_glue = "{Strategy}_{.value}"
    )
  cat("Wide format PSA data (first few rows):\n")
  print(head(df_psa_wide))
  
  # (2) Convert the PSA samples into a parameters data frame
  parameters_df <- data.frame(
    do.call(rbind, lapply(res$psa_samples, function(sim) {
      sapply(sim, function(x) {
        if (is.numeric(x) && length(x) > 1) {
          mean(x, na.rm = TRUE)
        } else {
          x
        }
      })
    }))
  )
  
  # (3) Create the PSA object using make_psa_obj()
  v_strategies <- c("MenABCWY", "MenACWY", "MenC")
  m_cost <- as.matrix(df_psa_wide[, paste0(v_strategies, "_Cost")])
  m_qalys <- as.matrix(df_psa_wide[, paste0(v_strategies, "_QALYs")])
  
  psa_obj <- make_psa_obj(
    cost = as.data.frame(m_cost),
    effectiveness = as.data.frame(m_qalys),
    parameters = parameters_df,
    strategies = v_strategies,
    currency = "$"
  )
  
  # Plot the PSA object (for a quick visual check)
  cat("Plotting PSA object for scenario:", res$scenario, "\n")
  print(plot(psa_obj))
  
  # (4) EVPI Calculation and Plot
  evpi_obj <- calc_evpi(psa_obj, wtp = seq(0, 100000, 5000))
  evpi_plot <- plot_evpi(evpi_obj) +
    theme_minimal() +
    labs(
      title = paste("EVPI Curve -", res$scenario),
      x = "Cost-Effectiveness Threshold",
      y = "EVPI ($)"
    )
  print(evpi_plot)
  
  # (5) One-Way Sensitivity Analysis (OWSA)
  # Helper function to automatically set ranges using the collapsed parameters
  auto_set_ranges <- function(psa_obj, params_of_interest) {
    psa_paramvals <- psa_obj$parameters
    param_ranges <- list()
    for (param in params_of_interest) {
      if (param %in% colnames(psa_paramvals)) {
        param_ranges[[param]] <- c(min = min(psa_paramvals[[param]], na.rm = TRUE),
                                   max = max(psa_paramvals[[param]], na.rm = TRUE))
      } else {
        stop(paste("Parameter", param, "not found in PSA samples."))
      }
    }
    return(param_ranges)
  }
  
  # Define parameters of interest (adjust as needed)
  params_of_interest <- c("c_MenABCWY", "c_MenACWY", "c_MenC", "p_bg_mort", "p_B", "p_C", "p_B_DeadIMD", "ve_MenABCWY_forACWY")
  param_ranges <- auto_set_ranges(psa_obj, params_of_interest)
  
  owsa_results <- owsa(
    sa_obj    = psa_obj,
    params    = params_of_interest,
    ranges    = param_ranges,
    nsamp     = n_samp,
    outcome   = "nmb",
    wtp       = 50000,
    strategies = v_strategies,
    poly.order = 2
  )
  
  tornado_plot <- owsa_tornado(owsa_results) +
    labs(
      title = paste("One-Way Sensitivity Analysis -", res$scenario),
      x = "Parameter Value",
      y = "Net Monetary Benefit (NMB)"
    ) +
    theme_minimal()
  print(tornado_plot)
  
  # (6) Two-Way Sensitivity Analysis (TWSA)
  # Define candidate parameters (or a subset)
  candidate_params <- c("c_MenABCWY", "c_MenACWY", "c_MenC", "p_bg_mort", "p_B", "p_C", "p_B_DeadIMD", "ve_MenABCWY_forACWY")
  available_params <- colnames(psa_obj$parameters)
  uncertain_params <- intersect(candidate_params, available_params)
  cat("Total number of uncertain parameters available for TWSA in", res$scenario, ":", length(uncertain_params), "\n")
  
  param_pairs <- combn(uncertain_params, 2, simplify = FALSE)
  cat("Total number of unique parameter pairs for TWSA in", res$scenario, ":", length(param_pairs), "\n")
  
  auto_set_twsa_ranges <- function(psa_obj, param1, param2) {
    psa_paramvals <- psa_obj$parameters
    range_pair <- list(
      c(min = min(psa_paramvals[[param1]], na.rm = TRUE),
        max = max(psa_paramvals[[param1]], na.rm = TRUE)),
      c(min = min(psa_paramvals[[param2]], na.rm = TRUE),
        max = max(psa_paramvals[[param2]], na.rm = TRUE))
    )
    names(range_pair) <- c(param1, param2)
    return(range_pair)
  }
  
  nsamp_twsa <- 10  # Number of samples for the TWSA metamodel (adjust as needed)
  twsa_results <- lapply(param_pairs, function(pair) {
    param1 <- pair[1]
    param2 <- pair[2]
    range_pair <- auto_set_twsa_ranges(psa_obj, param1, param2)
    twsa_result <- tryCatch({
      twsa(
        sa_obj     = psa_obj,
        param1     = param1,
        param2     = param2,
        ranges     = range_pair,
        nsamp      = nsamp_twsa,
        outcome    = "nmb",
        wtp        = 50000,
        strategies = v_strategies,
        poly.order = 2
      )
    }, error = function(e) {
      cat("Error for pair:", param1, "vs", param2, "\n", e$message, "\n")
      return(NULL)
    })
    list(pair = pair, result = twsa_result)
  })
  
  # Remove any pairs for which TWSA failed
  twsa_results <- twsa_results[!sapply(twsa_results, function(x) is.null(x$result))]
  cat("Number of successful TWSA pairs in", res$scenario, ":", length(twsa_results), "\n")
  
  # For example, plot the result for the first successful parameter pair
  if (length(twsa_results) > 0) {
    twsa_plot <- plot(twsa_results[[1]]$result) +
      labs(
        title = paste("Two-Way Sensitivity Analysis:", twsa_results[[1]]$pair[1], "vs", twsa_results[[1]]$pair[2], "-", res$scenario),
        x = paste(twsa_results[[1]]$pair[1], "Value"),
        y = paste(twsa_results[[1]]$pair[2], "Value"),
        fill = "Net Monetary Benefit (NMB)"
      ) +
      theme_minimal()
    print(twsa_plot)
  } else {
    cat("No successful TWSA results found in", res$scenario, "\n")
  }
  
  cat("========================================================\n")
}))
