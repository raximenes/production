# --------------------------------------------------------
# Deterministic Analysis for All Scenarios of the IMD Economic Model
# --------------------------------------------------------

# 1) Load necessary scripts (ensure that file paths are correct)
source("R/00_global_setup.R")
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
    survival_probability = survival_prob,
    params_det = params_det
    
  )
})
names(results) <- scenarios_to_run
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

# 8) Plot overall survival for "MenABCWY" in each scenario
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

# --------------------------------------------------------
# Probabilistic Sensitivity Analysis (PSA) for All Scenarios of the IMD Economic Model
# --------------------------------------------------------
n_samp <- 100
set.seed(123)

psa_results <- lapply(scenarios_to_run, function(scenario_name) {
  cat("Running PSA for scenario:", scenario_name, "\n")
  
  # Load base-case parameters for the current scenario and tag with scenario name
  params_det <- get_basecase_params(scenario_name)
  params_det$scenario_name <- scenario_name
  assign("scenario_name", scenario_name, envir = .GlobalEnv)
  
  # Generate PSA samples (make sure generate_psa_samples() is scenario‐aware if needed)
  samples_list <- generate_psa_samples(n_samp, scenario_name)
  
  # Run the PSA using the scenario-specific parameters
  df_psa_results <- run_psa_custom(
    psa_samp = samples_list,
    params_basecase = params_det,
    strategies = c("MenABCWY", "MenACWY", "MenC"),
    wtp = 50000,
    parallel = TRUE,
    progress = FALSE
  )
  
  # Summarise PSA outcomes by strategy
  df_psa_summary <- df_psa_results %>%
    group_by(Strategy) %>%
    summarise(
      Cost = mean(Cost),
      QALYs = mean(QALYs)
    ) %>%
    arrange(Cost)
  
  # Calculate the PSA CEA based on the summary (if desired)
  df_psa_cea <- calculate_icers(
    cost = df_psa_summary$Cost,
    effect = df_psa_summary$QALYs,
    strategies = df_psa_summary$Strategy
  )
  
  # Tag with scenario name
  df_psa_summary$Scenario <- scenario_name
  df_psa_cea$Scenario <- scenario_name
  
  list(
    scenario = scenario_name,
    psa_results = df_psa_results,   # raw PSA simulation results (with a "sim" column)
    psa_summary = df_psa_summary,
    psa_cea = df_psa_cea,
    samples_list = samples_list
  )
})

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
########################


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
  threshold_values <- seq(0, 3500000, by = 1000)
  
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
###################



