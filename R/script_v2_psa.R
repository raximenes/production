###########################################
# Deterministic Analysis for All Scenarios 
# of the IMD Economic Model
###########################################

# 1) Load necessary scripts (ensure that file paths are correct)
source("R/00_global_setup.R")
source("R/01_parameters_v2.R")
source("R/02_transition_array_v2.R")
source("R/03_markov_model_v2.R")
source("R/04_analysis_v2.R")  # Ensure the updated run_scenario() and related functions are used

# 2) Define the scenarios to run
scenarios_to_run <- c("Scenario A", "Scenario B", "Scenario C")

# 3) Loop through each scenario and run the analysis
results <- lapply(scenarios_to_run, function(scenario_name) {
  cat("Running scenario:", scenario_name, "\n")
  
  # Load base-case parameters for the current scenario
  params_det <- get_basecase_params(scenario_name)
  params_det$scenario_name <- scenario_name
  
  # Select strategies: For scenarios with separate initial/booster (Scenario C or D), 
  # use only "MenABCWY" and "MenACWY"; otherwise include all three strategies.
  if (scenario_name %in% c("Scenario C", "Scenario D")) {
    strategy_list <- c("MenABCWY", "MenACWY")
  } else {
    strategy_list <- c("MenABCWY", "MenACWY", "MenC")
  }
  
  # Run the deterministic simulation to evaluate cost and QALYs for each strategy
  df_det <- simulate_strategies(params_det, strategies = strategy_list, wtp = 50000)
  
  # Run the complete scenario analysis (including ICER calculation)
  res_det <- run_scenario(params_det, strategies = strategy_list, wtp = 50000)
  
  # Format the cost-effectiveness analysis (CEA) table
  table_cea <- format_table_cea(res_det$cea)
  
  # Generate the Markov trace for one strategy (here using "MenABCWY")
  transition_matrix <- build_transition_array(params_det, strategy = "MenABCWY")
  markov_results <- run_markov_model(transition_matrix, params_det)
  
  # Calculate overall survival for "MenABCWY"
  # (Overall survival = 1 - (Dead_IMD + Background_mortality))
  survival_prob <- 1 - rowSums(markov_results[, c("Dead_IMD", "Background_mortality")])
  
  # Return a list with all scenario results
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
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Scenario:", res$scenario, "\n\n")
  cat("Analysis - CEA Results:\n")
  print(res$analysis$cea)
  cat("\n--------------------------------------------------------\n")
}))

# 5) Plot the cost-effectiveness frontier for each scenario with a title including the scenario name
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting CE frontier for scenario:", res$scenario, "\n")
  frontier_plot <- plot(res$analysis$cea, label = "all", txtsize = 12) +
    theme(legend.position = c(0.8, 0.3)) +
    ggtitle(paste("Cost-effectiveness Frontier -", res$scenario))
  print(frontier_plot)
}))

# 6) Plot Markov traces for each scenario
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting Markov trace for scenario:", res$scenario, "\n")
  
  markov_results <- res$markov_trace
  n_cycles <- nrow(markov_results) - 1
  assign("n_cycles", n_cycles, envir = .GlobalEnv)
  
  v_names_states <- colnames(markov_results)[c(1,14,15,2,3,4,5,6,7,8,9,10,11,12,13)]
  assign("v_names_states", v_names_states, envir = .GlobalEnv)
  
  trace_plot <- plot_trace(markov_results) +
    ggtitle(paste("Markov Trace for MenABCWY -", res$scenario))
  
  print(trace_plot)
}))

# 7) Plot overall survival for "MenABCWY" in each scenario
invisible(lapply(results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting overall survival for scenario:", res$scenario, "\n")
  
  plot(res$survival_probability, type = 'l', ylim = c(0, 1),
       xlab = "Cycle (Year)", ylab = "Survival Probability",
       main = paste("Overall Survival - MenABCWY -", res$scenario))
  grid()
  
  cat("Approximate life expectancy for", res$scenario, ":", sum(res$survival_probability), "\n")
  cat("--------------------------------------------------------\n")
}))

###########################################
# Probabilistic Sensitivity Analysis (PSA) 
# for All Scenarios of the IMD Economic Model
###########################################
n_samp <- 100
set.seed(123)

psa_results <- lapply(scenarios_to_run, function(scenario_name) {
  cat("Running PSA for scenario:", scenario_name, "\n")
  
  # Load base-case parameters for the current scenario and tag with scenario name
  params_det <- get_basecase_params(scenario_name)
  params_det$scenario_name <- scenario_name
  assign("scenario_name", scenario_name, envir = .GlobalEnv)
  
  # Generate PSA samples. The generate_psa_samples() function is scenario‐aware.
  samples_list <- generate_psa_samples(n_samp, scenario_name)
  
  # Select strategies based on scenario
  if (scenario_name %in% c("Scenario C", "Scenario D")) {
    strategy_list <- c("MenABCWY", "MenACWY")
  } else {
    strategy_list <- c("MenABCWY", "MenACWY", "MenC")
  }
  
  # Run the PSA using the scenario-specific parameters in parallel mode
  df_psa_results <- run_psa_custom(
    psa_samp = samples_list,
    params_basecase = params_det,
    strategies = strategy_list,
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
  
  # Calculate the PSA CEA based on the summary outcomes
  df_psa_cea <- calculate_icers(
    cost       = df_psa_summary$Cost,
    effect     = df_psa_summary$QALYs,
    strategies = df_psa_summary$Strategy
  )
  
  df_psa_summary$Scenario <- scenario_name
  df_psa_cea$Scenario <- scenario_name
  
  list(
    scenario = scenario_name,
    psa_results = df_psa_results,   # raw PSA simulation results with a "sim" column
    psa_summary = df_psa_summary,
    psa_cea = df_psa_cea,
    samples_list = samples_list
  )
})
names(psa_results) <- scenarios_to_run

# 8) Print PSA CEA results for each scenario
invisible(lapply(psa_results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("PSA Results for scenario:", res$scenario, "\n")
  cat("PSA CEA Results:\n")
  print(res$psa_cea)
  cat("--------------------------------------------------------\n")
}))

# 9) Plot the cost-effectiveness frontier (PSA) for each scenario with title
invisible(lapply(psa_results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting PSA CE frontier for scenario:", res$scenario, "\n")
  frontier_plot <- plot(res$psa_cea, label = "all", txtsize = 12) +
    theme(legend.position = c(0.8, 0.3)) +
    ggtitle(paste("PSA Cost-effectiveness Frontier -", res$scenario))
  print(frontier_plot)
}))

# 10) Generate the Incremental Cost-Effectiveness Plane (ICEP) plot for each scenario
invisible(lapply(psa_results, function(res) {
  cat("--------------------------------------------------------\n")
  cat("Plotting ICEP for scenario:", res$scenario, "\n")
  
  df_psa_results <- res$psa_results
  ce_threshold <- 50000  # $50,000 per QALY
  reference_strategy <- df_psa_results %>%
    group_by(Strategy) %>%
    summarise(Cost = mean(Cost)) %>%
    arrange(Cost) %>%
    slice(1)
  
  df_psa_incremental <- df_psa_results %>%
    group_by(sim) %>%
    mutate(
      Incremental_Cost = Cost - Cost[Strategy == reference_strategy$Strategy],
      Incremental_QALYs = QALYs - QALYs[Strategy == reference_strategy$Strategy]
    ) %>%
    ungroup()
  
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
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
  df_psa_means <- df_psa_incremental %>%
    group_by(Strategy) %>%
    summarise(
      Mean_Incremental_Cost = mean(Incremental_Cost),
      Mean_Incremental_QALYs = mean(Incremental_QALYs)
    )
  
  icep_plot <- icep_plot +
    geom_point(data = df_psa_means, aes(x = Mean_Incremental_QALYs, y = Mean_Incremental_Cost, color = Strategy),
               size = 4, shape = 18) +
    geom_text(data = df_psa_means, aes(x = Mean_Incremental_QALYs, y = Mean_Incremental_Cost, label = Strategy),
              vjust = -1, hjust = 1, size = 4, color = "black") +
    annotate("text", label = paste("Cost-Effectiveness Threshold = $", 
                                   format(ce_threshold, big.mark = ","), "/QALY"),
             color = "red", size = 4, hjust = 0.5, vjust = -0.5, angle = 9)
  
  print(icep_plot)
}))

# 11) Additional PSA Plots for Each Scenario
invisible(lapply(psa_results, function(res) {
  df_psa_results <- res$psa_results
  
  # Plot 2: Cost Distribution by Strategy
  cost_density_plot <- ggplot(df_psa_results, aes(x = Cost, fill = Strategy)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Cost Distribution Across Strategies -", res$scenario),
         x = "Cost ($)", y = "Density", fill = "Strategy") +
    theme_minimal()
  print(cost_density_plot)
  
  # Plot 3: QALY Distribution by Strategy
  qaly_density_plot <- ggplot(df_psa_results, aes(x = QALYs, fill = Strategy)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("QALY Distribution Across Strategies -", res$scenario),
         x = "QALYs", y = "Density", fill = "Strategy") +
    theme_minimal()
  print(qaly_density_plot)
  
  # Plot 4: Cost-Effectiveness Acceptability Curve (CEAC)
  ce_threshold <- 50000
  threshold_values <- seq(0, 3500000, by = 1000)
  ceac_data <- data.frame()
  
  for (threshold in threshold_values) {
    df_nmb <- df_psa_results %>% mutate(NMB = QALYs * threshold - Cost)
    df_cost_effective <- df_nmb %>% group_by(sim) %>%
      summarise(BestStrategy = Strategy[which.max(NMB)]) %>% ungroup()
    df_prob <- df_cost_effective %>% group_by(BestStrategy) %>%
      summarise(Probability = n() / nrow(df_cost_effective)) %>% ungroup() %>%
      complete(BestStrategy = unique(df_psa_results$Strategy), fill = list(Probability = 0)) %>%
      mutate(CostEffectivenessThreshold = threshold)
    ceac_data <- bind_rows(ceac_data, df_prob)
  }
  
  ceac_data <- ceac_data %>% rename(Strategy = BestStrategy)
  
  ceac_plot <- ggplot(ceac_data, aes(x = CostEffectivenessThreshold, y = Probability, color = Strategy)) +
    geom_line(size = 1.2) +
    labs(title = paste("Cost-Effectiveness Acceptability Curve (CEAC) -", res$scenario),
         x = "Cost-Effectiveness Threshold per QALY", y = "Probability of Being Cost-Effective") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
  print(ceac_plot)
}))
