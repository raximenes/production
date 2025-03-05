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



# 8) Plot overall survival for "MenABCWY"
survival_prob <- 1 - rowSums(markov_results[, c("Dead_IMD", "Background_mortality")])
plot(survival_prob, type='l', ylim=c(0,1),
     xlab="Cycle (Year)", ylab="Survival Probability",
     main="Overall Survival - MenABCWY")
grid()
print(sum(survival_prob))  # Approximate life expectancy

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




# 9) Generate PSA samples -------------
n_samp <- 150
set.seed(123)
samples_list <- generate_psa_samples(n_samp)

# 10) Run Probabilistic Sensitivity Analysis (PSA)
df_psa_results <- run_psa_custom(
  psa_samp = samples_list,
  params_basecase = params_det,
  strategies = c("MenABCWY", "MenACWY", "MenC"),
  wtp = 50000,
  parallel = TRUE,
  progress = FALSE
)

# 11) Format PSA results and compute cost-effectiveness analysis
df_psa_summary <- df_psa_results %>%
  group_by(Strategy) %>%
  summarise(
    Cost = mean(Cost),
    QALYs = mean(QALYs)
  )

df_psa_summary <- df_psa_summary %>%
  arrange(Cost)

df_psa_cea <- calculate_icers(
  cost = df_psa_summary$Cost,
  effect = df_psa_summary$QALYs,
  strategies = df_psa_summary$Strategy
)

print(df_psa_cea)
format_table_cea(df_psa_cea)

# 11 a) Plot the cost-effectiveness frontier:
plot(df_psa_cea, label = "all", txtsize = 12) +
  theme(legend.position = c(0.8, 0.3))
############################
# Incremental Cost-Effectiveness Plane (ICEP) with PSA Results -----------
# Define the Willingness-To-Pay (WTP) threshold
ce_threshold <- 50000  # Example: 50,000 per QALY

# Choose the reference strategy (the least expensive one)
reference_strategy <- df_psa_results %>%
  group_by(Strategy) %>%
  summarise(Cost = mean(Cost)) %>%
  arrange(Cost) %>%
  slice(1)  # Selects the strategy with the lowest average cost

# Calculate incremental costs and QALYs relative to the reference strategy for all simulations
df_psa_incremental <- df_psa_results %>%
  group_by(sim) %>%
  mutate(
    Incremental_Cost = Cost - Cost[Strategy == reference_strategy$Strategy],
    Incremental_QALYs = QALYs - QALYs[Strategy == reference_strategy$Strategy]
  ) %>%
  ungroup()

# Plot the Incremental Cost-Effectiveness Plane (ICEP) with all PSA points
icep_plot <- ggplot(df_psa_incremental, aes(x = Incremental_QALYs, y = Incremental_Cost, color = Strategy)) +
  geom_point(alpha = 0.5, size = 1.5) +  # Plot all points with transparency
  geom_abline(slope = ce_threshold, intercept = 0, linetype = "dashed", color = "red") +  # WTP threshold line
  labs(
    title = "Incremental Cost-Effectiveness Plane (ICEP) - PSA Results",
    x = "Incremental QALYs",
    y = "Incremental Cost",
    color = "Strategy"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

# Add strategy means as highlighted points
df_psa_means <- df_psa_incremental %>%
  group_by(Strategy) %>%
  summarise(
    Mean_Incremental_Cost = mean(Incremental_Cost),
    Mean_Incremental_QALYs = mean(Incremental_QALYs)
  )

icep_plot <- icep_plot +
  geom_point(data = df_psa_means, aes(x = Mean_Incremental_QALYs, y = Mean_Incremental_Cost, color = Strategy),
             size = 4, shape = 18) +  # Highlighted points for means
  geom_text(data = df_psa_means, aes(x = Mean_Incremental_QALYs, y = Mean_Incremental_Cost, label = Strategy),
            vjust = -1, hjust = 1, size = 4, color = "black") +  # Labels for means
  annotate(
    "text",
    label = paste("Cost-Effectiveness Threshold = $", format(ce_threshold, big.mark = ","), "/QALY"),
    color = "red",
    size = 4,
    hjust = 0.5,
    vjust = -0.5,
    angle =  9
  ) 

# Display the plot
print(icep_plot)




# ====== Plot 2: Cost Distribution by Strategy ======
cost_density_plot <- ggplot(df_psa_results, aes(x = Cost, fill = Strategy)) +
  geom_density(alpha = 0.5) +  # Density plot to show cost variability
  labs(
    title = "Cost Distribution Across Strategies",
    x = "Cost ($)",
    y = "Density",
    fill = "Strategy"
  ) +
  theme_minimal()
# Display the plots
print(cost_density_plot)

# ====== Plot 3: QALY Distribution by Strategy ======
qaly_density_plot <- ggplot(df_psa_results, aes(x = QALYs, fill = Strategy)) +
  geom_density(alpha = 0.5) +  # Density plot to show QALY variability
  labs(
    title = "QALY Distribution Across Strategies",
    x = "QALYs",
    y = "Density",
    fill = "Strategy"
  ) +
  theme_minimal()

# Display the plots
print(qaly_density_plot)



# Cost-effectiveness Acceptability Curve (CEAC) ------------------------------

# Define a range of cost-effectiveness thresholds
threshold_values <- seq(0, 200000, by = 1000)

# Initialize a data frame to store the CEAC (Cost-Effectiveness Acceptability Curve) results
ceac_data <- data.frame()

# Loop over each cost-effectiveness threshold
for (threshold in threshold_values) {
  
  # Calculate the Net Monetary Benefit (NMB) for each strategy in each simulation
  df_nmb <- df_psa_results %>%
    mutate(NMB = QALYs * threshold - Cost)
  
  # Identify the cost-effective strategy (highest NMB) in each simulation
  df_cost_effective <- df_nmb %>%
    group_by(sim) %>%
    summarise(
      BestStrategy = Strategy[which.max(NMB)]
    ) %>%
    ungroup()
  
  # Calculate the probability of each strategy being cost-effective
  df_prob <- df_cost_effective %>%
    group_by(BestStrategy) %>%
    summarise(
      Probability = n() / nrow(df_cost_effective)
    ) %>%
    ungroup() %>%
    
    # Ensure all strategies appear, even if their probability is 0 at this threshold
    complete(
      BestStrategy = unique(df_psa_results$Strategy),
      fill = list(Probability = 0)
    ) %>%
    
    # Add the current threshold value
    mutate(CostEffectivenessThreshold = threshold)
  
  # Append to the final CEAC dataset
  ceac_data <- bind_rows(ceac_data, df_prob)
}

# Rename the column to "Strategy" (optional, for clarity)
ceac_data <- ceac_data %>%
  rename(Strategy = BestStrategy)

# Plot the CEAC
ceac_plot <- ggplot(ceac_data, 
                    aes(x = CostEffectivenessThreshold, 
                        y = Probability, 
                        color = Strategy)) +
  geom_line(size = 1.2) +
  labs(
    title = "Cost-Effectiveness Acceptability Curve (CEAC)",
    x = "Cost-Effectiveness Threshold per QALY",
    y = "Probability of Being Cost-Effective"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",        # Place legend at the bottom
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  )

# Display the plot
print(ceac_plot)










# 12) Wide data
df_psa_wide <- df_psa_results %>%
  select(sim, Strategy, Cost, QALYs) %>%
  pivot_wider(
    id_cols     = sim,
    names_from  = Strategy,
    values_from = c(Cost, QALYs),
    names_glue  = "{Strategy}_{.value}"
  )

head(df_psa_wide)

# Convert samples_list (a list of parameter lists) into a data frame.
# For any parameter that is a vector (length > 1), collapse it by taking the mean.
parameters_df <- data.frame(
  do.call(rbind, lapply(samples_list, function(sim) {
    sapply(sim, function(x) {
      if (is.numeric(x) && length(x) > 1) {
        mean(x, na.rm = TRUE)
      } else {
        x
      }
    })
  }))
)

# --- Create the PSA Object ---
# Here, df_psa_wide is assumed to have been created earlier with cost and QALYs for each simulation.
v_strategies <- c("MenABCWY", "MenACWY", "MenC")
m_cost <- as.matrix(df_psa_wide[, paste0(v_strategies, "_Cost")])
m_qalys <- as.matrix(df_psa_wide[, paste0(v_strategies, "_QALYs")])

psa_obj <- make_psa_obj(
  cost = as.data.frame(m_cost),               # Cost data frame
  effectiveness = as.data.frame(m_qalys),     # Effectiveness data frame (e.g., QALYs)
  parameters = parameters_df,                 # Parameters data frame (now all scalars)
  strategies = v_strategies,                  # Strategy names
  currency = "$"                              # Currency symbol
)
plot(psa_obj)

# # 13) Cost-effectiveness acceptability curve (CEAC)
#ceac_obj <- ceac(psa_obj, wtp = seq(0, 100000, 5000))
#plot(ceac_obj)

# 14) Expected Value of Perfect Information (EVPI) calculation
evpi_obj <- calc_evpi(psa_obj, wtp = seq(0, 100000, 5000))
plot_evpi(evpi_obj) +
  theme_minimal() +
  labs(
    title = "EVPI Curve",
    x = "Cost-Effectiveness",
    y = "EVPI ($)"
  )

###############################################################################
# OWSA Code with Time-Dependent Parameters Collapsed to Scalars
#
# Description:
#   In this code, we first convert the PSA samples (samples_list) into a 
#   parameters data frame. For parameters that are vectors (i.e. time-dependent
#   parameters such as "p_B", "p_bg_mort", etc.), we take the mean value across
#   cycles so that each parameter is represented by a single scalar value per 
#   simulation. Then we create the PSA object using make_psa_obj() and update the
#   auto_set_ranges() function to derive ranges based on these scalar values.
#   Finally, we run the one-way sensitivity analysis (OWSA) and generate a tornado
#   plot.
###############################################################################

# --- Step 1: Automatically Set Ranges for OWSA, Using the Collapsed Parameters ---
auto_set_ranges <- function(psa_obj, params_of_interest) {
  # Extract parameter values from the PSA object's parameters data frame.
  psa_paramvals <- psa_obj$parameters
  
  param_ranges <- list()
  for (param in params_of_interest) {
    # Because we collapsed time-dependent parameters to scalars, an exact match is expected.
    if (param %in% colnames(psa_paramvals)) {
      param_ranges[[param]] <- c(min = min(psa_paramvals[[param]], na.rm = TRUE),
                                 max = max(psa_paramvals[[param]], na.rm = TRUE))
    } else {
      stop(paste("Parameter", param, "not found in PSA samples."))
    }
  }
  return(param_ranges)
}

# Define the parameters of interest (using the base names)
params_of_interest <- c(
  "c_MenABCWY", "c_MenACWY", "c_MenC", "c_admin", 
  "p_bg_mort", "p_B", "p_C", "p_W", "p_Y", 
  "p_B_DeadIMD", "p_C_DeadIMD", "p_W_DeadIMD", "p_Y_DeadIMD", 
  "p_IMD_Scarring", "p_IMD_Single_Amput", "p_IMD_Multiple_Amput", 
  "p_IMD_Neuro_Disability", "p_IMD_Hearing_Loss", "p_IMD_Renal_Failure", 
  "p_IMD_Seizure", "p_IMD_Paralysis", 
  "ve_MenABCWY_forACWY", "ve_MenABCWY_forB", "ve_MenACWY", "ve_MenC", 
  "c_IMD_infection", "c_Scarring", "c_Single_Amput", "c_Multiple_Amput", 
  "c_Neuro_Disab", "c_Hearing_Loss", "c_Renal_Failure", "c_Seizure", 
  "c_Paralysis", "c_Dead", "c_Healthy", 
  "u_Healthy", "u_IMD", "u_Scarring", "u_Single_Amput", "u_Multiple_Amput", 
  "u_Neuro_Disability", "u_Hearing_Loss", "u_Renal_Failure", "u_Seizure", 
  "u_Paralysis", "u_Dead"
)

param_ranges <- auto_set_ranges(psa_obj, params_of_interest)

# --- Step 2: Run the One-Way Sensitivity Analysis (OWSA) ---
owsa_results <- owsa(
  sa_obj    = psa_obj,         # PSA object containing cost and effectiveness samples
  params    = params_of_interest,  # Parameter names of interest (base names)
  ranges    = param_ranges,        # Ranges automatically set from the PSA samples
  nsamp     = n_samp,              # Number of samples for the OWSA (e.g., same as PSA sample size)
  outcome   = "nmb",               # Outcome of interest (e.g., Net Monetary Benefit)
  wtp       = 50000,               # Willingness-to-pay threshold
  strategies = v_strategies,       # Strategies to consider
  poly.order = 2                   # Polynomial order for the metamodel
)

# --- Step 3: Plot the Tornado Plot for OWSA Results ---
owsa_tornado(owsa_results) +
  labs(
    title = "One-Way Sensitivity Analysis",
    x = "Parameter Value",
    y = "Net Monetary Benefit (NMB)"
  ) +
  theme_minimal()

# Plot a tornado plot displaying only parameters with at least a 10% relative impact,
# with full color and base text size 12.
tornado_plot <- owsa_tornado(
  owsa = owsa_results,         # owsa object
  return = "plot",
  txtsize = 12,
  min_rel_diff = 0.91,         # filter: only parameters with >=10% relative impact
  col = "full"                 # full color plot
)

# Display the plot
print(tornado_plot)
###########################
# Extract the tornado data as a data frame.
tornado_data <- owsa_tornado(
  owsa = owsa_results,
  return = "data",
  txtsize = 12,
  min_rel_diff = 0,  # initially extract all parameters
  col = "full",
  n_y_ticks = 8
)

# Check the structure of the tornado data.
str(tornado_data)
# (Assume the data frame contains columns, for example, "parameter" and "rel_diff", 
# where "rel_diff" is the relative change in outcome due to changes in that parameter.)

# Order the data by the absolute relative difference (largest first).
ordered_tornado_data <- tornado_data %>%
  arrange(desc(abs(rel_diff)))

# Optionally, select only the top N most relevant parameters (e.g., top 10).
top_n <- 15
top_tornado_data <- head(ordered_tornado_data, n = top_n)

# Now create a tornado plot for these top parameters using ggplot2.
library(ggplot2)
tornado_plot_manual <- ggplot(top_tornado_data, aes(x = reorder(parameter, abs(rel_diff)), y = rel_diff, fill = parameter)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  labs(
    title = paste("Tornado Plot: Top", top_n, "Most Relevant Parameters"),
    x = "Parameter",
    y = "Relative Change in Outcome (NMB)"
  ) +
  theme_minimal() +
  guides(fill = FALSE)

# Display the plot
print(tornado_plot_manual)




##############################################################################

###############################################################################
# Multi-Parameter Two-Way Sensitivity Analysis (TWSA) for Many Parameter Pairs
#
# Description:
#   This script performs a two-way sensitivity analysis on all unique pairs
#   of uncertain parameters from the PSA object (psa_obj). It assumes that the
#   PSA object's parameters data frame contains only scalar values (i.e. any 
#   time-dependent parameters have been collapsed by taking their mean). 
#
#   For each unique pair of parameters, the script:
#     1. Automatically computes the range (minimum and maximum) from the PSA
#        samples using the actual parameter names.
#     2. Runs the TWSA using the twsa() function.
#     3. Stores the result along with the parameter pair in a list.
#
#   An example plot is then produced for the first parameter pair.
#
# Note:
#   The PSA object (psa_obj) must already have been created (e.g., using make_psa_obj())
#   and should have a parameters data frame with the following columns (example):
#     c_MenABCWY, c_MenACWY, c_MenC, c_admin, p_bg_mort, p_B, p_C, p_W, p_Y, 
#     p_B_DeadIMD, p_C_DeadIMD, p_W_DeadIMD, p_Y_DeadIMD, p_IMD_Scarring, ...,
#     ve_MenABCWY_forACWY, ve_MenABCWY_forB, ve_MenACWY, ve_MenC, c_IMD_infection, 
#     n_cycles, cycle_length, d_c, d_e, coverage, c_Dead, c_Healthy, u_Healthy, u_Dead, etc.
###############################################################################

# --- Step 1: Define the Set of Uncertain Parameters for TWSA ---
# (Exclude deterministic parameters such as n_cycles, cycle_length, d_c, d_e, coverage, c_Dead, c_Healthy, u_Healthy, and u_Dead.)
candidate_params <- c(
  "c_MenABCWY", "c_MenACWY", "c_MenC", "c_admin", 
  "p_bg_mort", "p_B", "p_C", "p_W", "p_Y", 
  "p_B_DeadIMD", "p_C_DeadIMD", "p_W_DeadIMD", "p_Y_DeadIMD", 
  "p_IMD_Scarring", "p_IMD_Single_Amput", "p_IMD_Multiple_Amput", 
  "p_IMD_Neuro_Disability", "p_IMD_Hearing_Loss", "p_IMD_Renal_Failure", 
  "p_IMD_Seizure", "p_IMD_Paralysis", 
  "ve_MenABCWY_forACWY", "ve_MenABCWY_forB", "ve_MenACWY", "ve_MenC", 
  "c_IMD_infection", "c_Scarring", "c_Single_Amput", "c_Multiple_Amput", 
  "c_Neuro_Disab", "c_Hearing_Loss", "c_Renal_Failure", "c_Seizure", 
  "c_Paralysis", 
  "u_IMD", "u_Scarring", "u_Single_Amput", "u_Multiple_Amput", 
  "u_Neuro_Disability", "u_Hearing_Loss", "u_Renal_Failure", "u_Seizure", 
  "u_Paralysis"
)
candidate_params <- c(
  "c_MenABCWY", "c_MenACWY", "c_MenC", "p_bg_mort", "p_B", "p_C", 
  "p_B_DeadIMD","ve_MenABCWY_forACWY"
)
available_params <- colnames(psa_obj$parameters)
uncertain_params <- intersect(candidate_params, available_params)
cat("Total number of uncertain parameters available:", length(uncertain_params), "\n")

# --- Step 2: Generate All Unique Parameter Pairs ---
param_pairs <- combn(uncertain_params, 2, simplify = FALSE)
cat("Total number of unique parameter pairs for TWSA:", length(param_pairs), "\n")

# --- Step 3: Helper Function to Automatically Set Ranges for Two Parameters ---
# IMPORTANT CORRECTION: Use the actual parameter names as list names.
auto_set_twsa_ranges <- function(psa_obj, param1, param2) {
  psa_paramvals <- psa_obj$parameters  # Data frame with one row per simulation
  range_pair <- list(
    c(min = min(psa_paramvals[[param1]], na.rm = TRUE),
      max = max(psa_paramvals[[param1]], na.rm = TRUE)),
    c(min = min(psa_paramvals[[param2]], na.rm = TRUE),
      max = max(psa_paramvals[[param2]], na.rm = TRUE))
  )
  names(range_pair) <- c(param1, param2)
  return(range_pair)
}

# --- Step 4: Set TWSA Settings ---
v_strategies <- c("MenABCWY", "MenACWY", "MenC")
nsamp_twsa <- 10  # Number of samples for the TWSA metamodel (adjust as needed)

# --- Step 5: Loop Over Each Parameter Pair and Run TWSA ---
results_list <- lapply(param_pairs, function(pair) {
  param1 <- pair[1]
  param2 <- pair[2]
  
  # Get the range for these two parameters using the corrected helper function.
  range_pair <- auto_set_twsa_ranges(psa_obj, param1, param2)
  
  # Run the two-way sensitivity analysis.
  twsa_result <- tryCatch({
    twsa(
      sa_obj     = psa_obj,       # PSA object with cost/effectiveness samples
      param1     = param1,        # First parameter of the pair
      param2     = param2,        # Second parameter of the pair
      ranges     = range_pair,    # Ranges for these parameters (with proper names)
      nsamp      = nsamp_twsa,    # Number of samples for the metamodel
      outcome    = "nmb",         # Outcome of interest (e.g., Net Monetary Benefit)
      wtp        = 50000,         # Willingness-to-pay threshold
      strategies = v_strategies,  # Strategies to consider
      poly.order = 2            # Polynomial order for the metamodel
    )
  }, error = function(e) {
    cat("Error for pair:", param1, "vs", param2, "\n", e$message, "\n")
    return(NULL)
  })
  
  # Return a list element containing the parameter pair and its TWSA result.
  list(pair = pair, result = twsa_result)
})

# Remove any pairs for which TWSA failed.
results_list <- results_list[!sapply(results_list, function(x) is.null(x$result))]
cat("Number of successful TWSA pairs:", length(results_list), "\n")

# create a data frame summarizing the pairs:
pair_summary <- data.frame(
  Pair_Index = seq_along(results_list),
  Parameter1 = sapply(results_list, function(x) x$pair[1]),
  Parameter2 = sapply(results_list, function(x) x$pair[2]),
  stringsAsFactors = FALSE
)
print(pair_summary)
pair_summary[1,,]

# --- Step 6: Example Plot for the First Parameter Pair ---
if (length(results_list) > 0) {
  plot(results_list[[1]]$result) +
    labs(
      title = paste("Two-Way Sensitivity Analysis:", 
                    results_list[[1]]$pair[1], "vs", results_list[[1]]$pair[2]),
      x = paste(results_list[[1]]$pair[1], "Value"),
      y = paste(results_list[[1]]$pair[2], "Value"),
      fill = "Net Monetary Benefit (NMB)"
    ) +
    theme_minimal()
} else {
  cat("No successful TWSA results found.\n")
}



# (Optional) You can now further process results_list (e.g., save plots, summarize results, etc.)
