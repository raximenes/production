###############################################################################
# 03_markov_model.R
#
# Purpose:
#   This script defines the core Markov modeling functions:
#     1) run_markov_model() -> computes the state distribution over cycles,
#     2) make_cost_matrix() -> constructs the cost matrix for the model states,
#        correctly using time‐dependent (age‐specific) cost for IMD infection states,
#     3) make_utility_vector() -> creates the state utility vector, and
#     4) compute_costs_and_qalys() -> computes the total discounted cost and QALYs.
#
# Note on Corrections:
#   The original make_cost_matrix() function assumed that all cost parameters were 
#   scalars. Since c_IMD_infection is now a vector (age‐specific cost for infection),
#   the function has been modified so that for each cycle (except cycle 0, which is 
#   set to 0), the cost for infection states ("SeroB_Infect", "SeroC_Infect",
#   "SeroW_Infect", and "SeroY_Infect") is taken from the corresponding element 
#   of params$c_IMD_infection.
###############################################################################

# 1) Run the Markov chain
run_markov_model <- function(a_P, params) {
  n_cycles <- params$n_cycles
  
  v_m_init <- c(
    "Healthy" = 1,
    "SeroB_Infect" = 0, "SeroC_Infect" = 0, "SeroW_Infect" = 0, "SeroY_Infect" = 0,
    "Scarring" = 0, "Single_Amput" = 0, "Multiple_Amput" = 0, "Neuro_Disability" = 0,
    "Hearing_Loss" = 0, "Renal_Failure" = 0, "Seizure" = 0, "Paralysis" = 0,
    "Dead_IMD" = 0, "Background_mortality" = 0
  )
  
  m_M <- matrix(0, nrow = n_cycles + 1, ncol = length(v_m_init),
                dimnames = list(paste0("cycle_", 0:n_cycles), names(v_m_init)))
  m_M[1, ] <- v_m_init
  
  for (t in 1:n_cycles) {
    m_M[t + 1, ] <- m_M[t, ] %*% a_P[, , t]
  }
  
  return(m_M)
}

# 2) Corrected Cost & Utility Structures

# Corrected make_cost_matrix() function:
# This version creates a cost matrix with one row per cycle.
# For states with constant cost, the same value is assigned to every cycle.
# For the infection states (SeroB_Infect, SeroC_Infect, SeroW_Infect, SeroY_Infect),
# the cost is age-specific. We assume that cycle 0 (the initial state) has zero infection cost.
make_cost_matrix <- function(params, strategy) {
  n_cycles <- params$n_cycles
  state_names <- c("Healthy", "SeroB_Infect", "SeroC_Infect", "SeroW_Infect", "SeroY_Infect",
                   "Scarring", "Single_Amput", "Multiple_Amput", "Neuro_Disability",
                   "Hearing_Loss", "Renal_Failure", "Seizure", "Paralysis",
                   "Dead_IMD", "Background_mortality")
  
  # Initialize an empty matrix with (n_cycles+1) rows.
  m_cost <- matrix(NA, nrow = n_cycles + 1, ncol = length(state_names),
                   dimnames = list(paste0("cycle_", 0:n_cycles), state_names))
  
  # For states with constant costs, assign the same value to every cycle.
  m_cost[, "Healthy"]              <- params$c_Healthy
  m_cost[, "Scarring"]             <- params$c_Scarring
  m_cost[, "Single_Amput"]         <- params$c_Single_Amput
  m_cost[, "Multiple_Amput"]       <- params$c_Multiple_Amput
  m_cost[, "Neuro_Disability"]     <- params$c_Neuro_Disab
  m_cost[, "Hearing_Loss"]         <- params$c_Hearing_Loss
  m_cost[, "Renal_Failure"]        <- params$c_Renal_Failure
  m_cost[, "Seizure"]              <- params$c_Seizure
  m_cost[, "Paralysis"]            <- params$c_Paralysis
  m_cost[, "Dead_IMD"]             <- params$c_Dead
  m_cost[, "Background_mortality"] <- params$c_Dead
 
  # For the infection states, assign age-specific cost.
  # Assume that cycle 0 has no infection cost.
  m_cost[1, c("SeroB_Infect", "SeroC_Infect", "SeroW_Infect", "SeroY_Infect")] <- 0
  # Check if the vector length of c_IMD_infection matches the number of cycles.
  if (length(params$c_IMD_infection) >= n_cycles) {
    # For cycles 1 to n_cycles, assign the corresponding cost.
    for (t in 2:(n_cycles + 1)) {
      m_cost[t, "SeroB_Infect"] <- params$c_IMD_infection[t - 1]
      m_cost[t, "SeroC_Infect"] <- params$c_IMD_infection[t - 1]
      m_cost[t, "SeroW_Infect"] <- params$c_IMD_infection[t - 1]
      m_cost[t, "SeroY_Infect"] <- params$c_IMD_infection[t - 1]
    }
  } else {
    # If c_IMD_infection is not of expected length, fall back to a constant cost.
    m_cost[, c("SeroB_Infect", "SeroC_Infect", "SeroW_Infect", "SeroY_Infect")] <- params$c_IMD_infection[1]
  }
  
  # Add the vaccine cost in the specified cycle to the Healthy state.
  if(strategy == "MenABCWY") {
    m_cost[params$vaccination_time + 1, "Healthy"] <- m_cost[params$vaccination_time + 1, "Healthy"] + params$coverage * (params$c_MenABCWY + params$c_admin) * params$doses
  } else if(strategy == "MenACWY") {
    m_cost[params$vaccination_time + 1, "Healthy"] <- m_cost[params$vaccination_time + 1, "Healthy"] + params$coverage * (params$c_MenACWY + params$c_admin) * params$doses
  } else if(strategy == "MenC") {
    m_cost[params$vaccination_time + 1, "Healthy"] <- m_cost[params$vaccination_time + 1, "Healthy"] + params$coverage * (params$c_MenC + params$c_admin) * params$doses
  }
  return(m_cost)
}

# Construct the utility vector.
make_utility_vector <- function(params) {
  v_util <- c(
    "Healthy"              = params$u_Healthy,
    "SeroB_Infect"         = params$u_IMD,
    "SeroC_Infect"         = params$u_IMD,
    "SeroW_Infect"         = params$u_IMD,
    "SeroY_Infect"         = params$u_IMD,
    "Scarring"             = params$u_Scarring,
    "Single_Amput"         = params$u_Single_Amput,
    "Multiple_Amput"       = params$u_Multiple_Amput,
    "Neuro_Disability"     = params$u_Neuro_Disability,
    "Hearing_Loss"         = params$u_Hearing_Loss,
    "Renal_Failure"        = params$u_Renal_Failure,
    "Seizure"              = params$u_Seizure,
    "Paralysis"            = params$u_Paralysis,
    "Dead_IMD"             = params$u_Dead,
    "Background_mortality" = params$u_Dead
  )
  return(v_util)
}

# 3) Compute total discounted cost and QALYs.
compute_costs_and_qalys <- function(m_trace, params, strategy) {
  n_cycles <- params$n_cycles
  
  m_cost <- make_cost_matrix(params, strategy)
  v_util <- make_utility_vector(params)
  
  v_dwc <- 1 / ((1 + params$d_c)^(0:n_cycles))  # cost discount factors
  v_dwe <- 1 / ((1 + params$d_e)^(0:n_cycles))  # QALY discount factors
  
  # Weighted cycle correction factor (using Simpson's 1/3 method if available)
  v_wcc <- gen_wcc(n_cycles, method = "Simpson1/3")
  # Alternatively, for simplicity: v_wcc <- rep(1, n_cycles + 1)
  
  cycle_costs <- numeric(n_cycles + 1)
  cycle_qalys <- numeric(n_cycles + 1)
  
  for (t in 1:(n_cycles + 1)) {
    cycle_costs[t] <- sum(m_trace[t, ] * m_cost[t, ])
    cycle_qalys[t] <- sum(m_trace[t, ] * v_util)
  }
  
  discounted_costs <- cycle_costs * v_dwc * v_wcc
  discounted_qalys <- cycle_qalys * v_dwe * v_wcc
  
  total_cost  <- sum(discounted_costs)
  total_qalys <- sum(discounted_qalys)
  
  return(list(cost = total_cost, qalys = total_qalys))
}
################################################################################