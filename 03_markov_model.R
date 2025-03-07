###############################################################################
# 03_markov_model.R (Updated)
#
# Purpose:
#   This script defines the core Markov modeling functions for a model that
#   explicitly separates states by vaccination status, including:
#     1) run_markov_model() -> computes the distribution of the cohort across
#        states over the model cycles.
#     2) make_cost_matrix() -> constructs the (n_cycles+1) x n_states cost matrix,
#        with age-specific costs for all IMD infection states, for both vaccinated
#        and unvaccinated.
#     3) make_utility_vector() -> creates the vector of utilities for each state.
#     4) compute_costs_and_qalys() -> multiplies the Markov trace by cost & utility,
#        applies discounting, and returns total discounted costs & QALYs.
#
# Important:
#   - We now have separate healthy states: "Healthy_vaccinated" & "Healthy_unvaccinated".
#   - We also have 8 infection states total, e.g.: "SeroB_Infect_vacc"/"SeroB_Infect_unvacc", etc.
#   - Sequela states remain the same, and we have two death states.
#
# Make sure these are consistent with the states in build_transition_array().
###############################################################################

# (1) run_markov_model():
#     Evolve the distribution of a cohort across the states for n_cycles.
run_markov_model <- function(a_P, params) {
  n_cycles <- params$n_cycles
  
  # Initial distribution: coverage fraction starts in Healthy_vaccinated,
  # the remainder in Healthy_unvaccinated, everything else zero.
  v_m_init <- c(
    "Healthy_vaccinated"   = params$coverage,
    "Healthy_unvaccinated" = 1 - params$coverage,
    
    # Infection states
    "SeroB_Infect_vacc"   = 0,
    "SeroB_Infect_unvacc" = 0,
    "SeroC_Infect_vacc"   = 0,
    "SeroC_Infect_unvacc" = 0,
    "SeroW_Infect_vacc"   = 0,
    "SeroW_Infect_unvacc" = 0,
    "SeroY_Infect_vacc"   = 0,
    "SeroY_Infect_unvacc" = 0,
    
    # Sequela states
    "Scarring" = 0,
    "Single_Amput" = 0,
    "Multiple_Amput" = 0,
    "Neuro_Disability" = 0,
    "Hearing_Loss" = 0,
    "Renal_Failure" = 0,
    "Seizure" = 0,
    "Paralysis" = 0,
    
    # Death states
    "Dead_IMD" = 0,
    "Background_mortality" = 0
  )
  
  # Storage matrix for Markov trace (one row per cycle: cycle_0 to cycle_n)
  m_M <- matrix(0, nrow = n_cycles + 1, ncol = length(v_m_init),
                dimnames = list(paste0("cycle_", 0:n_cycles), names(v_m_init)))
  m_M[1, ] <- v_m_init
  
  # For each cycle, multiply the previous distribution by the transition array
  for (t in 1:n_cycles) {
    m_M[t + 1, ] <- m_M[t, ] %*% a_P[, , t]
  }
  
  return(m_M)
}

# (2) make_cost_matrix():
#     Creates an (n_cycles+1) x n_states matrix of costs, allowing for age/time-
#     dependent infection cost (c_IMD_infection), and adding vaccine cost at
#     cycle 0 to Healthy_vaccinated only.
make_cost_matrix <- function(params, strategy) {
  n_cycles <- params$n_cycles
  
  # State vector must match build_transition_array()
  state_names <- c(
    "Healthy_vaccinated",
    "Healthy_unvaccinated",
    
    "SeroB_Infect_vacc", "SeroB_Infect_unvacc",
    "SeroC_Infect_vacc", "SeroC_Infect_unvacc",
    "SeroW_Infect_vacc", "SeroW_Infect_unvacc",
    "SeroY_Infect_vacc", "SeroY_Infect_unvacc",
    
    "Scarring", "Single_Amput", "Multiple_Amput", "Neuro_Disability",
    "Hearing_Loss", "Renal_Failure", "Seizure", "Paralysis",
    
    "Dead_IMD", "Background_mortality"
  )
  
  # Initialize cost matrix
  m_cost <- matrix(NA, nrow = n_cycles + 1, ncol = length(state_names),
                   dimnames = list(paste0("cycle_", 0:n_cycles), state_names))
  
  # (A) Assign constant costs to Healthy, Sequela, and Death states.
  # If any of these are time-varying, adapt similarly to infection states.
  m_cost[, "Healthy_vaccinated"]   <- params$c_Healthy
  m_cost[, "Healthy_unvaccinated"] <- params$c_Healthy
  
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
  
  # (B) Assign age-specific cost to infection states.
  # We have 8 infection states: Sero{B,C,W,Y}_Infect_{vacc, unvacc}.
  # By assumption, cycle 0 has zero cost for infection states.
  inf_states <- c(
    "SeroB_Infect_vacc", "SeroB_Infect_unvacc",
    "SeroC_Infect_vacc", "SeroC_Infect_unvacc",
    "SeroW_Infect_vacc", "SeroW_Infect_unvacc",
    "SeroY_Infect_vacc", "SeroY_Infect_unvacc"
  )
  
  # Set cycle 0 cost to 0 for infection states
  m_cost[1, inf_states] <- 0
  
  # If c_IMD_infection has length >= n_cycles, it implies we have a time-specific cost.
  if (length(params$c_IMD_infection) >= n_cycles) {
    for (t in 2:(n_cycles + 1)) {
      # t-1 indexes the cost vector for cycle t.
      # The same cost is assigned to all infection states.
      m_cost[t, inf_states] <- params$c_IMD_infection[t - 1]
    }
  } else {
    # Fallback: if c_IMD_infection is not as long as needed, assume the first element.
    m_cost[, inf_states] <- params$c_IMD_infection[1]
    m_cost[1, inf_states] <- 0  # Ensure cycle 0 is zero.
  }
  
  # (C) Add vaccine cost at cycle 0 to "Healthy_vaccinated".
  # The coverage aspect is in run_markov_model, but the cost is only paid by those starting vaccinated.
  # So we do not multiply by coverage here (the state distribution accounts for coverage initially).
  if (strategy == "MenABCWY") {
    m_cost[1, "Healthy_vaccinated"] <- m_cost[1, "Healthy_vaccinated"] + (params$c_MenABCWY + params$c_admin)
  } else if (strategy == "MenACWY") {
    m_cost[1, "Healthy_vaccinated"] <- m_cost[1, "Healthy_vaccinated"] + (params$c_MenACWY + params$c_admin)
  } else if (strategy == "MenC") {
    m_cost[1, "Healthy_vaccinated"] <- m_cost[1, "Healthy_vaccinated"] + (params$c_MenC + params$c_admin)
  }
  
  return(m_cost)
}

# (3) make_utility_vector():
#     Creates a named vector of utilities for each state.
make_utility_vector <- function(params) {
  
  # Must match the states used above.
  v_util <- c(
    "Healthy_vaccinated"   = params$u_Healthy,
    "Healthy_unvaccinated" = params$u_Healthy,
    
    # 8 infection states => all have utility of IMD infection.
    "SeroB_Infect_vacc"   = params$u_IMD,
    "SeroB_Infect_unvacc" = params$u_IMD,
    "SeroC_Infect_vacc"   = params$u_IMD,
    "SeroC_Infect_unvacc" = params$u_IMD,
    "SeroW_Infect_vacc"   = params$u_IMD,
    "SeroW_Infect_unvacc" = params$u_IMD,
    "SeroY_Infect_vacc"   = params$u_IMD,
    "SeroY_Infect_unvacc" = params$u_IMD,
    
    # Sequela states
    "Scarring"         = params$u_Scarring,
    "Single_Amput"     = params$u_Single_Amput,
    "Multiple_Amput"   = params$u_Multiple_Amput,
    "Neuro_Disability" = params$u_Neuro_Disability,
    "Hearing_Loss"     = params$u_Hearing_Loss,
    "Renal_Failure"    = params$u_Renal_Failure,
    "Seizure"          = params$u_Seizure,
    "Paralysis"        = params$u_Paralysis,
    
    # Death states
    "Dead_IMD"             = params$u_Dead,
    "Background_mortality" = params$u_Dead
  )
  
  return(v_util)
}

# (4) compute_costs_and_qalys():
#     Multiplies Markov trace by the cost matrix and utility vector, then applies
#     discounting and returns total cost & QALYs.
compute_costs_and_qalys <- function(m_trace, params, strategy) {
  n_cycles <- params$n_cycles
  
  # Build the cost matrix and utility vector
  m_cost <- make_cost_matrix(params, strategy)
  v_util <- make_utility_vector(params)
  
  # Discounting vectors for costs & QALYs
  v_dwc <- 1 / ((1 + params$d_c)^(0:n_cycles))  # cost discount factors
  v_dwe <- 1 / ((1 + params$d_e)^(0:n_cycles))  # QALY discount factors
  
  # Weighted cycle correction factor (from darthtools)
  # e.g.: v_wcc <- gen_wcc(n_cycles, method = "Simpson1/3")
  # or, if not available, simply:
  v_wcc <- rep(1, n_cycles + 1)
  
  cycle_costs <- numeric(n_cycles + 1)
  cycle_qalys <- numeric(n_cycles + 1)
  
  # For each cycle, cost = sum of state probabilities × state cost.
  # Similarly for QALYs.
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
