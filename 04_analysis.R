###############################################################################
# 04_analysis.R (English version)
#
# Purpose:
#   1) simulate_strategies() -> runs the model for each strategy and calculates
#      cost, QALYs, and NMB (Net Monetary Benefit).
#   2) run_scenario() -> uses simulate_strategies() to generate cost-effectiveness
#      results, then computes ICERs using dampack.
#   3) run_psa_custom():
#       - Runs the model for each PSA sample (in parallel or sequentially).
#       - Allows measuring runtime.
#       - Displays a text progress bar in sequential mode (optional).
#
# Important:
#   - Make sure fill_sero_transition(), build_transition_array(), etc.,
#     are already defined and sourced before using these functions.
#   - When running in parallel, we export user-defined functions and load
#     necessary packages on each worker.
###############################################################################

# 1) simulate_strategies():
#    Runs the Markov model for each of the given strategies, returning a data
#    frame of cost, QALYs, and net monetary benefit (NMB) at a chosen WTP.
simulate_strategies <- function(l_params,
                                strategies = c("MenABCWY", "MenACWY", "MenC"),
                                wtp = 50000) {
  n_strat <- length(strategies)
  df_ce <- data.frame(
    Strategy = strategies,
    Cost     = numeric(n_strat),
    QALYs    = numeric(n_strat),
    NMB      = numeric(n_strat),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(strategies)) {
    strat_i <- strategies[i]
    # Build the transition array for this strategy
    a_P     <- build_transition_array(l_params, strat_i)
    # Run the Markov model
    m_trace <- run_markov_model(a_P, l_params)
    # Calculate total cost and QALYs
    res     <- compute_costs_and_qalys(m_trace, l_params, strat_i)
    
    # Store results
    df_ce$Cost[i]  <- res$cost
    df_ce$QALYs[i] <- res$qalys
    df_ce$NMB[i]   <- df_ce$QALYs[i] * wtp - df_ce$Cost[i]
  }
  return(df_ce)
}


# 2) run_scenario():
#    Runs simulate_strategies() for the specified strategies, then calculates ICERs
#    using the 'dampack' package. Returns a list containing:
#      - 'raw': the data frame of cost and QALYs
#      - 'cea': the cost-effectiveness analysis table (ICERs, etc.)
run_scenario <- function(l_params,
                         strategies = c("MenABCWY", "MenACWY", "MenC"),
                         wtp = 50000) {
  # 1) Run the strategies
  df_raw <- simulate_strategies(l_params, strategies, wtp)
  
  # 2) Load dampack to compute ICERs
  library(dampack)
  
  # 3) Calculate ICERs
  df_cea <- calculate_icers(
    cost       = df_raw$Cost,
    effect     = df_raw$QALYs,
    strategies = df_raw$Strategy
  )
  
  # 4) Return results
  out <- list(
    raw = df_raw[, c("Strategy", "Cost", "QALYs")],
    cea = df_cea
  )
  return(out)
}


# 3) run_psa_custom():
#    Performs a probabilistic sensitivity analysis (PSA) over multiple samples
#    of parameters. It can run either sequentially (with optional progress bar)
#    or in parallel (without progress messages to avoid conflicts).
#
#    Arguments:
#      - psa_samp: list of parameter sets (one element per simulation).
#      - params_basecase: base-case parameters (not strictly necessary here,
#                         but sometimes used for reference).
#      - strategies: vector of strategies to simulate.
#      - wtp: willingness-to-pay threshold for NMB.
#      - parallel: whether to run in parallel or sequential mode.
#      - n_cores: number of cores (defaults to detectCores()).
#      - progress: whether to show a progress bar (only in sequential mode).
#
#    Returns:
#      - df_out: a data frame of PSA results, with one row per simulation
#        per strategy (including cost, QALYs, and NMB).
run_psa_custom <- function(psa_samp,
                           params_basecase,
                           strategies = c("MenABCWY", "MenACWY", "MenC"),
                           wtp = 50000,
                           parallel = FALSE,
                           n_cores = parallel::detectCores(),
                           progress = TRUE) {
  #--------------------------------------
  # 0) Preliminary checks
  #--------------------------------------
  if (!is.list(psa_samp) || length(psa_samp) == 0) {
    stop("No PSA samples provided. Check the 'psa_samp' object.")
  }
  
  n_sim <- length(psa_samp)
  message("Starting PSA with ", n_sim, " simulations...")
  if (parallel) {
    message("Parallel mode is ON. Using ", n_cores, " cores.")
  } else {
    message("Parallel mode is OFF (sequential).")
  }
  
  # Start timer
  start_time <- Sys.time()
  
  #--------------------------------------
  # (A) SEQUENTIAL MODE
  #--------------------------------------
  if (!parallel) {
    if (progress) {
      pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
    }
    
    df_out <- data.frame(
      sim      = integer(0),
      Strategy = character(0),
      Cost     = numeric(0),
      QALYs    = numeric(0),
      NMB      = numeric(0),
      stringsAsFactors = FALSE
    )
    
    # Loop through each PSA sample
    for (i in seq_len(n_sim)) {
      params_i <- psa_samp[[i]]
      
      # Run all strategies for this sample
      df_ce_i <- simulate_strategies(params_i, strategies, wtp)
      df_ce_i$sim <- i
      
      # Append to results
      df_out <- rbind(df_out, df_ce_i[, c("sim","Strategy","Cost","QALYs","NMB")])
      
      # Update progress bar
      if (progress) setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
    
    # End timer
    end_time <- Sys.time()
    total_time_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
    message(sprintf("PSA completed in %.2f seconds (%.2f minutes).", 
                    total_time_sec, total_time_sec / 60))
    
    return(df_out)
    
  } else {
    #--------------------------------------
    # (B) PARALLEL MODE
    #--------------------------------------
    # No progress bar or messages here to avoid conflicts.
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # We need 'darthtools' and other functions for each worker:
    df_out <- tryCatch({
      foreach(
        i = seq_len(n_sim),
        .combine  = rbind,
        .packages = c("dplyr", "darthtools"), 
        .export   = c(
          "simulate_strategies",
          "build_transition_array",
          "fill_sero_transition",
          "run_markov_model",
          "compute_costs_and_qalys",
          "make_cost_matrix",
          "make_utility_vector"
        )
      ) %dopar% {
        params_i <- psa_samp[[i]]
        df_ce_i  <- simulate_strategies(params_i, strategies, wtp)
        df_ce_i$sim <- i
        
        # Return result for this simulation
        df_ce_i[, c("sim","Strategy","Cost","QALYs","NMB")]
      }
    }, finally = {
      # Stop cluster and revert to sequential
      stopCluster(cl)
      registerDoSEQ()
    })
    
    # End timer
    end_time <- Sys.time()
    total_time_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
    message(sprintf("PSA completed in %.2f seconds (%.2f minutes).", 
                    total_time_sec, total_time_sec / 60))
    
    return(df_out)
  }
}
