###############################################################################
# 04_analysis.R
#
# Purpose:
#   1) simulate_strategies()
#   2) run_scenario()
#   3) run_psa_custom():
#       - Parallel or sequential
#       - Automatic core detection
#       - Start/end timer to measure run time
#       - Text progress bar in sequential mode (if progress=TRUE)
#       - NO progress bar/messages in parallel mode (to avoid errors)
#
# NOTE on 'gen_wcc()' from darthtools:
#   We must load 'darthtools' in each worker. Hence we do .packages=c("dplyr","darthtools").
#
# IMPORTANT:
#   - Make sure that fill_sero_transition(), build_transition_array(), etc.
#     are defined in your other scripts and those scripts are already sourced
#     before running 'run_psa_custom()'.
###############################################################################

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
    a_P     <- build_transition_array(l_params, strat_i)
    m_trace <- run_markov_model(a_P, l_params)
    res     <- compute_costs_and_qalys(m_trace, l_params, strat_i)
    df_ce$Cost[i]  <- res$cost
    df_ce$QALYs[i] <- res$qalys
    df_ce$NMB[i]   <- df_ce$QALYs[i] * wtp - df_ce$Cost[i]
  }
  return(df_ce)
}


run_scenario <- function(l_params,
                         strategies = c("MenABCWY", "MenACWY", "MenC"),
                         wtp = 50000) {
  # Get the base-case parameters for the scenario
  #l_params <- get_basecase_params(scenario_name)
  
  # Simulate the strategies
  df_raw <- simulate_strategies(l_params, strategies, wtp)
 # library(dampack)  # ensure dampack is loaded (for calculate_icers)
  df_raw$Scenario <- l_params$scenario_name
  
  df_cea <- calculate_icers(
    cost       = df_raw$Cost,
    effect     = df_raw$QALYs,
    strategies = df_raw$Strategy
  )
  
  df_cea$Scenario <- l_params$scenario_name
  
  out <- list(
    raw = df_raw[, c("Scenario", "Strategy", "Cost", "QALYs")],
    cea = df_cea
  )
  return(out)
}


run_psa_custom <- function(psa_samp,
                           params_basecase,
                           strategies = c("MenABCWY", "MenACWY", "MenC"),
                           wtp = 50000,
                           parallel = FALSE,
                           n_cores = parallel::detectCores(),
                           progress = TRUE) {
  # Get the base-case parameters for the scenario
  #l_params <- get_basecase_params(scenario_name)
  l_params <- params_basecase
  #--------------------------------------
  # 0) Setup
  #--------------------------------------
  # Check if psa_samp is a list and has at least one sample
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
  
  # Timer: start
  start_time <- Sys.time()
  
  #--------------------------------------
  # (A) SEQUENTIAL MODE
  #--------------------------------------
  if (!parallel) {
    # Create a text progress bar in sequential mode if progress=TRUE
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
    
    for (i in seq_len(n_sim)) {
      params_i <- psa_samp[[i]]  # Use the i-th PSA sample directly
      
      df_ce_i <- simulate_strategies(params_i, strategies, wtp)
      df_ce_i$sim <- i
      
      df_out <- rbind(df_out, df_ce_i[, c("sim","Strategy","Cost","QALYs","NMB")])
      
      # Update progress bar in sequential mode
      if (progress) setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
    
    # Timer: end
    end_time <- Sys.time()
    total_time_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
    message(sprintf("PSA completed in %.2f seconds (%.2f minutes).", 
                    total_time_sec, total_time_sec/60))
    
    return(df_out)
    
  } else {
    #--------------------------------------
    # (B) PARALLEL MODE
    #--------------------------------------
    # No progress bar/messages here to avoid errors/conflicts.
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Because we need 'gen_wcc()' from 'darthtools' inside compute_costs_and_qalys(),
    # we must include "darthtools" in the .packages argument. 
    # Also, we export all user-defined functions used in the Markov logic.
    
    df_out <- tryCatch({
      foreach(
        i = seq_len(n_sim),
        .combine  = rbind,
        .packages = c("dplyr","darthtools"), 
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
        params_i <- psa_samp[[i]]  # Use the i-th PSA sample directly
        
        df_ce_i <- simulate_strategies(params_i, strategies, wtp)
        df_ce_i$sim <- i
        
        # Return result
        df_ce_i[, c("sim","Strategy","Cost","QALYs","NMB")]
      }
    }, finally = {
      # Ensure the cluster is stopped
      stopCluster(cl)
      registerDoSEQ()  # Restore sequential processing
    })
    
    # Timer: end
    end_time <- Sys.time()
    total_time_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
    message(sprintf("PSA completed in %.2f seconds (%.2f minutes).", 
                    total_time_sec, total_time_sec/60))
    
    return(df_out)
  }
}
################################################################################