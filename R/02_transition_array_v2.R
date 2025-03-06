#######################################
# 02_transition_array.R
#
# Purpose:
#   - Define the function build_transition_array() that creates the transition
#     probability array given the set of parameters (or a single draw from the PSA).
#   - Potentially define helper functions like fill_sero_transition().
###############################################################################

fill_sero_transition <- function(a_P,
                                 from_state,
                                 p_bg_mort,
                                 p_dead_imd,
                                 p_seq_list,
                                 sum_sq,
                                 t) {
  # Probability of surviving IMD (conditional on not background mortality)
  surv <- (1 - p_bg_mort) * (1 - p_dead_imd)
  
  a_P[from_state, "Background_mortality", t] <- p_bg_mort
  a_P[from_state, "Dead_IMD", t]             <- (1 - p_bg_mort) * p_dead_imd
  
  a_P[from_state, "Scarring", t]             <- surv * p_seq_list["Scarring"]
  a_P[from_state, "Single_Amput", t]         <- surv * p_seq_list["Single_Amput"]
  a_P[from_state, "Multiple_Amput", t]       <- surv * p_seq_list["Multiple_Amput"]
  a_P[from_state, "Neuro_Disability", t]     <- surv * p_seq_list["Neuro_Disability"]
  a_P[from_state, "Hearing_Loss", t]         <- surv * p_seq_list["Hearing_Loss"]
  a_P[from_state, "Renal_Failure", t]        <- surv * p_seq_list["Renal_Failure"]
  a_P[from_state, "Seizure", t]              <- surv * p_seq_list["Seizure"]
  a_P[from_state, "Paralysis", t]            <- surv * p_seq_list["Paralysis"]
  
  # If no sequela, return to Healthy
  a_P[from_state, "Healthy", t] <- surv * (1 - sum_sq)
  
  return(a_P)
}

build_transition_array <- function(params, strategy) {
  # This version is similar to your original code but adapted as a function
  # that uses a single "params" list (which can be for base-case or a single PSA draw).
  
  n_cycles <- params$n_cycles
  v_names_states <- c(
    "Healthy",
    "SeroB_Infect","SeroC_Infect","SeroW_Infect","SeroY_Infect",
    "Scarring","Single_Amput","Multiple_Amput","Neuro_Disability",
    "Hearing_Loss","Renal_Failure","Seizure","Paralysis",
    "Dead_IMD","Background_mortality"
  )
  n_states <- length(v_names_states)
  
  a_P <- array(0, dim=c(n_states, n_states, n_cycles),
               dimnames = list(v_names_states, v_names_states,
                               paste0("cycle_", 1:n_cycles)))
  
  # Extract needed param vectors
  p_bg_mort   <- params$p_bg_mort
  p_B_DeadIMD <- params$p_B_DeadIMD
  p_C_DeadIMD <- params$p_C_DeadIMD
  p_W_DeadIMD <- params$p_W_DeadIMD
  p_Y_DeadIMD <- params$p_Y_DeadIMD
  
  sum_sq <- (params$p_IMD_Scarring +
               params$p_IMD_Single_Amput +
               params$p_IMD_Multiple_Amput +
               params$p_IMD_Neuro_Disability +
               params$p_IMD_Hearing_Loss +
               params$p_IMD_Renal_Failure +
               params$p_IMD_Seizure +
               params$p_IMD_Paralysis)
  
  p_seq_list <- c(
    Scarring         = params$p_IMD_Scarring,
    Single_Amput     = params$p_IMD_Single_Amput,
    Multiple_Amput   = params$p_IMD_Multiple_Amput,
    Neuro_Disability = params$p_IMD_Neuro_Disability,
    Hearing_Loss     = params$p_IMD_Hearing_Loss,
    Renal_Failure    = params$p_IMD_Renal_Failure,
    Seizure          = params$p_IMD_Seizure,
    Paralysis        = params$p_IMD_Paralysis
  )
  
  # Helper: get infection prob depending on strategy
  # get_infection_prob <- function(t, group) {
  #   # This logic matches your original approach:
  #   if (strategy == "MenABCWY") {
  #     if (group == "B") {
  #       return(params$p_B[t] * (1 - params$ve_MenABCWY_forB[t]))
  #     } else {
  #       # For serogroups C, W, Y => uses ve_MenABCWY_forACWY
  #       return(params[[paste0("p_", group)]][t] * (1 - params$ve_MenABCWY_forACWY[t]))
  #     }
  #   } else if (strategy == "MenACWY") {
  #     if (group == "B") {
  #       return(params$p_B[t])
  #     } else {
  #       return(params[[paste0("p_", group)]][t] *
  #                (1 - params$ve_MenACWY[t]))
  #     }
  #   } else if (strategy == "MenC") {
  #     if (group == "B") {
  #       return(params$p_B[t])
  #     } else if (group == "C") {
  #       return(params$p_C[t] * (1 - params$ve_MenC[t]))
  #     } else if (group == "W") {
  #       return(params$p_W[t])
  #     } else {
  #       return(params$p_Y[t])
  #     }
  #   } else {
  #     stop("Unknown strategy. Must be one of MenABCWY, MenACWY, MenC.")
  #   }
  # }
  get_infection_prob <- function(t, group, strategy, params) {
    if (strategy == "MenABCWY") {
      if (group == "B") {
        # If ve_MenABCWY_forB is not defined, you can default it (e.g., to 0)
        ve_val <- if (!is.null(params$ve_MenABCWY_forB)) params$ve_MenABCWY_forB[t] else 0
        return(params$p_B[t] * (1 - ve_val))
      } else {
        return(params[[paste0("p_", group)]][t] * (1 - params$ve_MenABCWY_forACWY[t]))
      }
    } else if (strategy == "MenACWY") {
      if (group == "B") {
        return(params$p_B[t])
      } else {
        return(params[[paste0("p_", group)]][t] * (1 - params$ve_MenACWY[t]))
      }
    } else if (strategy == "MenC") {
      if (group == "B") {
        return(params$p_B[t])
      } else if (group == "C") {
        return(params$p_C[t] * (1 - params$ve_MenC[t]))
      } else if (group == "W") {
        return(params$p_W[t])
      } else {
        return(params$p_Y[t])
      }
    } else {
      stop("Unknown strategy. Must be one of MenABCWY, MenACWY, MenC.")
    }
  }
  
  # Fill transitions
  for (t in 1:n_cycles) {
    bgm <- p_bg_mort[t]
    pB  <- get_infection_prob(t, "B", strategy, params)
    pC  <- get_infection_prob(t, "C", strategy, params)
    pW  <- get_infection_prob(t, "W", strategy, params)
    pY  <- get_infection_prob(t, "Y", strategy, params)
    total_inf <- pB + pC + pW + pY
    
    # FROM Healthy
    a_P["Healthy","Background_mortality", t] <- bgm
    a_P["Healthy","SeroB_Infect", t]         <- (1 - bgm) * pB
    a_P["Healthy","SeroC_Infect", t]         <- (1 - bgm) * pC
    a_P["Healthy","SeroW_Infect", t]         <- (1 - bgm) * pW
    a_P["Healthy","SeroY_Infect", t]         <- (1 - bgm) * pY
    a_P["Healthy","Healthy", t]              <- (1 - bgm) * (1 - total_inf)
    
    # Sero B
    a_P <- fill_sero_transition(a_P, "SeroB_Infect", bgm, p_B_DeadIMD[t],
                                p_seq_list, sum_sq, t)
    # Sero C
    a_P <- fill_sero_transition(a_P, "SeroC_Infect", bgm, p_C_DeadIMD[t],
                                p_seq_list, sum_sq, t)
    # Sero W
    a_P <- fill_sero_transition(a_P, "SeroW_Infect", bgm, p_W_DeadIMD[t],
                                p_seq_list, sum_sq, t)
    # Sero Y
    a_P <- fill_sero_transition(a_P, "SeroY_Infect", bgm, p_Y_DeadIMD[t],
                                p_seq_list, sum_sq, t)
    
    # Sequela states
    sequelae_states <- c("Scarring","Single_Amput","Multiple_Amput","Neuro_Disability",
                         "Hearing_Loss","Renal_Failure","Seizure","Paralysis")
    for (st in sequelae_states) {
      a_P[st, "Background_mortality", t] <- bgm
      a_P[st, st, t] <- 1 - bgm
    }
    
    # Dead states => absorbing
    a_P["Dead_IMD","Dead_IMD", t]                     <- 1
    a_P["Background_mortality","Background_mortality", t] <- 1
  }
  
  # (Optional) checks: row sums, etc. 
  # Validates all probabilities in the transition array a_P are between 0 and 1 (inclusive) and 
  #   that each row sums to 1 (within a small tolerance, typically around 1e-6 or 1e-7).
  # If you have darthtools installed, you could do:Checks the entire 3D array (n_states × n_states × n_cycles) at once, covering all cycles and states.
  # check_transition_probability(a_P, verbose = TRUE)
  ### OR ###
  # Validate that each row sums to 1, warning if any deviate significantly.
  #row_sums <- apply(a_P[,,t], 1, sum)
 # if (any(abs(row_sums - 1) > 1e-6)) {
   # warning("Transition probabilities in cycle ", t, " do not sum to 1 for states: ",
   #         paste(names(row_sums[abs(row_sums - 1) > 1e-6]), collapse = ", "))
 # }
  
  return(a_P)
}
################################################################################
###############################################################################