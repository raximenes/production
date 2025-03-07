###############################################
# 02_transition_array.R
#
# Purpose:
#   - Define functions to build the transition probability array for a Markov model
#     that tracks vaccination status explicitly with separate "Healthy_vaccinated"
#     and "Healthy_unvaccinated" states, as well as separate infection states for
#     vaccinated vs. unvaccinated individuals.
#
# The key functions are:
#   1) fill_sero_transition()
#       - Fills transitions for an infection state (e.g. "SeroB_Infect_vacc").
#   2) get_infection_prob()
#       - Returns the infection probability for a given serogroup, strategy, and
#         vaccination status in cycle t.
#   3) build_transition_array()
#       - Builds the 3D transition probability array (n_states x n_states x n_cycles).
#
# Explanation:
#   - Instead of a single "Healthy" state, we now have "Healthy_vaccinated" and
#     "Healthy_unvaccinated" so que possamos associar riscos de infecção distintos.
#   - Para manter a consistência do status vacinal, as infecções também são
#     duplicadas em: "SeroB_Infect_vacc" vs. "SeroB_Infect_unvacc", e assim por
#     diante para os demais sorogrupos.
#   - Com isso, quando o paciente sobrevive à infecção sem sequela, ele retorna
#     ao estado "Healthy_vaccinated" ou "Healthy_unvaccinated" correto.
###############################################

# (1) fill_sero_transition():
#  Fills the row of the transition matrix (a_P) for a single infection state,
#  por exemplo "SeroB_Infect_vacc". Precisamos saber para qual estado "Healthy_*"
#  retornar se o indivíduo sobrevive sem sequela.
#
#  Args:
#   - a_P:                3D array de transição
#   - from_state:         string, e.g. "SeroB_Infect_vacc"
#   - to_healthy_state:   string, e.g. "Healthy_vaccinated"
#   - p_bg_mort:          probabilidade de mortalidade de base no ciclo t
#   - p_dead_imd:         probabilidade de morrer por IMD (condicional a não ser morte de base)
#   - p_seq_list:         vetor nomeado com as probabilidades de sequela
#   - sum_sq:             soma das probabilidades de sequela
#   - t:                  índice do ciclo
fill_sero_transition <- function(a_P,
                                 from_state,
                                 to_healthy_state,
                                 p_bg_mort,
                                 p_dead_imd,
                                 p_seq_list,
                                 sum_sq,
                                 t) {
  # Prob de sobreviver à infecção por IMD (condicional a não ser mortalidade de base)
  surv <- (1 - p_bg_mort) * (1 - p_dead_imd)
  
  # Transição para 'Background_mortality' (morte de base)
  a_P[from_state, "Background_mortality", t] <- p_bg_mort
  
  # Transição para 'Dead_IMD' (morte por IMD, se não foi morte de base)
  a_P[from_state, "Dead_IMD", t] <- (1 - p_bg_mort) * p_dead_imd
  
  # Transições para cada sequela
  a_P[from_state, "Scarring", t]         <- surv * p_seq_list["Scarring"]
  a_P[from_state, "Single_Amput", t]     <- surv * p_seq_list["Single_Amput"]
  a_P[from_state, "Multiple_Amput", t]   <- surv * p_seq_list["Multiple_Amput"]
  a_P[from_state, "Neuro_Disability", t] <- surv * p_seq_list["Neuro_Disability"]
  a_P[from_state, "Hearing_Loss", t]     <- surv * p_seq_list["Hearing_Loss"]
  a_P[from_state, "Renal_Failure", t]    <- surv * p_seq_list["Renal_Failure"]
  a_P[from_state, "Seizure", t]          <- surv * p_seq_list["Seizure"]
  a_P[from_state, "Paralysis", t]        <- surv * p_seq_list["Paralysis"]
  
  # Se não houver sequela, o indivíduo volta para o estado saudável
  a_P[from_state, to_healthy_state, t] <- surv * (1 - sum_sq)
  
  return(a_P)
}

# (2) get_infection_prob():
#   Retorna a probabilidade de infecção para um sorogrupo (B, C, W ou Y),
#   considerando a estratégia (MenABCWY, MenACWY ou MenC), o ciclo t, e se o
#   indivíduo é vacinado ou não (vaccinated = TRUE/FALSE).
get_infection_prob <- function(params, strategy, group, t, vaccinated) {
  # Probabilidade de infecção de base para esse sorogrupo
  p_base <- params[[paste0("p_", group)]][t]
  
  # Se o indivíduo é vacinado, aplicamos a eficácia vacinal apropriada
  if (vaccinated) {
    if (strategy == "MenABCWY") {
      if (group == "B") {
        # Eficácia específica para B em MenABCWY
        return(p_base * (1 - params$ve_MenABCWY_forB[t]))
      } else {
        # Serogrupos C, W, Y usam ve_MenABCWY_forACWY
        return(p_base * (1 - params$ve_MenABCWY_forACWY[t]))
      }
    } else if (strategy == "MenACWY") {
      if (group == "B") {
        # MenACWY não cobre B, então prob = p_base
        return(p_base)
      } else {
        # Para C, W, Y => aplica ve_MenACWY
        return(p_base * (1 - params$ve_MenACWY[t]))
      }
    } else if (strategy == "MenC") {
      if (group == "B") {
        # MenC não cobre B
        return(p_base)
      } else if (group == "C") {
        # Redução pela eficácia de MenC
        return(p_base * (1 - params$ve_MenC[t]))
      } else {
        # W ou Y => não há cobertura
        return(p_base)
      }
    } else {
      stop("Estratégia desconhecida. Use: 'MenABCWY', 'MenACWY' ou 'MenC'.")
    }
    
  } else {
    # Não vacinado => não há redução alguma (probabilidade = p_base)
    return(p_base)
  }
}

# (3) build_transition_array():
#   Constrói o array 3D de probabilidade de transição com states x states x cycles,
#   diferenciando estados saudáveis e infectados por status vacinal.
build_transition_array <- function(params, strategy) {
  n_cycles <- params$n_cycles
  
  # Definimos o vetor de nomes de estados:
  #   - 2 estados saudáveis (vac/unvac)
  #   - 8 estados de infecção (4 sorogrupos × 2 status vacinais)
  #   - 8 estados de sequela
  #   - 2 estados de morte (IMD e background)
  
  v_names_states <- c(
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
  n_states <- length(v_names_states)
  
  # Cria o array de transição (tudo começa em zero)
  a_P <- array(0,
               dim = c(n_states, n_states, n_cycles),
               dimnames = list(v_names_states, v_names_states,
                               paste0("cycle_", 1:n_cycles)))
  
  # Extrair vetores de parâmetros relevantes
  p_bg_mort   <- params$p_bg_mort   # mortalidade de base (vetor)
  p_B_DeadIMD <- params$p_B_DeadIMD
  p_C_DeadIMD <- params$p_C_DeadIMD
  p_W_DeadIMD <- params$p_W_DeadIMD
  p_Y_DeadIMD <- params$p_Y_DeadIMD
  
  # Soma das probabilidades de sequela
  sum_sq <- (params$p_IMD_Scarring +
               params$p_IMD_Single_Amput +
               params$p_IMD_Multiple_Amput +
               params$p_IMD_Neuro_Disability +
               params$p_IMD_Hearing_Loss +
               params$p_IMD_Renal_Failure +
               params$p_IMD_Seizure +
               params$p_IMD_Paralysis)
  
  # Vetor nomeado de probabilidades de sequela (para facilitar o loop)
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
  
  # Grupos de soros para facilitar os loops
  sero_groups <- c("B", "C", "W", "Y")
  
  # Lista de mortalidade por IMD (por sorogrupo)
  p_deadIMD_list <- list(
    B = p_B_DeadIMD,
    C = p_C_DeadIMD,
    W = p_W_DeadIMD,
    Y = p_Y_DeadIMD
  )
  
  # Preencher transições para cada ciclo
  for (t in 1:n_cycles) {
    bgm <- p_bg_mort[t]  # mortalidade de base no ciclo t
    
    # (A) Transições a partir de estados saudáveis (vac/unvac)
    for (vacc_status in c(TRUE, FALSE)) {
      
      from_healthy <- if (vacc_status) {
        "Healthy_vaccinated"
      } else {
        "Healthy_unvaccinated"
      }
      
      # Morte de base
      a_P[from_healthy, "Background_mortality", t] <- bgm
      
      # Probabilidade de infecção por cada sorogrupo
      pB <- get_infection_prob(params, strategy, "B", t, vaccinated = vacc_status)
      pC <- get_infection_prob(params, strategy, "C", t, vaccinated = vacc_status)
      pW <- get_infection_prob(params, strategy, "W", t, vaccinated = vacc_status)
      pY <- get_infection_prob(params, strategy, "Y", t, vaccinated = vacc_status)
      
      total_inf <- pB + pC + pW + pY
      
      # Alocar transições para os estados de infecção correspondentes
      if (vacc_status) {
        # Se é vacinado
        a_P[from_healthy, "SeroB_Infect_vacc", t] <- (1 - bgm) * pB
        a_P[from_healthy, "SeroC_Infect_vacc", t] <- (1 - bgm) * pC
        a_P[from_healthy, "SeroW_Infect_vacc", t] <- (1 - bgm) * pW
        a_P[from_healthy, "SeroY_Infect_vacc", t] <- (1 - bgm) * pY
      } else {
        # Se é não vacinado
        a_P[from_healthy, "SeroB_Infect_unvacc", t] <- (1 - bgm) * pB
        a_P[from_healthy, "SeroC_Infect_unvacc", t] <- (1 - bgm) * pC
        a_P[from_healthy, "SeroW_Infect_unvacc", t] <- (1 - bgm) * pW
        a_P[from_healthy, "SeroY_Infect_unvacc", t] <- (1 - bgm) * pY
      }
      
      # Permanecer saudável caso não ocorra infecção nem morte de base
      a_P[from_healthy, from_healthy, t] <- (1 - bgm) * (1 - total_inf)
    }
    
    # (B) Transições a partir dos estados de infecção (vac/unvac)
    for (group in sero_groups) {
      # Probabilidade de morrer por IMD para esse sorogrupo e ciclo
      this_p_dead_imd <- p_deadIMD_list[[group]][t]
      
      # from_state VAC
      vac_infect_name <- paste0("Sero", group, "_Infect_vacc")
      a_P <- fill_sero_transition(
        a_P,
        from_state       = vac_infect_name,
        to_healthy_state = "Healthy_vaccinated",
        p_bg_mort        = bgm,
        p_dead_imd       = this_p_dead_imd,
        p_seq_list       = p_seq_list,
        sum_sq           = sum_sq,
        t                = t
      )
      
      # from_state UNVAC
      unvac_infect_name <- paste0("Sero", group, "_Infect_unvacc")
      a_P <- fill_sero_transition(
        a_P,
        from_state       = unvac_infect_name,
        to_healthy_state = "Healthy_unvaccinated",
        p_bg_mort        = bgm,
        p_dead_imd       = this_p_dead_imd,
        p_seq_list       = p_seq_list,
        sum_sq           = sum_sq,
        t                = t
      )
    }
    
    # (C) Transições a partir dos estados de sequela
    #     (a prob de mortalidade de base se aplica e, se não morrer, fica onde está).
    
    sequelae_states <- c(
      "Scarring","Single_Amput","Multiple_Amput","Neuro_Disability",
      "Hearing_Loss","Renal_Failure","Seizure","Paralysis"
    )
    
    for (st in sequelae_states) {
      a_P[st, "Background_mortality", t] <- bgm
      a_P[st, st, t] <- 1 - bgm
    }
    
    # (D) Estados de morte são absorventes
    a_P["Dead_IMD","Dead_IMD", t]                     <- 1
    a_P["Background_mortality","Background_mortality", t] <- 1
  }
  
  # Retorna o array de transição completamente preenchido
  return(a_P)
}
