#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write functions to run linear and linear mixed effects models for development & stress analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### LM model (no random effects) ###

f_lm_model <- function(dataframe,
                       independentVar,
                       dependentVar,
                       covariates
) {
  
  # write formula based on variables
  lm_formula <- paste0("scale(", dependentVar, ") ~ ", independentVar, " + ", paste0(covariates, collapse = " + "))
  
  # run LM model
  lm_model <- lm(as.formula(lm_formula), data = dataframe)
  
  # extract results
  lm_results <- tidy(lm_model) %>% 
    clean_names %>% 
    rename("t_value" = "statistic")
  
  # return model & results as list
  return(list("model" = lm_model, "results" = lm_results))
  
}

### LME model (random effects) ###

f_lme_model <- function(dataframe, 
                        independentVar,
                        dependentVar,
                        covariates,
                        randomEffects
) {
  
  # define random effect term based on number of effects specified
  randomEffects_formula <- map(randomEffects, ~ paste0("( 1 | ", .x, " )")) %>% paste0(collapse = " + ")
  
  # write formula based on variables
  lme_formula <- paste0("scale(", dependentVar, ") ~ ", independentVar, " + ", paste0(covariates, collapse = " + "), " + ", randomEffects_formula)
  
  # run LME model
  lme_model <- lme4::lmer(as.formula(lme_formula), 
                          data = dataframe,
                          REML = FALSE,
                          control = lmerControl(optimizer ="Nelder_Mead")
  )
  
  # extract results
  lme_results <- lme_model %>% 
    summary %>% 
    coefficients %>% 
    as.data.frame() %>% 
    rownames_to_column("term") %>% 
    clean_names()
  
  # return model & results as list
  return(list("model" = lme_model, "results" = lme_results))
  
}