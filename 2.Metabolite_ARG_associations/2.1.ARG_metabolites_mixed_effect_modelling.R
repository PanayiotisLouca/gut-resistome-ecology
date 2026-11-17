## Author: 
  # Panayiotis Louca 

## Purpose of script: 
  #  

## Purpose of script: 
  #  Perform mixed effect modelling to test associations between ARGs and bacterial metabolites 

## Clear environment 
  rm(list = ls()) 

## Set seed 
  set.seed(1234)

## Set functions: 
      
## load up packages: 

  ### core 
  library(tidyverse)
  library(readxl)
    
  # mixed effect modelling 
  library(glmmTMB)
  library(lmerTest)
  library(broom.mixed)
  
# -------------------------------------------------------------------------- # 

# ----------------- # 
#   Import data  ---- 
# ----------------- # 
    
  df <- readRDS(file.path("/scratch/users/k2480753/ARG_analysis/data/ARG_species_DATA.rds"))
    
  # create vector of ARG and microbial metab column names 
  arg_col_names = grep("_RPKM$", names(df), value = TRUE)
  metab_cols <- names(dat_stool_metabs)[grepl("^m[0-9]+", names(dat_stool_metabs))]
  stool_SCFA_columns = c("SCFA_Methylbutyrate", "SCFA_Acetate", "SCFA_Butyrate", "SCFA_Hexanoate",
                       "SCFA_Isobutyrate", "SCFA_Isovalerate", "SCFA_Propionate", "SCFA_Valerate")

  
# -------------------------------------------------------------------------- # 
  
# ----------------- # 
#   SLURM setup  ---- 
# ----------------- # 

  args <- commandArgs(trailingOnly = TRUE)
  index <- as.integer(args[1])
  
  n_jobs <- 10
  
  combos <- expand.grid(
    arg = arg_col_names,
    metabs = c(metab_cols,stool_SCFA_columns),
    stringsAsFactors = FALSE
  )
  
  total <- nrow(combos)
  per_job <- ceiling(total / n_jobs)
  
  start <- (index - 1) * per_job + 1
  end <- min(index * per_job, total)
  
  selected_combos <- combos[start:end, ]
  
# -------------------------------------------------------------------------- # 
  
# ------------------------------ # 
#   Error & warning function  ---- 
# ------------------------------ # 
  
  catch <- function(expr) {
    warn <- err <- NULL
    
    value <- withCallingHandlers(
      tryCatch(expr, error = function(e) {
        err <<- e$message
        return(NULL)
      }),
      warning = function(w) {
        warn <<- w$message
        invokeRestart("muffleWarning")
      }
    )
    
    list(value = value, warning = warn, error = err)
  }
  
  
# -------------------------------------------------------------------------- # 
    
# ---------------------------------------------- # 
#   Test ARG ~ stool metabolite associations  ---- 
# ---------------------------------------------- # 

  results <- list()
  
  for (i in 1:nrow(selected_combos)) {
      
    arg <- selected_combos$arg[i]
    metab  <- selected_combos$metab[i]
    
    model_formula = as.formula(paste0(arg, " ~ ", metab, " + age_z + sex + bmi_z + gm_batch + (1 | fid)"))
    
    
    model_result <- catch({
      
      model <- glmmTMB(
        model_formula,
        data = df,
        family = gaussian()
      )
      
      tidy_out <- broom.mixed::tidy(model, conf.int = TRUE)
      
      effect <- tidy_out %>% filter(term == metab)
      
      if (nrow(effect) == 0) {
        stop("No coefficient found")
      }
      
      data.frame(
        ARG = arg,
        metabolite_annot = metab,
        beta = effect$estimate,
        lower_CI = effect$conf.low,
        upper_CI = effect$conf.high,
        p_value = effect$p.value,
        se = effect$std.error,
        n  = nobs(model),
        formula = deparse1(formula(model))
      )
    })
    
    if (is.null(model_result$value)) {
      
      results[[length(results) + 1]] <- data.frame(
        ARG      = arg,
        metabolite_annot = metab,
        beta     = NA_real_,
        lower_CI = NA_real_,
        upper_CI = NA_real_,
        p_value        = NA_real_,
        se       = NA_real_,
        n        = NA_real_,
        formula  = NA_real_,
        error    = ifelse(is.null(model_result$error),   NA_character_, model_result$error),
        warning  = ifelse(is.null(model_result$warning), NA_character_, model_result$warning)
      )
      
    } else {
      
      results[[length(results) + 1]] <- data.frame(
        ARG        = model_result$value$ARG,
        metabolite_annot = metab,
        beta       = model_result$value$beta,
        lower_CI   = model_result$value$lower_CI,
        upper_CI   = model_result$value$upper_CI,
        p_value    = model_result$value$p_value,
        se         = model_result$value$se,
        n          = nobs(model),
        formula    = deparse1(formula(model)),
        error      = ifelse(is.null(model_result$error),   NA_character_, model_result$error),
        warning    = ifelse(is.null(model_result$warning), NA_character_, model_result$warning)
      )
    }
    }
  
  results_df <- bind_rows(results)
  
# -------------------------------------------------------------------------- # 
  
# ----------------- # 
#   Save output  ---- 
# ----------------- # 
  
  out_path <- paste0(
    "/scratch/users/k2480753/ARG_analysis/ARG_stool_metabs/results/ARG_stool_metab_results_",
    index,
    ".csv"
  )
  
  write.csv(results_df, out_path, row.names = FALSE)

# ************************************************************************ # 
#                           END OF SCRIPT                                  # 
# ************************************************************************ # 