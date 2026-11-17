## Author: 
  # Panayiotis Louca 

## Purpose of script: 
  #  Perform mixed effect modelling to test associations between ARGs and bacterial species 

## Notes: 
  #  

## Clear environment 
  rm(list = ls()) 

## Set seed 
  set.seed(1234)
      
## load up packages: 

  ### core 
  library(tidyverse)
    
  # mixed effect modelling 
  library(glmmTMB)
  library(lmerTest)
  library(broom.mixed)
  
  # parallelisation 
  library(parallel)
  library(future.apply)
  
  # Increase to 2 GiB 
  options(future.globals.maxSize = 2000 * 1024^2)  

# -------------------------------------------------------------------------- # 
  
# ----------------- # 
#   Import data  ---- 
# ----------------- # 
  
  df <- readRDS("/scratch/users/k2480753/ARG_analysis/data/ARG_species_DATA.rds")
  
  # Vectorise column names 
  arg_columns <- grep("_RPKM$", names(df), value = TRUE)
  species_columns <- grep("species_", names(df), value = TRUE)

# -------------------------------------------------------------------------- # 
  
# ----------------- # 
#   SLURM setup  ---- 
# ----------------- # 
    
  args <- commandArgs(trailingOnly = TRUE)
  index <- as.integer(args[1])
  
  n_jobs <- 25
  
  combos <- expand.grid(
    arg = arg_columns,
    species = species_columns,
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
  
# ---------------- # 
#   Model loop  ---- 
# ---------------- # 
  
  results <- list()
  
  for (i in 1:nrow(selected_combos)) {
    
    arg <- selected_combos$arg[i]
    sp  <- selected_combos$species[i]
    
    model_result <- catch({
      
      model <- glmmTMB(
        as.formula(paste(arg, "~", sp, "+ age_z + sex + bmi_z + gm_batch + (1 | fid))")),
        data = df,
        family = gaussian()
      )
      
      tidy_out <- broom.mixed::tidy(model, conf.int = TRUE)
      
      effect <- tidy_out %>% filter(term == sp)
      
      if (nrow(effect) == 0) {
        stop("No coefficient found")
      }
      
      data.frame(
        ARG = arg,
        species = sp,
        beta = effect$estimate,
        se = effect$std.error,
        p = effect$p.value
      )
    })
    
    if (is.null(model_result$value)) {
      
      results[[length(results) + 1]] <- data.frame(
        ARG     = arg,
        species = sp,
        beta    = NA_real_,
        se      = NA_real_,
        p       = NA_real_,
        error   = ifelse(is.null(model_result$error),   NA_character_, model_result$error),
        warning = ifelse(is.null(model_result$warning), NA_character_, model_result$warning)
      )
      
    } else {
      
      results[[length(results) + 1]] <- data.frame(
        ARG        = model_result$value$ARG,
        species    = model_result$value$species,
        beta       = model_result$value$beta,
        se         = model_result$value$se,
        p          = model_result$value$p,
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
    "/scratch/users/k2480753/ARG_analysis/ARG_species_associations/results/ARG_species_results_",
    index,
    ".csv"
  )
  
  write.csv(results_df, out_path, row.names = FALSE)
  