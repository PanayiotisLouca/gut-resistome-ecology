## Author: 
  # Panayiotis Louca 

## Purpose of script: 
  #  Perform mixed effect modelling to test associations between ARGs and host metabolic health parameters 

## Clear environment 
  rm(list = ls()) 

## Set seed 
  set.seed(1234)
      
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

  df <- readRDS(file.path("/Users/panayiotislouca/Documents/KCL_files/ARG_characterisation/Analysis/Data/ARG_species_DATA.rds"))
  
  # create vector of ARG col names 
  arg_col_names = grep("_RPKM$", names(df), value = TRUE)
  
# -------------------------------------------------------------------------- # 
  
# ---------------------------------------------- # 
#   Import high interest ARG classifications  ---- 
# ---------------------------------------------- # 
  
  high_interest_classifications <- read.csv(file.path("/Users/panayiotislouca/Documents/KCL_files/ARG_characterisation/Analysis/High_interest_ARG_carriage/High_interest_ARG_classifications.csv"))
  
# -------------------------------------------------------------------------- # 
  
# ----------- # 
#   Merge  ---- 
# ----------- # 
  
  df <- merge(df, high_interest_classifications,
                       by = "iid",
                       all.x = TRUE)
  
# -------------------------------------------------------------------------- # 

# --------------------------------- # 
#   Classify high-interest ARGs  ---- 
# --------------------------------- # 
  
  df <- df %>%
    mutate(
      high_interest_binary = factor(
        if_else(n_high_interest == 0, "None", "Any"),
        levels = c("None", "Any")
      ))
  
  
# -------------------------------------------------------------------------- # 
  
# -------------------------------------------- # 
#   Model TyG ~ High interest ARG carriage  ---- 
# -------------------------------------------- # 
  
  model_carriage <- glmmTMB(
    TyG_index ~ high_interest_binary +
      age_z + sex + bmi_z + gm_batch + taxa_shannon +
      (1 | fid),
    data = df,
    family = gaussian())
  
  model_carriage_output <- broom.mixed::tidy(
    model_carriage,
    effects = "fixed",
    conf.int = TRUE) %>%
    filter(term == "high_interest_binaryAny")
  
  model_carriage_output
  
# -------------------------------------------------------------------------- # 
  
# ------------------------------ # 
#   Model TyG ~ ARG richness  ---- 
# ------------------------------ # 
  
  model_richness <- glmmTMB(
    TyG_index ~ amr_gene_richness +
      age_z + sex + bmi_z + gm_batch + taxa_shannon +
      (1 | fid),
    data = df,
    family = gaussian())
  
  model_richness_output <- broom.mixed::tidy(
    model_richness,
    effects = "fixed",
    conf.int = TRUE) %>%
    filter(term == "amr_gene_richness")
  
  model_richness_output
  
# -------------------------------------------------------------------------- # 

# ----------------------------- # 
#   Model TyG ~ ARG Shannon  ---- 
# ----------------------------- # 
  
  model_shannon <- glmmTMB(
    TyG_index ~ amr_gene_shannon +
      age_z + sex + bmi_z + gm_batch + taxa_shannon +
      (1 | fid),
    data = df,
    family = gaussian())
  
  model_shannon_output <- broom.mixed::tidy(
    model_shannon,
    effects = "fixed",
    conf.int = TRUE) %>%
    filter(term == "amr_gene_shannon")
  
  model_shannon_output
  
  
# -------------------------------------------------------------------------- # 

# -------------------------------------- # 
#   Model TyG ~ n high interest ARGs  ---- 
# -------------------------------------- # 
  
  model_nARG <- glmmTMB(
    TyG_index ~ n_high_interest +
      age_z + sex + bmi_z + gm_batch + taxa_shannon +
      (1 | fid),
    data = df,
    family = gaussian()
  )
  
  model_nARG_output <- broom.mixed::tidy(
    model_nARG,
    effects = "fixed",
    conf.int = TRUE) %>%
    filter(term == "n_high_interest")
  
  model_nARG_output

# -------------------------------------------------------------------------- # 

# --------------------- # 
#   Combine results  ---- 
# --------------------- # 
  
  results_df <- bind_rows(model_carriage_output,
                          model_richness_output,
                          model_shannon_output)
  
  # adjust for multiple testing 
  results_df <- results_df %>%
    mutate(fdr = p.adjust(p.value, method = "BH"))
  
  results_df
  
# -------------------------------------------------------------------------- # 

# ------------------------------------------------------------------------------------------- # 
#   Model clinically prioritised ARG carriage as the outcome against metabolic parameters  ---- 
# ------------------------------------------------------------------------------------------- # 
  
  model_high_interest_TyG <- glmmTMB(
    high_interest_binary ~ TyG_index +
      age_z + sex + bmi_z + gm_batch + taxa_shannon +
      (1 | fid),
    data = df,
    family = binomial(link = "logit"))
  
  model_high_interest_TyG_output <-
    broom.mixed::tidy(
      model_high_interest_TyG,
      effects = "fixed",
      conf.int = TRUE,
      exponentiate = TRUE
    ) %>%
    filter(term == "TyG_index")
  
  model_high_interest_TyG_output
  
  
# -------------------------------------------------------------------------- # 

# -------------------------- # 
#   Model Triglycerides   ---- 
# -------------------------- # 
  
  model_high_interest_TG <- glmmTMB(
    high_interest_binary ~ serum_TG +
      age_z + sex + bmi_z + gm_batch + taxa_shannon +
      (1 | fid),
    data = df,
    family = binomial(link = "logit"))
  
  model_high_interest_TG_output <-
    broom.mixed::tidy(
      model_high_interest_TG,
      effects = "fixed",
      conf.int = TRUE,
      exponentiate = TRUE
    ) %>%
    filter(term == "serum_TG")
  
  model_high_interest_TG_output
  
# -------------------------------------------------------------------------- # 

# ------------------- # 
#   Model Glucose  ---- 
# ------------------- # 
  
  model_high_interest_Gluc <- glmmTMB(
    high_interest_binary ~ serum_Gluc +
      age_z + sex + bmi_z + gm_batch + taxa_shannon +
      (1 | fid),
    data = df,
    family = binomial(link = "logit"))
  
  model_high_interest_Gluc_output <-
    broom.mixed::tidy(
      model_high_interest_Gluc,
      effects = "fixed",
      conf.int = TRUE,
      exponentiate = TRUE
    ) %>%
    filter(term == "serum_Gluc")
  
  model_high_interest_Gluc_output
  
# -------------------------------------------------------------------------- # 

# --------------------- # 
#   Combine results  ---- 
# --------------------- # 
  
  results_df <- bind_rows(
    model_high_interest_TyG_output,
    model_high_interest_TG_output,
    model_high_interest_Gluc_output)
  
  results_df <- results_df %>%
    mutate(fdr = p.adjust(p.value, "BH"))
  
  results_df
  
  
  
  results_TyG <- bind_rows(
    model_carriage_output  %>% mutate(model = "high_interest_binary"),
    model_richness_output  %>% mutate(model = "amr_gene_richness"),
    model_shannon_output   %>% mutate(model = "amr_gene_shannon"),
    model_nARG_output      %>% mutate(model = "n_high_interest")   # optional but useful
  ) %>%
    mutate(
      fdr = p.adjust(p.value, method = "BH"),
      outcome = "TyG_index",
      direction = "TyG ~ resistome"
    ) %>%
    select(model, term, estimate, std.error, conf.low, conf.high, p.value, fdr, everything())
  
# -------------------------------------------------------------------------- # 

# ------------------------ # 
#   Write results file  ---- 
# ------------------------ # 
  
  write.csv(results_TyG,
            file.path("/Users/panayiotislouca/Documents/KCL_files/ARG_characterisation/Analysis/ARG_host_phenotype_associations/ARG_metabolic_health_mixed_model_RESULTS.csv"),
            row.names = FALSE)
  
  
# ************************************************************************ # 
#                           END OF SCRIPT                                  # 
# ************************************************************************ # 