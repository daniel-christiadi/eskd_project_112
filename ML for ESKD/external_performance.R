# setup -------------------------------------------------------------------

libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "randomForestSRC", "pec", 
               "purrr", "boot")
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]

if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)

lib.ok <- sapply(libraries, require, character.only=TRUE)
predict_horizon <- 5
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)

# longi_cov_full <- c("albumin", "alkphos", "bicarb", "calcium", "chloride",
#                     "egfr", "glucose", "haemoglobin", "phosphate", "platelet",
#                     "potassium", "sodium", "wcc")
# base_cov_full <- c("sex", "age_init", "gn_fct")

longi_cov_top5 <- c("albumin", "bicarb", "chloride", "egfr", "haemoglobin")
base_cov_reduced <- c("age_init")

external_data_file_name <- "experiment_test.rds" #please input file name here!
analysis_models <- c("top5")
boot_n <- 5
source("utility_functions_for_WA.R")

# performance -------------------------------------------------------------

external_dt <- tibble(
    dt_path = list(external_data_file_name),
    base_cov = list(base_cov_reduced),
    longi_cov = list(longi_cov_top5),
    analysis_models = analysis_models
)

external_dt <- external_dt %>% 
    mutate(dt = map(dt_path, prepare_dt))

external_dt <- external_dt %>% 
    mutate(landmark = list(landmark)) %>% 
    unnest(landmark) %>% 
    select(dt, base_cov, longi_cov, landmark, analysis_models)

external_dt <- external_dt %>% 
    mutate(lmx_dt = pmap(list(dt, base_cov, longi_cov, landmark), 
                         extract_landmark_dt, predict_horizon))

external_dt <- external_dt %>% 
    mutate(error = pmap(list(lmx_dt, landmark, base_cov, longi_cov, analysis_models), 
                        performance_landmark_dt))

external_dt <- external_dt %>% 
    mutate(brier = pmap(list(lmx_dt, landmark, base_cov, longi_cov, analysis_models),
                                 brier_landmark_dt))

external_dt <- external_dt %>% 
    select(landmark, analysis_models, error, brier)

write_rds(external_dt, "external_performance.rds")

# end ---------------------------------------------------------------------