# setup -------------------------------------------------------------------

libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "randomForestSRC", "pec", 
               "gtsummary","purrr", "patchwork", "gt")
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]
if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)
lib.ok <- sapply(libraries, require, character.only=TRUE)
cross_val <- 5
seed <- 7
predict_horizon <- 5
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)

longi_cov_full <- c("albumin", "alkphos", "bicarb", "calcium", "chloride",
                    "egfr", "glucose", "haemoglobin", "phosphate", "platelet",
                    "potassium", "sodium", "wcc")

longi_cov_top5 <- c("albumin", "bicarb", "chloride", "egfr", "haemoglobin")
base_cov_full <- c("sex", "age_init", "gn_fct")
base_cov_reduced <- c("age_init")
experiment_fold <- "fold1"

source("utility_functions.R")

# model -------------------------------------------------------------------

experiment_dt <- tibble(
    dt_path = "main_dense3_dt.rds",
    base_cov = list(base_cov_full, base_cov_reduced),
    longi_cov = list(longi_cov_full, longi_cov_top5)
)

experiment_dt <- experiment_dt %>% 
    mutate(dt = map(dt_path, prepare_dt),
           folds = map(dt, create_folds, seed)) %>% 
    unnest_longer(folds)

experiment_dt <- experiment_dt %>% 
    select(folds) %>% 
    map(~ .x) %>% 
    bind_rows() %>% 
    bind_cols(experiment_dt) %>% 
    mutate(id_name = rep(c("experiment_full", "experiment_top5"), each = 5),
           id = str_to_lower(id)) %>% 
    select(id_name, id, base_cov, longi_cov, dt, splits) %>% 
    filter(id == experiment_fold) %>% 
    select(-id) %>% 
    mutate(train_dt = pmap(list(dt, splits), extract_train_dt),
           test_dt = pmap(list(dt, splits), extract_test_dt)) %>% 
    select(id_name:longi_cov, train_dt, test_dt)

pwalk(list(experiment_dt$train_dt, experiment_dt$longi_cov, 
                  experiment_dt$base_cov, experiment_dt$id_name),
                 final_models, landmark, predict_horizon)

# write_rds(experiment_dt$test_dt[[1]], "experiment_test.rds")
external_data_file_name <- "experiment_test.rds"
analysis_models <- c("experiment_full", "experiment_top5")

external_dt <- tibble(
    dt_path = rep(external_data_file_name, times = 2),
    base_cov = list(base_cov_full, base_cov_reduced),
    longi_cov = list(longi_cov_full, longi_cov_top5),
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
    mutate(performance = pmap(list(lmx_dt, landmark, base_cov, longi_cov, analysis_models), 
                              performance_landmark_dt))

# test cif ----------------------------------------------------------------

dt <- read_rds("experiment_test.rds")

random_test_id <- dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    slice_sample(n = 1,) %>% 
    pull(id)

random_test_dt <- tibble(dt = list(filter(dt, id == random_test_id)),
                         base_cov = list(base_cov_reduced),
                         longi_cov = list(longi_cov_top5),
                         analysis_name = c("experiment_top5"))

random_test_dt <- random_test_dt %>% 
    mutate(dt = map(dt, prepare_dt_for_dynamic_plot)) %>% 
    mutate(landmark = map(dt, extract_landmark, landmark)) %>% 
    unnest_longer(landmark) %>% 
    mutate(lmx_dt = pmap(list(dt, base_cov, longi_cov, landmark), 
                         extract_landmark_dt, predict_horizon)) %>% 
    mutate(pred_cif = pmap(list(lmx_dt, landmark, base_cov, longi_cov, 
                                analysis_name), predicted_cif, predict_horizon))

random_test_dt %>% 
    select(pred_cif) %>% 
    unnest_longer(pred_cif) %>% 
    map(~ .x) %>% 
    bind_rows() %>% 
    group_by(pred_relyear) %>% 
    summarise(mean_eskd_cif = mean(eskd_cif),
              mean_death_cif = mean(death_cif)) %>% View()

summarise_pts_info(test_dt$dt, 
                   base_cov = test_dt$base_cov,
                   longi_cov = test_dt$longi_cov)

# dynamic plot ------------------------------------------------------------







    
