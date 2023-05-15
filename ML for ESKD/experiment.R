# setup -------------------------------------------------------------------
libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "randomForestSRC", "pec", 
               "gtsummary","purrr", "patchwork", "gt","boot")
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

model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                  str_c(c(external_dt$base_cov[[8]], external_dt$longi_cov[[8]]), collapse = " + "), 
                                  sep = "~"))

external_dt

landmark_dt <- external_dt$lmx_dt[[8]] %>% 
    drop_na()

landmark_dt_nodeath <- landmark_dt %>% 
    filter(status_compete_num != 2)

dummy_dt <- landmark_dt_nodeath %>% 
    filter(status_compete_num == 1) %>% head(n=1)

dummy_dt$status_compete_num <- 2

landmark_dt_fakedeath <- landmark_dt_nodeath %>% 
    add_row(dummy_dt)

model_name <- str_c(external_dt$analysis_models[[8]], 1, "locf_models.gz", sep = "_")
model_path <- file.path(current_dir, model_name)
rsf_model <- readRDS(model_name)

brier_eskd <- pec(rsf_model, formula = model_formula, landmark_dt_fakedeath, cause = 1)

error_eskd <- predict.rfsrc(rsf_model, data = landmark_dt_fakedeath)

brier_eskd_fakedeath <- pec(rsf_model, data = landmark_dt_fakedeath, cause = 1)


error <- boot::boot(data = landmark_dt, statistic = error_bootstrap,
                    R = boot_n, rsf_model = rsf_model)

error_name <- str_c(analysis_models, landmark, "error.rds", sep = "_")
error_path <- file.path(current_dir, error_name)
write_rds(error, error_path)

performance <- data.table(LM = landmark,
                          N = nrow(landmark_dt))
performance



# test cif ----------------------------------------------------------------

dt <- read_rds("experiment_test.rds")

eskd_dt <- dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    filter(status_compete == "eskd")

death_dt <- dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    filter(status_compete == "death")

censored_dt <- dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    filter(status_compete == "no")

random_test_id <- censored_dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    slice_sample(n = 1,) %>% 
    pull(id)

random_test_id <- 13651

random_test_dt <- tibble(dt = list(filter(dt, id == random_test_id)),
                         base_cov = list(base_cov_reduced),
                         longi_cov = list(longi_cov_top5),
                         analysis_name = c("experiment_top5"))

random_test_dt <- random_test_dt %>% 
    mutate(dt = map(dt, prepare_dt_for_dynamic_plot)) %>% 
    mutate(landmark = map(dt, extract_landmark_for_dynamic_plot, landmark)) %>% 
    unnest_longer(landmark) %>% 
    mutate(lmx_dt = pmap(list(dt, base_cov, longi_cov, landmark), 
                         extract_landmark_dt, predict_horizon)) %>% 
    mutate(pred_cif = pmap(list(lmx_dt, landmark, base_cov, longi_cov, 
                                analysis_name), predicted_cif, predict_horizon))

summarise_pts_info(random_test_dt$dt, 
                   base_cov = random_test_dt$base_cov,
                   longi_cov = random_test_dt$longi_cov)

create_dynamic_plot(dt = random_test_dt$dt,
                    base_cov = random_test_dt$base_cov,
                    longi_cov = random_test_dt$longi_cov,
                    landmark = random_test_dt$landmark,
                    pred_cif = random_test_dt$pred_cif)



# just_in_case ------------------------------------------------------------

error_bootstrap <- function(rsf_model, data, index){
    d <- data[index, ]
    rsf_prediction <- predict.rfsrc(object = rsf_model, 
                                    newdata = d, na.action = "na.omit")
    error_temp <- unname(apply(rsf_prediction$err.rate, 2, mean, na.rm = TRUE))
    return(error_temp)
}

eskd_brier_bootstrap <- function(rsf_model, model_formula, data, index){
    d <- data[index, ]
    brier_eskd <- pec(rsf_model, formula = model_formula, data = d,
                      cause = 1)
    return(crps(brier_eskd))
}

death_brier_bootstrap <- function(rsf_model, model_formula, data, index){
    d <- data[index, ]
    brier_death <- pec(rsf_model, formula = model_formula, data = d,
                       cause = 2)
    return(crps(brier_death))
}

error <- boot::boot(data = test_lm, statistic = error_bootstrap,
                    R = 10, rsf_model = rsf_model)

boot::boot.ci(error, type = "norm", index = 2)

results_eskd_brier <- boot::boot(data = test_lm, statistic = eskd_brier_bootstrap,
                                 R = 10, rsf_model = rsf_model, 
                                 model_formula = model_formula)

results_death_brier <- boot::boot(data = test_lm, statistic = death_brier_bootstrap,
                                  R = 10, rsf_model = rsf_model, 
                                  model_formula = model_formula)



features <- read_rds("main_features.rds")
outcomes <- read_rds("main_outcome.rds")

outcomes %>% arrange(age_init)    

dense3 <- read_rds("main_dense3_dt.rds")

dense3 %>% distinct(id, .keep_all = TRUE) %>% arrange(age_init)
