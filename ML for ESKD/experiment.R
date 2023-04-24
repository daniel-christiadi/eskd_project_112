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


# experiment_starts -------------------------------------------------------

dt1 <- experiment_dt$train_dt[[1]]
dt1$type <- "train"
dt2 <- experiment_dt$test_dt[[1]]
dt2$type <- "test"

combined_dt <- rbind(dt1, dt2)

combined_dt %>% 
    filter(relyear == 0) %>% 
    select(sex:gn_cat, time_compete:status_compete, !!longi_cov_full, type) %>%
    tbl_summary(by = type) %>% add_p()

dt <- experiment_dt$train_dt[[1]]
test_dt <- experiment_dt$test_dt[[1]]

long_dt <- dt %>%
    select(id, relyear, !!longi_cov_full)

impute_dt <- long_dt %>%
    filter(relyear == 0)

temp_dt <- mice(impute_dt, m = 1, maxit = 50, method = "pmm", seed = seed)
impute_dt <- complete(temp_dt)

dt <- impute_dt %>%
    bind_rows(long_dt) %>%
    arrange(id) %>%
    distinct(id, relyear, .keep_all = TRUE) %>%
    group_by(id) %>%
    fill(!!longi_cov_full) %>%
    ungroup() %>%
    left_join(select(dt, id:status_compete_num),
              by = c("id", "relyear")) %>%
    arrange(id, relyear) %>%
    as.data.table()

super_dt <- NULL

for (i in seq_along(landmark)) {
    temp <- cutLM(data = dt, outcome = list(time = "time_compete",
                                            status = "status_compete_num"),
                  LM = landmark[i], horizon = landmark[i] + predict_horizon,
                  covs = list(fixed = base_cov_full, varying = longi_cov_full),
                  format = "long", id = "id", rtime = "relyear",
                  right = FALSE)
    super_dt <- rbind(super_dt, temp)
}

test_dt <- test_dt %>% 
    group_by(id) %>% 
    fill(!!longi_cov_full) %>% 
    ungroup() %>% 
    as.data.table()

test_super_dt <- NULL

for (i in seq_along(landmark)) {
    temp <- cutLM(data = test_dt, outcome = list(time = "time_compete",
                                                 status = "status_compete_num"),
                  LM = landmark[i], horizon = landmark[i] + predict_horizon,
                  covs = list(fixed = base_cov_full, varying = longi_cov_full),
                  format = "long", id = "id", rtime = "relyear", 
                  right = FALSE) 
    test_super_dt <- rbind(test_super_dt, temp)
}


lmx <- 1
dt_lm <- super_dt %>%
    filter(LM == lmx)

test_lm <- test_super_dt %>% 
    filter(LM == lmx) %>% 
    drop_na()

test_lm_nodeath <- test_lm %>% 
    filter(status_compete_num != 2)

model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                  str_c(c(base_cov_full, longi_cov_full), collapse = " + "),
                                  sep = "~"))

tune_rsf <- tune.rfsrc(model_formula, data = dt_lm,
                       ntreeTry = 500, nodesizeTry = c(1:9, seq(10, 100, 5)))

tune_matrix <- rbind(tune_matrix, tune_rsf$optimal)

rsf_model <- rfsrc(formula = model_formula, data = dt_lm, ntree = 1000,
                   splitrule = "logrankCR", importance = FALSE, statistics = FALSE,
                   mtry = tune_rsf$optimal[[2]], nodesize = tune_rsf$optimal[[1]], 
                   save.memory = TRUE)

rsf_prediction <- predict.rfsrc(object = rsf_model,
                                newdata = test_lm, na.action = "na.omit")

brier_eskd <- pec(rsf_model, formula = model_formula, data = test_lm_nodeath,
                  cause = 1)

brier_death <- pec(rsf_model, formula = model_formula, data = test_lm,
                   cause = 2)
    

output_tibble <- tibble(
    landmark = landmark,
    nodesize = tune_matrix[,1],
    mtry = tune_matrix[,2],
)    
tune_tibble_name <- str_c(model_name, "locf_tune.rds", sep = "_")
write_rds(output_tibble, file = tune_tibble_name)



results_previous <- read_rds("main_dense3_vs_acr3_performance.rds")

error_eskd <- c(10.1596, 13.82450, 13.578645, 15.20753, 13.36848)

confidence_interval <- function(vector, interval) {
    vec_sd <- sd(vector)
    n <- length(vector)
    vec_mean <- mean(vector)
    error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
    result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
    return(result)
}

confidence_interval(error_eskd, 0.95)


vimp <- read_rds("dense3_locf_vimp.rds")



subsample.rfsrc(rsf_model, B = 10)





# original_script ---------------------------------------------------------

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





    
