# setup -------------------------------------------------------------------

libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "mice", "randomForestSRC", "pec", 
               "gtsummary", "patchwork", "purrr")
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]
if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)
lib.ok <- sapply(libraries, require, character.only=TRUE)

cross_val <- 5
seed <- 7
predict_horizon <- 5
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)

longi_cov <- c("albumin", "alkphos", "bicarb", "calcium", "chloride",
                "egfr", "glucose", "haemoglobin", "phosphate", "platelet",
                "potassium", "sodium", "wcc")

base_cov <- c("sex", "age_init", "gn_fct")
approach <- list("main", "sens")
source("utility_functions.R")

# performance cross-validation --------------------------------------------

main_sens_dt <- tibble(
    method = c("main", "sens"),
    file_name = rep("dense3_dt.rds", times=2),
    base_cov = list(base_cov, base_cov),
    longi_cov = list(longi_cov, longi_cov)
)

main_sens_dt <- main_sens_dt %>% 
    unite(method, file_name, sep = "_", col = "dt_path", remove = FALSE)

main_sens_dt <- main_sens_dt %>% 
    mutate(dt = map(dt_path, prepare_dt))

main_sens_dt <- main_sens_dt %>% 
    mutate(folds = map(dt, create_folds, seed)) %>% 
    unnest_longer(folds)

main_sens_dt <- main_sens_dt %>% 
    select(folds) %>% 
    map(~.x) %>% 
    bind_rows() %>% 
    bind_cols(main_sens_dt) %>% 
    mutate(id = str_to_lower(id)) %>% 
    select(method, id, base_cov, longi_cov, dt, splits)

main_sens_dt <- main_sens_dt %>% 
    mutate(train_dt = pmap(list(dt, splits), extract_train_dt),
           test_dt = pmap(list(dt, splits), extract_test_dt)) %>% 
    select(method, id, base_cov, longi_cov, train_dt, test_dt)

main_sens_dt <- main_sens_dt %>% 
    mutate(locf_perf = pmap(list(train_dt, test_dt, longi_cov, base_cov), 
                                     cross_val_performance_locf_map, landmark,
                                     predict_horizon))

write_rds(main_sens_dt, "main_vs_sens_performance.rds")

# visualisation -----------------------------------------------------------

main_sens_dt <- read_rds("main_vs_sens_performance.rds")

cv_performance <- main_sens_dt %>% 
    unnest(locf_perf) %>% 
    select(method, LM:brier_death)

cv_performance <- cv_performance %>% 
    summarise(across(error_eskd:brier_death, \(x) confidence_interval (x, 0.95)),
              .by = c("method", "LM")) %>% 
    unnest_wider(error_eskd:brier_death, names_sep = "_")

cv_performance <- cv_performance %>% 
    pivot_longer(cols = error_eskd_mean:brier_death_upper,
                 names_to = "performance", 
                 values_to = "values") %>% 
    separate_wider_delim(cols = performance, delim = "_",
                         names = c("metrics", "event", "stats")) %>% 
    pivot_wider(names_from = stats, 
                values_from = values)

error <- cv_performance %>% 
    filter(metrics == "error")

error %>% 
    group_by(event) %>% gt()

error1 <- create_error_graph(error, "main",
                             "Main Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years")

error2 <- create_error_graph(error, "sens",
                             "Sensitivity Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years")

error1 + error2 + 
    plot_annotation(title = "Comparison Main vs Sensitivity Analysis Approach", 
                    theme = theme(plot.title = element_text(hjust = 0.5))) 

brier <- cv_performance %>% 
    filter(metrics == "brier")

brier %>% 
    group_by(event) %>% gt()

brier1 <- create_brier_graph(brier, "main",
                             "Main Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years")

brier2 <- create_brier_graph(brier, "sens",
                             "Sensitivity Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years")

brier1 + brier2 + 
    plot_annotation(title = "Comparison Main vs Sensitivity Analysis Approach", 
                    theme = theme(plot.title = element_text(hjust = 0.5))) 

# end ---------------------------------------------------------------------