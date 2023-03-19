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

cv_performance <- main_sens_dt %>% 
    unnest(locf_perf) %>% 
    select(method, LM:brier_death)

error1 <- cv_performance %>% 
    filter(method == "main") %>% 
    select(-(contains("brier"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    filter(method == "main") %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>% 
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Main Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30)) 

error2 <- cv_performance %>% 
    filter(method == "sens") %>% 
    select(-(contains("brier"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>%
    filter(method == "sens") %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>%  
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Sensitivity Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30))

error1 + error2 + plot_annotation(
    title = "Comparison Main vs Sensitivity Analysis Approach"
)

brier1 <- cv_performance %>% 
    filter(method == "main") %>% 
    select(-(contains("error"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    filter(method == "main") %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Main Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

brier2 <- cv_performance %>% 
    filter(method == "sens") %>% 
    select(-(contains("error"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    filter(method == "sens") %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Sensitivity Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

brier1 + brier2 + plot_annotation(
    title = "Comparison Main vs Sensitivity Analysis Approach"
)

# end ---------------------------------------------------------------------