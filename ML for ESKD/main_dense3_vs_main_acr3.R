libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "mice", "randomForestSRC", "pec", 
               "patchwork", "purrr")
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]
if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)
lib.ok <- sapply(libraries, require, character.only=TRUE)
source("utility_functions.R")

cross_val <- 5
seed <- 7
predict_horizon <- 5 
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)

longi_cov <- c("albumin", "alkphos", "bicarb", "calcium", "chloride",
                "egfr", "glucose", "haemoglobin", "phosphate", "platelet",
                "potassium", "sodium", "wcc")

longi_cov_acr <- c("albumin", "albcreat_ratio", "alkphos", "bicarb", "calcium", 
                   "chloride",  "egfr", "glucose", "haemoglobin", "phosphate", 
                   "platelet", "potassium", "sodium", "wcc")

base_cov <- c("sex", "age_init", "gn_fct")

# dense3 vs acr3 performance ----------------------------------------------

main_dt <- tibble(
    dt_path = c("main_dense3_dt.rds", "main_acr3_dt.rds"),
    base_cov = list(base_cov, base_cov),
    longi_cov = list(longi_cov, longi_cov_acr)
)

main_dt <- main_dt %>% 
    mutate(dt = map(dt_path, prepare_dt),
           folds = map(dt, create_folds, seed)) %>% 
    unnest_longer(folds)

main_dt <- main_dt %>% 
    select(folds) %>% 
    map(~.x) %>% 
    bind_rows() %>% 
    bind_cols(main_dt) %>% 
    select(id, base_cov, longi_cov, dt, splits) %>% 
    mutate(id = str_to_lower(id), 
           id_name = rep(c("dense3", "acr3"), each = 5)) %>% 
    relocate(id_name, .before = id)

main_dt <- main_dt %>% 
    mutate(train_dt = pmap(list(dt, splits), extract_train_dt),
           test_dt = pmap(list(dt, splits), extract_test_dt)) %>% 
    select(id_name, id, base_cov, longi_cov, train_dt, test_dt)

main_dt <- main_dt %>% 
    mutate(locf_perf = pmap(list(train_dt, test_dt, longi_cov, base_cov), 
                                     cross_val_performance_locf_map, landmark,
                                     predict_horizon),
           lme_perf = pmap(list(train_dt, test_dt, longi_cov, base_cov), 
                           cross_val_performance_lme_map, landmark,
                           predict_horizon),
           lmepoly_perf = pmap(list(train_dt, test_dt, longi_cov, base_cov), 
                               cross_val_performance_lmepoly_map, landmark,
                               predict_horizon))

write_rds(main_dt, "main_dense3_vs_acr3_performance.rds")


# dense3 vs acr3 LOCF visualisation --------------------------------------------

main_dt <- read_rds("main_dense3_vs_acr3_performance.rds")

locf_perf <- main_dt %>% 
    unnest_longer(locf_perf) %>% 
    select(locf_perf) %>% 
    map(~.x) %>% 
    bind_rows() 
    
dense3_acr3_locf_performance <- main_dt %>% 
    unnest_longer(locf_perf) %>% 
    select(id_name) %>% 
    bind_cols(locf_perf)

error1 <- dense3_acr3_locf_performance %>% 
    filter(id_name == "dense3") %>% 
    select(-(contains("brier"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>% 
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Dense3 Error Rate (1 - Harrell's C-index) per Landmark with LOCF 
         Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30)) 

error2 <- dense3_acr3_locf_performance %>% 
    filter(id_name == "acr3") %>% 
    select(-(contains("brier"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>% 
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "acr3 Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30)) 

error1 + error2 + plot_annotation(
    title = "Comparison Dense3 vs Acr3 Main LOCF"
)

brier1 <- dense3_acr3_locf_performance %>% 
    filter(id_name == "dense3") %>% 
    select(-(contains("error"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Dense3 Integrated Brier Score per Landmark with LOCF 
         Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

brier2 <- dense3_acr3_locf_performance %>%
    filter(id_name == "acr3") %>% 
    select(-(contains("error"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "acr3 Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

brier1 + brier2 + plot_annotation(
    title = "Comparison Dense3 vs Acr3 Main LOCF"
)

# main_dense3 exploration -------------------------------------------------

dense3_performance <- main_dt %>% 
    filter(id_name == "dense3") %>% 
    select(id_name, lme_perf:lmepoly_perf)

dense3_lme_performance <- dense3_performance %>% 
    select(lme_perf) %>% 
    unnest_longer(lme_perf) %>% 
    map(~.x) %>% 
    bind_rows()

dense3_lmepoly_performance <- dense3_performance %>% 
    select(lmepoly_perf) %>% 
    unnest_longer(lmepoly_perf) %>% 
    map(~.x) %>% 
    bind_rows()
    
error2 <- dense3_lme_performance %>% 
    select(-(contains("brier"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>% 
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Dense3 Error Rate (1 - Harrell's C-index) per Landmark with LME 
         Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30)) 

error3 <- dense3_lmepoly_performance %>% 
    select(-(contains("brier"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>% 
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Dense3 Error Rate (1 - Harrell's C-index) per Landmark with LME Poly 
         Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30)) 

error1 + error2 + error3 + plot_annotation(
    title = "Comparison Dense3 Main Exploration"
)

brier2 <- dense3_lme_performance %>%
    select(-(contains("error"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Dense3 Integrated Brier Score per Landmark with LME 
         Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

brier3 <- dense3_lmepoly_performance %>%
    select(-(contains("error"))) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Dense3 Integrated Brier Score per Landmark with LME Poly
         Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

brier1 + brier2 + brier3 + plot_annotation(
    title = "Comparison Dense3 Main Exploration"
)

# end ---------------------------------------------------------------------