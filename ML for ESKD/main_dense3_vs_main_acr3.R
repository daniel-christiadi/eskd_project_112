# setup -------------------------------------------------------------------

libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "mice", "randomForestSRC", "pec", 
               "patchwork", "purrr", "gt")
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

longi_cov_acr <- c("albumin", "albcreat_ratio", "alkphos", "bicarb", "calcium", 
                   "chloride",  "egfr", "glucose", "haemoglobin", "phosphate", 
                   "platelet", "potassium", "sodium", "wcc")

base_cov <- c("sex", "age_init", "gn_fct")
source("utility_functions.R")

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

main_dt_performance <- main_dt %>% 
    mutate(locf_perf = pmap(list(train_dt, test_dt, longi_cov, base_cov), 
                                     cross_val_performance_locf_map, landmark,
                                     predict_horizon),
           lme_perf = pmap(list(train_dt, test_dt, longi_cov, base_cov), 
                           cross_val_performance_lme_map, landmark,
                           predict_horizon),
           lmepoly_perf = pmap(list(train_dt, test_dt, longi_cov, base_cov), 
                               cross_val_performance_lmepoly_map, landmark,
                               predict_horizon))

write_rds(main_dt_performance, "main_dense3_vs_acr3_performance.rds")

# dense3 vs acr3 LOCF visualisation --------------------------------------------

main_dt_performance <- read_rds("main_dense3_vs_acr3_performance.rds")

locf_perf <- main_dt_performance %>% 
    unnest_longer(locf_perf) %>% 
    select(locf_perf) %>% 
    map(~.x) %>% 
    bind_rows() 
    
dense3_acr3_locf_performance <- main_dt_performance %>% 
    unnest_longer(locf_perf) %>% 
    select(id_name) %>% 
    bind_cols(locf_perf)

dense3_acr3_locf_performance <- dense3_acr3_locf_performance %>% 
    summarise(across(error_eskd:brier_death, \(x) confidence_interval (x, 0.95)),
              .by = c("id_name", "LM")) %>% 
    unnest_wider(error_eskd:brier_death, names_sep = "_")

dense3_acr3_locf_performance <- dense3_acr3_locf_performance %>% 
    pivot_longer(cols = error_eskd_mean:brier_death_upper,
                 names_to = "performance", 
                 values_to = "values") %>% 
    separate_wider_delim(cols = performance, delim = "_",
                         names = c("metrics", "event", "stats")) %>% 
    pivot_wider(names_from = stats, 
                values_from = values)

dense3_acr3_locf_performance <- dense3_acr3_locf_performance %>% 
    rename("method" = "id_name")

dense3_acr3_error <- dense3_acr3_locf_performance %>% 
    filter(metrics == "error")

dense3_acr3_error %>% 
    group_by(event) %>% gt()

error1 <- create_error_graph(dense3_acr3_error, method = "dense3",
                             "Dense3 Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years")

error2 <- create_error_graph(dense3_acr3_error, method = "acr3",
                             "acr3 Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years")

error1 + error2 + plot_annotation(
    title = "Comparison Dense3 vs Acr3 Main LOCF",
    theme = theme(plot.title = element_text(hjust = 0.5))
)

dense3_acr3_brier <- dense3_acr3_locf_performance %>% 
    filter(metrics == "brier")

dense3_acr3_brier %>% 
    group_by(event) %>% gt()

brier1 <- create_brier_graph(dense3_acr3_brier,
                             method = "dense3",
                             "Dense3 Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years")

brier2 <- create_brier_graph(dense3_acr3_brier,
                             method = "acr3",
                             "acr3 Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years")

brier1 + brier2 + plot_annotation(
    title = "Comparison Dense3 vs Acr3 Main LOCF",
    theme = theme(plot.title = element_text(hjust = 0.5)))

# main_dense3 exploration -------------------------------------------------

dense3_performance <- main_dt_performance %>% 
    filter(id_name == "dense3") %>% 
    select(id_name, lme_perf:lmepoly_perf)

dense3_lme_performance <- dense3_performance %>% 
    select(lme_perf) %>% 
    unnest_longer(lme_perf) %>% 
    map(~.x) %>% 
    bind_rows()

dense3_lme_performance <- dense3_lme_performance %>% 
    summarise(across(error_eskd:brier_death, \(x) confidence_interval (x, 0.95)),
              .by = c("LM")) %>% 
    unnest_wider(error_eskd:brier_death, names_sep = "_")

dense3_lme_performance <- dense3_lme_performance %>% 
    pivot_longer(cols = error_eskd_mean:brier_death_upper,
                 names_to = "performance", 
                 values_to = "values") %>% 
    separate_wider_delim(cols = performance, delim = "_",
                         names = c("metrics", "event", "stats")) %>% 
    pivot_wider(names_from = stats, 
                values_from = values)

dense3_lme_performance <- dense3_lme_performance %>% 
    mutate(method = "lme")

dense3_lmepoly_performance <- dense3_performance %>% 
    select(lmepoly_perf) %>% 
    unnest_longer(lmepoly_perf) %>% 
    map(~.x) %>% 
    bind_rows()

dense3_lmepoly_performance <- dense3_lmepoly_performance %>% 
    summarise(across(error_eskd:brier_death, \(x) confidence_interval (x, 0.95)),
              .by = c("LM")) %>% 
    unnest_wider(error_eskd:brier_death, names_sep = "_")

dense3_lmepoly_performance <- dense3_lmepoly_performance %>% 
    pivot_longer(cols = error_eskd_mean:brier_death_upper,
                 names_to = "performance", 
                 values_to = "values") %>% 
    separate_wider_delim(cols = performance, delim = "_",
                         names = c("metrics", "event", "stats")) %>% 
    pivot_wider(names_from = stats, 
                values_from = values)

dense3_lmepoly_performance <- dense3_lmepoly_performance %>% 
    mutate(method = "lme_poly")

dense3_lmes_performance <- rbind(dense3_lme_performance, dense3_lmepoly_performance)
    
dense3_lmes_error <- dense3_lmes_performance %>% 
    filter(metrics == "error")

dense3_lmes_error %>% 
    group_by(event) %>% gt()

error1 <- create_error_graph(dense3_acr3_error, method = "dense3",
                             "Dense3 Error Rate (1 - Harrell's C-index) per Landmark 
                             LOCF Prediction Horizon 5 years")

error2 <- create_error_graph(dense3_lmes_error,
                             method = "lme",
                             "Dense3 Error Rate (1 - Harrell's C-index) per Landmark 
                             LME Prediction Horizon 5 years")

error3 <- create_error_graph(dense3_lmes_error,
                             method = "lme_poly",
                             "Dense3 Error Rate (1 - Harrell's C-index) per Landmark 
                             LME Poly Prediction Horizon 5 years")

error1 + error2 + error3 + plot_annotation(
    title = "Comparison Dense3 Main Exploration",
    theme = theme(plot.title = element_text(hjust = 0.5))
)

dense3_lmes_brier <- dense3_lmes_performance %>% 
    filter(metrics == "brier")

dense3_lmes_brier %>% 
    group_by(event) %>% gt()

brier1 <- create_brier_graph(dense3_acr3_brier,
                             method = "dense3",
                             "Dense3 Integrated Brier Score per Landmark LOCF Prediction Horizon 5 years")

brier2 <- create_brier_graph(dense3_lmes_brier,
                             method = "lme",
                             "Dense3 Integrated Brier Score per Landmark LME Prediction Horizon 5 years")

brier3 <- create_brier_graph(dense3_lmes_brier,
                             method = "lme_poly",
                             "Dense3 Integrated Brier Score per Landmark LME Poly Prediction Horizon 5 years")

brier1 + brier2 + brier3 + plot_annotation(
    title = "Comparison Dense3 Main Exploration",
    theme = theme(plot.title = element_text(hjust = 0.5))
)

# dense3_vimp -------------------------------------------------------------

dense3_vimp <- main_dt %>% 
    filter(id_name == "dense3") %>% 
    mutate(vimp = pmap(list(train_dt, test_dt, longi_cov, base_cov),
                       cross_val_vimp_locf_map, landmark, predict_horizon))

write_rds(dense3_vimp, "dense3_locf_vimp_new.rds")

eskd_vimp <- dense3_vimp %>% 
    select(vimp) %>% 
    unnest(vimp) %>% 
    summarise(across(sex:wcc, \(x) confidence_interval (x, 0.95)),
              .by = c("event", "LM")) %>% 
    unnest_wider(sex:wcc, names_sep = "_") %>%
    filter(event == "eskd")


eskd_vimp <- eskd_vimp %>% 
    pivot_longer(cols = sex_mean:wcc_upper, names_to = "test", 
                 values_to = "values") %>% 
    separate_wider_delim(cols = test, delim = "_", 
                         names = c("test", "stats", "stats1"),
                         too_few = "align_start") %>% 
    mutate(stats1 = if_else(is.na(stats1), stats, stats1)) %>% 
    mutate(new_test = case_match(test,
                                 "age" ~ "age_init", .default = test)) %>% 
    select(LM, new_test, stats1, values) %>% 
    rename("stats" = "stats1", "test" = "new_test") %>% 
    pivot_wider(id_cols = c("LM", "test"),
                names_from = stats,
                values_from = values)

eskd_vimp %>% 
    arrange(LM) %>%  
    mutate(LM = factor(LM, levels = c(0.5, 1, 1.5, 2, 2.5, 3),
                       labels = c("0.5 years", "1 year", "1.5 years", "2 years",
                                  "2.5 years", "3 years"))) %>%
    ggplot(mapping = aes(x = fct_reorder(test, mean), y = mean, fill = test)) +
    geom_bar(position = "dodge", stat = "identity") + 
    geom_errorbar(mapping = aes(ymin = lower, ymax = upper), stat = "identity") + 
    facet_wrap(~ LM) + scale_y_continuous(expand = c(0,0)) +
    coord_flip() + 
    labs(title = "VIMP for ESKD by Landmark Times", 
         x = element_blank(), y = "VIMP") +
    theme(legend.position = "none", 
          plot.title = element_text(size = 17, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 15)) 

death_vimp <- dense3_vimp %>% 
    select(vimp) %>% 
    unnest(vimp) %>% 
    summarise(across(sex:wcc, \(x) confidence_interval (x, 0.95)),
              .by = c("event", "LM")) %>% 
    unnest_wider(sex:wcc, names_sep = "_") %>% 
    filter(event == "death")

death_vimp <- death_vimp %>% 
    pivot_longer(cols = sex_mean:wcc_upper, names_to = "test", 
                 values_to = "values") %>% 
    separate_wider_delim(cols = test, delim = "_", 
                         names = c("test", "stats", "stats1"),
                         too_few = "align_start") %>% 
    mutate(stats1 = if_else(is.na(stats1), stats, stats1)) %>% 
    mutate(new_test = case_match(test,
                                 "age" ~ "age_init",
                                 "gn" ~ "gn_fct", .default = test)) %>% 
    select(LM, new_test, stats1, values) %>% 
    rename("stats" = "stats1", "test" = "new_test") %>% 
    pivot_wider(id_cols = c("LM", "test"),
                names_from = stats,
                values_from = values)

death_vimp %>% 
    arrange(LM) %>%  
    mutate(LM = factor(LM, levels = c(0.5, 1, 1.5, 2, 2.5, 3),
                       labels = c("0.5 years", "1 year", "1.5 years", "2 years",
                                  "2.5 years", "3 years"))) %>%
    ggplot(mapping = aes(x = fct_reorder(test, mean), y = mean, fill = test)) +
    geom_bar(position = "dodge", stat = "identity") + 
    geom_errorbar(mapping = aes(ymin = lower, ymax = upper), stat = "identity") + 
    facet_wrap(~ LM) + scale_y_continuous(expand = c(0,0)) +
    coord_flip() + 
    labs(title = "VIMP for Death by Landmark Times",
         x = element_blank(), y = "VIMP") +
    theme(legend.position = "none", 
          plot.title = element_text(size = 17, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 15))

vimp_rank <- rbind(eskd_vimp, death_vimp)

# age_init, egfr, haemoglobin, sodium, and wcc pvalue are not normally distributed
vimp_rank %>% 
    select(LM:mean) %>% 
    rename("vimp_mean" = "mean") %>% 
    group_by(test) %>% 
    summarise(normality = rstatix::shapiro_test(vimp_mean)$p.value)

vimp_rank %>% 
    select(test:mean) %>% 
    summarise(median_vimp = median(mean),.by = test) %>% 
    arrange(-median_vimp) %>% 
    gt()

# end ---------------------------------------------------------------------