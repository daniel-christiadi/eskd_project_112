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
    dataset = c("main_dense3_dt.rds", "main_acr3_dt.rds"),
    base_cov = list(base_cov, base_cov),
    longi_cov = list(longi_cov, longi_cov_acr)
)

main_dt <- main_dt %>% 
    mutate(locf_perf = pmap(list(dataset, longi_cov, base_cov), cross_val_performance_locf))

main_dt <- main_dt %>% 
    mutate(lme_perf = pmap(list(dataset, longi_cov, base_cov), cross_val_performance_lme))

main_dt <- main_dt %>% 
    mutate(lmepoly_perf = pmap(list(dataset, longi_cov, base_cov), cross_val_performance_lmepoly))

# write_rds(main_dt, "main_dense3_vs_acr3_performance.rds")

# dense3 vs acr3 LOCF visualisation --------------------------------------------

main_dt <- read_rds("main_dense3_vs_acr3_performance.rds")

# main_dt %>% 
#     select(locf_perf) %>% 
#     unnest_longer(locf_perf) %>% 
#     map(~ .x)
#     bind_rows()
#     

locf_performance <- main_dt %>% 
    select(locf_perf) %>% 
    unnest_longer(locf_perf)
    
dense3_acr3_performance <- locf_performance$locf_perf

dense3_acr3_performance <- dense3_acr3_performance %>% 
    separate(col = analysis, into = c("main", "data", "dt", "method")) %>% 
    select(LM:fold, data)

error1 <- dense3_acr3_performance %>% 
    filter(data == "dense3") %>% 
    select(-(brier_eskd:data)) %>% 
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

error2 <- dense3_acr3_performance %>% 
    filter(data == "acr3") %>% 
    select(-(brier_eskd:data)) %>% 
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

brier1 <- dense3_acr3_performance %>% 
    filter(data == "dense3") %>% 
    select(-(error_eskd:error_death), -(fold:data)) %>% 
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

brier2 <- dense3_acr3_performance %>%
    filter(data == "acr3") %>% 
    select(-(error_eskd:error_death), -(fold:data)) %>% 
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

lme_performance <- main_dt %>% 
    filter(dataset == "main_dense3_dt.rds") %>% 
    select(lme_perf) %>% 
    unnest_longer(lme_perf)

lme_performance <- lme_performance$lme_perf

lme_performance <- lme_performance %>% 
    separate(col = analysis, into = c("main", "data", "dt", "method")) %>% 
    select(LM:fold, data)

lmepoly_performance <- main_dt %>% 
    filter(dataset == "main_dense3_dt.rds") %>% 
    select(lmepoly_perf) %>% 
    unnest_longer(lmepoly_perf)

lmepoly_performance <- lmepoly_performance$lmepoly_perf

lmepoly_performance <- lmepoly_performance %>% 
    separate(col = analysis, into = c("main", "data", "dt", "method")) %>% 
    select(LM:fold, data)

error2 <- lme_performance %>% 
    select(-(brier_eskd:data)) %>% 
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

error3 <- lmepoly_performance %>% 
    select(-(brier_eskd:data)) %>% 
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

brier2 <- lme_performance %>%
    select(-(error_eskd:error_death), -(fold:data)) %>% 
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

brier3 <- lmepoly_performance %>%
    select(-(error_eskd:error_death), -(fold:data)) %>% 
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


