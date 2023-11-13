library(tidyverse)
library(skimr)
library(gt)
library(gtsummary)
library(patchwork)
library(rstatix)
library(data.table)

# setup -------------------------------------------------------------------

dense3_dt <- read_rds("main_dense3_dt.rds")
acr3_dt <- read_rds("main_acr3_dt.rds")
sens3_dt <- read_rds("sens_dense3_dt.rds")

longi_cov_full <- c("albumin", "alkphos", "bicarb", "calcium", "chloride",
                    "egfr", "glucose", "haemoglobin", "phosphate", "platelet",
                    "potassium", "sodium", "wcc")
base_cov_full <- c("sex", "age_init", "gn_fct")

longi_cov_acr <- c("albumin", "albcreat_ratio", "alkphos", "bicarb", "calcium", 
                   "chloride",  "egfr", "glucose", "haemoglobin", "phosphate", 
                   "platelet", "potassium", "sodium", "wcc")

predict_horizon <- 5
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)
seed <- 7
cross_val <- 5

longi_cov_top5 <- c("albumin", "bicarb", "chloride", "egfr", "haemoglobin")
base_cov_reduced <- c("age_init")

kfre_full <- c("albcreat_ratio", "egfr", "albumin", "phosphate", "bicarb", 
               "calcium")

test_dataset <- "experiment_test.rds"
source("utility_functions.R")


# baseline_summary --------------------------------------------------------

dense3_missing_relyear0 <- dense3_dt %>% 
    filter(relyear == 0, test %in% longi_cov_full) %>% 
    select(id, relyear, test, value) %>% 
    pivot_wider(id_cols = c("id", "relyear"),
                names_from = test,
                values_from = value) %>% 
    skim()

dense3_missing_relyear0 %>% 
    select(skim_variable, n_missing, complete_rate) %>% gt()

acr3_missing_relyear0 <- acr3_dt %>% 
    filter(relyear == 0, test %in% longi_cov_acr) %>% 
    select(id, relyear, test, value) %>% 
    pivot_wider(id_cols = c("id", "relyear"),
                names_from = test,
                values_from = value) %>% 
    skim()

acr3_missing_relyear0 %>% 
    select(skim_variable, n_missing, complete_rate) %>% gt()

sens3_summary <- sens3_dt %>% 
    filter(relyear == 0, test %in% longi_cov_full) %>% 
    pivot_wider(names_from = test, 
                values_from = value)

dense3_summary <- dense3_dt %>% 
    filter(relyear == 0, test %in% longi_cov_full) %>% 
    pivot_wider(names_from = test,
                values_from = value)

dense3_summary %>% tbl_summary()

acr3_summary <- acr3_dt %>% 
    filter(relyear == 0, test %in% longi_cov_acr) %>% 
    pivot_wider(names_from = test,
                values_from = value)

acr3_summary %>% tbl_summary()

sens3_summary <- sens3_dt %>% 
    filter(relyear == 0, test %in% longi_cov_full) %>% 
    pivot_wider(names_from = test, 
                values_from = value)

sens3_summary %>% tbl_summary()

dense3_dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    filter(age_init >= 70) %>% 
    count(status_compete)

# compiling internal performance ------------------------------------------

main_sens_dt <- read_rds("main_vs_sens_performance.rds")
main_acr3_dt <- read_rds("main_dense3_vs_acr3_performance.rds")
main_top_dt <- read_rds("main_dense3_locf_exploration_performance.rds")

temp <- main_sens_dt %>% 
    select(locf_perf) %>% 
    unnest_longer(locf_perf) %>% 
    map(~.x) %>% 
    bind_rows()

sens_dt <- main_sens_dt %>% 
    unnest_longer(locf_perf) %>% 
    select(method) %>% 
    bind_cols(temp) %>% 
    filter(method == "sens")

temp <- main_acr3_dt %>% 
    select(locf_perf) %>% 
    unnest_longer(locf_perf) %>% 
    map(~ .x) %>% 
    bind_rows()

dense3_locf <- main_acr3_dt %>% 
    unnest_longer(locf_perf) %>% 
    select(id_name) %>% 
    bind_cols(temp) %>% 
    filter(id_name == "dense3")

dense3_locf <- dense3_locf %>% 
    rename("method" = "id_name") 

dense3_locf$method <- "locf"

acr3_dt <- main_acr3_dt %>% 
    unnest_longer(locf_perf) %>% 
    select(id_name) %>% 
    bind_cols(temp) %>% 
    filter(id_name == "acr3")

acr3_dt <- acr3_dt %>% 
    rename("method" = "id_name")

dense3_lmes <- main_acr3_dt %>% 
    filter(id_name == "dense3")

temp <- dense3_lmes %>% 
    select(lme_perf) %>% 
    unnest_longer(lme_perf) %>% 
    map(~ .x) %>% 
    bind_rows()

dense3_lme <- dense3_lmes %>% 
    unnest_longer(lme_perf) %>% 
    select(id_name) %>% 
    bind_cols(temp) %>% 
    rename("method" = "id_name")

dense3_lme$method <- "lme"

temp <- dense3_lmes %>% 
    select(lmepoly_perf) %>% 
    unnest_longer(lmepoly_perf) %>% 
    map(~ .x) %>% 
    bind_rows()

lme_poly <- dense3_lmes %>% 
    unnest_longer(lmepoly_perf) %>% 
    select(id_name) %>% 
    bind_cols(temp) %>% 
    rename("method" = "id_name")

lme_poly$method <- "lmepoly"

temp <- main_top_dt %>% 
    select(performance) %>% 
    unnest_longer(performance) %>% 
    map(~ .x) %>% 
    bind_rows()

main_top_dt <- main_top_dt %>% 
    unnest_longer(performance) %>% 
    select(id_name) %>% 
    bind_cols(temp) %>% 
    rename("method" = "id_name")

topten <- main_top_dt %>% 
    filter(method == "top10")

topfive <- main_top_dt %>% 
    filter(method == "top5")

performance <- data.table::rbindlist(list(dense3_locf, sens_dt, acr3_dt, 
                                          dense3_lme, lme_poly, topten, topfive))

# write_rds(performance, "internal_performance.rds")

# internal performance analysis -------------------------------------------

performance <- performance %>% 
    mutate(method = factor(method, 
                           levels = c("locf", "lme", "lmepoly", "top10", "top5", "sens", "acr3"),
                           labels = c("LOCF", "LME", "LME Poly", "Top 10", "Top 5", "sens", "acr3")))

error_eskd <- performance %>% 
    select(method, LM, contains("error")) %>% 
    pivot_longer(cols = contains("error"),
                 names_to = "event",
                 values_to = "error") %>% 
    separate_wider_delim(event, delim = "_", names = c("not_used", "event")) %>% 
    select(-not_used) %>% 
    filter(event == "eskd") %>% 
    arrange(LM)

error_eskd <- error_eskd %>% 
    mutate(c_index = 100-error)

error_death <- performance %>% 
    select(method, LM, contains("error")) %>% 
    pivot_longer(cols = contains("error"),
                 names_to = "event",
                 values_to = "error") %>% 
    separate_wider_delim(event, delim = "_", names = c("not_used", "event")) %>% 
    select(-not_used) %>% 
    filter(event == "death") %>% 
    arrange(LM)

error_death <- error_death %>% 
    mutate(c_index = 100-error)

brier_eskd <- performance %>% 
    select(method, LM, contains("brier")) %>% 
    pivot_longer(cols = contains("brier"),
                 names_to = "event",
                 values_to = "brier") %>% 
    separate_wider_delim(event, delim = "_", names = c("not_used", "event")) %>% 
    select(-not_used) %>% 
    filter(event == "eskd") %>% 
    arrange(LM)

brier_death <- performance %>% 
    select(method, LM, contains("brier")) %>% 
    pivot_longer(cols = contains("brier"),
                 names_to = "event",
                 values_to = "brier") %>% 
    separate_wider_delim(event, delim = "_", names = c("not_used", "event")) %>% 
    select(-not_used) %>% 
    filter(event == "death") %>% 
    arrange(LM)

# passed Levene Test and normality of residuals
error_eskd <- error_eskd %>% 
    group_by(LM) %>% 
    nest() %>% 
    mutate(levene = map(data, \(data) levene_test(data, c_index ~ method))) %>% 
    mutate(lin_reg = map(data, \(data) lm(c_index ~ method, data = data))) %>% 
    mutate(residual = map(lin_reg, \(lin_model) residuals(lin_model))) %>% 
    mutate(normality = map(residual, \(resid) shapiro_test(resid)))
    
anova_eskd <- error_eskd %>% 
    select(LM, data) %>% 
    mutate(anova_testing = map(data, \(data) anova_test(data, error ~ method))) 

anova_eskd$anova_testing[[2]] #LM 1

#LME and acr3 with 95% CI [-9.0939009, -0.5587974] P 0.0189
anova_eskd$data[[2]] %>% 
    tukey_hsd(error ~ method)

# all ns
krus_eskd <- error_eskd %>% 
    select(LM, data) %>% 
    mutate(kruskal = map(data, \(data) kruskal_test(data, c_index ~ method)))

table_error_eskd <- error_eskd %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    group_by(LM, method) %>% 
    get_summary_stats(c_index, type = "quantile") %>% 
    select(LM, method, "25%", "50%", "75%") %>% 
    mutate(event = "eskd")

# LM 0.5 and 3 did not pass Levene test, LM 2 and 3 failed normality of residuals
error_death <- error_death %>% 
    group_by(LM) %>% 
    nest() %>% 
    mutate(levene = map(data, \(data) levene_test(data, c_index ~ method))) %>%
    mutate(lin_reg = map(data, \(data) lm(c_index ~ method, data = data))) %>% 
    mutate(residual = map(lin_reg, \(lin_model) residuals(lin_model))) %>% 
    mutate(normality = map(residual, \(resid) shapiro_test(resid)))

# all ns
krus_death <- error_death %>% 
    select(LM, data) %>% 
    mutate(kruskal = map(data, \(data) kruskal_test(data, c_index ~ method)))

table_error_death <- error_death %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    group_by(LM, method) %>% 
    get_summary_stats(c_index, type = "quantile") %>% 
    select(LM, method, "25%", "50%", "75%") %>% 
    mutate(event = "death")

table_error <- table_error_eskd %>% 
    bind_rows(table_error_death) 

# write_csv(table_error, "error_iqr.csv")

# all but LM 3 did not pass Levene Test and none pass the normality of residual
brier_eskd <- brier_eskd %>% 
    group_by(LM) %>% 
    nest() %>% 
    mutate(levene = map(data, \(data) levene_test(data, brier ~ method))) %>% 
    mutate(lin_reg = map(data, \(data) lm(brier ~ method, data = data))) %>% 
    mutate(residual = map(lin_reg, \(lin_model) residuals(lin_model))) %>% 
    mutate(normality = map(residual, \(resid) shapiro_test(resid))) 

krus_brier_eskd <- brier_eskd %>% 
    select(LM, data) %>% 
    mutate(kruskal = map(data, \(data) kruskal_test(data, brier ~ method)))

#all ns
krus_brier_eskd %>% 
    mutate(man_whit = map(data, \(data) wilcox_test(data, brier ~ method, p.adjust.method = "bonferroni"))) %>% 
    unnest(man_whit) %>% View()

table_brier_eskd <- brier_eskd %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    group_by(LM, method) %>% 
    get_summary_stats(brier, type = "quantile") %>% 
    select(LM, method, "25%", "50%", "75%") %>% 
    mutate(event = "eskd")

#all pass Levene test, but all did not pass the normality of residuals
brier_death <- brier_death %>% 
    group_by(LM) %>% 
    nest() %>% 
    mutate(levene = map(data, \(data) levene_test(data, brier ~ method))) %>% 
    mutate(lin_reg = map(data, \(data) lm(brier ~ method, data = data))) %>% 
    mutate(residual = map(lin_reg, \(lin_model) residuals(lin_model))) %>% 
    mutate(normality = map(residual, \(resid) shapiro_test(resid))) 

krus_brier_death <- brier_death %>% 
    select(LM, data) %>% 
    mutate(kruskal = map(data, \(data) kruskal_test(data, brier ~ method)))

# all ns
krus_brier_death %>% 
    mutate(man_whit = map(data, \(data) wilcox_test(data, brier ~ method, p.adjust.method = "bonferroni"))) %>% 
    unnest(man_whit) %>% View()

table_brier_death <- brier_death %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    group_by(LM, method) %>% 
    get_summary_stats(brier, type = "quantile") %>% 
    select(LM, method, "25%", "50%", "75%") %>% 
    mutate(event = "death")

table_brier <- table_brier_eskd %>% 
    bind_rows(table_brier_death)

# write_csv(table_brier, "brier_iqr.csv")

# visualisation of performance --------------------------------------------

main_graph <- c("LOCF", "LME", "LME Poly", "Top 10", "Top 5")
sens_graph <- c("LOCF", "sens", "acr3")

error_eskd %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    filter(method %in% sens_graph) %>% 
    mutate(LM = factor(LM)) %>% 
    ggplot(mapping = aes(x = LM, y = c_index, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "ESKD Harrell's C-index with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "C-index (%)",
         fill = "Method") + 
    theme(plot.title = element_text(size = 17, hjust = 0.5), 
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(50, 100, 5), limits = c(70, 95))

error_death %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    filter(method %in% sens_graph) %>% 
    mutate(LM = factor(LM)) %>% 
    ggplot(mapping = aes(x = LM, y = c_index, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "Death Harrell's C-index with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "C-index (%)",
         fill = "Method") + 
    theme(plot.title = element_text(size = 17, hjust = 0.5), 
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(50, 100, 5), limits = c(70, 95))

brier_eskd %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    filter(method %in% sens_graph) %>% 
    mutate(LM = factor(LM)) %>% 
    ggplot(mapping = aes(x = LM, y = brier, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "ESKD Integrated Brier Score with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "Brier score",
         fill = "Method") +
    theme(plot.title = element_text(size = 17, hjust = 0.5),
          axis.title.x = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(0, 0.07, 0.005), limits = c(0, 0.07))

brier_death %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    filter(method %in% sens_graph) %>% 
    mutate(LM = factor(LM)) %>% 
    ggplot(mapping = aes(x = LM, y = brier, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "Death Integrated Brier Score with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "Brier score",
         fill = "Method") +
    theme(plot.title = element_text(size = 17, hjust = 0.5),
          axis.title.x = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(0, 0.07, 0.005), limits = c(0, 0.07))

# external performance ----------------------------------------------------
error_temp0.5 <- read_rds("wa_performance/top5_0.5_error.rds")
error_temp1 <- read_rds("wa_performance/top5_1_error.rds")
error_temp1.5 <- read_rds("wa_performance/top5_1.5_error.rds")
error_temp2 <- read_rds("wa_performance/top5_2_error.rds")
error_temp2.5 <- read_rds("wa_performance/top5_2.5_error.rds")
error_temp3 <- read_rds("wa_performance/top5_3_error.rds")

error_temp0.5 <- data.table(error_temp0.5$t)
error_temp0.5[, LM := 0.5]
error_temp1 <- data.table(error_temp1$t)
error_temp1[, LM := 1]
error_temp1.5 <- data.table(error_temp1.5$t)
error_temp1.5[, LM := 1.5]
error_temp2 <- data.table(error_temp2$t)
error_temp2[, LM := 2]
error_temp2.5 <- data.table(error_temp2.5$t)
error_temp2.5[, LM := 2.5]
error_temp3 <- data.table(error_temp3$t)
error_temp3[, LM := 3]

ext_error <- rbindlist(list(error_temp0.5, error_temp1, 
                            error_temp1.5, error_temp2,
                            error_temp2.5, error_temp3))

ext_error <- ext_error %>% 
    mutate(method = "external") %>% 
    mutate(method = factor(method, levels = c("external"), labels = c("External"))) %>% 
    mutate(event = "eskd") %>% 
    rename("error" = "V2") %>% 
    select(LM, method, event, error) %>% 
    mutate(error = 100*error)

ext_error <- ext_error %>% 
    mutate(c_index = 100-error)

int_ext_error <- error_eskd %>% 
    select(LM, data) %>% 
    unnest(cols = c(data)) %>% 
    filter(method == "Top 5") %>% 
    bind_rows(ext_error) 

# all LM pass Levene and Shapiro test
int_ext_error %>% 
    group_by(LM) %>% 
    nest() %>% 
    mutate(levene = map(data, \(data) levene_test(data, c_index ~ method))) %>% 
    mutate(shapiro = map(data, \(data) shapiro_test(data$c_index))) %>% 
    mutate(wilcox = map(data, \(data) wilcox_test(formula = data$c_index ~ data$method)))

# all ns
int_ext_error %>% 
    filter(LM == 3) %>% 
    wilcox_test(c_index ~ method)

int_ext_error %>% 
    group_by(LM, method) %>% 
    get_summary_stats(c_index, type = "quantile") 

int_ext_error %>% 
    mutate(LM = factor(LM)) %>% 
    ggplot(mapping = aes(x = LM, y = c_index, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "ESKD Harrell's C-index with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "C-index (%)",
         fill = "Method") + 
    theme(plot.title = element_text(size = 17, hjust = 0.5), 
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(50, 100, 5), limits = c(70, 95))


# compare KFRE ------------------------------------------------------------

main_dt <- tibble(
    dt_path = c("main_dense3_dt.rds"),
    base_cov = list(base_cov_reduced),
    longi_cov = list(longi_cov_top5)
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
           id_name = rep(c("dense3"), each = 5)) %>% 
    relocate(id_name, .before = id)

main_dt <- main_dt %>% 
    mutate(train_dt = pmap(list(dt, splits), extract_train_dt),
           test_dt = pmap(list(dt, splits), extract_test_dt)) %>% 
    select(id_name, id, base_cov, longi_cov, train_dt, test_dt)


# training model for cross-validation -------------------------------------

final_models_modified <- function(dt, longi_cov, base_cov, model_name, ...){
    # creating trained models with all data for external validation #
    long_dt <- dt %>%
        select(id, relyear, !!longi_cov)
    
    impute_dt <- long_dt %>%
        filter(relyear == 0)
    
    temp_dt <- mice(impute_dt, m = 1, maxit = 50, method = "pmm", seed = seed)
    impute_dt <- complete(temp_dt)
    
    dt <- impute_dt %>%
        bind_rows(long_dt) %>%
        arrange(id) %>%
        distinct(id, relyear, .keep_all = TRUE) %>%
        group_by(id) %>%
        fill(!!longi_cov) %>%
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
                      covs = list(fixed = base_cov, varying = longi_cov),
                      format = "long", id = "id", rtime = "relyear",
                      right = FALSE)
        super_dt <- rbind(super_dt, temp)
    }
    
    tune_matrix = NULL
    
    for (lmx in landmark){
        print(glue::glue("Processing landmark: {lmx}"))
        
        dt_lm <- super_dt %>%
            filter(LM == lmx)
        
        model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                          str_c(c(base_cov, longi_cov), collapse = " + "),
                                          sep = "~"))
        
        tune_rsf <- tune.rfsrc(model_formula, data = dt_lm,
                               ntreeTry = 500, nodesizeTry = c(1:9, seq(10, 100, 5)))
        
        tune_matrix <- rbind(tune_matrix, tune_rsf$optimal)
        
        rsf_model <- rfsrc(formula = model_formula, data = dt_lm, ntree = 1000,
                           splitrule = "logrankCR", importance = FALSE, statistics = FALSE,
                           mtry = tune_rsf$optimal[[2]], nodesize = tune_rsf$optimal[[1]], 
                           save.memory = TRUE)
        
        rsf_model_name <- str_c(model_name, lmx, "locf_models.rds", sep = "_")
        saveRDS(rsf_model, file = rsf_model_name)
    }
}

main_dt %>% 
    mutate(model_name = str_c("exper_top5_fold", seq(1, 5))) %>% 
    mutate(training = pmap(list(train_dt, longi_cov, base_cov, model_name),
                           final_models_modified, landmark, predict_horizon))

compare_dt <- main_dt %>% 
    select(id, base_cov, longi_cov, test_dt)

variables_extraction <- unique(c(longi_cov_top5, kfre_full))
base_variables_extraction <- c(base_cov_reduced, "sex")

compare_dt <- compare_dt %>% 
    rename("fold" = "id") %>% 
    select(fold, test_dt) %>% 
    unnest(test_dt) %>% 
    nest(data = everything(), .by = c("id", "fold")) %>% 
    mutate(landmark = map(data, extract_landmark_for_dynamic_plot, landmark))
    
compare_dt <- compare_dt %>% 
    mutate(base_cov = list(base_variables_extraction),
           longi_cov = list(variables_extraction)) %>% 
    mutate(data_points = map_int(landmark, length)) %>% 
    filter(data_points >= 3) %>% 
    select(-data_points) %>% 
    unnest_longer(landmark)

extract_landmark_dt_modified <- function(dt, base_cov, longi_cov, landmark, ...){
    # external or test data only #
    # extracting dataset per landmark #
    base_cov <- unlist(base_cov)
    
    longi_cov <- unlist(longi_cov)

    dt <- cutLM(data = dt, outcome = list(time = "time_compete",
                                          status = "status_compete_num"),
                LM = landmark, horizon = landmark + predict_horizon,
                covs = list(fixed = base_cov, varying = longi_cov),
                format = "long", id = "id", rtime = "relyear", 
                right = FALSE) 
    dt
}

compare_dt <- compare_dt %>% 
    mutate(lmx_dt = pmap(list(data, base_cov, longi_cov, landmark), 
                         extract_landmark_dt_modified, predict_horizon))

landmark_dt <- compare_dt %>% 
    select(fold, lmx_dt) %>% 
    unnest(lmx_dt) %>% 
    group_by(id, fold) %>% 
    fill(!!variables_extraction, .direction = "downup") %>% 
    ungroup(fold, id) %>% 
    drop_na()

kfre_dt <- landmark_dt %>% 
    group_by(id, fold) %>% 
    slice_tail() %>% 
    ungroup()

#total 1396 test ids with 137 eskd and 137 death
kfre_dt %>% 
    count(status_compete_num)

kfre_dt <- kfre_dt %>% 
    mutate(corrected_age = age_init + LM,
           log_acr = log(albcreat_ratio/0.113),
           alb = albumin/10,
           phos = phosphate/0.3229,
           ca = calcium/0.2495,
           sex_num = as.numeric(sex)-1,
           kfre2y_pred = unlist(pmap(list(corrected_age, sex_num, egfr, log_acr, alb, phos, bicarb, ca), kfre8_model_2y)),
           kfre5y_pred = unlist(pmap(list(corrected_age, sex_num, egfr, log_acr, alb, phos, bicarb, ca), kfre8_model_5y)))

predicted_cif_modified <- function(landmark_dt, landmark, analysis_models, ...){
    # measuring CIF for ESKD and Death event #
    model_name <- str_c(analysis_models, landmark, "locf_models.gz", sep = "_")
    rsf_model <- readRDS(model_name)
    
    rsf_prediction <- predict.rfsrc(object = rsf_model,
                                    newdata = landmark_dt)
    
    column_cif <- sindex(jump.times = rsf_prediction$time.interest,
                         eval.times = seq(landmark, 
                                          landmark + predict_horizon, 
                                          0.25))
    
    column_cif[1] <- 1
    
    cif_dt <- data.table(
        pred_relyear = seq(landmark, landmark + predict_horizon, 0.25),
        eskd_cif = rsf_prediction$cif[, column_cif, "CIF.1"],
        death_cif = rsf_prediction$cif[, column_cif, "CIF.2"]
    )
    cif_dt
}

folds <- str_c("fold", seq(1, 5))
compile_dt <- NULL

for (fold_roi in folds) {
    for (lmx in landmark){
        model_name <- str_c("exper_top5", fold_roi, lmx, "locf_models.rds", sep = "_")
        
        print(model_name)
        
        model <- read_rds(model_name)

        loop_dt <- landmark_dt %>%
            filter(fold == fold_roi & LM == lmx)

        rsf_prediction <- predict.rfsrc(object = model,
                                        newdata = loop_dt)

        column_cif <- sindex(jump.times = rsf_prediction$time.interest,
                             eval.times = seq(lmx,
                                              lmx + predict_horizon,
                                              0.25))

        column_cif[1] <- 1

        column_name <- seq(lmx, lmx + predict_horizon, 0.25)

        eskd_matrix <- as.matrix(rsf_prediction$cif[, column_cif, "CIF.1"])
        colnames(eskd_matrix) <- column_name
        cif_dt <- data.table(eskd_matrix)

        result_dt <- loop_dt %>%
            select(fold:id) %>%
            bind_cols(cif_dt) %>%
            pivot_longer(cols = !c("fold", "id"),
                         names_to = "pred_relyear",
                         values_to = "eskd_cif") %>%
            nest(data = pred_relyear:eskd_cif, .by = c("fold", "id")) %>%
            mutate(LM = lmx)

        compile_dt <- rbind(compile_dt, result_dt)
    }
}

kfre_dt <- kfre_dt %>% 
    rename("max_LM" = "LM")

compile_dt <- compile_dt %>% 
    left_join(y = select(kfre_dt, fold, id, max_LM), by = c("fold", "id"))

summarising_cif_pred <- function(dt, max_lmx, year_interest){
    predict_horizon <- max_lmx + year_interest
    
    eskd_cif <- dt %>% 
        filter(pred_relyear == predict_horizon) %>% 
        pull(eskd_cif)
    
    eskd_cif
}

compile_dt <- compile_dt %>% 
    mutate(rsf2y_pred = pmap_dbl(list(data, max_LM, 2), 
                                 summarising_cif_pred),
           rsf5y_pred = pmap(list(data, max_LM, 5), 
                             summarising_cif_pred))

rsf2y_pred <- compile_dt %>% 
    select(fold, id, rsf2y_pred) %>% 
    group_by(fold, id) %>% 
    summarise(mean_rsf2y_pred = mean(rsf2y_pred))

rsf2y_pred <- rsf2y_pred %>% 
    rename("rsf2y_pred" = "mean_rsf2y_pred")

rsf5y_pred <- compile_dt %>% 
    select(fold, id, rsf5y_pred) %>% 
    unnest(rsf5y_pred)

kfre_dt <- kfre_dt %>% 
    arrange(fold, id) %>% 
    select(fold:status_compete_num, max_LM, kfre2y_pred, kfre5y_pred) %>% 
    left_join(rsf2y_pred, by = c("fold", "id")) %>% 
    left_join(rsf5y_pred, by = c("fold", "id"))

kfre_dt <- kfre_dt %>% 
    mutate(event2y = if_else((max_LM+2)>=time_compete, status_compete_num, 0),
           deathcens_2y = if_else(event2y == 2, 0, event2y),
           event5y = if_else((max_LM+5)>=time_compete, status_compete_num, 0),
           deathcens_5y = if_else(event5y == 2, 0, event5y))

final_dt <- kfre_dt %>% 
    # filter(status_compete_num !=2) %>% # remove if death censored
    mutate(kfre2y_pred_class = if_else(kfre2y_pred >= 0.15, "eskd", "no"),
           kfre5y_pred_class = if_else(kfre5y_pred >= 0.15, "eskd", "no"),
           rsf2y_pred_class = if_else(rsf2y_pred >= 0.15, "eskd", "no"),
           rsf5y_pred_class = if_else(rsf5y_pred >= 0.15, "eskd", "no")) %>% 
    mutate(kfre2y_pred_class = factor(kfre2y_pred_class, levels = c("no", "eskd")),
           kfre5y_pred_class = factor(kfre5y_pred_class, levels = c("no", "eskd")),
           rsf2y_pred_class = factor(rsf2y_pred_class, levels = c("no", "eskd")),
           rsf5y_pred_class = factor(rsf5y_pred_class, levels = c("no", "eskd"))) %>% 
    mutate(deathcens_2y_class = if_else(deathcens_2y == 1, "eskd", "no"),
           deathcens_5y_class = if_else(deathcens_5y == 1, "eskd", "no")) %>% 
    mutate(deathcens_2y_class = factor(deathcens_2y_class, levels = c("no", "eskd")),
           deathcens_5y_class = factor(deathcens_5y_class, levels = c("no", "eskd")))


death_cens_dt <- kfre_dt %>% 
    mutate(death_cens_status = if_else(status_compete_num == 2, 0, status_compete_num)) 

metrics_compiled <- NULL

for (fold_interest in folds){
    dt_interest <- death_cens_dt %>% 
        filter(fold == fold_interest)
    
    kfre2y_c <- survcomp::concordance.index(dt_interest$kfre2y_pred, 
                                            dt_interest$time_compete, 
                                            dt_interest$status_compete_num)
    
    rsf2y_c <- survcomp::concordance.index(dt_interest$rsf2y_pred, 
                                           dt_interest$time_compete, 
                                           dt_interest$status_compete_num)
    
    kfre5y_c <- survcomp::concordance.index(dt_interest$kfre5y_pred, 
                                            dt_interest$time_compete, 
                                            dt_interest$status_compete_num)
    
    rsf5y_c <- survcomp::concordance.index(dt_interest$rsf5y_pred, 
                                            dt_interest$time_compete, 
                                            dt_interest$status_compete_num)
    
    
    metrics_dt <- data.table(
        fold = fold_interest,
        metrics_var = "c_index",
        rsf2y = rsf2y_c$c.index,
        kfre2y = kfre2y_c$c.index,
        rsf5y = rsf5y_c$c.index,
        kfre5y = kfre5y_c$c.index
    )
    metrics_compiled <- rbind(metrics_compiled, metrics_dt)
}

# library(yardstick)
# metrics <- c("accuracy", "sens", "spec", "ppv")

# metrics_compiled <- NULL
# 
# for (fold_interest in folds){
#     dt_interest <- final_dt %>% 
#         filter(fold == fold_interest)
#     
#     rsf2y_metrics <- conf_mat(dt_interest, truth = deathcens_2y_class, estimate = rsf2y_pred_class) %>% 
#         summary(event_level = "second") %>% 
#         filter(.metric %in% metrics) %>% 
#         pull(.estimate)
#     
#     kfre2y_metrics <- conf_mat(dt_interest, truth = deathcens_2y_class, estimate = kfre2y_pred_class) %>% 
#         summary(event_level = "second") %>% 
#         filter(.metric %in% metrics) %>% 
#         pull(.estimate)
#     
#     rsf5y_metrics <- conf_mat(dt_interest, truth = deathcens_5y_class, estimate = rsf5y_pred_class) %>% 
#         summary(event_level = "second") %>% 
#         filter(.metric %in% metrics) %>% 
#         pull(.estimate)
#     
#     kfre5y_metrics <- conf_mat(dt_interest, truth = deathcens_5y_class, estimate = kfre5y_pred_class) %>% 
#         summary(event_level = "second") %>% 
#         filter(.metric %in% metrics) %>% 
#         pull(.estimate)
#     
#     metrics_dt <- data.table(
#         fold = fold_interest,
#         metrics_var = metrics,
#         rsf2y = rsf2y_metrics,
#         kfre2y = kfre2y_metrics,
#         rsf5y = rsf5y_metrics,
#         kfre5y = kfre5y_metrics
#     )
#     metrics_compiled <- rbind(metrics_compiled, metrics_dt)
# }

#accuracy, sens, spec, ppv are normally distributed
metrics_compiled %>% 
    filter(metrics_var == "c_index") %>% 
    select(rsf2y, kfre2y) %>% 
    pivot_longer(cols = everything(), names_to = "group", values_to = "value") %>% 
    t_test(value ~ group)

metrics_compiled %>% 
    filter(metrics_var == "c_index") %>% 
    select(rsf5y, kfre5y) %>% 
    pivot_longer(cols = everything(), names_to = "group", values_to = "value") %>% 
    t_test(value ~ group)

metrics_compiled %>% 
    filter(metrics_var == "c_index") %>% 
    select(rsf2y:kfre2y) %>% 
    map(\(mu) confidence_interval(mu, 0.95))

metrics_compiled %>% 
    filter(metrics_var == "c_index") %>% 
    select(rsf5y:kfre5y) %>% 
    map(\(mu) confidence_interval(mu, 0.95))
 
# end ---------------------------------------------------------------------