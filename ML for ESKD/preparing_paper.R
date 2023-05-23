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

longi_cov_acr <- c("albumin", "albcreat_ratio", "alkphos", "bicarb", "calcium", 
                   "chloride",  "egfr", "glucose", "haemoglobin", "phosphate", 
                   "platelet", "potassium", "sodium", "wcc")

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

error_death <- performance %>% 
    select(method, LM, contains("error")) %>% 
    pivot_longer(cols = contains("error"),
                 names_to = "event",
                 values_to = "error") %>% 
    separate_wider_delim(event, delim = "_", names = c("not_used", "event")) %>% 
    select(-not_used) %>% 
    filter(event == "death") %>% 
    arrange(LM)

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
    mutate(levene = map(data, \(data) levene_test(data, error ~ method))) %>% 
    mutate(lin_reg = map(data, \(data) lm(error ~ method, data = data))) %>% 
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
    mutate(kruskal = map(data, \(data) kruskal_test(data, error ~ method)))

table_error_eskd <- error_eskd %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    group_by(LM, method) %>% 
    get_summary_stats(error, type = "quantile") %>% 
    select(LM, method, "25%", "50%", "75%") %>% 
    mutate(event = "eskd")

# LM 0.5 and 3 did not pass Levene test, LM 2 and 3 failed normality of residuals
error_death <- error_death %>% 
    group_by(LM) %>% 
    nest() %>% 
    mutate(levene = map(data, \(data) levene_test(data, error ~ method))) %>%
    mutate(lin_reg = map(data, \(data) lm(error ~ method, data = data))) %>% 
    mutate(residual = map(lin_reg, \(lin_model) residuals(lin_model))) %>% 
    mutate(normality = map(residual, \(resid) shapiro_test(resid)))

# all ns
krus_death <- error_death %>% 
    select(LM, data) %>% 
    mutate(kruskal = map(data, \(data) kruskal_test(data, error ~ method)))

table_error_death <- error_death %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    group_by(LM, method) %>% 
    get_summary_stats(error, type = "quantile") %>% 
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
    ggplot(mapping = aes(x = LM, y = error, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "ESKD Error Rate (1-Harrell's C-index) with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "Error Rate (%)",
         fill = "Method") + 
    theme(plot.title = element_text(size = 17, hjust = 0.5), 
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0,27))

error_death %>% 
    select(LM, data) %>% 
    unnest(data) %>% 
    filter(method %in% sens_graph) %>% 
    mutate(LM = factor(LM)) %>% 
    ggplot(mapping = aes(x = LM, y = error, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "Death Error Rate (1-Harrell's C-index) with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "Error Rate (%)",
         fill = "Method") + 
    theme(plot.title = element_text(size = 17, hjust = 0.5), 
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0,27))

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

int_ext_error <- error_eskd %>% 
    select(LM, data) %>% 
    unnest(cols = c(data)) %>% 
    filter(method == "Top 5") %>% 
    bind_rows(ext_error)

# all LM pass Levene and Shapiro test
int_ext_error %>% 
    group_by(LM) %>% 
    nest() %>% 
    mutate(levene = map(data, \(data) levene_test(data, error ~ method))) %>% 
    mutate(shapiro = map(data, \(data) shapiro_test(data$error))) %>% 
    mutate(wilcox = map(data, \(data) wilcox_test(formula = data$error ~ data$method)))

# all ns
int_ext_error %>% 
    filter(LM == 3) %>% 
    wilcox_test(error ~ method)

int_ext_error %>% 
    group_by(LM, method) %>% 
    get_summary_stats(error, type = "quantile")

int_ext_error_eskd %>% 
    mutate(LM = factor(LM)) %>% 
    ggplot(mapping = aes(x = LM, y = error, fill = method)) + 
    geom_boxplot() + theme_bw() + 
    ggsci::scale_fill_lancet() + 
    labs(title = "ESKD Error Rate (1-Harrell's C-index) with Prediction Horizon 5 years",
         x = "Landmark Time (years)", y = "Error Rate (%)",
         fill = "Method") + 
    theme(plot.title = element_text(size = 17, hjust = 0.5), 
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 15)) + 
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0,27))

# end ---------------------------------------------------------------------