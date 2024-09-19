libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "randomForestSRC", "pec", 
               "gtsummary", "future", "furrr", "patchwork")
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

longi_cov_top10 <- c("albumin", "alkphos", "bicarb", "chloride", "egfr",
                     "haemoglobin", "platelet",
                     "potassium", "sodium", "wcc")

longi_cov_top5 <- c("albumin", "bicarb", "chloride", "egfr", "haemoglobin")
base_cov_full <- c("sex", "age_init", "gn_fct")
base_cov_reduced <- c("age_init")

source("utility_functions.R")

# dense3 locf exploration -------------------------------------------------

dense3_dt <- tibble(
    dt_path = "main_dense3_dt.rds",
    base_cov = list(base_cov_reduced),
    longi_cov = list(longi_cov_top5)
)

dense3_dt <- dense3_dt %>% 
    mutate(dt = map(dt_path, prepare_dt),
           folds = map(dt, create_folds, seed)) %>% 
    unnest_longer(folds)

dense3_dt <- dense3_dt %>% 
    select(folds) %>% 
    map(~.x) %>% 
    bind_rows() %>% 
    bind_cols(dense3_dt) %>% 
    select(id, base_cov, longi_cov, dt, splits) %>% 
    mutate(id_name = rep(c("top5"), each = 5)) %>% 
    mutate(id = str_to_lower(id))

dense3_dt <- dense3_dt %>% 
    mutate(train_dt = pmap(list(dt, splits), extract_train_dt),
           test_dt = pmap(list(dt, splits), extract_test_dt)) %>% 
    select(id_name, id, base_cov, longi_cov, train_dt, test_dt)

dense3_dt <- dense3_dt %>% 
    mutate(lmx = list(landmark)) %>% 
    select(id_name:longi_cov, lmx, test_dt) %>% 
    unnest(cols = lmx) %>% 
    mutate(lmx_dt = pmap(list(test_dt, base_cov, longi_cov, lmx), extract_landmark_dt, predict_horizon))
    
dense3_dt <- dense3_dt %>% 
    mutate(analysis_models = "top5") %>% 
    mutate(auc_eskd = pmap(list(lmx_dt, lmx, base_cov, longi_cov, id, analysis_models), tvroc_landmark_dt_eskd, predict_horizon),
           auc_death = pmap(list(lmx_dt, lmx, base_cov, longi_cov, id, analysis_models), tvroc_landmark_dt_death, predict_horizon))

# write_rds(dense3_dt, "internal_roc.rds")


tvroc_landmark_dt_eskd <- function(landmark_dt, landmark, base_cov, 
                              longi_cov, id, analysis_models, ...){
    # external or test data only #
    # measuring error rate per landmark #
    landmark_dt <- landmark_dt %>% 
        drop_na()
    
    current_dir <- getwd()
    model_name <- str_c("exper", analysis_models, id, landmark, "locf_models.rds", sep = "_")
    model_path <- file.path(current_dir, model_name)
    rsf_model <- readRDS(model_path)
    
    rsf_prediction <- predict.rfsrc(object = rsf_model, 
                                    newdata = landmark_dt, na.action = "na.omit")
    
    roc <- timeROC (T = landmark_dt$time_compete,
                    delta = landmark_dt$status_compete_num,
                    weighting = "aalen",
                    marker = rsf_prediction$predicted[,1], cause = 1,
                    times = predict_horizon)
    
    tvroc <- unname(roc$AUC_2)
}

tvroc_landmark_dt_death <- function(landmark_dt, landmark, base_cov, 
                              longi_cov, id, analysis_models, ...){
    # external or test data only #
    # measuring error rate per landmark #
    landmark_dt <- landmark_dt %>% 
        drop_na()
    
    current_dir <- getwd()
    model_name <- str_c("exper", analysis_models, id, landmark, "locf_models.rds", sep = "_")
    model_path <- file.path(current_dir, model_name)
    rsf_model <- readRDS(model_path)
    
    rsf_prediction <- predict.rfsrc(object = rsf_model, 
                                    newdata = landmark_dt, na.action = "na.omit")
    
    roc <- timeROC (T = landmark_dt$time_compete,
                    delta = landmark_dt$status_compete_num,
                    weighting = "aalen",
                    marker = rsf_prediction$predicted[,2], cause = 2,
                    times = predict_horizon)
    
    tvroc <- unname(roc$AUC_2)
}

auc_eskd <- dense3_dt %>% 
    select(id, lmx, auc_eskd) %>% 
    unnest_wider(col = auc_eskd, names_sep = "t0") %>% 
    rename("auc_eskd" = "auc_eskdt02") %>% 
    select(!auc_eskdt01)

roc <- dense3_dt %>% 
    select(id, auc_death) %>% 
    unnest_wider(col = auc_death, names_sep = "t0") %>% 
    rename("auc_death" = "auc_deatht02") %>% 
    select(auc_death) %>% 
    bind_cols(auc_eskd) %>% 
    relocate(auc_death, .after = auc_eskd)


roc <- roc %>% 
    pivot_longer(cols = auc_eskd:auc_death,
                 names_to = "event",
                 values_to = "auc") %>% 
    separate(col = event,
             into = c("nauc", "event"), sep = "_") %>% 
    select(!nauc)

roc %>% 
    mutate(LM = factor(lmx),
           event = factor(event, levels = c("eskd", "death"),
                          labels = c("ESKD", "Death"))) %>% 
    ggplot(mapping = aes(x = LM, y = auc, color = event, fill = event)) + geom_boxplot() +
    theme_bw() + 
    ggsci::scale_fill_lancet() +
    scale_color_manual(values = c("blue", "red")) +
    coord_cartesian(ylim = c(0.79,0.9)) +
    scale_y_continuous(labels = scales::label_percent(), breaks = c(0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9)) +
    labs(title = "Time-dependent ROC curve Internal Cohort",
         y = "Area Under the Curve") +
    theme(plot.title = element_text(size = 17, hjust = 0.5), 
          axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          legend.title = element_blank(),
          legend.text = element_text(size = 15))

roc %>% 
    group_by(lmx, event) %>% 
    summarise(median_auc = median(auc),
              iqr_25 = quantile(auc, 0.25),
              iqr_75 = quantile(auc, 0.75)) %>% 
    ungroup() %>% 
    group_by(event) %>% 
    summarise(median = median(median_auc))







