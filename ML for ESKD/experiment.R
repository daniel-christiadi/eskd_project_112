# setup -------------------------------------------------------------------

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

future_pwalk(list(experiment_dt$train_dt, experiment_dt$longi_cov, 
                  experiment_dt$base_cov, experiment_dt$id_name),
                 final_models, landmark, predict_horizon)

# write_rds(experiment_dt$test_dt[[1]], "experiment_test.rds")


# test cif ----------------------------------------------------------------

experiment_top5_model <- vector("list", length = length(landmark))

model_name <- "experiment_top5"

# i <- 1

for (i in seq_along(landmark)){
    file_name <- str_c(model_name, landmark[i], "locf_models.gz", 
                       sep = "_")
    # file_name <- str_c(c("top5", "model", landmark[i], "landmark_5.gz"), 
    #                    collapse = "_")
    # file_name <- str_c(c("top5", "model", landmark[i], "landmark_2.gz"), 
    #                    collapse = "_")
    print(file_name)
    experiment_top5_model[[i]] <- readRDS(file = file_name)
    # model_2y[[i]] <- readRDS(file = file_name)
} 
    
test_dt <- experiment_dt$test_dt[[1]]

test_dt <- test_dt %>% 
    select(id, !!base_cov_reduced, relyear, !!longi_cov_top5, 
           time_compete:status_compete_num) %>% 
    group_by(id) %>% 
    fill(!!longi_cov_top5) %>% 
    ungroup() %>% 
    as.data.table()

lmx_test_dt <- NULL

for (i in seq_along(landmark)) {
    temp <- cutLM(data = test_dt, outcome = list(time = "time_compete",
                                                 status = "status_compete_num"),
                  LM = landmark[i], horizon = landmark[i] + predict_horizon,
                  covs = list(fixed = base_cov_reduced, varying = longi_cov_top5),
                  format = "long", id = "id", rtime = "relyear", 
                  right = FALSE) 
    lmx_test_dt <- rbind(lmx_test_dt, temp)
}

lmx_test_dt <- lmx_test_dt %>% 
    arrange(id, LM)

lmx_test_dt <- lmx_test_dt %>% 
    select(-(time_compete:status_compete_num))

super_eskd_cif <- NULL
super_death_cif <- NULL

index <- 1
for (index in seq_along(landmark)){
    temp  <- lmx_test_dt %>% 
        filter(LM == landmark[index])
    
    # rsf_prediction <- predict.rfsrc(object = model_2y[[index]],
    #                                 newdata = temp, na.action = "na.impute")
    
    rsf_prediction <- predict.rfsrc(object = experiment_top5_model[[index]],
                                    newdata = temp, na.action = "na.impute")
    
    column_cif <- sindex(jump.times = rsf_prediction$time.interest,
                         eval.times = seq(landmark[index], 
                                          landmark[index] + predict_horizon, 
                                          0.25))
    
    column_cif[1] <- 1
    eskd_cif <- data.table(id = temp$id)
    
    eskd_cif[, c(as.character(seq(landmark[index], landmark[index] + predict_horizon, 0.25))) := as.data.table(rsf_prediction$cif[, column_cif, "CIF.1"])]
    
    eskd_cif <- melt(eskd_cif, id.vars = c("id"), 
                     variable.name = "relyear",
                     value.name = "cif")
    
    super_eskd_cif <- rbind(super_eskd_cif, eskd_cif)
    
    death_cif <- data.table(id = temp$id)
    
    death_cif[, c(as.character(seq(landmark[index], landmark[index] + predict_horizon, 0.25))) := as.data.table(rsf_prediction$cif[, column_cif, "CIF.2"])]
    
    death_cif <- melt(death_cif, id.vars = c("id"), 
                      variable.name = "relyear",
                      value.name = "cif")
    
    super_death_cif <- rbind(super_death_cif, death_cif)
}



    
