# setup -------------------------------------------------------------------

libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "randomForestSRC", "pec", 
               "gtsummary", "future", "furrr", "patchwork")
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]

if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)

lib.ok <- sapply(libraries, require, character.only=TRUE)
predict_horizon <- 5
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)

longi_cov_full <- c("albumin", "alkphos", "bicarb", "calcium", "chloride",
                    "egfr", "glucose", "haemoglobin", "phosphate", "platelet",
                    "potassium", "sodium", "wcc")
base_cov_full <- c("sex", "age_init", "gn_fct")

longi_cov_top5 <- c("albumin", "bicarb", "chloride", "egfr", "haemoglobin")
base_cov_reduced <- c("age_init")

external_data_file_name <- "experiment_test.rds" #please input file name here!
analysis_models <- c("full", "top5")
source("utility_functions_for_WA.R")

# performance -------------------------------------------------------------

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

# debug(performance_landmark_dt)
# performance_landmark_dt(landmark_dt = external_dt$lmx_dt[[1]],
#                         base_cov = external_dt$base_cov[[1]],
#                         longi_cov = external_dt$longi_cov[[1]],
#                         landmark = external_dt$landmark[1],
#                         analysis_models = external_dt$analysis_models[1])
# 
landmark_dt <- external_dt$lmx_dt[[1]]
base_cov <- external_dt$base_cov[[1]]
longi_cov <- external_dt$longi_cov[[1]]
landmark <- external_dt$landmark[1]
analysis_models <- external_dt$analysis_models[1]


landmark_dt <- landmark_dt %>%
    drop_na()

landmark_dt <- landmark_dt %>%
    filter(status_compete_num != 2)

model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                  str_c(c(base_cov, longi_cov), collapse = " + "),
                                  sep = "~"))
current_dir <- getwd()
model_name <- str_c(analysis_models, landmark, "locf_models.gz", sep = "_")
model_path <- file.path(current_dir, model_name)
rsf_model <- readRDS(model_path)

debug(pec)
predictSurvProb(rsf_model, landmark_dt, times = landmark + predict_horizon, na.action = "na.omit")

rsf_prediction <- predict(object = rsf_model, outcome = "test",
                          newdata = landmark_dt, na.action = "na.omit",)

brier_eskd <- pec(rsf_model, formula = model_formula, data = landmark_dt)

brier_death <- pec(rsf_model, formula = model_formula, data = landmark_dt,
                   cause = 2)
error <- unname(apply(rsf_prediction$err.rate, 2, mean, na.rm = TRUE))

performance <- data.table(LM = landmark,
                          N = rsf_prediction$n,
                          error_eskd = 100 * error[1],
                          error_death = 100 * error[2])
performance

external_dt <- external_dt %>% 
    mutate(performance = pmap(list(lmx_dt, landmark, base_cov, longi_cov, analysis_models), 
                              performance_landmark_dt))

external_dt <- external_dt %>% 
    select(landmark, analysis_models, performance)

write_rds(external_dt, "external_performance.rds")

# end ---------------------------------------------------------------------