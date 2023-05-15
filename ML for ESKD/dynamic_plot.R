# setup -------------------------------------------------------------------

libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "randomForestSRC", "pec", 
               "gtsummary","purrr", "patchwork", "gt","boot")
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]
if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)
lib.ok <- sapply(libraries, require, character.only=TRUE)
cross_val <- 5
seed <- 7
predict_horizon <- 5
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)

longi_cov_top5 <- c("albumin", "bicarb", "chloride", "egfr", "haemoglobin")
base_cov_reduced <- c("age_init")

source("utility_functions.R")

external_filename <- "experiment_test.rds"

# dynamic_plot ------------------------------------------------------------

dt <- read_rds(external_filename)

censored_dt <- dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    filter(status_compete == "no")


random_test_id <- eskd_dt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    slice_sample(n = 1,) %>% 
    pull(id)

random_test_id <- 10719

random_test_dt <- tibble(dt = list(filter(dt, id == random_test_id)),
                         base_cov = list(base_cov_reduced),
                         longi_cov = list(longi_cov_top5),
                         analysis_name = c("experiment_top5"))

random_test_dt <- random_test_dt %>% 
    mutate(dt = map(dt, prepare_dt_for_dynamic_plot)) %>% 
    mutate(landmark = map(dt, extract_landmark_for_dynamic_plot, landmark)) %>% 
    unnest_longer(landmark) %>% 
    mutate(lmx_dt = pmap(list(dt, base_cov, longi_cov, landmark), 
                         extract_landmark_dt, predict_horizon)) %>% 
    mutate(pred_cif = pmap(list(lmx_dt, landmark, base_cov, longi_cov, 
                                analysis_name), predicted_cif, predict_horizon))

summarise_pts_info(random_test_dt$dt, 
                   base_cov = random_test_dt$base_cov,
                   longi_cov = random_test_dt$longi_cov)

create_dynamic_plot(dt = random_test_dt$dt,
                    base_cov = random_test_dt$base_cov,
                    longi_cov = random_test_dt$longi_cov,
                    landmark = random_test_dt$landmark,
                    pred_cif = random_test_dt$pred_cif)

# end ---------------------------------------------------------------------