# setup -------------------------------------------------------------------

libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "randomForestSRC", "pec", 
               "gtsummary","purrr", "patchwork", "gt", "assertthat")
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

external_filename <- "experiment_test.rds" #change the file name to the file name that contains the patient for analysis

# dynamic_plot ------------------------------------------------------------

dt <- read_rds(external_filename)

id <- dt %>% 
    distinct(id) %>% 
    pull(id)

assert_that(length(id) == 1) #ensure only one patient per analysis

graph_dt <- tibble(dt = list(filter(dt, id == id)),
                         base_cov = list(base_cov_reduced),
                         longi_cov = list(longi_cov_top5),
                         analysis_name = c("top5"))

graph_dt <- graph_dt %>% 
    mutate(dt = map(dt, prepare_dt_for_dynamic_plot)) %>% 
    mutate(landmark = map(dt, extract_landmark_for_dynamic_plot, landmark)) %>% 
    unnest_longer(landmark) %>% 
    mutate(lmx_dt = pmap(list(dt, base_cov, longi_cov, landmark), 
                         extract_landmark_dt, predict_horizon)) %>% 
    mutate(pred_cif = pmap(list(lmx_dt, landmark, base_cov, longi_cov, 
                                analysis_name), predicted_cif, predict_horizon))

summarise_pts_info(graph_dt$dt, 
                   base_cov = graph_dt$base_cov,
                   longi_cov = graph_dt$longi_cov)

create_dynamic_plot(dt = graph_dt$dt,
                    base_cov = graph_dt$base_cov,
                    longi_cov = graph_dt$longi_cov,
                    landmark = graph_dt$landmark,
                    pred_cif = graph_dt$pred_cif)

# end ---------------------------------------------------------------------