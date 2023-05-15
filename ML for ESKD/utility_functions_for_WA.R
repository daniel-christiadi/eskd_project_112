library(tidyverse)
library(data.table)
library(dynpred)
library(prodlim)
library(randomForestSRC)
library(pec)
library(rsample)
library(nlme)
library(mice)
library(boot)

prepare_dt <- function(dt_path){
    if (str_detect(dt_path, pattern = regex(pattern = "RDS", 
                                            ignore_case = TRUE))){
        dt <- read_rds(dt_path)
    } else{
        dt <- fread(dt_path)
    }
    
    dt <- as.data.table(dt)
    
    dt <- dt %>% 
        pivot_wider(names_from = test, 
                    values_from = value, 
                    values_fn = mean)
    
    dt <- dt %>% 
        mutate(status_compete = replace(status_compete, 
                                        status_compete == "no", "censored")) %>% 
        mutate(status_compete_num = case_when(
            status_compete == "censored" ~ 0,
            status_compete == "eskd" ~ 1,
            .default = 2
        )) %>% 
        relocate(status_compete_num, .after = status_compete) %>% 
        mutate(sex = factor(sex, levels = c("Female", "Male")),
               gn_fct = if_else(is.na(gn_cat), "no", "yes")) %>% 
        relocate(gn_fct, .after = gn_cat) %>% 
        mutate(gn_fct = factor(gn_fct, levels = c("no", "yes"))) %>% 
        as.data.table()
    
    dt
}

extract_landmark_dt <- function(dt, base_cov, longi_cov, landmark, ...){
    ## external or test data only
    dt <- dt %>% 
        group_by(id) %>% 
        fill(!!longi_cov) %>% 
        ungroup() %>% 
        as.data.table()
    
    dt <- cutLM(data = dt, outcome = list(time = "time_compete",
                                          status = "status_compete_num"),
                LM = landmark, horizon = landmark + predict_horizon,
                covs = list(fixed = base_cov, varying = longi_cov),
                format = "long", id = "id", rtime = "relyear", 
                right = FALSE) 
    
    dummy_dt <- dt %>% 
        filter(status_compete_num == 1) %>% 
        head(n = 1)
    
    dummy_dt$status_compete_num <- 2
    
    dt <- dt %>% 
        add_row(dummy_dt, .before = 2)
    
    dt
}

error_bootstrap <- function(rsf_model, data, index){
    d <- data[index, ]
    rsf_prediction <- predict.rfsrc(object = rsf_model, 
                                    newdata = d, na.action = "na.omit")
    error_temp <- unname(apply(rsf_prediction$err.rate, 2, mean, na.rm = TRUE))
    return(error_temp)
}


performance_landmark_dt <- function(landmark_dt, landmark, base_cov, 
                                    longi_cov, analysis_models, ...){
    ## external or test data only
    landmark_dt <- landmark_dt %>% 
        drop_na()
    
    current_dir <- getwd()
    model_name <- str_c(analysis_models, landmark, "locf_models.gz", sep = "_")
    model_path <- file.path(current_dir, model_name)
    rsf_model <- readRDS(model_path)
    
    error <- boot::boot(data = landmark_dt, statistic = error_bootstrap,
                        R = boot_n, rsf_model = rsf_model)
    
    performance <- data.table(LM = landmark,
                              error_eskd = 100 * error$t[,1],
                              error_death = 100 * error$t[,2])
    performance
}

brier_bootstrap <- function(rsf_model, data, index, formula, cause){
    d <- data[index, ]
    output <- pec(object = rsf_model, 
                  data = d,
                  formula = formula,
                  cause = cause)
    ibs <- crps(output)[[2]]
    return(ibs)
}

brier_landmark_dt <- function(landmark_dt, landmark, base_cov, 
                              longi_cov, analysis_models, ...){
    ## external or test data only
    landmark_dt <- landmark_dt %>% 
        drop_na()
    
    model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                      str_c(c(base_cov, longi_cov), collapse = " + "), 
                                      sep = "~"))
    current_dir <- getwd()
    model_name <- str_c(analysis_models, landmark, "locf_models.gz", sep = "_")
    model_path <- file.path(current_dir, model_name)
    rsf_model <- readRDS(model_path)
    
    brier_eskd <- boot::boot(data = landmark_dt, statistic = brier_bootstrap,
                             R = boot_n, rsf_model = rsf_model, cause = 1, 
                             formula = model_formula)
    
    performance <- data.table(brier = brier_eskd$t[,1])
    performance
}