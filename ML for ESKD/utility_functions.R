library(tidyverse)
library(data.table)
library(dynpred)
library(prodlim)
library(randomForestSRC)
library(pec)
library(rsample)
library(nlme)
library(mice)
library(timeROC)
library(boot)

recovery_check <- function(nested_df){
    # check whether eGFR < 15 is consistent with ESKD or severe AKI using modality status #
    modality <- nested_df$modality
    # modality <- modality[!duplicated(modality)]
    eskd <- c('Haemodiafiltration', 'Haemofiltration', 'HD', 'PD', 'Nocturnal HD')
    recovered <- c('CKD', 'CKD (4-5)', 'Outpatient')
    result <- list()
    
    for (index in seq_along(modality)){
        if ((modality[index] %in% eskd) & (modality[index+1] %in% recovered)){
            result <- append(result, 1)
        } else {
            result <- append(result, 0)
        }
    }
    as_vector(result)
}

extracting_freq <- function(meta_dt){
    # count the frequency of the clinicopathological test performed during observation #
    test_names <- meta_dt %>% 
        count(test) %>% 
        pull(test)
    
    id_dt <- meta_dt %>% 
        distinct(id) %>% 
        data.table()
    
    setkey(id_dt, id)
    
    dt <- meta_dt %>% 
        select(id, test, value) %>% 
        data.table()
    
    setkey(dt, id)
    
    for (test_name in test_names){
        m <- dt[test == test_name]
        n <- m[, .N, by = id]
        setnames(n, old = 'N', new = test_name)
        setkey(n, id)
        id_dt <- n[id_dt]
    }
    id_dt
}

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

prepare_dt_for_dynamic_plot <- function(dt){
    dt <- as.data.table(dt)
    
    dt <- dt %>% 
        pivot_wider(names_from = test, 
                    values_from = value)
    
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

create_folds <- function(dt, seed){
    # create individual fold for K-fold performance analysis #
    set.seed(seed)
    
    folds <- dt %>% 
        select(id, sex:status_compete_num) %>% 
        distinct(id, .keep_all = TRUE) %>% 
        vfold_cv(v = cross_val, strata = "status_compete")
    folds
}

extract_train_dt <- function(dt, folds_splits){
    temp_dt <- analysis(folds_splits)
    temp_id <- temp_dt$id
    
    train_dt <- dt %>% 
        filter(id %in% temp_id)
    
    train_dt
}

extract_test_dt <- function(dt, folds_splits){
    temp_dt <- assessment(folds_splits)
    temp_id <- temp_dt$id
    
    test_dt <- dt %>% 
        filter(id %in% temp_id)
    
    test_dt
}

extract_landmark_for_dynamic_plot <- function(prepared_dt, ...){
    max_test_time <- prepared_dt %>% 
        filter(relyear <= max(landmark)) %>% 
        slice_max(relyear) %>% 
        mutate(relyear = plyr::round_any(relyear, 0.5, f = ceiling)) %>% 
        pull(relyear)
    
    if(max_test_time >= max(landmark)){
        landmark <- seq(from = min(landmark), to = max(landmark), by = 0.5)
    } else {
        landmark <- seq(from = min(landmark), to = max_test_time + 0.5, 0.5)
    }
    landmark
}

extract_landmark_dt <- function(dt, base_cov, longi_cov, landmark, ...){
    # external or test data only #
    # extracting dataset per landmark #
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
    dt
}

performance_landmark_dt <- function(landmark_dt, landmark, base_cov, 
                                    longi_cov, analysis_models, ...){
    # external or test data only #
    # measuring error rate per landmark #
    landmark_dt <- landmark_dt %>% 
        drop_na()

    current_dir <- getwd()
    model_name <- str_c(analysis_models, landmark, "locf_models.gz", sep = "_")
    model_path <- file.path(current_dir, model_name)
    rsf_model <- readRDS(model_path)
    
    error <- boot::boot(data = landmark_dt, statistic = error_bootstrap,
                        R = boot_n, rsf_model = rsf_model)
    
    error_name <- str_c(analysis_models, landmark, "error.rds", sep = "_")
    error_path <- file.path(current_dir, error_name)
    write_rds(error, error_path)
    
    performance <- data.table(LM = landmark,
                              N = nrow(landmark_dt))
    performance
}

tvroc_landmark_dt <- function(landmark_dt, landmark, base_cov, 
                              longi_cov, analysis_models, ...){
    # external or test data only #
    # measuring error rate per landmark #
    landmark_dt <- landmark_dt %>% 
        drop_na()
    
    current_dir <- getwd()
    model_name <- str_c(analysis_models, landmark, "locf_models.gz", sep = "_")
    model_path <- file.path(current_dir, model_name)
    rsf_model <- readRDS(model_path)
    
    d_roc <- boot::boot(data = landmark_dt, statistic = tvroc_bootstrap,
                        R = boot_n, rsf_model = rsf_model, 
                        predict_horizon = predict_horizon)
    
    roc_name <- str_c(analysis_models, landmark, "roc.rds", sep = "_")
    roc_path <- file.path(current_dir, roc_name)
    write_rds(d_roc, roc_path)
    
    performance <- data.table(LM = landmark,
                              N = nrow(landmark_dt))
    performance
}

brier_landmark_dt <- function(landmark_dt, landmark, base_cov, 
                              longi_cov, analysis_models, ...){
    # external or test data only #
    # measuring integrated brier score per landmark #
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
    
    brier_name <- str_c(analysis_models, landmark, "brier.rds", sep = "_")
    brier_path <- file.path(current_dir, brier_name)
    write_rds(brier_eskd, brier_path)
    
    performance <- data.table(LM = landmark,
                              N = nrow(landmark_dt))
    performance
}


cross_val_performance_locf_map <- function(train_dt, test_dt, longi_cov, 
                                           base_cov, ...){
    # measuring error rate and integrated brier score for LOCF approach #
    long_train <- train_dt %>%
        select(id, relyear, !!longi_cov)
    
    ## imputing covariates at relyear == 0   
    impute_dt <- long_train %>% 
        filter(relyear == 0)
    
    temp_dt <- mice(impute_dt, m = 1, maxit = 50, method = "pmm", seed = seed)
    impute_dt <- complete(temp_dt)
    
    train_dt <- impute_dt %>% 
        bind_rows(long_train) %>% 
        arrange(id) %>% 
        distinct(id, relyear, .keep_all = TRUE) %>%
        group_by(id) %>%
        fill(!!longi_cov) %>% 
        ungroup() %>% 
        left_join(select(train_dt, id:status_compete_num),
                  by = c("id", "relyear")) %>% 
        arrange(id, relyear) %>% 
        as.data.table()
    
    test_dt <- test_dt %>% 
        group_by(id) %>% 
        fill(!!longi_cov) %>% 
        ungroup() %>% 
        as.data.table()
    
    train_super_dt <- NULL
    
    for (i in seq_along(landmark)) {
        temp <- cutLM(data = train_dt, outcome = list(time = "time_compete",
                                                      status = "status_compete_num"),
                      LM = landmark[i], horizon = landmark[i] + predict_horizon,
                      covs = list(fixed = base_cov, varying = longi_cov),
                      format = "long", id = "id", rtime = "relyear", 
                      right = FALSE) 
        train_super_dt <- rbind(train_super_dt, temp)
    }
    
    test_super_dt <- NULL
    
    for (i in seq_along(landmark)) {
        temp <- cutLM(data = test_dt, outcome = list(time = "time_compete",
                                                     status = "status_compete_num"),
                      LM = landmark[i], horizon = landmark[i] + predict_horizon,
                      covs = list(fixed = base_cov, varying = longi_cov),
                      format = "long", id = "id", rtime = "relyear", 
                      right = FALSE) 
        test_super_dt <- rbind(test_super_dt, temp)
    }
    
    model_performance <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Processing landmark: {lmx}"))
        
        train_lm <- train_super_dt %>% 
            filter(LM == lmx)
        
        test_lm <- test_super_dt %>% 
            filter(LM == lmx) %>% 
            drop_na()
        
        model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                          str_c(c(base_cov, longi_cov), collapse = " + "), 
                                          sep = "~"))
        
        tune_rsf <- tune.rfsrc(model_formula, data = train_lm, 
                               ntreeTry = 500, nodesizeTry = c(1:9, seq(10, 100, 5)))
        
        rsf_model <- rfsrc(formula = model_formula, data = train_lm, ntree = 1000, 
                           splitrule = "logrankCR", importance = TRUE, statistics = TRUE, 
                           mtry = tune_rsf$optimal[[2]], nodesize = tune_rsf$optimal[[1]])
        
        rsf_prediction <- predict.rfsrc(object = rsf_model,
                                        newdata = test_lm, na.action = "na.omit")
        
        brier_eskd <- pec(rsf_model, formula = model_formula, data = test_lm,
                          cause = 1)
        
        brier_death <- pec(rsf_model, formula = model_formula, data = test_lm,
                           cause = 2)
        
        error_temp <- unname(apply(rsf_prediction$err.rate, 2, mean, na.rm = TRUE))
        
        performance_temp <- data.table(LM = lmx,
                                       N = rsf_prediction$n,
                                       error_eskd = 100 * error_temp[1],
                                       error_death = 100 * error_temp[2],
                                       brier_eskd = crps(brier_eskd)[2],
                                       brier_death = crps(brier_death)[2])
        
        model_performance <- rbind(model_performance, performance_temp)
    }
    model_performance
}

cross_val_vimp_locf_map <- function(train_dt, test_dt, longi_cov, 
                                    base_cov, ...){
    # measuring VIMP for LOCF approach #
    long_train <- train_dt %>%
        select(id, relyear, !!longi_cov)
    
    ## imputing covariates at relyear == 0   
    impute_dt <- long_train %>% 
        filter(relyear == 0)
    
    temp_dt <- mice(impute_dt, m = 1, maxit = 50, method = "pmm", seed = seed)
    impute_dt <- complete(temp_dt)
    
    train_dt <- impute_dt %>% 
        bind_rows(long_train) %>% 
        arrange(id) %>% 
        distinct(id, relyear, .keep_all = TRUE) %>%
        group_by(id) %>%
        fill(!!longi_cov) %>% 
        ungroup() %>% 
        left_join(select(train_dt, id:status_compete_num),
                  by = c("id", "relyear")) %>% 
        arrange(id, relyear) %>% 
        as.data.table()
    
    test_dt <- test_dt %>% 
        group_by(id) %>% 
        fill(!!longi_cov) %>% 
        ungroup() %>% 
        as.data.table()
    
    train_super_dt <- NULL
    
    for (i in seq_along(landmark)) {
        temp <- cutLM(data = train_dt, outcome = list(time = "time_compete",
                                                      status = "status_compete_num"),
                      LM = landmark[i], horizon = landmark[i] + predict_horizon,
                      covs = list(fixed = base_cov, varying = longi_cov),
                      format = "long", id = "id", rtime = "relyear", 
                      right = FALSE) 
        train_super_dt <- rbind(train_super_dt, temp)
    }
    
    test_super_dt <- NULL
    
    for (i in seq_along(landmark)) {
        temp <- cutLM(data = test_dt, outcome = list(time = "time_compete",
                                                     status = "status_compete_num"),
                      LM = landmark[i], horizon = landmark[i] + predict_horizon,
                      covs = list(fixed = base_cov, varying = longi_cov),
                      format = "long", id = "id", rtime = "relyear", 
                      right = FALSE) 
        test_super_dt <- rbind(test_super_dt, temp)
    }
    
    model_vimp <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Processing landmark: {lmx}"))
        
        train_lm <- train_super_dt %>% 
            filter(LM == lmx)
        
        test_lm <- test_super_dt %>% 
            filter(LM == lmx) %>% 
            drop_na()
        
        model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                          str_c(c(base_cov, longi_cov), collapse = " + "), 
                                          sep = "~"))
        
        tune_rsf <- tune.rfsrc(model_formula, data = train_lm, 
                               ntreeTry = 500, nodesizeTry = c(1:9, seq(10, 100, 5)))
        
        rsf_model <- rfsrc(formula = model_formula, data = train_lm, ntree = 1000, 
                           splitrule = "logrankCR", importance = TRUE, statistics = TRUE, 
                           mtry = tune_rsf$optimal[[2]], nodesize = tune_rsf$optimal[[1]])
        
        vimp_temp <- data.table(LM = lmx,
                                event = c("eskd", "death"))
        
        vimp_temp <- cbind(vimp_temp, 100*t(rsf_model$importance))
        model_vimp <- rbind(model_vimp, vimp_temp)
    }
    model_vimp
}

create_lme_formula <- function(y_name, baseline_covariates){
    lme_formula <- as.formula(str_c({{y_name}},
                                    str_c(c("relyear", baseline_covariates), collapse = " + "),
                                    sep = "~"))
}

create_lme_model <- function(lme_formula, dataset){
    lme(lme_formula,
        random = ~ 1|id,
        data = dataset, na.action = na.exclude,
        control = list(msMaxIter = 1000, msMaxEval = 1000, maxIter = 1000, opt = "optim"))
}

create_lme_model_complex <- function(y_name, baseline_covariates, dataset){
    X_terms <- str_c(c("relyear", "I(relyear^2)", baseline_covariates), 
                     collapse = "+")
    blueprint <- str_c(c(y_name, X_terms), collapse = " ~ ")
    lme_formula <- as.formula(blueprint)
    model <- lme(lme_formula, 
                 random = ~ relyear|id,
                 data = dataset, na.action = na.exclude,
                 control = list(msMaxIter = 1000, msMaxEval = 1000, maxIter = 1000, opt = "optim"))
    
}

cross_val_performance_lme_map <- function(train_dt, test_dt, longi_cov, 
                                          base_cov, ...){
    # measuring error rate and integrated brier score for simple LME approach #
    model_list <- vector(mode = "list", length = length(longi_cov))
    
    for (index in seq_along(longi_cov)){
        model_formula <- create_lme_formula(longi_cov[index],
                                            baseline_covariates = base_cov)
        model_list[[index]] <- create_lme_model(model_formula, train_dt)
    }
    
    names(model_list) <- longi_cov
    train_super_dt <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Processing train dataset for landmark: {lmx}"))
        
        temp_dt <- train_dt %>% 
            distinct(id, .keep_all = TRUE) %>% 
            select(id, !!base_cov, time_compete:status_compete_num) %>% 
            mutate(relyear = lmx)
        
        for (model_index in seq_along(model_list)){
            response_name <- names(model_list)[[model_index]]
            b <- ranef(model_list[[model_index]])
            sigma <- model_list[[model_index]]$sigma
            betas <- fixef(model_list[[model_index]])
            
            terms_X <- as.formula(str_c(c("~relyear", base_cov),
                                        collapse = "+")) # for simple linear mixed model
            
            terms_Z <- formula(model_list[[model_index]]$modelStruct$reStruct[[1]])
            all_vars <- unique(c(all.vars(terms_X), all.vars(terms_Z)))
            dataset <- model_list[[model_index]]$data
            
            dataset <- dataset %>%
                distinct(id, .keep_all = TRUE) %>% 
                select(id, !!all_vars) %>% 
                mutate(relyear = lmx)
            
            X_matrix <- model.matrix(terms_X, data = dataset)
            Z_matrix <- model.matrix(terms_Z, data = dataset)
            level <- as.vector(c(X_matrix %*% betas) + rowSums(Z_matrix * b))
            
            temp_dt[[response_name]] <- level
        }
        
        cut_temp_dt <- cutLM(data = temp_dt, 
                             outcome = list(time = "time_compete",
                                            status = "status_compete_num"),
                             LM = lmx, 
                             horizon = lmx + predict_horizon,
                             covs = list(fixed = base_cov, 
                                         varying = longi_cov),
                             format = "long", id = "id", rtime = "relyear", 
                             right = FALSE)
        
        train_super_dt <- rbind(train_super_dt, cut_temp_dt)
    }
    
    train_super_dt <- train_super_dt %>% 
        select(-relyear) %>% 
        relocate(LM, .after = id) %>% 
        arrange(id, LM)
    
    test_pred_dt <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Processing test dataset for landmark: {lmx}"))
        
        test_temp_long_dt <- test_dt %>% 
            filter(relyear <= lmx)
        
        test_dt_y_pred <- NULL
        
        for (model_index in seq_along(model_list)){
            response_name <- names(model_list)[[model_index]]
            
            dataset <- test_temp_long_dt %>% 
                drop_na(!!response_name)
            
            dataset_outcome <- dataset %>% 
                distinct(id, .keep_all = TRUE) %>% 
                select(id, time_compete:status_compete_num, !!base_cov) %>% 
                mutate(relyear = lmx)
            
            sigma <- model_list[[model_index]]$sigma
            D <- getVarCov(model_list[[model_index]])
            attr(D, "group.levels") <- NULL
            
            dataset_X_split <- split(data.frame(
                model.matrix(terms_X, data = dataset)), dataset$id)
            
            dataset_Z_split <- split(data.frame(
                model.matrix(terms_Z, data = dataset)), dataset$id)
            
            dataset_Y_split <- split(dataset[[response_name]], dataset$id)
            
            dataset_V <- lapply(dataset_Z_split, function(x){
                as.matrix(x) %*% tcrossprod(D, as.matrix(x)) + diag(sigma^2, nrow(as.matrix(x)))
            })
            
            test_b <- lapply(seq_len(length(dataset_X_split)), 
                             function(i, x, y, z, v) { 
                                 tcrossprod(D, as.matrix(z[[i]])) %*%
                                     solve(v[[i]]) %*% (y[[i]] - as.matrix(x[[i]]) %*% fixef(model_list[[model_index]]))},
                             x = dataset_X_split, y = dataset_Y_split, z = dataset_Z_split, v = dataset_V)
            
            test_b <- matrix(unlist(test_b), ncol = dim(b)[2], byrow = TRUE)
            dataset_X_pred <- model.matrix(terms_X, data = dataset_outcome)
            dataset_Z_pred <- model.matrix(terms_Z, data = dataset_outcome)
            
            dataset_outcome$test <- response_name
            dataset_outcome$value <- as.vector(c(dataset_X_pred %*% fixef(model_list[[model_index]])) + 
                                                   rowSums(dataset_Z_pred * test_b))
            
            test_dt_y_pred <- rbind(test_dt_y_pred, dataset_outcome)
        }
        test_dt_y_pred <- test_dt_y_pred %>% 
            pivot_wider(names_from = test, values_from = value)
        
        test_pred_dt <- rbind(test_pred_dt, test_dt_y_pred)
    }
    
    test_super_dt <- NULL
    
    for(lmx in landmark){
        test_super_temp <- cutLM(data = test_pred_dt, outcome = list(time = "time_compete",
                                                                     status = "status_compete_num"),
                                 LM = lmx,
                                 horizon = lmx + predict_horizon,
                                 covs = list(fixed = base_cov,
                                             varying = longi_cov),
                                 format = "long", id = "id", rtime = "relyear", right = FALSE)
        test_super_dt <- rbind(test_super_dt, test_super_temp)
    }
    
    test_super_dt <- test_super_dt %>% 
        select(-relyear) %>% 
        relocate(LM, .after = id) %>% 
        arrange(id, LM)
    
    model_performance <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Training model for landmark: {lmx}"))
        train_lm <- train_super_dt %>% 
            filter(LM == lmx)
        
        test_lm <- test_super_dt %>% 
            filter(LM == lmx)
        
        model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                          str_c(c(base_cov, longi_cov), collapse = " + "), 
                                          sep = "~"))
        
        tune_rsf <- tune.rfsrc(model_formula, data = train_lm, 
                               ntreeTry = 500, nodesizeTry = c(1:9, seq(10, 100, 5)))
        
        rsf_model <- rfsrc(formula = model_formula, data = train_lm, ntree = 1000, 
                           splitrule = "logrankCR", importance = TRUE, statistics = TRUE, 
                           mtry = tune_rsf$optimal[[2]], nodesize = tune_rsf$optimal[[1]])
        
        rsf_prediction <- predict.rfsrc(object = rsf_model,
                                        newdata = test_lm, na.action = "na.omit")
        
        test_no_missing <- test_lm %>%
            drop_na()
        
        brier_eskd <- pec(rsf_model, formula = model_formula, data = test_no_missing,
                          cause = 1)
        
        brier_death <- pec(rsf_model, formula = model_formula, data = test_no_missing,
                           cause = 2)
        
        error_temp <- unname(apply(rsf_prediction$err.rate, 2, mean, na.rm = TRUE))
        
        performance_temp <- data.table(LM = lmx,
                                       N = rsf_prediction$n,
                                       error_eskd = 100 * error_temp[1],
                                       error_death = 100 * error_temp[2],
                                       brier_eskd = crps(brier_eskd)[2],
                                       brier_death = crps(brier_death)[2])
        
        model_performance <- rbind(model_performance, performance_temp)
    }
    model_performance
}

cross_val_performance_lmepoly_map <- function(train_dt, test_dt, longi_cov,
                                              base_cov, ...){
    # measuring error rate and integrated brier score for polynomial LME approach #
    model_list <- vector(mode = "list", length = length(longi_cov))
    
    for (index in seq_along(longi_cov)){
        model <- create_lme_model_complex(longi_cov[index], base_cov,
                                          dataset = train_dt)
        model_list[[index]] <- model
    }
    
    names(model_list) <- longi_cov
    train_super_dt <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Processing train dataset for landmark: {lmx}"))
        
        temp_dt <- train_dt %>% 
            distinct(id, .keep_all = TRUE) %>% 
            select(id, !!base_cov, time_compete:status_compete_num) %>% 
            mutate(relyear = lmx)
        
        for (model_index in seq_along(model_list)){
            response_name <- names(model_list)[[model_index]]
            b <- ranef(model_list[[model_index]])
            sigma <- model_list[[model_index]]$sigma
            betas <- fixef(model_list[[model_index]])
            
            terms_X <- as.formula(str_c(c("~relyear", "I(relyear^2)", base_cov),
                                        collapse = "+")) # for complex mixed model
            
            terms_Z <- formula(model_list[[model_index]]$modelStruct$reStruct[[1]])
            all_vars <- unique(c(all.vars(terms_X), all.vars(terms_Z)))
            dataset <- model_list[[model_index]]$data
            
            dataset <- dataset %>%
                distinct(id, .keep_all = TRUE) %>% 
                select(id, !!all_vars) %>% 
                mutate(relyear = lmx)
            
            X_matrix <- model.matrix(terms_X, data = dataset)
            Z_matrix <- model.matrix(terms_Z, data = dataset)
            level <- as.vector(c(X_matrix %*% betas) + rowSums(Z_matrix * b))
            
            temp_dt[[response_name]] <- level
        }
        
        cut_temp_dt <- cutLM(data = temp_dt, 
                             outcome = list(time = "time_compete",
                                            status = "status_compete_num"),
                             LM = lmx, 
                             horizon = lmx + predict_horizon,
                             covs = list(fixed = base_cov, 
                                         varying = longi_cov),
                             format = "long", id = "id", rtime = "relyear", 
                             right = FALSE)
        
        train_super_dt <- rbind(train_super_dt, cut_temp_dt)
    }
    
    train_super_dt <- train_super_dt %>% 
        select(-relyear) %>% 
        relocate(LM, .after = id) %>% 
        arrange(id, LM)
    
    test_pred_dt <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Processing test dataset for landmark: {lmx}"))
        
        test_temp_long_dt <- test_dt %>% 
            filter(relyear <= lmx)
        
        test_dt_y_pred <- NULL
        
        for (model_index in seq_along(model_list)){
            response_name <- names(model_list)[[model_index]]
            
            dataset <- test_temp_long_dt %>% 
                drop_na(!!response_name)
            
            dataset_outcome <- dataset %>% 
                distinct(id, .keep_all = TRUE) %>% 
                select(id, time_compete:status_compete_num, !!base_cov) %>% 
                mutate(relyear = lmx)
            
            sigma <- model_list[[model_index]]$sigma
            D <- getVarCov(model_list[[model_index]])
            attr(D, "group.levels") <- NULL
            
            dataset_X_split <- split(data.frame(
                model.matrix(terms_X, data = dataset)), dataset$id)
            
            dataset_Z_split <- split(data.frame(
                model.matrix(terms_Z, data = dataset)), dataset$id)
            
            dataset_Y_split <- split(dataset[[response_name]], dataset$id)
            
            dataset_V <- lapply(dataset_Z_split, function(x){
                as.matrix(x) %*% tcrossprod(D, as.matrix(x)) + diag(sigma^2, nrow(as.matrix(x)))
            })
            
            test_b <- lapply(seq_len(length(dataset_X_split)), 
                             function(i, x, y, z, v) { 
                                 tcrossprod(D, as.matrix(z[[i]])) %*%
                                     solve(v[[i]]) %*% (y[[i]] - as.matrix(x[[i]]) %*% fixef(model_list[[model_index]]))},
                             x = dataset_X_split, y = dataset_Y_split, z = dataset_Z_split, v = dataset_V)
            
            test_b <- matrix(unlist(test_b), ncol = dim(b)[2], byrow = TRUE)
            dataset_X_pred <- model.matrix(terms_X, data = dataset_outcome)
            dataset_Z_pred <- model.matrix(terms_Z, data = dataset_outcome)
            
            dataset_outcome$test <- response_name
            dataset_outcome$value <- as.vector(c(dataset_X_pred %*% fixef(model_list[[model_index]])) + 
                                                   rowSums(dataset_Z_pred * test_b))
            
            test_dt_y_pred <- rbind(test_dt_y_pred, dataset_outcome)
        }
        
        test_dt_y_pred <- test_dt_y_pred %>% 
            pivot_wider(names_from = test, values_from = value)
        
        test_pred_dt <- rbind(test_pred_dt, test_dt_y_pred)
    }
    
    test_super_dt <- NULL
    
    for(lmx in landmark){
        test_super_temp <- cutLM(data = test_pred_dt, outcome = list(time = "time_compete",
                                                                     status = "status_compete_num"),
                                 LM = lmx,
                                 horizon = lmx + predict_horizon,
                                 covs = list(fixed = base_cov,
                                             varying = longi_cov),
                                 format = "long", id = "id", rtime = "relyear", right = FALSE)
        test_super_dt <- rbind(test_super_dt, test_super_temp)
    }
    
    test_super_dt <- test_super_dt %>% 
        select(-relyear) %>% 
        relocate(LM, .after = id) %>% 
        arrange(id, LM)
    
    model_performance <- NULL
    
    for (lmx in landmark){
        print(glue::glue("Training model for landmark: {lmx}"))
        train_lm <- train_super_dt %>% 
            filter(LM == lmx)
        
        test_lm <- test_super_dt %>% 
            filter(LM == lmx)
        
        model_formula <- as.formula(str_c("Surv(time_compete, status_compete_num)",
                                          str_c(c(base_cov, longi_cov), collapse = " + "), 
                                          sep = "~"))
        
        tune_rsf <- tune.rfsrc(model_formula, data = train_lm, 
                               ntreeTry = 500, nodesizeTry = c(1:9, seq(10, 100, 5)))
        
        rsf_model <- rfsrc(formula = model_formula, data = train_lm, ntree = 1000, 
                           splitrule = "logrankCR", importance = TRUE, statistics = TRUE, 
                           mtry = tune_rsf$optimal[[2]], nodesize = tune_rsf$optimal[[1]])
        
        rsf_prediction <- predict.rfsrc(object = rsf_model,
                                        newdata = test_lm, na.action = "na.omit")
        
        test_no_missing <- test_lm %>%
            drop_na()
        
        brier_eskd <- pec(rsf_model, formula = model_formula, data = test_no_missing,
                          cause = 1)
        
        brier_death <- pec(rsf_model, formula = model_formula, data = test_no_missing,
                           cause = 2)
        
        error_temp <- unname(apply(rsf_prediction$err.rate, 2, mean, na.rm = TRUE))
        
        performance_temp <- data.table(LM = lmx,
                                       N = rsf_prediction$n,
                                       error_eskd = 100 * error_temp[1],
                                       error_death = 100 * error_temp[2],
                                       brier_eskd = crps(brier_eskd)[2],
                                       brier_death = crps(brier_death)[2])
        
        model_performance <- rbind(model_performance, performance_temp)
    }
    model_performance
}

confidence_interval <- function(vector, interval) {
    vec_sd <- sd(vector)
    n <- length(vector)
    vec_mean <- mean(vector)
    error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
    result <- c("mean" = vec_mean, 
                "lower" = vec_mean - error, 
                "upper" = vec_mean + error)
    final_result <- list(result)
}

create_error_graph <- function(dt, method, title){
    error <- dt %>% 
        filter(method == !!method) %>% 
        mutate(LM = factor(LM),
               event = factor(event, levels = c("eskd", "death"), 
                              labels = c("ESKD", "Death"))) %>% 
        ggplot(mapping = aes(x= event, y = mean, color = event)) + 
        geom_point() + 
        geom_pointrange(mapping = aes(ymin = lower, ymax = upper), stat = "identity") +
        facet_grid(~ LM) + 
        scale_color_manual(values = c("blue", "red")) + 
        labs(title = title,
             y = "Error Rate (%)") + 
        theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                           axis.text.x = element_text(size = 12),
                           axis.title.x = element_blank(),
                           axis.title.y = element_text(size = 15),
                           axis.text.y = element_text(size = 10),
                           legend.title = element_blank(),
                           legend.position = "none") + 
        scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0,27))
}

create_brier_graph <- function(dt, method, title){
    brier <- dt %>% 
        filter(method == !!method) %>% 
        mutate(LM = factor(LM),
               event = factor(event, levels = c("eskd", "death"), 
                              labels = c("ESKD", "Death"))) %>% 
        ggplot(mapping = aes(x= event, y = mean, color = event)) + 
        geom_point() + 
        geom_pointrange(mapping = aes(ymin = lower, ymax = upper), stat = "identity") +
        facet_grid(~ LM) + 
        scale_color_manual(values = c("blue", "red")) + 
        labs(title = title,
             y = "Brier Score") + 
        theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                           axis.title.x = element_blank(),
                           axis.text.x = element_text(size = 12),
                           axis.title.y = element_text(size = 15),
                           axis.text.y = element_text(size = 10),
                           legend.title = element_blank(),
                           legend.position = "none") + 
        scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06))
}

final_models <- function(dt, longi_cov, base_cov, model_name, ...){
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
    output_tibble <- tibble(
        landmark = landmark,
        nodesize = tune_matrix[,1],
        mtry = tune_matrix[,2],
    )    
    tune_tibble_name <- str_c(model_name, "locf_tune.rds", sep = "_")
    write_rds(output_tibble, file = tune_tibble_name)
}

brier_bootstrap <- function(rsf_model, data, index, formula, cause){
    d <- data[index, ]
    output <- modified_pec(object = rsf_model, 
                           data = d,
                           formula = formula,
                           cause = cause)
    ibs <- crps(output)[[2]]
    return(ibs)
}

error_bootstrap <- function(rsf_model, data, index){
    d <- data[index, ]
    rsf_prediction <- predict.rfsrc(object = rsf_model, 
                                    newdata = d, na.action = "na.omit")
    error_temp <- unname(apply(rsf_prediction$err.rate, 2, mean, na.rm = TRUE))
    return(error_temp)
}

tvroc_bootstrap <- function(rsf_model, data, index, predict_horizon){
    d <- data[index, ]
    rsf_prediction <- predict.rfsrc(object = rsf_model, 
                                    newdata = d, na.action = "na.omit")
    roc <- timeROC (T = d$time_compete,
                    delta = d$status_compete_num,
                    weighting = "aalen",
                    marker = rsf_prediction$predicted[,1], cause = 1,
                    times = predict_horizon)
    
    tvroc <- unname(roc$AUC)
    tvroc
}


summarise_pts_info <- function(dt, base_cov, longi_cov){
    base_cov <- base_cov %>% 
        pluck(1)
    
    longi_cov <- longi_cov %>% 
        pluck(1)
    
    dt %>%
        pluck(1) %>% 
        select(relyear, !!base_cov, time_compete:status_compete, all_of(longi_cov)) %>%
        rename("Intial_age" = age_init,
               "Test_time" = relyear,
               "Event_time" = time_compete,
               "Status" = status_compete) %>%
        gt() %>%
        tab_header(title = md("**Example Patient**")) %>%
        tab_spanner(label = "Pathology Test",
                    columns = !!longi_cov)
}

predicted_cif <- function(landmark_dt, landmark, base_cov, longi_cov,
                          analysis_models, ...){
    # measuring CIF for ESKD and Death event #
    model_name <- str_c(analysis_models, landmark, "locf_models.gz", sep = "_")
    rsf_model <- readRDS(model_name)
    
    rsf_prediction <- predict.rfsrc(object = rsf_model,
                                    newdata = landmark_dt, na.action = "na.impute")
    
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

create_dynamic_plot <- function(dt, base_cov, longi_cov, landmark, pred_cif){
    base_cov_temp <- base_cov %>% 
        pluck(1) 
    
    longi_cov_temp <- longi_cov %>% 
        pluck(1) 
    
    longi_plot <- dt %>% 
        pluck(1) %>% 
        filter(relyear <= max(landmark)) %>% 
        select(relyear, !!base_cov_temp, !!longi_cov_temp) %>% 
        # drop_na() %>% 
        pivot_longer(!!longi_cov_temp, 
                     names_to = "test",
                     values_to = "value") %>% 
        ggplot(aes(x = relyear, y = value)) + geom_point() + geom_line() + 
        theme_bw() + 
        facet_grid(test ~., switch = "y", scales = "free_y",
                   labeller = as_labeller(c(albumin = "sAlb", 
                                            bicarb = "Bicarbonate",
                                            chloride = "Chloride",
                                            egfr = "eGFR",
                                            haemoglobin = "Hb"))) + 
        coord_cartesian(xlim = c(-0.05,3.05)) +
        scale_x_continuous(expand = c(0, 0)) +
        labs(y = NULL, x = NULL) + 
        theme(strip.text.y.left = element_text(size = 15),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 10))
    
    surv_plot <- pred_cif %>% 
        bind_rows() %>% 
        group_by(pred_relyear) %>% 
        summarise(mean_eskd_cif = mean(eskd_cif),
                  mean_death_cif = mean(death_cif)) %>% 
        pivot_longer(cols = -pred_relyear,
                     names_to = "event",
                     values_to = "value") %>% 
        mutate(event = factor(event, 
                              levels = c("mean_eskd_cif", "mean_death_cif"),
                              labels = c("ESKD", "Death"))) %>%
        ggplot(aes(x = pred_relyear, y = value, colour = event)) + geom_line() + 
        scale_y_continuous(position = "right", labels = scales::label_percent(),
                           breaks = seq(0, 0.5, 0.1), limits = c(0, 0.55)) +
        scale_colour_manual(name = "Event CIF", values = c("blue", "red"),
                            labels = c("ESKD", "Death")) +
        labs(x = NULL, y = NULL) + 
        coord_cartesian(xlim = c(3, 8)) + theme_bw() + 
        scale_x_continuous(expand = c(0, 0), breaks = seq(4,8,1)) +
        theme(legend.position = c(0.1,0.9),
              axis.text.x = element_text(size = 12),
              axis.text.y.right = element_text(size = 15),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 12)) 
    
    relyear_lab <- "Time (years)"
    p_lab <- ggplot() + 
        annotate(geom = "text", x = 1, y = 1, label = relyear_lab, size = 5) + 
        coord_cartesian(clip = "off") + theme_void() 
    
    dynamic_plot <- (longi_plot | surv_plot) / p_lab + 
        plot_annotation(title = "Dynamic Plot Prediction") + 
        plot_layout(heights = c(1,0.01)) & 
        theme(plot.title = element_text(size = 20, hjust = 0.5))
    
    dynamic_plot
}

predicted_cif_model_loaded <- function(landmark_dt, landmark, base_cov, 
                                       longi_cov, ...){
    # need to all load model first as a list into model  
    index <- match(landmark, landmark)
    
    rsf_prediction <- predict.rfsrc(object = model[[index]],
                                    newdata = landmark_dt, na.action = "na.impute")
    
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

kfre8_model_2y <- function(age, sex, egfr, log_acr, albumin, phosphate, 
                          bicarb, calcium){
    1-0.9629^exp(-0.1848*(age/10-7.036)+0.2906*(sex-0.5642)-0.4156*(egfr/5-7.222)+0.3480*(log_acr-5.137)-0.3569*(albumin-3.997)+0.1582*(phosphate-3.916)-0.01199*(bicarb-25.57)-0.1581*(calcium-9.355))
}

kfre8_model_5y <- function(age, sex, egfr, log_acr, albumin, phosphate, 
                          bicarb, calcium){
    1-0.8636^exp(-0.1848*(age/10-7.036)+0.2906*(sex-0.5642)-0.4156*(egfr/5-7.222)+0.3480*(log_acr-5.137)-0.3569*(albumin-3.997)+0.1582*(phosphate-3.916)-0.01199*(bicarb-25.57)-0.1581*(calcium-9.355))
}