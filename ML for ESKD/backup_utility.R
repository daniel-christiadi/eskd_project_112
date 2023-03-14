cross_val_performance_locf <- function(dt_path, longi_cov, base_cov){
    dt <- read_rds(dt_path)
    
    analysis_name <- tools::file_path_sans_ext(dt_path)
    analysis_name <- str_c(analysis_name, "locf", sep = "_")
    
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
        mutate(gn_fct = factor(gn_fct, levels = c("no", "yes")))
    
    set.seed(seed)
    
    folds <- dt %>% 
        select(id, sex:status_compete_num) %>% 
        distinct(id, .keep_all = TRUE) %>% 
        rsample::vfold_cv(v = cross_val, strata = "status_compete")
    
    cv_performance <- NULL
    
    for (fold in 1:cross_val){
        print(glue::glue("Current fold: {fold}"))
        
        train_id <- folds$splits[[fold]] %>% 
            rsample::analysis() %>% 
            pull(id)
        
        test_id <- folds$splits[[fold]] %>% 
            rsample::assessment() %>% 
            pull(id)
        
        train_dt <- dt %>% 
            filter(id %in% train_id)
        
        test_dt <- dt %>% 
            filter(id %in% test_id)
        
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
        
        model_performance$fold <- fold
        cv_performance <- rbind(cv_performance, model_performance)
    }
    cv_performance$analysis <- analysis_name
    cv_performance
}

cross_val_performance_lme <- function(dt_path, longi_cov, base_cov){
    dt <- read_rds(dt_path)
    
    analysis_name <- tools::file_path_sans_ext(dt_path)
    analysis_name <- str_c(analysis_name, "lme", sep = "_")
    
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
        data.table()
    
    set.seed(seed)
    
    folds <- dt %>% 
        select(id, sex:status_compete_num) %>% 
        distinct(id, .keep_all = TRUE) %>% 
        rsample::vfold_cv(v = cross_val, strata = "status_compete")
    
    cv_performance <- NULL
    
    for (fold in 1:cross_val){
        print(glue::glue("Current fold: {fold}"))
        
        train_id <- folds$splits[[fold]] %>% 
            rsample::analysis() %>% 
            pull(id)
        
        train_dt <- dt[id %in% train_id]
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
        
        test_id <- folds$splits[[fold]] %>% 
            rsample::assessment() %>% 
            pull(id)
        
        test_dt <- dt[id %in% test_id]
        
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
        
        model_vimp <- NULL
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
        
        model_performance$fold <- fold
        cv_performance <- rbind(cv_performance, model_performance)
    }
    cv_performance$analysis <- analysis_name 
    cv_performance
}

cross_val_performance_lmepoly <- function(dt_path, longi_cov, base_cov){
    dt <- read_rds(dt_path)
    
    analysis_name <- tools::file_path_sans_ext(dt_path)
    analysis_name <- str_c(analysis_name, "lmepoly", sep = "_")
    
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
        data.table()
    
    set.seed(seed)
    
    folds <- dt %>% 
        select(id, sex:status_compete_num) %>% 
        distinct(id, .keep_all = TRUE) %>% 
        rsample::vfold_cv(v = cross_val, strata = "status_compete")
    
    cv_performance <- NULL
    
    for (fold in 1:cross_val){
        print(glue::glue("Current fold: {fold}"))
        
        train_id <- folds$splits[[fold]] %>% 
            rsample::analysis() %>% 
            pull(id)
        
        train_dt <- dt[id %in% train_id]
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
        
        test_id <- folds$splits[[fold]] %>% 
            rsample::assessment() %>% 
            pull(id)
        
        test_dt <- dt[id %in% test_id]
        
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
        
        model_vimp <- NULL
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
        
        model_performance$fold <- fold
        cv_performance <- rbind(cv_performance, model_performance)
    }
    cv_performance$analysis <- analysis_name
    cv_performance
}