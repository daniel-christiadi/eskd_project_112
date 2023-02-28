libraries <- c("tidyverse", "data.table", "skimr", "nlme", "rsample", 
               "dynpred", "prodlim", "mice", "randomForestSRC", "pec", 
               "gtsummary", "patchwork")
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]
if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)
lib.ok <- sapply(libraries, require, character.only=TRUE)

cross_val <- 5
seed <- 7
predict_horizon <- 5
landmark <- c(0.5, 1, 1.5, 2, 2.5, 3)

longi_cov <- c("albumin", "alkphos", "bicarb", "calcium", "chloride",
                "egfr", "glucose", "haemoglobin", "phosphate", "platelet",
                "potassium", "sodium", "wcc")

base_cov <- c("sex", "age_init", "gn_fct")
approach <- list("main", "sens")

# k-fold ------------------------------------------------------------------

cv_vimp <- NULL
cv_performance <- NULL

for (method in approach){
    print(glue::glue("Analysing: {method}"))
    
    dataset_name <- str_c(method, "dense3_dt.rds", sep = "_")
    
    dt <- read_rds(dataset_name)
    
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
        
        model_vimp <- NULL
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
            
            vimp_temp <- data.table(LM = lmx,
                                    event = c("eskd", "death"))
            
            vimp_temp <- cbind(vimp_temp, 100*t(rsf_model$importance))
            model_vimp <- rbind(model_vimp, vimp_temp)
            
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
        
        model_vimp$fold <- fold
        model_vimp$method <- method
        cv_vimp <- rbind(cv_vimp, model_vimp)
        
        model_performance$fold <- fold
        model_performance$method <- method
        cv_performance <- rbind(cv_performance, model_performance)
    }
}

write_rds(cv_performance, "main_vs_sens_performance.rds")
write_rds(cv_vimp, "main_vs_sens_vimp.rds")

# visualisation -----------------------------------------------------------

error1 <- cv_performance %>% 
    select(-(brier_eskd:fold)) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    filter(method == "main") %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>% 
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Original Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30)) 

error2 <- cv_performance %>% 
    select(-(brier_eskd:fold)) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>%
    filter(method == "sens") %>% 
    pivot_longer(cols = c("error_eskd", "error_death"),
                 names_to = "error") %>% 
    mutate(error = factor(error, levels = c("error_eskd", "error_death"),
                          labels = c("ESKD", "Death"))) %>%  
    ggplot(aes(x = error, y = value, fill = error)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Sensitivity Error Rate (1 - Harrell's C-index) per Landmark with LOCF Prediction Horizon 5 years", 
         y = "Error Rate (%)") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") + 
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0,30))

error1 + error2 + plot_annotation(
    title = "Comparison Original vs Sensitivity Analysis Approach"
)

brier1 <- cv_performance %>% 
    select(-(error_eskd:error_death), -fold) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    filter(method == "main") %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Original Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

ggsave("locf_brier5.png", plot1, width = 12, height = 7, scale = 0.5)


brier2 <- cv_performance %>% 
    select(-(error_eskd:error_death), -fold) %>% 
    arrange(LM) %>% 
    mutate(LM = as_factor(LM)) %>% 
    filter(method == "sens") %>% 
    pivot_longer(cols = c("brier_eskd", "brier_death"),
                 names_to = "brier") %>% 
    mutate(brier = factor(brier, levels = c("brier_eskd", "brier_death"),
                          labels = c("ESKD", "Death"))) %>%    
    ggplot(aes(x = brier, y = value, fill = brier)) + geom_boxplot() +
    facet_grid(~ LM) + scale_fill_manual(values = c("blue", "red")) + 
    labs(title = "Sensitivity Integrated Brier Score per Landmark with LOCF Prediction Horizon 5 years",
         y = "Brier Score") + 
    theme_bw() + theme(axis.title.x = element_blank(), 
                       legend.title = element_blank(),
                       legend.position = "none") +
    scale_y_continuous(breaks = seq(0, 0.06, 0.005), limits = c(0, 0.06)) #universal setting

brier1 + brier2 + plot_annotation(
    title = "Comparison Original vs Sensitivity Analysis Approach"
)

