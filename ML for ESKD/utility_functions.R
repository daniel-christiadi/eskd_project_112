library(tidyverse)
library(data.table)
library(nlme)


create_lme_formula <- function(y_name, baseline_covariates){
    lme_formula <- as.formula(str_c({{y_name}},
                                    str_c(c("relyear",baseline_cov), collapse = " + "),
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

extract_cif_feature <- function(outcome_dt, train_test_dt, outcome){
    temp <- train_test_dt %>% 
        arrange(id, relyear) %>% 
        group_by(id, relyear) %>%
        summarise(mean_cif = 100*mean(cif), .groups = "keep") %>% 
        ungroup() %>% 
        group_by(id) %>% 
        mutate(relyear_before = lag(relyear),
               mean_cif_before = lag(mean_cif),
               slope = (mean_cif - mean_cif_before)/(relyear - relyear_before),
               slope_before = lag(slope),
               delta_slope = (slope - slope_before)) %>% 
        mutate(max_slope = max(slope, na.rm = TRUE),
               max_delta = max(delta_slope, na.rm = TRUE),
               max_cif = max(mean_cif, na.rm = TRUE)) %>% 
        slice_max(order_by = slope) %>% 
        distinct(id, .keep_all = TRUE) %>% 
        rename(cif_atmaxslope = mean_cif) %>% 
        select(id, cif_atmaxslope, max_slope:max_cif) %>% 
        ungroup()
    
    cif_atmaxdelta <- train_test_dt %>% 
        arrange(id, relyear) %>% 
        group_by(id, relyear) %>%
        summarise(mean_cif = 100*mean(cif), .groups = "keep") %>% 
        ungroup() %>% 
        group_by(id) %>% 
        mutate(relyear_before = lag(relyear),
               mean_cif_before = lag(mean_cif),
               slope = (mean_cif - mean_cif_before)/(relyear - relyear_before),
               slope_before = lag(slope),
               delta_slope = (slope - slope_before)) %>% 
        slice_max(order_by = delta_slope) %>% 
        distinct(id, .keep_all = TRUE) %>% 
        rename(cif_atmaxdelta = mean_cif) %>% 
        pull(cif_atmaxdelta)
    
    temp <- temp %>% 
        mutate(cif_atmaxdelta = cif_atmaxdelta) %>% 
        left_join(y = select(outcome_dt, id, time_compete, status_compete)) %>% 
        distinct(id, .keep_all = TRUE) %>% 
        mutate(event = if_else(
            condition = time_compete <= 8 & status_compete == !!outcome,
            true = "yes",
            false = "no"),
        event = as.factor(event))
        # select(-(time_compete:status_compete))
    
    temp
}

# mean_cif_youden <- cutpointr(train_eskd_cif, mean_cif, event, direction = ">=",
#                              pos_class = "yes", neg_class = "no",
#                              method = maximize_metric, metric = youden, 
#                              boot_runs = 100)
# 
# mean_cif_youden %>% add_metric(list(ppv,npv))
# 
# mean_cif_mcc <- cutpointr(train_eskd_cif, mean_cif, event, direction = ">=",
#                           pos_class = "yes", neg_class = "no",
#                           method = maximize_metric, metric = mcc,
#                           boot_runs = 100)
# mean_cif_mcc %>% add_metric(list(ppv, npv))
# 
# slope_youden <- cutpointr(train_eskd_cif, slope, event, direction = ">=",
#                           pos_class = "yes", neg_class = "no",
#                           method = maximize_metric, metric = youden, 
#                           boot_runs = 100, na.rm = TRUE)
# 
# slope_mcc <- cutpointr(train_eskd_cif, slope, event, direction = ">=",
#                        pos_class = "yes", neg_class = "no",
#                        method = maximize_metric, metric = mcc, 
#                        boot_runs = 100, na.rm = TRUE)


# rocr_prediction <- prediction(model_prediction, train_eskd_cif$event)
# perf_prec <- performance(rocr_prediction, "prec", x.measure = "cutoff")
# perf_rec <- performance(rocr_prediction, "rec", x.measure = "cutoff")
# perf <- performance(rocr_prediction, "prec", "rec")
# best.sum <- which.max(perf_prec@y.values[[1]]+perf_rec@y.values[[1]])
# close.intersection <- which.min(abs(perf_prec@y.values[[1]]-perf_rec@y.values[[1]]))
# 
# plotdata <- data.frame(x = perf@x.values[[1]],
#                        y = perf@y.values[[1]], 
#                        p = perf@alpha.values[[1]])
# 
# perf@alpha.values[[1]][best.sum]
# perf@alpha.values[[1]][close.intersection]
# 
# plot1 <- ggplot(data = plotdata) +
#     geom_path(aes(x = x, y = y)) + xlab(perf@x.name) + ylab(perf@y.name) + 
#     theme_bw() + labs(title = "ESKD Precision-Recall for Mean CIF, Slope, and Delta Slope Model 5y only")
# 
# plot1 + 
#     geom_point(data = plotdata[close.intersection, ], aes(x = x, y = y), 
#                col = "red") + 
#     annotate("text", 
#              x = plotdata[close.intersection, ]$x + 0.1,
#              y = plotdata[close.intersection, ]$y,
#              label = paste("p =", round(plotdata[close.intersection, ]$p, 3))) + 
#     geom_point(data = plotdata[best.sum, ], aes(x = x, y = y), 
#                col = "red") + 
#     annotate("text", 
#              x = plotdata[best.sum, ]$x + 0.1,
#              y = plotdata[best.sum, ]$y,
#              label = paste("p =", round(plotdata[best.sum, ]$p, 3)))
