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


