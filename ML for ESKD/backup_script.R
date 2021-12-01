

# Checking recovery post AKI
hd_recovery <- function(vect){
    status <- vector(mode = 'logical', length = length(vect))
    category <- c('Haemodiafiltration', 'Haemofiltration', 'HD', 'PD', 'Nocturnal HD')
    recovered <- c('CKD', 'CKD (4-5)', 'Other/None', 'Outpatient')
    for (index in seq_along(vect)){
        if ((vect[index] %in% category) & (vect[index + 1] %in% recovered)){
        }
    }
    status
}


# Baseline kidney function 
baseline<-function(list_column){
    dt <- as.data.table(list_column)
    baseline_vect = vector(mode = 'numeric', length = nrow(dt))
    for (index in 1:nrow(list_column)){
        row <- dt[index,]
        r = row$intense_code
        
        if (row$intense == 'yes'){
            d = row$date
            prev.year <- dt[(date >= d - (3600*24*365.25)) & (intense_code == (r - 1))]
            baseline_vect[index] = median(prev.year$value)
        }
        else {
            baseline_vect[index] <- NA
        }
    }
    baseline_vect
}

# Engineering AKI
dt <- melted[test == 'creatinine']
setorder(dt, id, date)
dt <- dt %>% 
    group_by(id) %>% 
    mutate(lag1 = lag(date, 1L),
           lag2 = lag(date, 2L), 
           lead1 = lead(date, 1L))
dt <- dt %>% 
    group_by(id) %>% 
    mutate(intense = case_when(
        difftime(date, lag1, units = 'days') <= 5 ~ 'yes',
        difftime(date, lag2, units = 'days') <= 5 ~ 'yes',
        difftime(date, lead1, units = 'days') >= -5 & difftime(date, lead1, units = 'days') <= 0 ~ 'yes',
        TRUE ~ 'no'))

# identify increase frequency of measuring kidney function and calculating baseline kidney function
dt <- dt %>% 
    group_by(id) %>% 
    mutate(intense_code = rleid(intense))
dt <- as.data.table(dt)
dt[, baseline_function := baseline(.SD), by = id]
dt[, aki_stage := fcase(
    value / baseline_function >= 3, 3,
    value / baseline_function >= 2, 2,
    value / baseline_function >= 1.5, 1
)]

setkey(dt, id, date)

# incorporating egfr, aki and modality
setkey(egfr, id)
dt <- unique(egfr[dt, allow.cartesian=TRUE])
dt[, c('test', 'i.test', 'lag1', 'lag2', 'lead1') := NULL]
setnames(dt, old = 'value', new = 'egfr')
setnames(dt, old = 'i.value', new = 'creat')
dt[, date2 := date]
setcolorder(dt, c('id', 'date', 'date2', 'egfr', 'creat', 'baseline_function', 'intense', 'intense_code', 'aki_stage'))
setnames(dt, old = 'date2', new = 'end.date')
setkey(dt, id, date, end.date)

aki <- modality[dt, roll = TRUE]
aki[, aki_stage := fifelse(
    modality %in% c('Haemodiafiltration', 'Haemofiltration', 'HD', 'PD', 'Nocturnal HD'), 0, aki_stage
)]

aki[, recovery := hd_recovery(modality), by = id]
setkey(aki, id, date)
aki <- aki[(aki_stage > 0) & (modality != 'Transplant')]
aki <- aki[, .SD[which.max(aki_stage)], 
           by = c('id', 'intense_code')][, .(id, intense_code, modality, aki_stage, i.end.date)]
aki[, n := .N, by = c('id', 'intense_code')]
aki[, modality := NULL]
setnames(aki, old = 'i.end.date', new = 'date')
setcolorder(aki, neworder = c('id', 'date', 'intense_code', 'aki_stage', 'n'))

# AKI plot
aki_plot <- aki %>% 
    ggplot(aes(x = aki_stage)) + geom_bar() + xlab(NULL) + 
    scale_x_discrete(breaks = c('1', '2', '3'), 
                     labels = c('AKI Stage I', 'AKI Stage II', 'AKI Stage III')) +
    labs(title = 'AKI Incidence by Stage') 
aki_plot

# Relative year for AKI
aki$relyear <- (aki$date - outcomes_dt[.(aki$id), date_egfr_init])/365.25
setattr(aki$relyear, 'units', 'years')

aki <- aki %>% 
    filter(relyear > 0) %>% drop_na()

setkey(aki, id)
# fwrite(aki, 'aki_data.gz')

# As the dataset containing the clinico-pathological of ts_1y is too big for 
# tsfresh-featurisation, the dataset is divided into 2 files using disk.frame and fst
library(disk.frame) 
library(fst)

shard(ts_1y, 
      shardby = "id",
      nchunks = 2,outdir = file.path(getwd(),"temp"))

ts_1ya <- read_fst("temp/1.fst")
fwrite(ts_1ya, file = "ts_1ya.gz")

ts_1yb <- read_fst("temp/2.fst")
fwrite(ts_1yb, file = "ts_1yb.gz")

# The tsfresh-featurisation is performed in python (see ts_featurisation_script.ipynb)
# After tsfresh-featurisation, the output dataset (ts_1ya_feat.gz and ts_1yb_feat.gz) 
# are combined into ts_1y_feat.gz
ts_1ya_feat <- fread("ts_1ya_feat.gz")
ts_1yb_feat <- fread("ts_1yb_feat.gz")
ts_1y_feat <- rbind(ts_1ya_feat, ts_1yb_feat)
# fwrite(ts_1y_feat, file = "ts_1y_feat.gz")

# As the tsfresh-featurised data contains sum and length values which indicate 
# the frequency of test, we need to remove the features. Additionally, as eGFR, 
# which is used to define ESKD, is derived from creatinine, creatinine 
# featurisation is removed.
ts_1y_feat <- fread("ts_1y_feat.gz")
ts_2y_feat <- fread("ts_2y_feat.gz")
ts_2y_feat <- ts_2y_feat %>% 
    select(!ends_with("sum_values") & !ends_with("__length")) %>% 
    select(!contains("creatinine"))

ts_1y_feat <- ts_1y_feat %>% 
    select(!ends_with("sum_values") & !ends_with("__length")) %>% 
    select(!contains("creatinine"))

# As the tsfresh-featurised data contains infinite and nan values, 
# we convert these values to NA
for (j in seq_along(ts_1y_feat)){
    data.table::set(ts_1y_feat,i = which(is.infinite(ts_1y_feat[[j]])|is.nan(ts_1y_feat[[j]])), j = j, NA)
}

for (j in seq_along(ts_2y_feat)){
    data.table::set(ts_2y_feat,i = which(is.infinite(ts_2y_feat[[j]])|is.nan(ts_2y_feat[[j]])), j = j, NA)
}

setnames(ts_1y_feat, old = "V1", new = "id")
setnames(ts_2y_feat, old = "V1", new = "id")

# Creation of train-ready dataset
dt <- outcomes_dt %>% 
    select(-birth.year, -value_egfr_last, -value_eskd) %>% 
    mutate(egfr_y = (date_egfr_last - date_egfr_init)/31557600) %>% 
    mutate(eskd_y = (date_eskd - date_egfr_init)/31557600) %>% 
    mutate(death_y = (date.death - date_egfr_init)/365.25)

setattr(dt$egfr_y, "units", "years")
setattr(dt$eskd_y, "units", "years")
setattr(dt$death_y, "units", "years")

dt <- dt %>% 
    mutate(duration = as.numeric(case_when(
        !is.na(eskd_y) ~ eskd_y,
        TRUE ~ egfr_y))) %>%
    mutate(status = if_else(condition = duration == egfr_y,
                            true = "censored",
                            false = "eskd",
                            missing = "missing")) %>% 
    mutate(duration_comp = as.numeric(case_when(
        status == "eskd" ~ eskd_y,
        status == "censored" & !is.na(death_y) ~ death_y,
        TRUE ~ egfr_y))) %>% 
    mutate(status_comp = case_when(
        duration_comp == eskd_y ~ "eskd",
        duration_comp == death_y ~ "death",
        TRUE ~ "censored"
    )) %>% 
    select(-(date_egfr_init:death_y))

dt <- dt %>%
    mutate(htn = factor(htn, levels = c("no", "htn")),
           dkd = factor(dkd, levels = c("no", "dkd")),
           gn = factor(gn),
           gender = factor(gender),
           status = factor(status, levels = c("censored", "eskd")),
           status_comp = factor(status_comp, levels = c("censored", "eskd", 
                                                        "death")))

# Combining the outcome dataset with the tsfresh-featurised clinicopathological data. 
# For tsfresh-featurised using 1 year data - the duration and duration_comp will 
# need to be substracted by 1. Similarly, tsfresh-featurised using 2 years of 
# data, the duration and duration_comp will need o be substracted by 2
dt_1y <- dt %>% 
    right_join(ts_1y_feat, by = "id")

dt_1y <- dt_1y %>% 
    mutate(duration = duration-1,
           duration_comp = duration_comp-1)
fwrite(dt_1y, "dt_1y.gz")

dt_2y <- dt %>% 
    right_join(ts_2y_feat, by = "id")

dt_2y <- dt_2y %>% 
    mutate(duration = duration-2,
           duration_comp = duration_comp-2)
fwrite(dt_2y, "dt_2y.gz")