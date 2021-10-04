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