# setup -------------------------------------------------------------------

library(tidyverse)

external_performance_file_name <- "external_performance.rds"

external_dt <- read_rds(external_performance_file_name)

external_dt <- external_dt %>% 
    select(landmark, analysis_models, performance)

write_rds(external_dt, "external_performance_without_dataset.rds")

# end ---------------------------------------------------------------------