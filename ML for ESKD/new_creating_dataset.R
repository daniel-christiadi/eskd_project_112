libraries <- c('readxl', 'tidyverse', 'janitor', 'skimr', 'heatmaply', 
               'data.table', 'gtsummary', 'lubridate')
installed <- installed.packages()[,'Package']
new.libs <- libraries[!(libraries %in% installed)]
if(length(new.libs)) install.packages(new.libs,repos="http://cran.csiro.au",
                                      dependencies=TRUE)
lib.ok <- sapply(libraries, require, character.only=TRUE)
approach <- "main" # "sens" for sensitivity analysis
extraction_date <<- as_datetime("2022-01-14")
heatmap_breaks<-c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
source("utility_functions.R")

# text data1 --------------------------------------------------------------

# lab_results <- fread(file = 'parsed_labresults.tsv', 
#                      sep = '\t')
# 
# setnames(lab_results, 'PatientObjectID', 'id')
# setnames(lab_results, 'Date_Local', 'date')
# setnames(lab_results, 'ResultName', 'test1')
# setnames(lab_results, 'OtherResultName', 'test2')
# setnames(lab_results, 'ResultString', 'value1')
# setnames(lab_results, "Result", "value2")
# setnames(lab_results, 'Units', 'unit1')
# setnames(lab_results, 'OtherUnits', 'unit2')
# setnames(lab_results, "MappedUnits", 'unit3')
# 
# pathology <- lab_results %>% 
#     select(id, date, test1, test2, value1, value2, unit1, unit2, unit3) %>%
#     mutate(test = case_when(
#         test1 == "Other" ~ test2,
#         .default = test1
#     )) %>%
#     mutate(unit = case_when(
#         unit3 == "Other" & unit1 == "NULL" ~ unit2,
#         unit3 == "Other" & unit1 == "*" ~ unit2,
#         unit3 == "Other" ~ unit1,
#         .default = unit3
#     )) %>%
#     mutate(value1 = parse_number(value1),
#            value2 = parse_number(value2),
#            value = if_else(is.na(value1), value2, value1)) %>% 
#     select(id, date, test, value, unit)

# fwrite(pathology, "pathology.gz")

# text data2 --------------------------------------------------------------

pathology <- fread('pathology.gz')

pathology <- pathology %>% 
    drop_na(value)

# test name Glucose in mmol/L, serum and HgbA1c (Serum) in %
glucose <- pathology %>% 
    filter(test == 'Glucose, serum') %>% 
    select(id:value)

glucose$test <- 'glucose'

hba1c <- pathology %>% 
    filter(test %in% c('HgbA1c (Serum)', "Hb A1C")) %>% 
    filter(value < 30) %>% # removed unrealistic result with no unit n = 2 
    select(id:value)

hba1c$test <- 'hba1c'

# test name Haemoglobin (g/L)
hb <- pathology %>% 
    filter(test %in% c('Haemoglobin','HAEMOGLOBIN')) %>% 
    filter(value <= 300) %>% # removing unrealistic Hb = 1187 g/L
    select(id:value)

hb[id == 19924 & value == 12.3]$value <- 123
hb$test <- 'haemoglobin'

# test name Haematocrit (%) : Haematocrit can be calculated (RBC*MCV)/10
hct <- pathology %>% 
    filter(test == 'Haematocrit') %>% 
    select(id:value)

hct$test <- 'haematocrit'

# test name Red Blood Cell Count in x10^12/L
rbc <- pathology %>% 
    filter(test == 'Red Blood Cell Count') %>% 
    select(id:value)

rbc$test <- 'rbc'

# test name Red Cell Distribution Width. (%)
rdw <- pathology %>% 
    filter(test %in% c('Red Cell Distribution Width.', "RDW")) %>%
    select(id:value)

rdw$test <- 'rdw'

# test name MCH in pg - MCH is calculated Hb (g/L) divided by RBC(in millions/microL)
mch <- pathology %>% 
    filter(test %in% c('Mean Corp. Haem.', 'MCH')) %>% 
    select(id:value)

mch$test <- 'mch'

# test name MCHC in g/L -  MCHC is calculated by Hb / Hct
mchc <- pathology %>% 
    filter(test == 'MCHC') %>% 
    select(id:value)

mchc$test <- 'mchc'

# test name MCV if fL
mcv <- pathology %>% 
    filter(test %in% c('Mean Cell Volume', 'MCV')) %>% 
    select(id:value)

mcv$test <- 'mcv'

# test name White Blood Cell Count in 10^9/L
wcc <- pathology %>% 
    filter(test == 'White Blood Cell Count') %>% 
    select(id:value)

wcc$test <- 'wcc'

# test name Actual Neutrophils (10^9/L)
neut <- pathology %>% 
    filter(test %in% c('Total Neuts', 'Actual Neutrophils')) %>% 
    select(id:value)

neut$test <- 'neutrophils'

# test name Actual Lymphocytes (10^9/L)
lymph <- pathology %>% 
    filter(test == 'Actual Lymphocytes') %>% 
    mutate(value = case_when(unit == 'x10*6/L' ~ value/1000,
                             .default = value)) %>% 
    select(id:value)

lymph$test <- 'lymphocytes'

# test name Actual Eosinophils (10^9/L)
eosin <- pathology %>% 
    filter(test == 'Actual Eosinophils') %>% 
    select(id:value)

eosin$test <- 'eosinophils'

# test name Actual Monocytes (10^9/L)
monocyte <- pathology %>% 
    filter(test == 'Actual Monocytes') %>% 
    select(id:value)

monocyte$test <- 'monocytes'

# test name Actual Basophils (10^9/L)
baso <- pathology %>% 
    filter(test == 'Actual Basophils') %>% 
    select(id:value)

baso$test <- 'basophils'

# test name Platelet Count (10^9/L)
platelet <- pathology %>% 
    filter(test %in% c('Platelets', 'Platelet Count')) %>% 
    select(id:value)

platelet$test <- 'platelet'

# test name Ferritin (Serum) (microg/L)
ferritin <- pathology %>% 
    filter(test == 'Ferritin (Serum)') %>% 
    select(id:value)

ferritin$test <- 'ferritin'

# test name Iron (Serum) in umol/L
iron <- pathology %>% 
    filter(test == 'Iron (Serum)') %>% 
    select(id:value)

iron$test <- 'iron'

# test name Transferrin (Serum) in umol/ - conversion rate from g/L to umol/L is 12.5628
transferrin <- pathology %>% 
    filter(test == 'Transferrin (Serum)') %>% 
    mutate(value = case_when(unit %in% c('g/l', 'g/L') ~ value*12.5628,
                             .default =  value)) %>% 
    select(id:value)

transferrin$test <- 'transferrin'

# test name Transferrin Saturation (Serum) (%)
tsat <- pathology %>% 
    filter(test == 'Transferrin Saturation (Serum)') %>% 
    select(id:value)

tsat$test <- 'tsat'

# test name Ca, adj; Calcium, serum; (mmol/L)
ca_adj <- pathology %>% 
    filter(test == "Ca, adj") %>% 
    select(id:value)

ca <- pathology %>% 
    filter(test == "Calcium, serum") %>% 
    filter(value < 10, value != 0) %>% 
    select(id:value)

# test name Phosphorus, serum (mmol/L)
phosph <- pathology %>% 
    filter(test == 'Phosphorus, serum') %>% 
    filter(value < 10) %>% # removing likely false data 
    select(id:value)

phosph$test <- 'phosphate'

# test name PTH Intact in pmol/L
pth <- pathology %>% 
    filter(test %in% c('Intact PTH', 'PTH Intact')) %>% 
    mutate(value = case_when(unit == 'pg/ml' ~ value*0.106,
                             .default = value)) %>% 
    select(id:value)

pth$test <- 'pth'

# test name Alkaline Phosphatase (Serum) in U/L
alkphos <- pathology %>% 
    filter(test %in% c('Alkaline Phosphatase (Serum)', 'Alk. Phos.', 
                       'ALKP')) %>%
    select(id:value)

alkphos$test <- 'alkphos'

# test name Magnesium (Serum) in mmol/L
mg <- pathology %>% 
    filter(test == 'Magnesium (Serum)') %>%
    filter(value <= 10) %>% # removed n = 6 values with no unit 
    select(id:value)

mg$test <- 'magnesium'

# test name Sodium, serum in mmol/L
sodium <- pathology %>% 
    filter(test == 'Sodium, serum') %>%
    filter(value != 5.8, value != 1238) %>%  # removal n = 2 values with unrealistic Na level 5.8 and 1238
    select(id:value) 

sodium$test <- 'sodium'

# test name Potassium, serum in mmol/L
potassium <- pathology %>% 
    filter(test == 'Potassium, serum') %>%
    filter(value <= 10) %>% # removal n = 5 values with unrealistic K level > 10
    select(id:value) 

potassium$test <- 'potassium'

# test name Chloride, serum in mmol/L
chloride <- pathology %>% 
    filter(test == 'Chloride, serum') %>% 
    filter(value > 70) %>%  # removal n = 9 values with unrealistic level < 20
    select(id:value)

chloride$test <- 'chloride'

# test name Bicarbonate (Serum) in mmol/L
bicarb <- pathology %>% 
    filter(test == 'Bicarbonate (Serum)') %>% 
    select(id:value)

bicarb$test <- 'bicarbonate'

# test name Alanine Aminotransferase in U/L
alt <- pathology %>% 
    filter(test == 'Alanine Aminotransferase') %>% 
    select(id:value)

alt$test <- 'alt'

# test name Albumin (Serum) in g/L
alb <- pathology %>% 
    filter(test == 'Albumin (Serum)') %>% 
    select(id:value) %>% 
    filter(value %>% between(10, 80)) # removed further n = 3 unrealistic value alb < 5 g/L and > 100

alb$test <- 'albumin'

# test name Gamma Glutamyl Transferase (Serum) in U/L
ggt <- pathology %>% 
    filter(test == 'Gamma Glutamyl Transferase (Serum)') %>% 
    select(id:value)

ggt$test <- 'ggt'

# test name Bilirubin Total (Serum) in umol/L
bili <- pathology %>% 
    filter(test == 'Bilirubin Total (Serum)') %>% 
    select(id:value) 

bili$test <- 'bilirubin'

# test name Globulin (Serum) in g/L
globulin <- pathology %>% 
    filter(test == 'Globulin (Serum)') %>% 
    select(id:value)

globulin$test <- 'globulin'

# test name Protein (Serum) in g/L
dt <- pathology %>% 
    filter(test %in% c('Protein (Serum)','Total Protein')) %>% 
    pivot_wider(id_cols = c("id", "date"),
                names_from = test,
                values_from = value,
                values_fn = mean)
setnames(dt, old = "Total Protein", "total_protein")  
setnames(dt, old = "Protein (Serum)", "serum_protein")  

protein <- dt %>% 
    mutate(value = if_else(is.na(total_protein), 
                           serum_protein, total_protein),
           test = "protein") %>% 
    select(id, date, test, value)

# test name Creatinine (Serum) in umol/L
creatinine <- pathology %>% 
    filter(test == 'Creatinine (Serum)') %>% 
    mutate(value = case_when(
        (unit == 'mmol/L' & value <= 10) ~ value*1000,
        .default = value
    )) %>% 
    filter(value > 0) %>% 
    select(id:value)

creatinine$test <- 'creatinine'

# test name eGFR mL/min/1.73m2
egfr <- pathology %>% 
    filter(test == "eGFR") %>% 
    mutate(unit_factor = if_else(unit == "NULL", "from_null", 
                                 "from_nonnull")) %>% 
    pivot_wider(id_cols = c("id", "date"),
                names_from = unit_factor,
                values_from = value,
                values_fn = mean) %>% 
    mutate(value = if_else(is.na(from_nonnull), from_null, from_nonnull),
           test = "egfr") %>% 
    filter(value <= 90) %>% 
    select(id, date, test, value)

# test name Urate in umol/L
urate <- pathology %>% 
    filter(test %in% c('Urate', 'Uric Acid (Serum)')) %>% 
    mutate(value = if_else(unit %in% c('mmol/l', 'mmol/L'), value*1000, 
                           value),
           value = if_else(value < 10, value*1000, value)) %>% 
    select(id:value)

urate$test <- 'urate'

# test name Urea (serum) in mmol/L
urea <- pathology %>% 
    filter(test == 'Urea (serum)') %>% 
    mutate(unit_factor = if_else(unit == "NULL", "from_null", 
                                 "from_nonnull")) %>% 
    pivot_wider(id_cols = c("id", "date"),
                names_from = unit_factor,
                values_from = value,
                values_fn = mean) %>% 
    mutate(value = if_else(is.na(from_nonnull), from_null, from_nonnull),
           test = "urea") %>% 
    filter(value < 200) %>%
    select(id, date, test, value)

# test name C-Reactive Protein in mg/L
crp <- pathology %>% 
    filter(test %in% c('CRP', 'C-Reactive Protein')) %>%
    mutate(test = "crp") %>% 
    select(id, date, test, value)

# test name Erythrocyte Sedimentation Rate in mm/hr
esr <- pathology %>% 
    filter(test %in% c('ESR', 'Erythrocyte Sedimentation Rate', "E.S.R.")) %>% 
    mutate(test = "esr") %>% 
    select(id:value)

# test name Protein/Creatinine Ratio in mg/mmol
protcreat <- pathology %>% 
    filter(test %in% c('Protein/creatinine', 'Protein/Creatinine', 
                       'Protein/Creatinine Ratio',
                       'Protein Creatinine Ratio (Urine)')) %>% 
    mutate(value = if_else(unit == 'g/g', value*100, value),
           test = "protcreat_ratio") %>% 
    select(id, date, test, value)

protcreat$test <- 'protcreat_ratio'

# test name Albumin/Creatinine Ratio in mg/mmol
albcreat <- pathology %>% 
    filter(test %in% c('Albumin/creatinine', 'Albumin/Creatinine', 
                       'Albumin/Creatinine Ratio', 
                       'Albumin/Creatinine Ratio (Urine)')) %>% 
    mutate(test = "albcreat_ratio") %>% 
    select(id, date, test, value)

albcreat$test <- 'albcreat_ratio'


# new_pathology -----------------------------------------------------------

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "diabetes",na = "NULL", skip = 1,
                 col_names = c("id", "date", "glucose", "hba1c"),
                 col_types = c("numeric", "date", "numeric", "numeric"))
dt <- dt %>% 
    arrange(id, date) %>% 
    mutate(hba1c = replace(hba1c, 
                           id == 22476 & hba1c == 87,
                           values = 8.7)) %>% 
    pivot_longer(cols = c("glucose", "hba1c"),
                 names_to = "test",
                 values_to = "value") %>% 
    drop_na(value)

glucose <- dt %>% 
    filter(test == "glucose") %>%
    full_join(y = glucose, by = c("id", "date"), multiple = "all") %>% 
    rename(test_x = "test.x", value_x = "value.x",
           test_y = "test.y", value_y = "value.y") %>% 
    mutate(test_x = if_else(is.na(test_x), test_y, test_x),
           value_x = if_else(is.na(value_x), value_y, value_x)) %>% 
    select(id, date, test_x, value_x) %>% 
    rename(test = test_x, value = value_x) %>% 
    distinct() %>% 
    filter(value > 0)

hba1c <- dt %>% 
    filter(test == "hba1c") %>% 
    full_join(y = hba1c, by = c("id", "date"), multiple = "all") %>% 
    rename(test_x = "test.x", value_x = "value.x",
           test_y = "test.y", value_y = "value.y") %>% 
    mutate(test_x = if_else(is.na(test_x), test_y, test_x),
           value_x = if_else(is.na(value_x), value_y, value_x)) %>% 
    select(id, date, test_x, value_x) %>% 
    rename(test = test_x, value = value_x) %>% distinct()

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "PCR", na = "NULL")
dt <- dt %>% 
    rename(id = PatientObjectID, date = ResultDate_Local,
           test = ResultName, value_num = ResultNumeric, value = Result) %>% 
    select(id, date, test, value_num, value) %>% 
    mutate(value = parse_number(value),
           test = "protcreat_ratio") %>% 
    select(-value_num)

protcreat <- dt %>% 
    full_join(y = protcreat, by = c("id", "date"), multiple = "all") %>% 
    rename(test_x = "test.x", value_x = "value.x",
           test_y = "test.y", value_y = "value.y") %>% 
    mutate(test_x = if_else(is.na(test_x), test_y, test_x),
           value_x = if_else(is.na(value_x), value_y, value_x)) %>% 
    select(id, date, test_x, value_x) %>% 
    rename(test = test_x, value = value_x) %>% distinct()

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "ACR", na = "NULL")
dt <- dt %>% 
    rename(id = PatientObjectID, date = ResultDate_Local, test = ResultName,
           test2 = OtherResultName, value = Result, 
           value_num = ResultNumeric) %>% 
    mutate(test = if_else(test == "Other", test2, test),
           value = parse_number(value),
           test = "albcreat_ratio") %>% 
    select(id, date, test, value)

albcreat <- dt %>% 
    full_join(y = albcreat, by = c("id", "date"), multiple = "all") %>% 
    rename(test_x = "test.x", value_x = "value.x",
           test_y = "test.y", value_y = "value.y") %>% 
    mutate(test_x = if_else(is.na(test_x), test_y, test_x),
           value_x = if_else(is.na(value_x), value_y, value_x)) %>% 
    select(id, date, test_x, value_x) %>% 
    rename(test = test_x, value = value_x) %>% distinct()

# converting protein/creatinine ratio to albumin/creatinine ratio. Current unit of uPCR is in mg/mmol. To convert the unit to mg/g, the value is multiplied by 8.85. In contrast to convert from mg/g to mg/mmol, the value is multiplied by 0.113
albcreat <- protcreat %>% 
    mutate(value_mgg = value*8.85) %>% 
    mutate(log_acr_mgg = case_when(
        value_mgg < 40 ~ 0.9518+0.1264*log(value_mgg),
        between(value_mgg, 40, 60, 
                incbounds = TRUE) ~ -1.2568+0.7251*log(value_mgg),
        between(value_mgg, 60, 250, 
                incbounds = TRUE) ~ -6.7837+2.0751*log(value_mgg),
        between(value_mgg, 250, 1000, 
                incbounds = TRUE) ~ -2.9649+1.3834*log(value_mgg),
        .default = -0.0239+0.9577*log(value_mgg)
    )) %>% 
    mutate(acr_mgg = exp(log_acr_mgg)) %>% 
    mutate(acr_value = acr_mgg*0.113) %>% 
    select(id, date, acr_value) %>% 
    full_join(albcreat, by = c("id", "date"), multiple = "all") %>% 
    distinct() %>% 
    arrange(id, date) %>% 
    replace_na(list(test = "albcreat_ratio")) %>% 
    mutate(value = if_else(is.na(value),
                           acr_value,
                           value)) %>% 
    select(id, date, test, value)

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "haem1", na = "NULL") 
haem <- dt %>% 
    arrange(PatientObjectID, LabTestDate_Local) %>% 
    rename(id = PatientObjectID, date = LabTestDate_Local, 
           hct = HctNum, hb = HgbNum, rbc = RBCNum, rdw = RBWNum,
           mchc = MCHCNum, mch = MCHNum, mcv = MCVNum, wbc = WBCNum,
           plt = PlateletsNum)

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "haem2", na = "NULL") 
haem <- dt %>% 
    arrange(PatientObjectID, LabTestDate_Local) %>% 
    rename(id = PatientObjectID, date = LabTestDate_Local,
           neut = NeutrophilsNum, lymph = LymphsNum, eosin = EosinophilsNum,
           mono = MonosNum, baso = BasosNum) %>% 
    full_join(y = haem, by = c("id", "date"), multiple = "all")

neut <- haem %>% 
    select(id, date, neut) %>% 
    full_join(y = neut, by = c("id", "date"), multiple = "all") %>%
    mutate(neut = if_else(is.na(neut), value, neut)) %>% 
    select(id, date, neut) %>% 
    mutate(test = "neutrophils") %>% 
    rename(value = neut) %>% 
    relocate(test, .before = value) %>% 
    distinct() 

lymph <- haem %>% 
    select(id, date, lymph) %>% 
    full_join(y = lymph, by = c("id", "date"), multiple = "all") %>% 
    mutate(lymph = if_else(lymph > 500, value, lymph),
           lymph = if_else(is.na(lymph), value, lymph)) %>% 
    select(id, date, lymph) %>% 
    mutate(test = "lymphocytes") %>% 
    rename(value = lymph) %>% 
    relocate(test, .before = value) %>% 
    distinct()

eosin <- haem %>% 
    select(id, date, eosin) %>% 
    full_join(y = eosin, by = c("id", "date"), multiple = "all") %>% 
    mutate(eosin = if_else(is.na(eosin), value, eosin)) %>% 
    select(id, date, eosin) %>% 
    mutate(test = "eosinophils") %>% 
    rename(value = eosin) %>% 
    relocate(test, .before = value) %>% distinct()

monocyte <- haem %>% 
    select(id, date, mono) %>% 
    full_join(y = monocyte, by = c("id", "date"), multiple = "all") %>% 
    mutate(mono = if_else(is.na(mono), value, mono)) %>% 
    select(id, date, mono) %>% 
    mutate(test = "monocytes") %>% 
    rename(value = mono) %>% 
    relocate(test, .before = value) %>% distinct()

baso <- haem %>% 
    select(id, date, baso) %>% 
    full_join(y = baso, by = c("id", "date"), multiple = "all") %>% 
    mutate(baso = if_else(is.na(baso), value, baso)) %>% 
    select(id, date, baso) %>% 
    mutate(test = "basophils") %>% 
    rename(value = baso) %>% 
    relocate(test, .before = value) %>% distinct()

hct <- haem %>% 
    select(id, date, hct) %>% 
    full_join(y = hct, by = c("id", "date"), multiple = "all") %>% 
    mutate(hct = if_else(is.na(hct), value, hct)) %>% 
    select(id, date, hct) %>% 
    mutate(test = "haematocrit") %>% 
    rename(value = hct) %>% 
    relocate(test, .before = value) %>% distinct()

hb <- haem %>% 
    select(id, date, hb) %>% 
    full_join(y = hb, by = c("id", "date"), multiple = "all") %>%
    mutate(hb = if_else(is.na(hb), value, hb)) %>% 
    select(id, date, hb) %>% 
    mutate(test = "haemoglobin") %>%
    rename(value = hb) %>% 
    relocate(test, .before = value) %>% distinct()

rbc <- haem %>% 
    select(id, date, rbc) %>% 
    full_join(y = rbc, by = c("id", "date"), multiple = "all") %>% 
    mutate(rbc = if_else(is.na(rbc), value, rbc)) %>% 
    select(id, date, rbc) %>% 
    mutate(test = "rbc") %>% 
    rename(value = rbc) %>% 
    relocate(test, .before = value) %>% distinct()

rdw <- haem %>% 
    select(id, date, rdw) %>% 
    full_join(y = rdw, by = c("id", "date"), multiple = "all") %>% 
    mutate(rdw = if_else(is.na(rdw), value, rdw)) %>% 
    select(id, date, rdw) %>% 
    mutate(test = "rdw") %>% 
    rename(value = rdw) %>% 
    relocate(test, .before = value) %>% distinct()

mch <- haem %>% 
    select(id, date, mch) %>% 
    full_join(y = mch, by = c("id", "date"),multiple = "all") %>% 
    mutate(mch = if_else(is.na(mch), value, mch)) %>% 
    select(id, date, mch) %>% 
    mutate(test = "mch") %>% 
    rename(value = mch) %>% 
    relocate(test, .before = value) %>% distinct()

mchc <- haem %>% 
    select(id, date, mchc) %>% 
    full_join(y = mchc, by = c("id", "date"),multiple = "all") %>% 
    mutate(mchc = if_else(is.na(mchc), value, mchc)) %>% 
    select(id, date, mchc) %>% 
    mutate(test = "mchc") %>% 
    rename(value = mchc) %>% 
    relocate(test, .before = value) %>% distinct()

mcv <- haem %>% 
    select(id, date, mcv) %>% 
    full_join(y = mcv, by = c("id", "date"), multiple = "all") %>% 
    mutate(mcv = if_else(is.na(mcv), value, mcv)) %>% 
    select(id, date, mcv) %>% 
    mutate(test = "mcv") %>% 
    rename(value = mcv) %>% 
    relocate(test, .before = value) %>% distinct()

wcc <- haem %>% 
    select(id, date, wbc) %>% 
    full_join(y = wcc, by = c("id", "date"), multiple = "all") %>% 
    mutate(wbc = if_else(is.na(wbc), value, wbc)) %>% 
    select(id, date, wbc) %>% 
    mutate(test = "wcc") %>% 
    rename(value = wbc) %>% 
    relocate(test, .before = value) %>% distinct()

platelet <- haem %>% 
    select(id, date, plt) %>% 
    full_join(y = platelet, by = c("id", "date"), multiple = "all") %>% 
    mutate(plt = if_else(is.na(plt), value, plt)) %>% 
    select(id, date, plt) %>% 
    mutate(test = "platelet") %>% 
    rename(value = plt) %>% 
    relocate(test, .before = value) %>% distinct()

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "iron", na = "NULL", skip = 1,
                 col_names = c("id", "date", "ferritin", "iron",
                               "transferrin", "tsat"))
dt <- dt %>% 
    arrange(id, date)

ferritin <- dt %>% 
    select(id, date, ferritin) %>% 
    full_join(y = ferritin, by = c("id", "date"), multiple = "all") %>% 
    mutate(ferritin = if_else(is.na(ferritin), value, ferritin)) %>% 
    filter(ferritin < 60000) %>% 
    select(id, date, ferritin) %>% 
    mutate(test = "ferritin") %>% 
    rename(value = ferritin) %>% 
    relocate(test, .before = value) %>% distinct()

iron <- dt %>% 
    select(id, date, iron) %>% 
    full_join(y = iron, by = c("id", "date"), multiple = "all") %>% 
    mutate(iron = if_else(is.na(iron), value, iron)) %>% 
    select(id, date, iron) %>% 
    mutate(test = "iron") %>% 
    rename(value = iron) %>% 
    relocate(test, .before = value) %>% distinct()

transferrin <- dt %>% 
    select(id, date, transferrin) %>% 
    full_join(y = transferrin, by = c("id", "date"), multiple = "all") %>% 
    mutate(transferrin = if_else(is.na(transferrin), value, transferrin)) %>% 
    select(id, date, transferrin) %>% 
    mutate(test = "transferrin") %>%
    rename(value = transferrin) %>% 
    relocate(test, .before = value) %>% distinct()

tsat <- dt %>% 
    select(id, date, tsat) %>% 
    full_join(y = tsat, by = c("id", "date"), multiple = "all") %>% 
    mutate(tsat = if_else(is.na(tsat), value, tsat)) %>% 
    select(id, date, tsat) %>% 
    mutate(test = "tsat") %>% 
    rename(value = tsat) %>% 
    relocate(test, .before = value) %>% distinct()

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "calcium", na = "NULL", skip = 1,
                 col_names = c("id", "date", "calcium", "adjca",
                               "phosphate", "pth", "alkphos"))
dt <- dt %>% 
    arrange(id, date)

ca_adj <- dt %>% 
    select(id, date, adjca) %>% 
    full_join(y = ca_adj, by = c("id", "date"), multiple = "all") %>% 
    mutate(adjca = if_else(is.na(adjca), value, adjca)) %>% 
    select(id, date, adjca) %>% 
    mutate(test = "adjca") %>% 
    rename(value = adjca) %>% 
    relocate(test, .before = value) %>% distinct() 

ca <- dt %>% 
    select(id, date, calcium) %>% 
    full_join(y = ca, by = c("id", "date"), multiple = "all") %>% 
    mutate(value = if_else(is.na(value), calcium, value)) %>%
    filter(value < 10) %>%
    select(id, date, value) %>% 
    mutate(test = "calcium") %>% 
    relocate(test, .before = value) %>% distinct()

phosph <- dt %>% 
    select(id, date, phosphate) %>% 
    full_join(y = phosph, by = c("id", "date"), multiple = "all") %>% 
    mutate(phosphate = if_else(is.na(phosphate), value, phosphate)) %>% 
    filter(phosphate < 10) %>% 
    select(id, date, phosphate) %>% 
    mutate(test = "phosphate") %>% 
    rename(value = phosphate) %>% 
    relocate(test, .before = value) %>% distinct()

pth <- dt %>% 
    select(id, date, pth) %>% 
    full_join(y = pth, by = c("id","date"), multiple = "all") %>% 
    mutate(value = if_else(is.na(value), pth, value)) %>% 
    select(id, date, value) %>% 
    mutate(test = "pth") %>% 
    relocate(test, .before = value) %>% distinct() %>% drop_na(value)

alkphos <- dt %>% 
    select(id, date, alkphos) %>% 
    full_join(y = alkphos, by = c("id", "date"), multiple = "all") %>% 
    mutate(value = if_else(is.na(value), alkphos, value)) %>% 
    select(id, date, value) %>% 
    mutate(test = "alkphos") %>% 
    relocate(test, .before = value) %>% distinct() %>% drop_na(value)

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "mag", na = "NULL", skip = 1,
                 col_names = c("id", "date", "magnesium"))
dt <- dt %>% 
    arrange(id, date)

mg <- dt %>% 
    full_join(y = mg, by = c("id", "date"), multiple = "all") %>% 
    mutate(magnesium = if_else(is.na(magnesium), value, magnesium)) %>% 
    filter(magnesium < 10) %>% 
    select(id, date, magnesium) %>% 
    mutate(test = "magnesium") %>% 
    rename(value = magnesium) %>% 
    relocate(test, .before = value) %>% distinct()

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "electrolytes",na = "NULL", skip = 1,
                 col_names = c("id", "date", "sodium", "potassium",
                               "chloride", "bicarb", "creatinine", 
                               "urate", "urea"),
                 col_types = c("numeric", "date", "numeric", "numeric",
                               "numeric", "numeric", "numeric",
                               "numeric", "numeric"))
dt <- dt %>% 
    arrange(id, date)

sodium <- dt %>% 
    select(id, date, sodium) %>% 
    full_join(y = sodium, by = c("id", "date"), multiple = "all") %>% 
    mutate(sodium = if_else(is.na(sodium), value, sodium)) %>% 
    filter(sodium != 5.8, sodium != 1238) %>% 
    select(id, date, sodium) %>% 
    mutate(test = "sodium") %>% 
    rename(value = sodium) %>% 
    relocate(test, .before = value) %>% distinct()

potassium <- dt %>% 
    select(id, date, potassium) %>% 
    full_join(y = potassium, by = c("id", "date"), multiple = "all") %>% 
    mutate(potassium = if_else(is.na(potassium), value, potassium)) %>% 
    filter(potassium != 0, potassium < 10) %>% 
    select(id, date, potassium) %>% 
    mutate(test = "potassium") %>% 
    rename(value = potassium) %>% 
    relocate(test, .before = value) %>% distinct()

chloride <- dt %>% 
    select(id, date, chloride) %>% 
    full_join(y = chloride, by = c("id", "date"), multiple = "all") %>% 
    filter(chloride > 70) %>% 
    select(id, date, chloride) %>% 
    mutate(test = "chloride") %>% 
    rename(value = chloride) %>% 
    relocate(test, .before = value) %>% distinct()

bicarb <- dt %>% 
    select(id, date, bicarb) %>% 
    full_join(y = bicarb, by = c("id", "date"), multiple = "all") %>% 
    mutate(bicarb = if_else(is.na(bicarb), value, bicarb)) %>% 
    filter(bicarb < 50) %>% 
    select(id, date, bicarb) %>% 
    mutate(test = "bicarb") %>% 
    rename(value = bicarb) %>% 
    relocate(test, .before = value) %>% distinct() 

creatinine <- dt %>% 
    select(id, date, creatinine) %>% 
    full_join(y = creatinine, by = c("id", "date"), multiple = "all") %>%
    mutate(value = if_else(is.na(value), creatinine, value),
           value = if_else(value < 10, value*1000, value)) %>% 
    select(id, date, value) %>% 
    mutate(test = "creatinine") %>% 
    relocate(test, .before = value) %>% distinct()

urate <- dt %>% 
    select(id, date, urate) %>% 
    full_join(y = urate, by = c("id", "date"), multiple = "all") %>% 
    mutate(value = if_else(is.na(value), urate, value),
           value = if_else(value < 40, value*1000, value)) %>% 
    select(id, date, value) %>% 
    mutate(test = "urate") %>% 
    relocate(test, .before = value) %>% distinct() %>% drop_na(value)

urea <- dt %>% 
    select(id, date, urea) %>% 
    full_join(y = urea, by = c("id", "date"), multiple = "all") %>%
    mutate(urea = if_else(is.na(urea), value, urea)) %>% 
    select(id, date, urea) %>% 
    mutate(test = "urea") %>% 
    rename(value = urea) %>% 
    relocate(test, .before = value) %>% distinct()

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "lfts", na = "NULL", skip = 1,
                 col_names = c("id", "date", "alkphos", "alt",
                               "albumin", "ggt", "bilirubin", "protein"),
                 col_types = c("numeric", "date", "skip", "numeric", 
                               "numeric", "numeric", "numeric", "numeric"))
dt <- dt %>% 
    arrange(id, date)

alt <- dt %>% 
    select(id, date, alt) %>% 
    full_join(y = alt, by = c("id", "date"), multiple = "all") %>% 
    mutate(alt = if_else(is.na(alt), value, alt)) %>% 
    select(id, date, alt) %>% 
    mutate(test = "alt") %>% 
    rename(value = alt) %>% 
    relocate(test, .before = value) %>% distinct()

alb <- dt %>% 
    select(id, date, albumin) %>% 
    full_join(y = alb, by = c("id", "date"), multiple = "all") %>% 
    mutate(albumin = if_else(is.na(albumin), value, albumin)) %>% 
    filter(between(albumin, 10, 100)) %>% 
    select(id, date, albumin) %>% 
    mutate(test = "albumin") %>% 
    rename(value = "albumin") %>% 
    relocate(test, .before = value) %>% distinct()

ggt <- dt %>% 
    select(id, date, ggt) %>% 
    full_join(y = ggt, by = c("id", "date"), multiple = "all") %>% 
    mutate(value = if_else(is.na(value), ggt, value)) %>% 
    select(id, date, value) %>% 
    mutate(test = "ggt") %>% 
    relocate(test, .before = value) %>% distinct() 

bili <- dt %>% 
    select(id, date, bilirubin) %>% 
    full_join(y = bili, by = c("id", "date"), multiple = "all") %>% 
    mutate(value = if_else(is.na(value), bilirubin, value)) %>% 
    select(id, date, value) %>% 
    mutate(test = "bilirubin") %>% 
    relocate(test, .before = value) %>% distinct() 

protein <- dt %>% 
    select(id, date, protein) %>% 
    full_join(y = protein, by = c("id", "date"), multiple = "all") %>%
    mutate(value = if_else(is.na(value), protein, value)) %>% 
    select(id, date, value) %>% 
    mutate(test = "protein") %>% 
    relocate(test, .before = value) %>% distinct() %>% drop_na(value)

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "crp", skip = 1, na = "NULL",
                 col_names = c("id", "date", "crp"),
                 col_types = c("numeric", "date", "numeric"))
dt <- dt %>% 
    arrange(id, date)

crp <- dt %>% 
    full_join(y = crp, by = c("id", "date"), multiple = "all") %>% 
    mutate(crp = if_else(is.na(crp), value, crp)) %>% 
    select(id, date, crp) %>% 
    mutate(test = "crp") %>% 
    rename(value = crp) %>% 
    relocate(test, .before = value) %>% distinct() 

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "globulin", skip = 1, na = "NULL",
                 col_names = c("id", "date", "globulin"))

globulin <- dt %>% 
    arrange(id, date) %>% 
    full_join(y = globulin, by = c("id", "date"),multiple = "all") %>% 
    mutate(globulin = if_else(is.na(globulin), value, globulin)) %>% 
    select(id, date, globulin) %>% 
    mutate(test = "globulin") %>% 
    rename(value = globulin) %>% 
    relocate(test, .before = value) %>% distinct() 

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "egfr", skip = 1, na = "NULL",
                 col_names = c("id", "date", "egfrnum", "egfr", "creatinine"),
                 col_types = c("numeric", "date", "numeric", "text",
                               "numeric"))

dt <- dt %>% 
    arrange(id, date) %>% 
    mutate(egfr = parse_number(egfr)) %>% 
    select(-egfrnum)

egfr <- dt %>% 
    select(id, date, egfr) %>% 
    full_join(y = egfr, by = c("id", "date"), multiple = "all") %>% 
    mutate(value = if_else(is.na(value), egfr, value)) %>% 
    filter(value <= 90) %>% 
    select(id, date, value) %>% 
    mutate(test = "egfr") %>% 
    relocate(test, .before = value) %>% 
    distinct() %>% 
    arrange(id, date)

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "esr", skip = 1, na = "NULL",
                 col_names = c("id", "date", "esrnum", "esr"),
                 col_types = c("numeric", "date", "numeric", "text"))

dt <- dt %>% 
    arrange(id, date) %>% 
    mutate(esr = parse_number(esr)) %>% 
    select(-esrnum)

esr <- dt %>% 
    full_join(y = esr, by = c("id", "date"), multiple = "all") %>% 
    mutate(esr = if_else(is.na(esr), value, esr)) %>% 
    select(id, date, esr) %>% 
    mutate(test = "esr") %>% 
    rename(value = esr) %>% 
    relocate(test, .before = value) %>% distinct()

# vital signs -------------------------------------------------------------

dt <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                 sheet = "Exams", na = "NULL", skip = 1,
                 col_types = c("numeric", "date", "numeric", "numeric",
                               "numeric", "numeric", "numeric", "numeric",
                               "numeric", "numeric"),
                 col_names = c("id", "date", "sitsys", "sitdias", "sithr",
                               "stasys", "stadias", "stahr", "weight",
                               "height"))

dt <- dt %>% 
    arrange(id, date)

sitsys <- dt %>% 
    select(id, date, sitsys) %>% 
    filter(between(sitsys, 60, 250)) %>% 
    mutate(test = "sitsys") %>% 
    rename(value = sitsys) %>% 
    relocate(test, .before = value) %>% distinct()

sitdias <- dt %>% 
    select(id, date, sitdias) %>% 
    filter(sitdias > 20) %>% 
    mutate(test = "sitdias") %>% 
    rename(value = sitdias) %>% 
    relocate(test, .before = value) %>% distinct()

sithr <- dt %>% 
    select(id, date, sithr) %>%
    mutate(test = "sithr") %>% 
    rename(value = sithr) %>% 
    relocate(test, .before = value) %>% distinct()

stasys <- dt %>% 
    select(id, date, stasys) %>% 
    filter(between(stasys, 50, 250)) %>% 
    mutate(test = "stasys") %>% 
    rename(value = stasys) %>% 
    relocate(test, .before = value) %>% distinct()

stadias <- dt %>% 
    select(id, date, stadias) %>% 
    filter(stadias > 20) %>% 
    mutate(test = "stadias") %>% 
    rename(value = stadias) %>% 
    relocate(test, .before = value) %>% distinct()

stahr <- dt %>% 
    select(id, date, stahr) %>% 
    filter(between(stahr, 30, 160)) %>%
    mutate(test = "stahr") %>% 
    rename(value = stahr) %>% 
    relocate(test, .before = value) %>% distinct()

weight <- dt %>% 
    select(id, date, weight) %>% 
    filter(between(weight, 20, 310)) %>% 
    mutate(test = "weight") %>% 
    rename(value = weight) %>% 
    relocate(test, .before = value) %>% distinct()

height <- dt %>% 
    select(id, date, height) %>% 
    filter(between(height, 50, 250)) %>% 
    mutate(test = "height") %>% 
    rename(value = height) %>% 
    relocate(test, .before = value) %>% distinct()

vital_signs <- na.omit(unique(rbindlist(
    list(sitsys, sitdias, sithr, stasys, stadias, stahr, weight, height)
)))

dt <- vital_signs %>% 
    pivot_wider(id_cols = c("id", "date"),
                names_from = test,
                values_from = value,
                values_fn = mean) %>% 
    mutate(height_m = height / 100,
           bmi = weight/(height_m^2),
           postural_sys = stasys - sitsys,
           postural_dias = stadias - sitdias,
           postural_hr = stahr - sithr)

bmi_posturals <- dt %>% 
    select(id, date, bmi:postural_hr) %>% 
    pivot_longer(cols = bmi:postural_hr,
                 names_to = "test",
                 values_to = "value")

# clinicopathology compilation --------------------------------------------

melted <- na.omit(unique(rbindlist(
    list(alb, albcreat, alkphos, alt, baso, bicarb, bili, ca, ca_adj, chloride, 
         creatinine, crp, egfr, eosin, esr, ferritin, ggt, globulin, 
         glucose, hb, hba1c, hct, iron, lymph, mch, mchc, mcv, mg, monocyte,
         neut, phosph, platelet, potassium, protcreat, protein, pth, rbc, rdw, 
         sodium, transferrin, tsat, urate, urea, wcc, vital_signs,
         bmi_posturals))))

# write_rds(melted, "melted.rds")

# demographic and modality ------------------------------------------------

# Extracting from lab_results (raw text), there were 10577 distinct IDs. 
# While, from pathology, there were 10,564 distinct IDs. 
# From the eGFR perspective (new extraction), there were 9,700 distinct IDs.
# Extracting patients only, we had extracted 15,481 distinct IDs. 
# While, on melted had 10198 distinct IDs (based on test == "egfr"). 

patient <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                      sheet = "IDs", na = "NULL",
                      col_names = c("id", "birth_year", "sex"),
                      col_types = c("numeric", "numeric", "text"),
                      skip = 1)

patient <- patient %>% distinct(id, .keep_all = TRUE)

dt <- read_excel("~/PhD Project 1/ML for ESKD/need_gender_clarification.xlsx",
                 sheet = "Sheet1", na = "NULL",
                 col_names = c("id", "sex"),
                 col_types = c("numeric", "text"),
                 skip = 1)

patient <- patient %>% 
    left_join(y = dt, by = "id") %>% 
    mutate(sex.x = case_when(
        sex.x == "Unknown" ~ sex.y,
        sex.x == "Indeterminate" ~ sex.y,
        .default = sex.x
    )) %>% 
    select(id, birth_year, sex.x) %>% 
    rename(sex = sex.x) %>% 
    arrange(id)

dt <- read_excel("~/PhD Project 1/ML for ESKD/Diagnosis.xlsx",
                 sheet = "GN Diagnosis", skip = 1, na = "NULL",
                 col_names = c("id", "diagnosis", "diag_date"),
                 col_types = c("numeric", "text", "date"))

lupus <- tibble(unified_diag = "lupus", 
                diagnosis = c('Systemic lupus erythematosus glomerulonephritis syndrome', 
                              'Systemic lupus erythematosus glomerulonephritis syndrome, World Health Organization (WHO) class III', 
                              'Systemic lupus erythematosus glomerulonephritis syndrome, World Health Organization (WHO) class IV', 
                              'Systemic lupus erythematosus glomerulonephritis syndrome, World Health Organization (WHO) class V', 
                              'Systemic lupus erythematosus glomerulonephritis syndrome, World Health Organization (WHO) class VI'))

unspecific <- tibble(unified_diag = "unspecific",
                     diagnosis = c('Acute glomerulonephritis', 
                                   'Chronic diffuse glomerulonephritis',
                                   'Chronic focal glomerulonephritis', 
                                   'Chronic glomerulonephritis',
                                   'Crescentic glomerulonephritis',
                                   'Focal AND segmental proliferative glomerulonephritis',
                                   'Glomerulonephritis', 
                                   'Immune-complex glomerulonephritis',
                                   'Necrotizing glomerulonephritis',
                                   'Proliferative glomerulonephritis', 
                                   'Renal vasculitis',
                                   'Vasculitis secondary to drug',
                                   'Acute post-streptococcal glomerulonephritis',
                                   'Post-infectious glomerulonephritis'))

mp_mcgn <- tibble(unified_diag = "mp_mcgn",
                  diagnosis = c('Acute nephritic syndrome, diffuse mesangial proliferative glomerulonephritis', 
                                'C3 Glomerulonephritis',
                                'Mesangial proliferative glomerulonephritis',
                                'Mesangial proliferative glomerulonephritis',
                                'Mesangiocapillary glomerulonephritis',
                                'Mesangiocapillary glomerulonephritis, type I',
                                'Mesangiocapillary glomerulonephritis, type II'))

anca <- tibble(unified_diag = "anca",
               diagnosis = c('Eosinophilc Granulomatosis with Polyangiitis',
                             'Granulomatous Polyangiitis', 
                             'Microscopic polyangiitis',
                             'Microscopic Polyangiitis',
                             'Primary pauci-immune necrotizing and crescentic glomerulonephritis',
                             "Wegener's granulomatosis"))
fg_itg <- tibble(unified_diag = "fg_itg",
                 diagnosis = c('Fibrillary glomerulonephritis', 
                               'Immunotactoid Glomerulopathy'))

fsgs <- tibble(unified_diag = "fsgs",
               diagnosis = c('Focal Segmental Glomerulosclerosis'))
igan <- tibble(unified_diag = "igan",
               diagnosis = c('IgA nephropathy', 'Primary IgA nephropathy'))
mn <- tibble(unified_diag = "mn",
             diagnosis = c('Membranous glomerulonephritis',
                           'Membranous glomerulonephritis - stage III',
                           'Nephrotic syndrome, diffuse membranous glomerulonephritis')) 
mcd <- tibble(unified_diag = "mcd",
              diagnosis = c('Minimal change disease',
                            'Steroid-dependent minimal change glomerulonephritis'))

gn_reference <- rbindlist(list(lupus, unspecific, mp_mcgn, anca, fg_itg, 
                               fsgs, igan, mn, mcd))

gn <- dt %>% 
    select(-diag_date) %>% 
    left_join(gn_reference, by = "diagnosis", multiple = "all") %>% 
    select(-diagnosis) %>% 
    rename(diagnosis = unified_diag) %>% 
    arrange(id) %>% drop_na()

gn <- gn %>% 
    mutate(diagnosis = replace(diagnosis, 
                               id == 4339 & diagnosis == "mcd", NA),
           diagnosis = replace(diagnosis, 
                               id == 4551 & diagnosis == "mn", NA),
           diagnosis = replace(diagnosis, 
                               id == 10143 & diagnosis == "igan", NA),
           diagnosis = replace(diagnosis, 
                               id == 15613 & diagnosis == "mcd", NA),
           diagnosis = replace(diagnosis, 
                               id == 18367 & diagnosis == "unspecific", 
                               NA)) %>%
    drop_na() %>% 
    distinct(id, .keep_all = TRUE) 

death <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                    sheet = "Mortality", na = "NULL",
                    col_names = c("id", "death_date"),
                    col_types = c("numeric", "date"),
                    skip = 1)
death <- death %>% drop_na()

modality <- read_excel("~/PhD Project 1/ML for ESKD/new_pathology.xlsx",
                       sheet = "Modality", na = "NULL",
                       col_names = c("id", "start", "stop", "modality"),
                       col_types = c("numeric", "date", "date", "text"),
                       skip = 1)
modality <- modality %>% 
    arrange(id, start) %>% 
    replace_na(replace = list(stop = extraction_date))

modality <- modality %>% 
    group_by(id) %>% 
    nest() %>% 
    mutate(recovery = map(data, recovery_check)) %>% 
    mutate(any_recovery = map_dbl(recovery, max)) %>% 
    unnest(cols = c("data", "recovery")) %>% 
    mutate(corr_recovery = lag(recovery)) %>% 
    ungroup()

modality <- modality %>% 
    mutate(dur_recovery = if_else(
        corr_recovery == 1, as.duration(stop-start)/dmonths(1), 0, missing = 0))

eskd_reference <- c("Haemodiafiltration", "Haemofiltration", "HD",
                    "Nocturnal HD", "PD", "Transplant")

eskd_atbegin <- modality %>% 
    arrange(id, start) %>% 
    group_by(id) %>% 
    slice_head() %>% 
    ungroup() %>% 
    mutate(eskd_atbegin = if_else(
        modality %in% eskd_reference, "yes", "no"
    )) %>%
    filter(eskd_atbegin == "yes") %>% 
    select(id, start, modality) 

# choosing pts with no recovery to ensure the validity of date of ESKD 
# by modality
eskd_event <- modality %>% 
    arrange(id, start) %>% 
    filter(modality %in% eskd_reference, !id %in% eskd_atbegin$id, 
           any_recovery == 0) %>% 
    group_by(id) %>% 
    slice_head() %>% 
    ungroup() %>% 
    select(id, start, modality)

eskd_bymodality <- rbindlist(list(eskd_atbegin, eskd_event))

eskd_bymodality <- eskd_bymodality %>% 
    rename("eskd_bymodality_date" = start) %>% 
    select(-modality)

# initial egfr, ESKD, and last egfr ---------------------------------------

# There were 15 IDs without baseline demographics (sex and birthyear) and 
# 495 IDs with pre-existing ESKD (Hx of ESKD by modality before initial eGFR). 
# Hence, a total of 9688 IDs for outcome_dt

initial_egfr <- melted[test == "egfr"] %>% 
    .[order(id, date)] %>% 
    .[, .SD[1], by = id] %>% 
    .[, .(id, date)]

initial_egfr <- initial_egfr %>% 
    left_join(patient, by = c("id")) %>% 
    drop_na(birth_year) %>% 
    left_join(death, by = c("id")) %>% 
    left_join(gn, by = c("id")) %>% 
    left_join(eskd_bymodality, by = c("id")) %>% 
    rename("initial_egfr_date" = "date") %>%
    mutate(preexisting_eskd = if_else(initial_egfr_date >= eskd_bymodality_date, 
                                      "yes", 
                                      "no", missing = "no")) %>% 
    filter(preexisting_eskd == "no") %>% 
    select(-preexisting_eskd)

last_egfr <- melted[test == "egfr"] %>% 
    .[order(id, date)] %>% 
    .[, .SD[.N], by = id] %>% 
    .[, .(id, date, value)]

setnames(last_egfr, old = c("date", "value"), 
         new = c("last_egfr_date", "last_egfr_value"))

eskd_egfr <- melted %>% 
    filter(test == "egfr", id %in% initial_egfr$id) %>% 
    select(-test) %>% 
    arrange(id, date) %>% 
    group_by(id) %>% 
    filter(value <= 15) %>% 
    slice_head() %>% 
    ungroup() %>% 
    rename("eskd_byegfr_date" = "date", "eskd_egfr_value" = "value")

# create outcome and feature dataset --------------------------------------

# When the pts do not experience event,
# Main approach : extraction date as the date of censoring
# Sensitivity analysis : the last pathology as the date of censoring

outcome_dt <- initial_egfr %>% 
    left_join(last_egfr, by = c("id")) %>% 
    left_join(eskd_egfr, by = c("id"))

outcome_dt <- outcome_dt %>% 
    mutate(eskd_date = case_when(
        !is.na(eskd_bymodality_date) & is.na(eskd_byegfr_date) ~ eskd_bymodality_date,
        !is.na(eskd_bymodality_date) ~ eskd_byegfr_date,
        is.na(eskd_bymodality_date) & !is.na(eskd_byegfr_date) & last_egfr_value < 17 ~ eskd_byegfr_date,
    )) %>%  
    mutate(eskd = if_else(!is.na(eskd_date), "eskd", "no"))

if (approach == "main"){
    outcome_dt[is.na(eskd_date), eskd_date := extraction_date]
    
    outcome_dt <- outcome_dt %>% 
        mutate(event_date = pmin(eskd_date, death_date, na.rm = TRUE),
               event = if_else(event_date == eskd_date, eskd, "death"))
    
    outcome_dt <- outcome_dt %>% 
        select(id, initial_egfr_date, birth_year, sex, diagnosis, 
               event_date, event) %>% 
        mutate(age_init = year(initial_egfr_date) - birth_year,
               time_compete = as.duration(event_date - initial_egfr_date)/dyears(1)) %>% 
        rename("status_compete" = "event", "gn_cat" = "diagnosis") %>% 
        filter(time_compete >= 0) %>% 
        select(id, initial_egfr_date, sex, age_init, gn_cat, time_compete, status_compete)
    
    outcome_dt <- outcome_dt %>% 
        mutate(initial_date = floor_date(initial_egfr_date, unit = "month")) %>% 
        select(-initial_egfr_date) %>% 
        relocate(initial_date, .before = sex)
        
    dt <- melted[id %chin% outcome_dt$id,]
    
    raw_features_freq <- extracting_freq(dt)
    
    raw_features_freq <- raw_features_freq %>% 
        select(-id)
    
    raw_features_freq[is.na(raw_features_freq)] <- 0
    
    heatmaply(x = t(as.matrix(raw_features_freq)), ylab = 'Features', xlab = 'Patients', 
              main = 'Data Density of Original Clinico-Pathological Data',
              dist_method = 'manhattan', hclust_method = 'ward.D2',
              Rowv = TRUE, Colv = TRUE, dendrogram = 'both', 
              showticklabels = c(FALSE, TRUE), 
              scale_fill_gradient_fun = scale_fill_viridis(trans='sqrt', 
                                                           breaks=heatmap_breaks,
                                                           labels=heatmap_breaks, 
                                                           option = 'B', 
                                                           direction = -1,
                                                           na.value='white'))
    
    features_dt <- dt %>% 
        mutate(resample_date = floor_date(date, unit = "month")) %>% 
        group_by(id, test, resample_date) %>% 
        mutate(monthly_value = mean(value)) %>% 
        ungroup() %>% 
        select(id, resample_date, test, monthly_value) %>% 
        distinct() %>% 
        rename("value" = "monthly_value", "date" = "resample_date")
} else{
    # "sens" definition of censoring
    outcome_dt[is.na(eskd_date), eskd_date := last_egfr_date]
    
    outcome_dt <- outcome_dt %>% 
        mutate(event_date = case_when(
            eskd == "eskd" ~ pmin(eskd_date, death_date, na.rm = TRUE),
            .default = pmax(eskd_date, death_date, na.rm = TRUE))) %>% 
        mutate(event = if_else(event_date == eskd_date, eskd, "death"))
    
    dt <- melted %>% 
        right_join(select(outcome_dt, id, event_date), by = c('id'))
    
    raw_features_freq <- extracting_freq(dt)
    
    raw_features_freq <- raw_features_freq %>% 
        select(-id)
    
    raw_features_freq[is.na(raw_features_freq)] <- 0
    
    heatmaply(x = t(as.matrix(raw_features_freq)), ylab = 'Features', xlab = 'Patients', 
              main = 'Data Density of Sensitivity Clinico-Pathological Data',
              dist_method = 'manhattan', hclust_method = 'ward.D2',
              Rowv = TRUE, Colv = TRUE, dendrogram = 'both', 
              showticklabels = c(FALSE, TRUE), 
              scale_fill_gradient_fun = scale_fill_viridis(trans='sqrt', 
                                                           breaks=heatmap_breaks,
                                                           labels=heatmap_breaks, 
                                                           option = 'B', 
                                                           direction = -1,
                                                           na.value='white'))
    
    features_dt <- dt %>% 
        mutate(resample_date = floor_date(date, unit = "month")) %>% 
        group_by(id, test, resample_date) %>% 
        mutate(monthly_value = mean(value)) %>% 
        ungroup() %>% 
        select(id, resample_date, test, monthly_value, event_date) %>% 
        distinct() %>% 
        rename("value" = "monthly_value") %>% 
        mutate(yr_frmevent = as.duration(event_date - resample_date)/dyears(1)) 
    
    # filtering clinico-pathologial test within 10 years of event_date
    features_dt <- features_dt %>% 
        filter(dplyr::between(yr_frmevent, 0, 10)) %>% 
        rename("date" = "resample_date") %>%
        select(id, date, test, value)
    
    outcome_dt <- features_dt %>% 
        arrange(id, date) %>% 
        filter(test == "egfr") %>% 
        group_by(id) %>% 
        slice_head() %>% 
        ungroup() %>% 
        select(id, date) %>% 
        rename("initial_date" = "date") %>%
        right_join(outcome_dt, by = "id") %>% 
        drop_na(initial_date) %>% 
        select(id, initial_date, birth_year, sex, diagnosis, event_date, event)
    
    outcome_dt <- outcome_dt %>% 
        mutate(age_init = year(initial_date) - birth_year,
               time_compete = as.duration(event_date-initial_date)/dyears(1)) %>% 
        rename("status_compete" = "event", "gn_cat" = "diagnosis") %>% 
        select(id, initial_date, sex, age_init, gn_cat, time_compete, status_compete)
}

features_file_name <- str_c(approach, "_features.rds")
write_rds(features_dt, features_file_name)
outcome_file_name <- str_c(approach, "_outcome.rds")
write_rds(outcome_dt, outcome_file_name)


# resampled dataset heatmap ------------------------------------------------

outcome_dt <- read_rds("main_outcome.rds")
features_dt <- read_rds("main_features.rds")

# outcome_dt <- read_rds("sens_outcome.rds")
# features_dt <- read_rds("sens_features.rds")
# approach <- "sens"

features_freq <- extracting_freq(features_dt)

features_freq <- features_freq %>% 
    select(-id)

features_freq[is.na(features_freq)] <- 0

heatmaply(x = t(as.matrix(features_freq)), ylab = 'Features', xlab = 'Patients', 
          main = 'Data Density of Resampled Clinico-Pathological Data',
          dist_method = 'manhattan', hclust_method = 'ward.D2',
          Rowv = TRUE, Colv = TRUE, dendrogram = 'both', 
          showticklabels = c(FALSE, TRUE), 
          scale_fill_gradient_fun = scale_fill_viridis(trans='sqrt', 
                                                       breaks=heatmap_breaks,
                                                       labels=heatmap_breaks, 
                                                       option = 'B', 
                                                       direction = -1,
                                                       na.value='white'))

# filter dataset based on test frequency ----------------------------------

meta_dt <- features_dt %>% 
    right_join(outcome_dt, by = c("id")) %>% 
    mutate(relyear = as.duration(date - initial_date)/dyears(1)) %>% 
    filter(relyear >= 0) %>% 
    select(id, test, relyear, value, sex, age_init, gn_cat, time_compete,
           status_compete) %>% 
    arrange(id, relyear)

meta_freq <- extracting_freq(meta_dt)

dense3_id <- meta_freq %>% 
    filter(egfr >= 3) %>% 
    select(id, albumin, alkphos, bicarb, calcium, chloride, glucose, 
           haemoglobin, phosphate, platelet, potassium, sodium, wcc) %>% 
    filter(if_all(.cols = !id, .fns = ~. >=3)) %>% 
    pull(id)

dense3_dt <- meta_dt %>% 
    filter(id %in% dense3_id)

file_name <- str_c(approach, "dense3_dt.rds", sep = "_")
write_rds(dense3_dt, file_name)

dense3_matrix <- meta_freq %>% 
    filter(id %in% dense3_id) %>% 
    select(-id)

dense3_matrix[is.na(dense3_matrix)] <- 0

heatmaply(x = t(as.matrix(dense3_matrix)), ylab = 'Features', xlab = 'Patients', 
          main = 'Data Density of Resampled Clinico-Pathological Data for Dense3',
          dist_method = 'manhattan', hclust_method = 'ward.D2',
          Rowv = TRUE, Colv = TRUE, dendrogram = 'both', 
          showticklabels = c(FALSE, TRUE), 
          scale_fill_gradient_fun = scale_fill_viridis(trans='sqrt', 
                                                       breaks=heatmap_breaks,
                                                       labels=heatmap_breaks, 
                                                       option = 'B', 
                                                       direction = -1,
                                                       na.value='white'))    


# sensitivity analysis incorporating albumin-creatinine ratio >= 3 times----

acr3_id <- meta_freq %>%
    filter(id %in% dense3_id, albcreat_ratio >= 3) %>%
    pull(id)

acr3_dt <- meta_dt %>% 
    filter(id %in% acr3_id)

file_name <- str_c(approach, "acr3_dt.rds", sep = "_")
write_rds(acr3_dt, file_name)








