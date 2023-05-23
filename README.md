
# Dynamic Survival Prediction of ESKD 

This repository provides scripts for the article *"Dynamic survival prediction of end-stage kidney disease using random survival forest for competing risk analysis"* by Christiadi et al.

## How to create dynamic prediction plot

![Dynamic_plot Screenshot]()

#### Download the models and scripts into in the same folder with dataset

```http
  xxx
```

#### Prepare the dataset in person-period format as below

|  id | test | relyear | value | sex | age_init | gn_cat | time_compete | status_compete |
| :-- | :----| :-------| ----- | :-- | :------- | :----- | :----------- | :------------- |
| 1 | albumin | 0 | 44 | Female | 47 | lupus | 12.7 | no |
| 1 | bicarb | 0 | 28 | Female | 47 | lupus | 12.7 | no |
| 1 | chloride | 0 | 98 | Female | 47 | lupus | 12.7 | no |
| 1 | egfr | 0 | 53.6 | Female | 47 | lupus | 12.7 | no |
| 1 | haemoglobin | 0 | 123 | Female | 47 | lupus | 12.7 | no |
| 1 | albumin | 0.33 | 42 | Female | 47 | lupus | 12.7 | no |
| 1 | bicarb | 0.33 | 28 | Female | 47 | lupus | 12.7 | no |
| 1 | chloride | 0.33 | 100 | Female | 47 | lupus | 12.7 | no |
| 1 | egfr | 0.33 | 62 | Female | 47 | lupus | 12.7 | no |
| 1 | haemoglobin | 0.33 | 127 | Female | 47 | lupus | 12.7 | no |
| 1 | albumin | 0.5 | 47 | Female | 47 | lupus | 12.7 | no |
| 1 | bicarb | 0.5 | 25 | Female | 47 | lupus | 12.7 | no |
| 1 | chloride | 0.5 | 101 | Female | 47 | lupus | 12.7 | no |
| 1 | egfr | 0.5 | 51 | Female | 47 | lupus | 12.7 | no |
| 1 | haemoglobin | 0.5 | 123 | Female | 47 | lupus | 12.7 | no |

#### Parameters 

| Parameter | Type     | Description                       |
| :-------- | :------- | :-------------------------------- |
| `id`      | `num` | **Required**. Patient's ID |
| `test`    | `chr` | **Required**. Clinicopathological Data: `albumin`, `bicarb`, `chloride`, `egfr`, `haemoglobin` |
| `relyear` | `num` | **Required**. Relative time in years from time 0 when the respective Clinicopathological tests performed |
| `value` | `num` | **Required**. Values for the respective Clinicopathological tests |
| `sex` | `chr` | **Required**. Patient's gender: `Male`, `Female` |
| `age_init` | `num` | **Required**. Age at time 0 |
| `gn_cat` | `chr` | **Required**. GN or vasculitis diagnosis: `anca`, `fsgs`, `igan`, `lupus`, `mcd`, `mn`, `mp_mcgn`, `unspecific`, `NA` |
| `time_compete` | `num` | **Required**. Time to event in years |
| `status_compete` | `chr` | **Required**. Event: `death`, `eskd`, `no` |

#### Run the `dynamic_plot.R` and change the name of external_filename into the file name that contains the patient for analysis. Ensure only one ID in the dataset as the script only run analysis for one patient only. 

## Dependencies

* data.table 1.14.8
* dynpred 0.1.2
* furrr 0.3.1
* gtsummary 1.6.2
* heatmaply 1.4.2
* janitor 2.2.0
* lubridate 1.9.2
* mice 3.14.0 
* nlme 3.1-162
* patchwork 1.1.2
* pec 2022.05.04
* prodlim 2019.11.13
* purrr 1.0.1
* randomForestSRC 3.2.1
* readxl 1.4.2
* rsample 1.1.1
* skimr 2.1.4
* tidyverse 2.0.0
* boot 1.3-28.1

## License

[MIT](https://choosealicense.com/licenses/mit/)

