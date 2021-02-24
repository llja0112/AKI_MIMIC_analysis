library(dplyr)

patients <- read.csv('data/patients_aki.csv')
admissions <- read.csv('data/admissions.csv')
intubations <- read.csv('data/intubation_events.csv')
dialysis <- read.csv('data/dialysis_events.csv')
labs <- read.csv('data/labs.csv')
icus <- read.csv('data/icus.csv')

icus <- distinct(icus, HADM_ID, .keep_all=TRUE)

patients <- select(patients, -ROW_ID)
admissions <- select(admissions, -ROW_ID)
icus <- select(icus, -ROW_ID)

demographics <- admissions %>% left_join(patients, by='SUBJECT_ID')
demographics <- demographics %>% left_join(icus, by=c('HADM_ID', 'SUBJECT_ID'))
demographics <- rename(demographics, patient_id = SUBJECT_ID)
# siteid
demographics <- rename(demographics, admission_date = ADMITTIME)
# days_since_admission
demographics <- rename(demographics, last_discharge_date = DISCHTIME)
demographics <- rename(demographics, severe_date = INTIME)
demographics <- rename(demographics, death_date = DOD)
# demographics <- rename(demographics, deceased = EXPIRE_FLAG)
demographics <- rename(demographics, sex = GENDER)
demographics <- rename(demographics, race = ETHNICITY)
# demographics <- rename(demographics, age = )

# still_in_hospital
demographics <- demographics %>% dplyr::mutate(severe = ifelse( (HADM_ID %in% icus$HADM_ID) , 1, 0))
demographics <- demographics %>% dplyr::mutate(severe = ifelse( (HADM_ID %in% icus$HADM_ID) , 1, 0))

demographics <- demographics %>% dplyr::select(patient_id,admission_date,last_discharge_date,severe_date,severe,death_date,deceased,sex,race)

demographics_filt <- demographics %>% dplyr::mutate(time_to_severe = ifelse(severe == 1, as.numeric(as.Date(severe_date) - as.Date(admission_date)),NA))
# demographics_filt <- demographics_filt %>% dplyr::mutate(death_date = ifelse(deceased == 1, gsub('T', ' ', death_date), NA))
# demographics_filt <- demographics_filt %>% dplyr::mutate(time_to_death = ifelse(deceased == 1, as.numeric(as.Date(gsub('T', ' ', death_date))) - as.Date(admission_date)),NA)
demographics_filt <- demographics_filt %>% dplyr::mutate(length_stay =as.numeric(as.Date(last_discharge_date) - as.Date(admission_date)))

# demographics <- demographics %>% dplyr::mutate(deceased = ifelse( death_date == '', 0, 1 ))
# demographics <- demographics %>% dplyr::mutate(severe_date = as.numeric( as.Date(INTIME) - as.Date(ADMITTIME) ) )

# dplyr::select(patient_id,siteid,admission_date,days_since_admission,last_discharge_date,still_in_hospital,severe_date,severe,death_date,deceased,sex,age_group,race,race_collected)
