library(dplyr)

patients <- read.csv('data/patients_aki.csv', stringsAsFactors = FALSE)
admissions <- read.csv('data/admissions.csv', stringsAsFactors = FALSE)
intubations <- read.csv('data/intubation_events.csv', stringsAsFactors = FALSE)
dialysis <- read.csv('data/dialysis_events.csv', stringsAsFactors = FALSE)
labs <- read.csv('data/labs.csv', stringsAsFactors = FALSE)
icus <- read.csv('data/icus.csv', stringsAsFactors = FALSE)

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
demographics <- rename(demographics, date_of_birth = DOB)
demographics <- rename(demographics, admission_id = HADM_ID)

# still_in_hospital exlcuded because all patinet either have been discharged or are deceased
demographics <- demographics %>% dplyr::mutate(severe = ifelse( (admission_id %in% icus$HADM_ID) , 1, 0))
demographics <- demographics %>% dplyr::mutate(deceased = ifelse( death_date == '' , 0, 1))

demographics <- demographics %>% dplyr::select(patient_id,admission_id,admission_date,last_discharge_date,severe_date,severe,death_date,deceased,sex,race,date_of_birth)

demographics_filt <- demographics %>% dplyr::mutate(time_to_severe = ifelse(severe == 1, as.numeric(as.Date(severe_date) - as.Date(admission_date)),NA))
demographics_filt <- demographics_filt %>% dplyr::mutate(death_date = ifelse(deceased == 1, gsub(' ', 'T', death_date), NA))
demographics_filt <- demographics_filt %>% dplyr::mutate(time_to_death = ifelse(deceased == 1, as.numeric(as.Date(death_date) - as.Date(admission_date)), NA) )
demographics_filt <- demographics_filt %>% dplyr::mutate(length_stay =as.numeric(as.Date(last_discharge_date, format='%Y-%m-%d') - as.Date(admission_date)))
demographics_filt <- demographics_filt %>% dplyr::mutate(age =as.numeric(as.Date(last_discharge_date, format='%Y-%m-%d') - as.Date(date_of_birth))/365)

# Reorder the columns to be more readable
demographics_filt <- demographics_filt %>% dplyr::select(patient_id,sex,age,race,length_stay,severe,time_to_severe,deceased,time_to_death)

# Cormorbids not included for now

message("Creating table for intubation...")

intubation <- intubations %>% left_join(admissions, by='HADM_ID')
intubation <- rename(intubation, patient_id = SUBJECT_ID.x)
intubation <- rename(intubation, admission_id = HADM_ID)
intubation <- rename(intubation, intubation_time = STARTTIME)
intubation <- rename(intubation, admission_time = ADMITTIME)

intubation <- intubation %>% dplyr::mutate(days_since_admission = as.numeric(as.Date(intubation_time) - as.Date(admission_time)))

intubation <- intubation %>% dplyr::select(patient_id, admission_id, days_since_admission)
intubation <- intubation[order(intubation$admission_id,intubation$days_since_admission),]
intubation <- intubation[!duplicated(intubation$admission_id),]


message("Creating table for RRT...")

rrt <- dialysis %>% left_join(admissions, by='HADM_ID')
rrt <- rename(rrt, patient_id = SUBJECT_ID.x)
rrt <- rename(rrt, admission_id = HADM_ID)
rrt <- rename(rrt, dialysis_time = STARTTIME)
rrt <- rename(rrt, admission_time = ADMITTIME)

rrt <- rrt %>% dplyr::mutate(days_since_admission = as.numeric(as.Date(dialysis_time) - as.Date(admission_time)))

rrt <- rrt %>% dplyr::select(patient_id, admission_id, days_since_admission)
rrt <- rrt[order(rrt$admission_id, rrt$days_since_admission),]
rrt <- rrt[!duplicated(rrt$admission_id),]

# ====================
# PART 2: AKI Detection Code
# ====================
# For the purposes of generating a table of AKI events, we will have to create a new data.table
# with the serum creatinine levels.
# We will then use this table, labs_cr_aki, to generate a summary table containing details of 
# each AKI event and the post-AKI recovery labs
message("Now proceeding to AKI detection code:")
message("Extracting serum creatinine values...")

labs_cr_aki <- labs %>% left_join(admissions, by='HADM_ID')
labs_cr_aki <- rename(labs_cr_aki, patient_id = SUBJECT_ID.x)
labs_cr_aki <- rename(labs_cr_aki, admission_id = HADM_ID)
labs_cr_aki <- rename(labs_cr_aki, lab_time = STARTTIME)
labs_cr_aki <- rename(labs_cr_aki, admission_time = ADMITTIME)

labs_cr_aki <- labs_cr_aki %>% dplyr::mutate(days_since_admission = as.numeric(as.Date(lab_time) - as.Date(admission_time)))

labs_cr_aki <- labs_cr_aki %>% dplyr::filter(days_since_admission >= -60)

# Generate separate demographics table for patients who do not have any sCr values fulfilling 
# the above (e.g. all the labs are before t= -90days or patient has no sCr value)
message("Removing patients who do not have any serum creatinine values...")
pts_valid_cr <- unique(labs_cr_aki$patient_id)
demog_no_cr <- demographics_filt[!(demographics_filt$patient_id %in% pts_valid_cr),]
demographics_filt <- demographics_filt[demographics_filt$patient_id %in% pts_valid_cr,]
# There are two possible scenarios which we have to consider when detecting each AKI event:
# (1) AKI occurs after admission
#   - easy to detect with the formal KDIGO definition as we only need to use older data points
#   - we can use an approach similar to what is used in the MIMIC-III data-set: use a rolling
#     time frame and detect the lowest Cr value to use as our baseline
# (2) Patient presents with an AKI at the time of admission 
#   - this makes it harder for us to determine what is the true baseline especially with limited
#     longitudinal data
#   - Hence one way to get around this is to generate a "retrospective" baseline (i.e. look into
#     future Cr values and find the minimum in a rolling timeframe) and use this as a surrogate
#     baseline
#   - Such a workaround is the only feasible way of dealing with missing retrospective data though
#     it is likely that we may miss quite a number of true AKIs using this method
