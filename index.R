library(dplyr)
# devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)
# source(utils.R)

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
# demographics_filt <- demographics_filt %>% dplyr::mutate(age_group = )

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
labs_cr_aki <- rename(labs_cr_aki, lab_time = CHARTTIME)
labs_cr_aki <- rename(labs_cr_aki, admission_time = ADMITTIME)
labs_cr_aki <- rename(labs_cr_aki, value = VALUE)

labs_cr_aki <- labs_cr_aki %>% dplyr::mutate(days_since_admission = as.numeric(as.Date(lab_time) - as.Date(admission_time)))
# labs_cr_aki <- labs_cr_aki %>% dplyr::mutate(value = ifelse(value != '', as.numeric(value), NA))
labs_cr_aki$value = as.numeric(labs_cr_aki$value)
labs_cr_aki <- labs_cr_aki[!is.na(labs_cr_aki$value),]

labs_cr_aki <- labs_cr_aki %>% dplyr::filter(days_since_admission >= -60)
# labs_cr_aki <- labs_cr_aki %>% select(patient_id, lab_time, days_since_admission, value)

# Generate separate demographics table for patients who do not have any sCr values fulfilling
# the above (e.g. all the labs are before t= -90days or patient has no sCr value)
# message("Removing patients who do not have any serum creatinine values...")
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

# Another outcome we are interested in is to look at acute kidney disease, AKD (in between AKI and CKD)
# We will use the definitions proposed for AKD as described by Chawla et. al. 2017 (ref (1))
# We are interested in renal recovery at the 7-day and 90-day timepoint
# We will use a cutoff of recovery to 1.25x baseline Cr as recovery, as used in ref (2)
# References:
# 1. Chawla, L., Bellomo, R., Bihorac, A. et al. Acute kidney disease and renal recovery: consensus report
#    of the Acute Disease Quality Initiative (ADQI) 16 Workgroup. Nat Rev Nephrol 13, 241-257 (2017).
#    https://doi.org/10.1038/nrneph.2017.2
# 2. Pannu, N., James, M., Hemmelgarn, B. & Klarenbach, S. Association between AKI, Recovery of Renal
#    Function, and Long-Term Outcomes after Hospital Discharge. Clinical Journal of the American Society of
#    Nephrology 8, 194-202 (2013). https://doi.org/10.2215/CJN.06480612

# Generate sCr levels at +7d (cr_7d) and +90d (cr_90d) timepoints (for determining post-AKI recovery, AKD)
labs_cr_aki <- data.table::setDT(labs_cr_aki)[,':='(cr_7d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+7,incbounds = TRUE)],1),cr_90d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+90,incbounds = TRUE)],1)),by=c('patient_id','days_since_admission')][]

# At this point, our table has these headers:
# patient_id  siteid  days_since_admission  value min_cr_90d min_cr_48h  min_cr_retro_7day min_cr_48h_retro  cr_7d cr_90d

# Now we have to start grading AKI severity at each time point
# This approach is similar to how the MIMIC-III dataset generates AKI severity
# Generate two columns using both the formal KDIGO AKI definition and the modified retrospective AKI definition
message("Generating KDIGO severity grades for each serum Cr value")
labs_cr_aki$aki_kdigo <- apply(labs_cr_aki,1,aki_kdigo_grade)
labs_cr_aki$aki_kdigo_retro <- apply(labs_cr_aki,1,aki_kdigo_grade_retro)
labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(aki_kdigo_final = max(aki_kdigo,aki_kdigo_retro))

# Generate two columns grading AKD severity at 7d and 90d (grade 0B/C is coded as 0.5)
labs_cr_aki$akd_7d <- apply(labs_cr_aki,1,akd_grade_7d)
labs_cr_aki$akd_90d <- apply(labs_cr_aki,1,akd_grade_90d)

# Now we are going to generate the start days of each AKI
labs_cr_aki_tmp <- labs_cr_aki
labs_cr_aki_tmp$valid = 1

# Find the day of the minimum Cr used for grading AKIs (taken as baseline)
message("Now finding the day at which the minimum serum Cr is achieved")
labs_cr_aki_tmp <- labs_cr_aki_tmp %>% dplyr::group_by(patient_id) %>% tidyr::complete(days_since_admission = tidyr::full_seq(days_since_admission,1)) %>% dplyr::mutate(value = zoo::na.fill(value,Inf))
labs_cr_aki_tmp2 <- labs_cr_aki_tmp
labs_cr_aki_tmp3 <- labs_cr_aki_tmp
labs_cr_aki_tmp2 <- labs_cr_aki_tmp2 %>% split(.$patient_id) %>% purrr::map(~pos_min(.$value,.$days_since_admission)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
colnames(labs_cr_aki_tmp2)[2] <- "day_min"
labs_cr_aki_tmp3 <- labs_cr_aki_tmp3 %>% split(.$patient_id) %>% purrr::map(~pos_min(.$value,.$days_since_admission,lag=FALSE)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
colnames(labs_cr_aki_tmp3)[2] <- "day_min_retro"
labs_cr_aki_tmp4 <- cbind(labs_cr_aki_tmp,"day_min" = labs_cr_aki_tmp2$day_min,"day_min_retro" = labs_cr_aki_tmp3$day_min_retro)
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[!is.na(labs_cr_aki_tmp4$valid),]

# Generate delta_cr
message("Now identifying maxima points of serum Cr")
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(min_cr_7d_final = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::mutate(delta_cr = value - min_cr_7d_final)

# Use the largest delta_cr to find the peak of each AKI
labs_cr_aki_delta_maxima <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id) %>% dplyr::filter(delta_cr %in% delta_cr[which.peaks(delta_cr,decreasing=FALSE)])
labs_cr_aki_delta_maxima$delta_is_max = 1
labs_cr_aki_delta_maxima <- labs_cr_aki_delta_maxima %>% dplyr::rename(delta_maxima = delta_cr) %>% dplyr::select(patient_id,days_since_admission,delta_maxima,delta_is_max)
labs_cr_aki_tmp4 <- merge(labs_cr_aki_tmp4,labs_cr_aki_delta_maxima,by=c("patient_id","days_since_admission"),all.x=TRUE)

# Generate a separate table (for reference) of all creatinine peaks not fulfilling KDIGO AKI criteria
message("Now generating a separate table for serum Cr peaks which do not reach AKI definitions")
labs_cr_nonaki <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final == 0,]
labs_cr_nonaki[is.na(labs_cr_nonaki)] <- 0
labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$delta_is_max > 0,]
labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)

# Filter for KDIGO grades > 0
message("Now generating tables of all AKI events")
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final > 0,]
labs_cr_aki_tmp4[is.na(labs_cr_aki_tmp4)] <- 0
# Filter for maxima of delta_cr (which should give us the peaks)
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$delta_is_max > 0,]

# Filter and reorder columns to generate our final table of all AKI events
labs_aki_summ <- labs_cr_aki_tmp4 %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)

labs_aki_summ <- labs_aki_summ %>% dplyr::distinct(patient_id,days_since_admission,.keep_all=TRUE)

# Final headers for labs_aki_summ:
# patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d
# days_since_admission - time at which peak Cr is achieved
# day_min - time at which Cr begins to rise

# Generate the highest Cr peak for non-AKI peaks detected
labs_nonaki_summ <- labs_cr_nonaki %>% dplyr::group_by(patient_id) %>% dplyr::slice(which.max(delta_cr))

# We also want to generate tables to determine (1) if AKIs occurred before/after severe disease onset
# (2) how long before/after disease severity
# These tables will help in segregating the populations for analysis later
message("Generating intermediate tables to determine if AKIs occured before or after severe COVID-19 onset")
severe_time <- demographics_filt %>% dplyr::select(patient_id,severe,time_to_severe)
labs_aki_severe <- merge(labs_aki_summ,severe_time,by="patient_id",all.x=TRUE)
labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = ifelse(!is.na(time_to_severe), time_to_severe - day_min,NA))
labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_before_aki = ifelse(severe_to_aki < 0,1,0))
# Final headers for labs_aki_severe:
# patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d  severe time_to_severe  severe_to_aki severe_before_aki

labs_nonaki_severe <- merge(labs_nonaki_summ,severe_time,by="patient_id",all.x=TRUE)
labs_nonaki_severe$severe_to_aki <- NA
labs_nonaki_severe$severe_before_aki <- 1

## Save the generated AKI tables for future reference / debugging (note: these will NOT be uploaded!!)
#write.csv(labs_aki_summ,"PatientAKIEvents.csv",row.names=FALSE)
#write.csv(labs_aki_severe,"PatientAKIEvents_Severe.csv",row.names=FALSE)


# =====================
# Demographics Table
# =====================
# message("Now generating the demographics table...")
# # demog_summ <- merge(demog_summ,comorbid,by="patient_id",all.x=TRUE)
# demog_summ <- demographics_filt %>% dplyr::select(patient_id,sex,age,race,severe,deceased,time_to_severe,time_to_death)
# demog_summ[is.na(demog_summ)] <- 0
# demog_summ$aki <- 0
# demog_summ$aki[demog_summ$patient_id %in% labs_aki_summ$patient_id] <- 1
# demog_summ$severe <- factor(demog_summ$severe,levels=c(0,1),labels=c("Non-severe","Severe"))
# demog_summ$deceased <- factor(demog_summ$deceased,levels=c(0,1),labels=c("Alive","Deceased"))
# demog_summ$aki <- factor(demog_summ$aki,levels=c(0,1),labels=c("No AKI","AKI"))
# # demog_summ[comorbid_list] <- lapply(demog_summ[comorbid_list],factor)
# # table_one_vars <- c("sex","age_group","race","severe","deceased","time_to_severe","time_to_death",comorbid_list)
# table_one_vars <- c("sex","age","race","severe","deceased","time_to_severe","time_to_death")
# table_one <- tableone::CreateTableOne(data=demog_summ,vars=table_one_vars,strata="aki")
# export_table_one <- print(table_one,showAllLevels=TRUE,formatOptions=list(big.mark=","))
# # write.csv(export_table_one,file=file.path('output', paste0(currSiteId, "_TableOne.csv")))
# write.csv(export_table_one,file=file.path('output', 'TableOne.csv'))
# # capture.output(summary(table_one),file=file.path('output', paste0(currSiteId, "_TableOne_Missingness.txt")))
# capture.output(summary(table_one),file=file.path('output', "TableOne_Missingness.txt"))
# message("TableOne with patient demographics should have been generated in CSV files at this point. Check for any errors.")
#
# demog_time_to_event_tmp <- demog_summ[,c("sex","age","race")]
# demog_time_to_event_tmp <- data.table::as.data.table(lapply(demog_time_to_event_tmp,factor))
# demog_time_to_event_tmp <- data.table::as.data.table(demog_time_to_event_tmp)[,sapply(demog_time_to_event_tmp,function(col) nlevels(col) > 1),with=FALSE]
# demog_list <- colnames(demog_time_to_event_tmp)
# demog_time_to_event <- demog_summ[,c("patient_id",demog_list)]
