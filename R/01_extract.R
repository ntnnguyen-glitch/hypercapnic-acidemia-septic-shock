# Extract the data used in the study of hypercapnic acidemia in septic shock from the raw MIMIC-IV v3.1 files.
source("R/00_setup.R")
con <- dbConnect(duckdb::duckdb(), file.path(WORK_DIR, "extract.duckdb"))
X <- function(sql) { t0 <- Sys.time(); dbExecute(con, sql)
  message(round(difftime(Sys.time(), t0, units = "mins"), 1), " min: ", substr(gsub("\\s+", " ", sql), 1, 60)) }
X("PRAGMA threads=8")
rm_ <- function(p) sprintf("read_csv('%s/%s.csv.gz', header=true, auto_detect=true, sample_size=-1)", MIMIC_DIR, p)

## MIMIC-IV: first ICU stay of each adult
X(paste0("create or replace table st as select i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, p.gender,
  p.anchor_age + (year(i.intime) - p.anchor_year) age, p.dod, a.deathtime, a.dischtime, a.hospital_expire_flag
  from ", rm_("icu/icustays"), " i join ", rm_("hosp/patients"), " p using(subject_id) join ", rm_("hosp/admissions"), " a using(hadm_id)
  qualify row_number() over(partition by i.subject_id order by i.intime) = 1"))
X("delete from st where age < 18")
X(paste0("create or replace table dx as select d.hadm_id, d.icd_code from ", rm_("hosp/diagnoses_icd"), " d join st using(hadm_id)"))
# sepsis (explicit ICD-9-CM 995.91, 995.92, 785.52; ICD-10-CM A40, A41, R65.2)
X("create or replace table sep as select distinct hadm_id from dx where regexp_matches(icd_code, '^(99591|99592|78552|A40|A41|R652)')")
X("create or replace table sta as select * from st where hadm_id in (select hadm_id from sep)")
X(paste0("create or replace table lab as select l.subject_id, s.stay_id, l.specimen_id, l.itemid, l.charttime, l.valuenum v, l.value from ",
  rm_("hosp/labevents"), " l join sta s on l.subject_id = s.subject_id where
  (l.itemid = 50912 and l.charttime between s.intime - interval 7 day and s.intime + interval 30 day) or
  (l.itemid in (50818,50820,50803,50882,50813,50862,51006,50971,50983,50902,50931,51222,51301,51265,50868,50821,50885,52033,50816)
   and l.charttime between s.intime - interval 24 hour and s.intime + interval 30 hour)"))
X(paste0("create or replace table ce as select c.stay_id, c.itemid, c.charttime, c.valuenum v from ", rm_("icu/chartevents"),
  " c join sta s using(stay_id) where c.itemid in (226512,224639,220045,220210,224690,220052,220181,225312,220277,223762,223761,220739,223900,223901,223835)
  and c.charttime between s.intime - interval 6 hour and s.intime + interval 24 hour"))
X(paste0("create or replace table pv as select p.stay_id, p.itemid, p.starttime, p.endtime from ", rm_("icu/procedureevents"),
  " p join sta using(stay_id) where p.itemid in (225792, 225794)"))
X(paste0("create or replace table inp as select n.stay_id, n.itemid, n.starttime, n.endtime from ", rm_("icu/inputevents"),
  " n join sta s using(stay_id) where n.itemid in (221906,221289,221662,222315,221749)"))

dbDisconnect(con, shutdown = TRUE)
