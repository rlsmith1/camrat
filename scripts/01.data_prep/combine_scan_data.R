
###
### Pull together 3D ROI stats output from AFNI and add metadata
###

# libraries -------------------------------------------------------------------

library(tidyverse)
library(janitor)

# load external data -------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
load(paste0(base_dir, "objects/25Aug2023_df_tbv.RDS")) # df_tbv (generated in calc_tbv_from_masks.R)
load(paste0(base_dir, "objects/25Aug2023_sigma_labels.RDS")) # df_sigma_labels (generated in sigma_atlas.R)
load(paste0(base_dir, "objects/ratdb_scans.Rdata")) # df_scans (generated in ratdb.R)

# subject metadata
df_eda_subject_info <- read_csv(paste0(base_dir, "data/subject_info.csv")) %>% 
  na.omit() %>% 
  clean_names %>% 
  dplyr::rename("sex" = "gender")

# info on if scan had gadolinium
df_gad <- read_xlsx(paste0(base_dir, "data/gadolinium_use_records.xlsx")) %>% 
  dplyr::select(subject, scan_date, gad_by_record)

# load MRI scan data (output from 3dROIstats) -------------------------------------------------------------------

data_dir <- paste0(base_dir, "data/3d_ROI_stats/")
data_files <- list.files(data_dir)

df_data_tmp <- map_dfr(
  
  .x = 1:length(data_files),
  .f = ~ read_table(paste0(data_dir, data_files[.x])) %>% 
    dplyr::select(contains("Mean"), contains("Volume")) %>% 
    pivot_longer(1:ncol(.), names_to = "metric", values_to = "value") %>% 
    separate(metric, into = c("metric", "idx"), sep = "_") %>% 
    mutate(subject = str_split(data_files[.x], pattern = "_") %>% 
             unlist %>% 
             .[[1]] %>% 
             str_remove("sub-"),
           timepoint = str_split(data_files[.x], pattern = "_") %>% 
             unlist %>% 
             .[[2]] %>% 
             str_remove("ses-PND") %>% 
             as.numeric %>% 
             as.factor,
           scan_type = str_split(data_files[.x], pattern = "_") %>% 
             unlist %>% 
             .[[3]],
           .before = 1
    )

) %>% 
  mutate(metric = ifelse(metric == "Mean", scan_type, paste0(scan_type, "_", metric))) %>% 
  dplyr::select(-scan_type)

  
# combine with metadata -------------------------------------------------------------------

df_data <- df_data_tmp %>% 
  
  # Add ROI information
  left_join(df_sigma_labels %>% mutate(idx = as.character(idx)), by = "idx") %>% 
  dplyr::select(-idx) %>% 
  
  # add subject information
  mutate(study = ifelse(substr(subject, start = 1, stop = 3) == "JWD", "MRC", "GSK"),
         id = str_remove(subject, "EDAA") %>% as.numeric) %>% 
  left_join(df_eda_subject_info) %>% 
  mutate(sex = ifelse(study == "MRC", "Male", sex),
         group = ifelse(study == "MRC", "external_control", group)) %>% 
  dplyr::select(-id) %>% 
  
  # add TBV
  left_join(df_tbv, by = join_by(subject, timepoint)) %>% 
  
  # add scan information
  left_join(df_scans %>% 
              mutate(pattern = ifelse(str_detect(ids, "EDA"), ".*_", ".*-"),
                     subject = str_remove(ids, pattern),
                     age = dages) %>% 
              dplyr::rename("scan_date" = "dates") %>% 
              dplyr::select(subject, age, scan_date) %>% 
              mutate(timepoint = cut(age, 
                                     breaks = c(0, 25, 35, 70, 310), 
                                     labels = c("20", "35", "63", "300")) 
              ),
            by = join_by(subject, timepoint)
  ) %>% 
  
  # reorder columns
  dplyr::select(study, subject, timepoint, sex, group, age, scan_date, tbv, hemisphere, everything()) %>% 
  distinct() %>% 
  arrange(study, subject, timepoint, metric, region_of_interest)

# format & save -----------------------------------------------------------

# LIST ROIS THAT DID NOT REGISTER WELL TO REMOVE FROM ANALYSIS
bad_reg_rois <- c("spinal_cord", "brainstem", "cerebell", "commissural_stria_terminalis", "central_canal")

# REMOVE SCANS AND REGIONS THAT WERE NOT SUCCESSFULLY REGISTERED 
df_data <- df_data %>% 
  
  # remove subjects that didn't register well
  filter(!(subject == "EDAA32" & timepoint == 300) &
           !(subject == "EDAA46" & timepoint == 63)) %>% 
  
  # remove regions that didn't register well
  filter(!str_detect(region_of_interest, paste0(bad_reg_rois, collapse = "|"))) %>% 
  
  # remove scans with gadolinium stain
  anti_join(df_gad %>% filter(gad_by_record == 1), 
            by = join_by(subject, scan_date)) %>% 
  
  # make timepoint a factor variable
  mutate(timepoint = factor(timepoint, levels = c(20, 35, 63, 300)))

# SAVE
save(df_data, file = paste0(base_dir, "objects/25Aug2023_df_data.RDS"))

