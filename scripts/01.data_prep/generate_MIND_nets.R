
### LIBRARIES
library(tidyverse)

### DATA
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
load(paste0(base_dir, "objects/21July2023_df_data.RDS")) # df_data, generated in combine_scan_data.R

### READ IN AND FORMAT ALL MIND CSVs
sessions <- c("ses-PND020", "ses-PND035", "ses-PND063", "ses-PND300")
studies <- c("JWD", "EDA")

df_mind <- tibble()

for (study in studies) {
  
  print(study)
  study_dir <- paste0(base_dir, "data/MIND_files/", study, "/ROI/")
  
  for (ses in sessions) {
    
    print(ses)
    ses_dir <- paste0(study_dir, ses, "_MT/output/")
    subjects <- list.files(ses_dir) %>% str_remove(".csv")
    
    for(sub in subjects) {
      
      df_tmp <- read_csv(paste0(ses_dir, sub, ".csv")) %>% 
        mutate(R1 = colnames(.), .before = 1) %>% 
        pivot_longer(2:ncol(.),
                     names_to = "R2",
                     values_to = "weight") %>% 
        mutate(subject = sub, 
               timepoint = ses,
               weight = as.numeric(weight), 
               .before = 1)
      
      df_mind <- df_mind %>% bind_rows(df_tmp)
      
    }
    
  }
  
}

### FORMAT
df_mind <- df_mind %>% 
  mutate(subject = str_remove(subject, "sub-"),
         timepoint = str_remove(timepoint, "ses-PND") %>% 
           as.numeric %>% 
           as.factor) %>% 
  
  # add metadata
  left_join(df_data %>% 
              ungroup %>% 
              dplyr::select(subject, timepoint, sex, group, age, scan_date) %>% 
              distinct(),
            by = c("subject", "timepoint")) %>% 
  
  # add TBV
  left_join(df_tbv, by = c("subject", "timepoint")) %>% 
  
  # finalize formatting
  filter(!(R1 == "commissural_stria_terminalis" | R2 == "commissural_stria_terminalis")) %>% 
  dplyr::select(subject, timepoint, sex, group, scan_date, age, tbv, everything())  %>% 
  mutate(study = ifelse(substr(subject, start = 1, stop = 3) == "JWD", "MRC", "GSK"), .before = 1) %>% 
  mutate(group = ifelse(is.na(group), "control", group),
         sex = ifelse(study == "MRC", "Male", sex) )

### SAVE
save(df_mind, file = paste0(base_dir, "objects/18July2023_df_mind.RDS"))

