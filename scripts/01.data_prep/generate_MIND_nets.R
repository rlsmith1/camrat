
### LIBRARIES
library(tidyverse)

### DATA
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
load(paste0(base_dir, "objects/25Aug2023_df_data.RDS")) # df_data, generated in combine_scan_data.R

### READ IN AND FORMAT ALL MIND CSVs
sessions <- c("ses-PND020", "ses-PND035", "ses-PND063", "ses-PND300")
studies <- c("MRC", "GSK")

df_mind_GM <- tibble()

for (study in studies) {
  
  print(study)
  study_dir <- paste0(base_dir, "data/MIND_files/", study, "_GM/ROI/")
  
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
        mutate(study = study,
               subject = sub, 
               timepoint = ses,
               weight = as.numeric(weight), 
               .before = 1)
      
      df_mind_GM <- df_mind_GM %>% bind_rows(df_tmp)
      
    }
    
  }
  
}

### FORMAT
df_mind_GM <- df_mind_GM %>% 
  mutate(subject = str_remove(subject, "sub-"),
         timepoint = str_remove(timepoint, "ses-PND") %>% 
           as.numeric %>% 
           as.factor) %>% 
  
  # add metadata
  left_join(df_data %>% 
              ungroup %>% 
              dplyr::select(subject, timepoint, sex, group, age, scan_date, tbv) %>% 
              distinct(),
            by = c("subject", "timepoint")) %>% 

  # finalize formatting
  dplyr::select(subject, timepoint, sex, group, scan_date, age, tbv, everything())  %>% 
  mutate(group = ifelse(is.na(group), "control", group),
         sex = ifelse(study == "MRC", "Male", sex) 
  )

### SAVE
save(df_mind_GM, file = paste0(base_dir, "objects/29Aug2023_df_mind_GM.RDS")) # df_mind_GM

