
###
### Calculate total brain volume for each scan by counting non-zero mask voxels. resolution 0.16 mm^3
###

# libraries ---------------------------------------------------------------

library(tidyverse)
library(RNifti)
library(janitor)

# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

# calculate eTIV ----------------------------------------------------------

df_tbv <- tibble()  
studies <- c("JWD", "EDA")
sessions <- c("ses-PND020", "ses-PND035", "ses-PND063", "ses-PND300")

for (study in studies) {
  
  print(study)
  
  for (ses in sessions) {
    
    print(paste0("session ", ses))
    mask_dir <- paste0(base_dir, "data/MT_scans/", study, "/", ses, "/mask/")
    
    # GET LIST OF SUBJECTS
    subjects <- list.files(mask_dir) %>% 
      str_remove("SIGMA_Anatomical_Brain_Atlas_in_") %>% 
      str_split("_") %>% 
      map(., 1) %>% 
      unlist()
    
    df_tmp <- tibble()
    
    for (sub in subjects) {
      
      print(paste0("converting ", sub))
      
      # read in atlas in subject space
      mask <- paste0(mask_dir, sub, "_", ses, "_MTR_deobliqued_RIP_shft_mask.nii")
      mask_data <- asNifti(mask)
      
      # combine and add labels
      df_combos <- expand_grid(
        
        x = 1:dim(mask_data)[1],
        y = 1:dim(mask_data)[2]
        
      )
      
      df_sub <- map2_dfr(
        
        .x = df_combos %>% pull(x),
        .y = df_combos %>% pull(y),
        .f = ~ tibble(
          
          subject = sub,
          session = ses,
          idx = mask_data[.x, .y, ]
          
        )
        
      ) %>% 
        group_by(subject, session) %>% 
        summarise(tbv = sum(idx) * (0.16^3))
      
      df_tmp <- df_tmp %>% bind_rows(df_sub)

    }
    
    df_tbv <- df_tbv %>% bind_rows(df_tmp) 
    
  }
  
  
}

df_tbv <- df_tbv %>%
  mutate(subject = str_remove(subject, "sub-"),
         timepoint = str_remove(session, "ses-PND") %>% 
           as.numeric %>% 
           as.factor
         ) %>% 
  dplyr::select(subject, timepoint, tbv)
  
save(df_tbv, file = paste0(base_dir, "objects/18July2023_df_tbv.RDS"))  
    
