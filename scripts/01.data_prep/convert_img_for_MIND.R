



# libraries ---------------------------------------------------------------

library(tidyverse)
library(RNifti)
library(janitor)

# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

# ATLAS LABELS AND HIERARCHY
load(paste0(base_dir, "objects/25Aug2023_sigma_labels.RDS")) # df_sigma_labels

# generate matrix for MIND and export to csv ------------------------------

# SCAN DIRECTORIES    
study <- "MRC" # GSK
ses <- "ses-PND020" # ses-PND035, ses-PND063, ses-PND300

mset_dir <- paste0(base_dir, "data/MT_scans/", study, "/", ses, "/mset/")
infile_dir <- paste0(base_dir, "data/MT_scans/", study, "/", ses, "/infile/")

# GET LIST OF SUBJECTS TO CONVERT IMAGES TO DFS FOR
    
subjects <- list.files(mset_dir) %>% 
  str_remove("SIGMA_Anatomical_Brain_Atlas_in_") %>% 
  str_split("_") %>% 
  map(., 1) %>% 
  unlist()

# Don't run EDAA 32 for now (PND300)
#subjects <- subjects[!str_detect(subjects, "32")]

# READ ATLAS FOR IMAGE, COMBINE WITH MT VALUES FOR IMAGE, EXPORT
for (sub in subjects) {
  
  print(paste0("converting ", sub))
  
  # read in atlas in subject space
  mset <- paste0(mset_dir, "SIGMA_Anatomical_Brain_Atlas_in_", sub, "_", ses, "_MTR_deobliqued_RIP_shft.nii")
  mset_data <- asNifti(mset)
  
  # read in MT subject scan
  infile <- paste0(infile_dir, sub, "_", ses, "_MTR_deobliqued_RIP_shft.nii")
  infile_data <- asNifti(infile)
  
  # combine and add labels
  df_combos <- expand_grid(
    
    x = 1:dim(mset_data)[1],
    y = 1:dim(mset_data)[2]
    
  )
  
  df_mind <- map2_dfr(
    
    .x = df_combos %>% pull(x),
    .y = df_combos %>% pull(y),
    .f = ~ tibble(idx = mset_data[.x, .y, ],
                  MT = infile_data[.x, .y, ])
    
  ) %>% 
    
    filter(idx != 0) %>% 
    left_join(df_sigma_labels %>% 
                mutate(idx = as.numeric(idx)), 
              by = "idx") %>% 
    mutate(Label = region_of_interest) %>% 
    dplyr::select(Label, MT) %>% 
    arrange(Label) %>% 
    
    # REMOVE REGIONS THAT DIDN'T REGISTER WELL
    filter(!str_detect(Label, "spinal_cord|brainstem|cerebellum|commissural_stria_terminalis"))
  
  # remove features that don't have at least 2 voxels
  low_voxel_feature <- df_mind %>% 
    count(Label) %>% 
    filter(n < 2) %>% 
    pull(Label)
  
  df_mind <- df_mind %>% filter(!(Label %in% low_voxel_feature))
  
  # export as csv
  write.csv(df_mind,
            file = paste0(base_dir, "data/MIND_files/", study, "/ROI/", ses, "_MT/input/", sub, ".csv"),
            row.names = FALSE)
  
}

    