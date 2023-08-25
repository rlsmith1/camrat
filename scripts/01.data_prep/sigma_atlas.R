
###
### Reads in SIGMA atlas as nifti file and generates polygons for plotting in ggplot
###

# libraries ---------------------------------------------------------------

library(tidyverse)
library(purrr)
library(RNifti)
library(janitor)
library(sf)
library(plotly)

# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
atlas_dir <- paste0(base_dir, "../atlases/SIGMA/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/")

# ATLAS IMAGE DATA
img_path <-  paste0(atlas_dir, "SIGMA_Anatomical_Brain_Atlas.nii")
m_img_data <- asNifti(img_path)

# ATLAS LABELS & HIERARCHY
structures_path <- paste0(atlas_dir, "SIGMA_Anatomical_Brain_Atlas_ListOfStructures.csv")
abbrev_path <- paste0(base_dir, "data/SIGMA_abbreviations.xlsx")

# format labels dataframe for future use ----------------------------------

df_sigma_labels <- read_csv(structures_path) %>% 
  clean_names() %>% 
  left_join(readxl::read_xlsx(abbrev_path) %>% 
              clean_names %>% 
              dplyr::rename("roi_abbreviation" = "abbreviation"),
            by = join_by(region_of_interest, system, territories, matter)
  ) %>% 
  filter(!is.na(roi_abbreviation)) %>% 
  mutate_at(c("matter", "territories", "system", "region_of_interest"), 
            ~tolower(.x) %>% 
              str_replace_all(" ", "_") %>% 
              str_replace_all("-", "_")) %>% 
  pivot_longer(2:3, names_to = "hemi", values_to = "idx") %>% 
  mutate(hemisphere = str_replace(hemi, "_.*", "")) %>% 
  select(-hemi, -original_atlas) %>% 

  # remove multiple instance of fi
  filter(!(region_of_interest == "fimbria_of_the_hippocampus" & idx %in% c(941, 942))) %>% 
  distinct() %>% 
  
  # correct misspellings
  dplyr::mutate(system = str_remove(system, "_system"),
                system = ifelse(system == "hippocampus_fomation", "hippocampal_formation", system),
                region_of_interest = str_replace(region_of_interest, "reto", "retro")
  )

# save for later use
save(df_sigma_labels, file = paste0(base_dir, "objects/25Aug2023_sigma_labels.RDS"))

    
# data format -------------------------------------------------------------

# CONVERT 3D MATRIX TO 3 COLUMN DATAFRAME
df_img_data <- m_img_data %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "x")

df_img_data <- m_img_data %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "x") %>% 
  pivot_longer(contains("V"), names_to = "y", values_to = "z") %>% 
  mutate(y = str_remove(y, "V"))

expand_grid(m_img_data)

# CONVERT 3D MATRIX INTO 4 COLUMN DATAFRAME, WHERE VALUE IS ROI IDX
doParallel::registerDoParallel()
l_img_data <- 1:dim(m_img_data)[1] %>% 
  map( ~ m_img_data[.x, , ] %>% 
         as_tibble %>% 
         mutate(x = .x, .before = 1)
  )

l_img_data <- 1:length(l_img_data) %>% 
  map( ~  l_img_data[[.x]] %>%     
         mutate(y = row_number(), .before = 2) %>% 
         pivot_longer(3:ncol(.), names_to = "z", values_to = "idx") %>% 
         mutate(z = str_remove(z, "V") %>% as.numeric()
         )
  )

df_img_data <- l_img_data %>% 
  bind_rows() %>% 
  filter(idx != 0) # removes voxels that are not part of the brain

# COMBINE IMAGE DATA WITH ATLAS LABELS BY IDX
df_sigma <- df_sigma_labels %>% 
  left_join(df_img_data) %>% 
  distinct() %>% 
  filter(!is.na(x)) # removes CSF

# Create polygon atlas for plotting -----------------------------------------------------

# FIND BORDERS FOR ALL REGIONS ACROSS ALL PLANES AND CREATE CONVEX HULL AS 3D OUTLINE

# coronal
y_slices <- df_sigma %>% filter(!is.na(y)) %>% arrange(y) %>% pull(y) %>% unique
df_coronal_atlas <- y_slices %>% 
  map_df( ~ df_sigma %>%
            filter(y == .x) %>% 
            st_as_sf(coords = c("x", "z")) %>% 
            group_by(matter, territories, system, sys_abbreviation, region_of_interest, roi_abbreviation, hemisphere) %>%
            summarize(geometry = st_union(geometry)) %>%
            st_convex_hull() %>% 
            dplyr::mutate(side = "coronal",
                          slice = .x,
                          .before = 1
            )
  )

# sagittal
x_slices <- df_sigma %>% filter(!is.na(x)) %>% arrange(x) %>% pull(x) %>% unique
df_sagittal_atlas <- x_slices %>% 
  map_df( ~ df_sigma %>%
            filter(x == .x) %>% 
            st_as_sf(coords = c("y", "z")) %>% 
            group_by(matter, territories, system, sys_abbreviation, region_of_interest, roi_abbreviation, hemisphere) %>%
            summarize(geometry = st_union(geometry)) %>%
            st_convex_hull() %>% 
            dplyr::mutate(side = "sagittal",
                          slice = .x,
                          .before = 1
            )
  )

# axial
z_slices <- df_sigma %>% filter(!is.na(z)) %>% arrange(z) %>% pull(z) %>% unique
df_axial_atlas <- z_slices %>% 
  map_df( ~ df_sigma %>%
            filter(z == .x) %>% 
            st_as_sf(coords = c("x", "y")) %>% 
            group_by(matter, territories, system, sys_abbreviation, region_of_interest, roi_abbreviation, hemisphere) %>%
            summarize(geometry = st_union(geometry)) %>%
            st_convex_hull() %>% 
            dplyr::mutate(side = "axial",
                          slice = .x,
                          .before = 1
            )
  )

# COMBINE
df_sigma_atlas <- df_axial_atlas %>% 
  bind_rows(df_sagittal_atlas) %>% 
  bind_rows(df_coronal_atlas) %>% 
  na.omit()

# SAVE
save(df_sigma_atlas, file = paste0(base_dir, "objects/25Aug2023_sigma_atlas_for_plotting.RDS"))
    
    