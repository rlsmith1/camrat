###
### Convert scan format as it is outputted from the scanner to BIDS
###

# INSTALL TIDYVERSE IF IT IS NOT INSTALLED
install.packages(setdiff(c("tidyverse"), rownames(installed.packages())))

# LOAD LIBRARIES
library(tidyverse)

# SET DIRECTORIES
orig_dir <- "/rds/project/rds-83zLpPieDV8/q10014/nii"
out_dir <- "/rds/project/rds-cAQcxgoLHoA/Livia"
data_dir <- paste0(out_dir, "/data")

# READ IN SCAN LIST AS TIBBLE
df_scans <- read_csv(paste0(out_dir, "/scan_list.txt"))

# change column names to id & date
df_scans <- df_scans %>% 
  bind_rows(tibble(colnames(.)[1], 
                   colnames(.)[2]) %>% 
              dplyr::rename_all(~colnames(df_scans)), 
            .) %>% 
  dplyr::rename_all(~c("id", "date"))

# ASSIGN SUBJECT AND SESSION NUMBER
df_scans <- df_scans %>% 
  group_by(id) %>% 
  mutate(sub = paste0("sub", id) %>% str_remove("B3526"),
         ses = paste0("ses-T", row_number())
  )

# LOOP OVER SUBJECTS AND SESSIONS (each line in df_scans) TO COPY SCANS TO OUR DIRECTORY IN BIDS FORMAT
for (i in 1:nrow(df_scans)) {
  
  # identify subject & session number
  id <- df_scans[i, ] %>% pull(id)
  date <- df_scans[i, ] %>% pull(date)
  sub <- df_scans[i, ] %>% pull(sub)
  ses <- df_scans[i, ] %>% pull(ses)
  print(paste0('Copying: ', sub, '; ' , ses))
  
  # set up BIDS folder structure in our data dir
  dir.create(paste0(data_dir, "/", sub, '/', ses, '/anat/'), recursive = T)  
  #dir.create(paste0(out_dir, sub, '/', ses, '/func/'), recursive = T) # no functional scans for ex-vivo
  
  current_data_dir <- paste0(data_dir, "/", sub, '/', ses, '/anat/')
  
  # copy over MT, T1, and PD files
  
    # MT
    print("Copying MT scan")
    mt_old <- paste0(orig_dir, "/", id, "/", date, '/Series_000_/mt.nii')
    mt_new <- paste0(current_data_dir, "/", sub, "_", ses, "_MTR.nii")
    file.copy(mt_old, mt_new) %>% try()
    
    # T1
    print("Copying T1 scan")
    t1_old <- paste0(orig_dir, "/", id, "/", date, '/Series_000_/t1.nii')
    t1_new <- paste0(current_data_dir, "/", sub, "_", ses, "_T1w.nii")
    file.copy(t1_old, t1_new) %>% try()
    
    # PD
    print("Copying PD scan")
    pd_old <- paste0(orig_dir, "/", id, "/", date, '/Series_000_/pd.nii')
    pd_new <- paste0(current_data_dir, "/", sub, "_", ses, "_PDw.nii")
    file.copy(pd_old, pd_new) %>% try()
  
}

print("finished!")

# CREATE SCAN TABLE FOR DATA PREP AND REGISTRATION
df_scan_table <- df_scans %>% 
  expand_grid(contrast = c("MTR", "PDw", "T1w")) %>% 
  mutate(init_scale = 1.0) %>% 
  dplyr::select(-id, -date)

write.table(df_scan_table, paste0(out_dir, "/scan_table_for_registration.txt"),
            sep = ",",  col.names = FALSE, row.names = FALSE, quote = FALSE)

