

# libraries ---------------------------------------------------------------

    library(tidyverse)
    library(readxl)



# original data (ratdb) --------------------------------------------------------------------

    base_dir <- "~/OneDrive - National Institutes of Health/projects/CamRat/camrat/"
    load(paste0(base_dir, "objects/ratdb_scans.Rdata"))

    

# create QC spreadsheet ---------------------------------------------------


    df_scans %>% 
        dplyr::select(ids, did, type, fmri, dages, dates) %>% 
        mutate(contrast = case_when(
            
            type == "m" ~ "MTR",
            type == "p" ~ "PDw",
            type == "t" ~ "T1w",
            fmri == 1 ~ "fMRI"
            
        ), .before = 3
        
        ) %>% 
        filter(!is.na(contrast)) %>% 
        dplyr::select(-type, -fmri) %>% 
        distinct() %>% 
        mutate(study = case_when(
            
            grepl("JWD", ids) ~ "JWD",
            grepl("EDA", ids) ~ "EDA"
            
        ), .before = 1
        
        ) %>% 
        dplyr::select(-ids) %>% 
        mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310))) %>% 
        mutate(timepoint = case_when(
            
            timepoint == 1 ~ 20,
            timepoint == 2 ~ 35,
            timepoint == 3 ~ 63,
            timepoint == 4 ~ 300
            
        )) %>% 
        mutate(timepoint = factor(timepoint, levels = c("20", "35", "63", "300"))) %>% 
        filter(!is.na(timepoint)) %>% 
        dplyr::rename("id" = "did", "age" = "dages", "scan_date" = "dates") %>% 
        select(study, id, timepoint, age, scan_date, contrast) %>% 
        arrange(study, id, timepoint, age, contrast) %>% 
        mutate(init_scale = case_when(
            
            timepoint == 20 ~ 0.75,
            timepoint == 35 ~ 0.90,
            timepoint == 63 ~ 1.00,
            timepoint == 300 ~ 1.00
            
        )) #%>% 
        
        #write.csv("outputs/20230124_scan_registration_qc.csv", row.names = FALSE)



# ROUND 1 -----------------------------------------------


    df_qc1 <- read_excel(paste0(base_dir, "scan_registration_qc (1).xlsx")) %>% 
        clean_names
    
    df_qc1 %>% 
        dplyr::select(-age, -scan_date) %>% 
        
        # fMRI doesn't go through animal warper
        filter(contrast != "fMRI") %>% 
        
        # worry about this later
        filter(!(timepoint == 20 & contrast == "T1w")) %>% 
        
        # don't rerun scans that have succeeded
        filter(registration_success != "y" | is.na(registration_success)) %>% # print(n=nrow(.))
        
        # change init_scale based on conditions from QC inspections
        mutate(init_scale = case_when(
            
            timepoint == 20 & contrast == "PDw" & !is.na(other_notes) ~ 0.72,
            timepoint == 20 & contrast == "PDw" & is.na(other_notes) ~ 0.78,
            timepoint == 35 & contrast == "PDw"~ 0.88,
            timepoint == 20 & contrast == "MTR" ~ 0.78,
            timepoint == 35 & contrast == "T1w" ~ 0.85,
            timepoint == 300 & init_scale == 1 ~ 0.95,
            TRUE ~ init_scale

        )) %>% 
        
        # format for upload to HPC
        mutate(ses = ifelse(timepoint < 100, 
                            paste0("ses-PND0", timepoint), 
                            paste0("ses-PND", timepoint)
        ),
        sub = case_when(
            
            study == "EDA" & id < 10 ~ paste0("sub-", study, "A0", id),
            study == "EDA" & id > 10 ~ paste0("sub-", study, "A", id),
            study == "JWD" & id < 10 ~ paste0("sub-", study, "0", id),
            study == "JWD" & id > 10 ~ paste0("sub-", study, id)
            
        )
        ) %>% 
        dplyr::select(sub, ses, contrast, init_scale) %>% 
        
        write.csv(paste0(base_dir, "scan_list2.csv"), row.names = FALSE)
    
    

# ROUND 2 -----------------------------------------------------------------

# forgot to take out ok-to-exist flag, rerun scans that had been run...
    
    df_qc2 <- read_csv(paste0(base_dir, "scan_list2.csv"))

    df_qc1 %>% 
        filter(!(registration_success %in% c("y", "?")) & !is.na(registration_success)) %>% 

        # format for upload to HPC
        mutate(ses = ifelse(timepoint < 100, 
                            paste0("ses-PND0", timepoint), 
                            paste0("ses-PND", timepoint)
        ),
        sub = case_when(
            
            study == "EDA" & id < 10 ~ paste0("sub-", study, "A0", id),
            study == "EDA" & id > 10 ~ paste0("sub-", study, "A", id),
            study == "JWD" & id < 10 ~ paste0("sub-", study, "0", id),
            study == "JWD" & id > 10 ~ paste0("sub-", study, id)
            
        )
        ) %>% 
        dplyr::select(sub, ses, contrast) %>% 
        
        inner_join(df_qc2) %>% 
        
        write.csv(paste0(base_dir, "scan_list3.csv"), row.names = FALSE)

        

# ROUND 3 -----------------------------------------------------------------

    df_qc3 <- read_excel(paste0(base_dir, "scan_registration_qc (2).xlsx")) %>% 
        clean_names
    
    
    df_qc3 %>% 
        dplyr::rename("contrast" = "coxtrast", "init_scale" = "ixit_scale") %>% 
        filter(contrast != "fMRI") %>% 
        filter(registration_success != "o" | is.na(registration_success)) %>% 
        mutate(init_scale = ifelse(
            
            grepl("[0-9]", other_notes),
            other_notes,
            init_scale
            
        ),
        
        init_scale = ifelse(
            
            is.na(registration_success),
            0.90,
            init_scale
            
        ),
        
        init_scale = ifelse(
            
            grepl("[a-z]", other_notes),
            0.8,
            init_scale
            
        ),
        
        init_scale = as.numeric(init_scale)
        
        ) %>% 
        
        # format for upload to HPC
        mutate(ses = ifelse(timepoint < 100, 
                            paste0("ses-PND0", timepoint), 
                            paste0("ses-PND", timepoint)
        ),
        sub = case_when(
            
            study == "EDA" & id < 10 ~ paste0("sub-", study, "A0", id),
            study == "EDA" & id > 10 ~ paste0("sub-", study, "A", id),
            study == "JWD" & id < 10 ~ paste0("sub-", study, "0", id),
            study == "JWD" & id > 10 ~ paste0("sub-", study, id)
            
        )
        ) %>% 
        dplyr::select(sub, ses, contrast, init_scale) %>% 
        
        write.csv(paste0(base_dir, "scan_list4.csv"), row.names = FALSE)
    
    


# ROUND 4 ------------------------------------------------------------------

    
    
    df_qc4 <- read_excel(paste0(base_dir, "scan_registration_qc (3).xlsx")) %>% 
        clean_names
    
    
    df_qc4 %>% 
        filter(contrast != "fMRI") %>% 
        filter(registration_success == "x") %>%
        mutate(init_scale = ifelse(
            
            grepl("[0-9]", other_notes),
            other_notes,
            init_scale
            
        ),
        
        init_scale = ifelse(
            
            grepl("[a-z]", other_notes),
            0.7,
            init_scale
            
        ),
        
        init_scale = as.numeric(init_scale)
        
        ) %>% 
        
        # format for upload to HPC
        mutate(ses = ifelse(timepoint < 100, 
                            paste0("ses-PND0", timepoint), 
                            paste0("ses-PND", timepoint)
        ),
        sub = case_when(
            
            study == "EDA" & id < 10 ~ paste0("sub-", study, "A0", id),
            study == "EDA" & id > 10 ~ paste0("sub-", study, "A", id),
            study == "JWD" & id < 10 ~ paste0("sub-", study, "0", id),
            study == "JWD" & id > 10 ~ paste0("sub-", study, id)
            
        )
        ) %>% 
        dplyr::select(sub, ses, contrast, init_scale) %>% 
        
        
        write.csv(paste0(base_dir, "scan_list5.csv"), row.names = FALSE)
    

    

# ROUND 5 -----------------------------------------------------------------

    df_qc5 <- read_excel(paste0(base_dir, "scan_registration_qc (4).xlsx")) %>% 
        clean_names
    
    df_qc5 %>% 
        filter(contrast != "fMRI" & registration_success != "o") %>% 
        count(timepoint, contrast) # look at ROI stats!
    
    # format for upload to HPC
    df_qc5 %>% 
        filter(contrast != "fMRI" & registration_success == "o") %>% 
        mutate(ses = ifelse(timepoint < 100, 
                            paste0("ses-PND0", timepoint), 
                            paste0("ses-PND", timepoint)
        ),
        sub = case_when(
            
            study == "EDA" & id < 10 ~ paste0("sub-", study, "A0", id),
            study == "EDA" & id > 10 ~ paste0("sub-", study, "A", id),
            study == "JWD" & id < 10 ~ paste0("sub-", study, "0", id),
            study == "JWD" & id > 10 ~ paste0("sub-", study, id)
            
        )
        ) %>% 
        dplyr::select(sub, ses, contrast, init_scale) %>% 
        
        write.csv(paste0(base_dir, "scan_list.csv"), row.names = FALSE)
    
    
    
    # try one more round on these PND020s....    
    df_qc5 %>% 
        filter(contrast != "fMRI") %>% 
        filter(registration_success == "x") %>%
        mutate(init_scale = ifelse(
            
            grepl("[0-9]", other_notes),
            other_notes,
            init_scale
            
        ),
        
        init_scale = ifelse(
            
            grepl("[a-z]", other_notes),
            0.7,
            init_scale
            
        ),
        
        init_scale = as.numeric(init_scale)
        
        ) %>% 
        
        # format for upload to HPC
        mutate(ses = ifelse(timepoint < 100, 
                            paste0("ses-PND0", timepoint), 
                            paste0("ses-PND", timepoint)
        ),
        sub = case_when(
            
            study == "EDA" & id < 10 ~ paste0("sub-", study, "A0", id),
            study == "EDA" & id > 10 ~ paste0("sub-", study, "A", id),
            study == "JWD" & id < 10 ~ paste0("sub-", study, "0", id),
            study == "JWD" & id > 10 ~ paste0("sub-", study, id)
            
        )
        ) %>% 
        dplyr::select(sub, ses, contrast, init_scale) %>% 
        
        
        write.csv(paste0(base_dir, "scan_list6.csv"), row.names = FALSE)
    
    

# ROUND 6 -----------------------------------------------------------------

    
    df_qc6 <- read_excel(paste0(base_dir, "scan_registration_qc (5).xlsx")) %>% 
        clean_names
    
    df_qc6 %>% 
        filter(contrast != "fMRI" & registration_success != "o") %>% 
        count(timepoint, contrast) # look at ROI stats!
    
    # format for upload to HPC
    df_qc6 %>% 
        filter(contrast != "fMRI" & registration_success == "o") %>% 
        mutate(ses = ifelse(timepoint < 100, 
                            paste0("ses-PND0", timepoint), 
                            paste0("ses-PND", timepoint)
        ),
        sub = case_when(
            
            study == "EDA" & id < 10 ~ paste0("sub-", study, "A0", id),
            study == "EDA" & id > 10 ~ paste0("sub-", study, "A", id),
            study == "JWD" & id < 10 ~ paste0("sub-", study, "0", id),
            study == "JWD" & id > 10 ~ paste0("sub-", study, id)
            
        )
        ) %>% 
        dplyr::select(sub, ses, contrast, init_scale) %>% 
        
        write.csv(paste0(base_dir, "scan_list.csv"), row.names = FALSE)
    
    