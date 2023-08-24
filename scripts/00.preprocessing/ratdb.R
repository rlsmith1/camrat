


# libraries ---------------------------------------------------------------


    library(R.matlab)
    library(tidyverse)
    library(purrr)
    library(data.table)
    library(lubridate)
    library(janitor)



# read in data ------------------------------------------------------------


    # from ratdb.mat
    setwd("/Users/rachelsmith/OneDrive - National Institutes of Health/CamRat/camrat/")
    l_ratdb <- readMat("data/ratdb.mat")

    ## note: cages uninteresting (cages that have been mislabeled that Steve's ratdb fixes automatically)
    
    # subject info files
    filenames <- paste0(path, list.files(path = "data/subject_data/"))
    
    l_subjs <- 1:length(filenames) %>% purrr::map(~read_csv(filenames[.x]) %>% as_tibble())
    
    ## subject file for every study - new subject for each day a rat is scanned
    ## note: EDAA in SUBJECT_id denotes GSK animals (Ethan's initials)

    # scan type (copied from ratdb in matlab: type = char(type))
    scan_type <- 'apmbt  ddapmt  dbdapmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt d dapmbt  dd pmt bapmbt  ddapmbt  dddd pmt bapmbt  ddapmbtd d ddddapmbt  dddd pmbbt db bdpmbdbtdb b    pmdtd bbd    bddbmt  pbd  tbddp  md d b  pmdtd b pmbbbdt db  bm   b dp t  dpmddbt  ddapmbtapmbt  dd pmddtbdd bapmbt  ddapmbt  dd pmt b pmt bapmbt  ddapmt  dbdapmbt  dd pm  dtd b tb m p pmdtd bapmbt  dd pmbbdt  b db pmbt  pmt b      apmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt  dd pmt bapmt  dbdapmbt d d pbbm tbdbd pt dm db  pmdtd bapmbt  ddapmbt  ddb  pm   bdd ttbd mpd   pmdtd b pbbmd tdbddb pb t md  d pmdtd bapmbt  ddapmt  dbdapmbt  ddapmbt  ddapmbt  ddapmbt  dd  ddapmbt  ddbapmbt d dapmbt  ddddapmbt  ddddapmbt  ddapmbt  ddapmbt  ddapmbt  dd pmtdbtbd mp d apmbt  ddapmdbdt  ddapmbt  apmbt  dddd b pmb tdb  b pmdtd b tmbdp d b pmdtd b pmbbdt  b dbtb m pd d pmt pdmbbtb bd  dd btd dbpm    pmbbt db bdpbdt  mdbp bbb pmt b pmbbt db bdapmbt  ddapmbt  ddddapmt  dbdapmbt d dapmbt  ddapmbt  ddapmdt  bdapmbt  ddapmbt  ddapmbt b ddapmbt  ddapmbt  ddapmbt d dapmbt  ddapmbt  ddddapmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt  ddbapmbt  ddapmbt  dd      apmbt d dapmbt  ddapmbt  ddapmbt d dapmbt  ddapmbt  ddapmbt  ddapmbt  ddtmdpbab fd apmbt  ddtm pba    apm  t      fbapmft db dbtmdbdpab f apm   t  bb apmbt  ddapm b t     b   pb  mat    apmt b    apmbt  ddapmbbt bdb fdap b  tm apm b  t     apmbt  ddapmbdtb fd p matbbpdbb d   maf t   apmbdtbb fd apmbdtb fd apmtb bapmbd   t   fdb apmdtbb fd apmtbapm  bd    b fdt apmbddtb dd apmtb bapm bd   t   fdb apmbdtdbd fd a  mpd bb fdt apmbt  ddapmt  dbdapmbt  ddapmdtb dapmbt  ddapmbt d dapmdbtb dapmbtb apmbt  d  dapmt  dbdapmbt  ddapmbtb f apmbtb f apmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt  ddapmbt  apmbt  ddapmbt  ddapmbt  ddapmbt  ddmbapmdbtb fd apmbt  ddapmbt  ddfapmdftb bd apmdbtb fd apmbt  ddapmdtbb fd apmt  dbdapmbt  ddapmbt d dapm dtbb fdf apmbtb  bapmtbtmbfbpadd  bapmbtb tmbfapb dd   pmfbtbadd   pmbtba  pmfbtbadd   pmfbtbadd   pmbtba  pmbta pmfbtbadd  pdbbbd m tf a pmfbtbadd   pm   t    fab       b pmfbtbadd   pm  ft  ab pmfbtbadd   pmddfbtba   pmbfbtadd   pmbta pmfbtadd pmfbtbadd   pm fbtbadd  pmdfbtba d  pmfbtbadd   pmfbtbadd   pmbfbtadd   pm fbtbadd  p bm bta  pmfbtbadd   pmfbtbadd  pdbbd m tf a pmfbtbadd   pmbfbtbadd   pm tbadf bd pm fbtbadd  pmbftbadd   pmfbtbadd  pmfbtbadd   pmfbtbadd   pmfbtbadd   pmfbtbadd   tnmdbb dpaadf dtnmdbb dpaadf dtnmdbbd ddpaaf tnmdbb dpaadf dtnmdbb dpaadf dtnmdbb dpaadf dtnmdbb dpaadf dtnmdbb dpaadf dndbdf mtpd dabatnmdbbdpdadaf  tnmdbb dpaadf dtnmdbb dpaadf dndbdf mdtdp abatnmdbb dpaadf dtnmdbbdpdadaf  tndmdbbn  dpaabfd mbdda dpfnbdat tnmdbbdpaadf d tnmdbb dpaadf dtnmdbb dpaadf dtnmdbb dpaadf dtnmdbbbdpaadf d tnmdbb dpaadf dtnmdbbdpdadaf  tnmdbb dpaadf dtnmdbb dpaadf dtnmdbb dpaadf dtnmdbb dpaadf dtnmdbb dpaadf dtnmdb dpdadabf tnmbb fpaaddd dtnmdbbdpdadaf  tmbddb papt fm abn na  dadtnmdbbb dpaadf dtmbb pa a      tmdbbb adf t  ddabt ddtndbb dpaaddtnmbbddpa pmbt ddndpadapmpm mf dd afa apmbt  apmtf b fapmtf b apm tfb apmtf b apmtf b apmtf b bamf t papmtf b apmtf b apmtf b apmtf b apm tfb apmtf b apmtf b apmtf b apmtf b apmtf b pmtf b apmtf b apmtf b apmtf b apm tfb apmtf b apmtf b apmtf b apmtf b apmtf b apmtf b apmtf b apmtf b apmtf b apmtf b apmtf b apm tfb apmtf b apmtf b apmtf b apmtf b tmpfabtmpfabtmpfabapmtfbtmpfa btm pfa btmpfab tmpfa  btmpfabtmpfabtmpfabtm  pf a  bapmtbftmpfabtmpfa bapmtf b apmtf b apmtf b apmtf b apmtf b '
    
    # to get asns, run in matlab: for i = 1:numel(asns);fprintf(1, '%s\n',asns{i});end;
    
    
    
# explore data ------------------------------------------------------------


    str(l_ratdb)
    length(l_ratdb)
    names(l_ratdb)

    

# format data into scans df and subject df --------------------------------


    # determine data type of each list element
    l_datatype <- 1:length(l_ratdb) %>% purrr::map(~l_ratdb[[.x]] %>% class()) %>% 
      unlist() %>%
      subset(. %in% c("list", "matrix"))
      
    # select only matrices
    l_ratdb_mats <- l_ratdb[which(l_datatype == "matrix")]
    
    # convert matrices to tibbles
    l_ratdb_dfs <- 1:length(l_ratdb_mats) %>% purrr::map(~l_ratdb_mats[[.x]] %>% as_tibble())
    
    # add variable names
    names(l_ratdb_dfs) <- names(l_ratdb)[which(l_datatype == "matrix")]
    
    # convert scans to df
    df_scans <- l_ratdb_dfs[-c(1, 4)] %>% 
      rbindlist() %>% 
      as_tibble() #%>% 
      mutate(variable = names(l_ratdb_dfs) %>% subset(!(. %in% c("ages", "matrices"))), .before = "V1")
    
    # add matrices back
    df_matrices <- l_ratdb$matrices %>% 
      t() %>% 
      as_tibble() #%>% 
      mutate(variable = c("matrix_x", "matrix_y", "matrix_z"), .before = "V1")
      
    names <- c(names(l_ratdb_dfs) %>% subset(!(. %in% c("ages", "matrices"))), "matrix_x", "matrix_y", "matrix_z")
    
    df_scans <- df_scans %>% bind_rows(df_matrices)
    df_scans <- t(df_scans) %>% as_tibble()
    colnames(df_scans) <- names

    # add variables that are in list form, but contain scan info
    
        # select scan-related variables in list form (already concatenated scan-specific variables in mat form)
        scans <- c("acqps", "asns", "dates", "dfiles", "ids", "meths", "sfroots", "visus") # scan-specific variables
        subj <- c("ages", "cages", "dobs", "sdates", "subjs") # subject-specific variables
        
        scan_lists <- l_ratdb[which(l_datatype == "list")] %>% names() %>%
            subset(. %in% scans)
    
        # select these variables from ratdb
        l_ratdb_lists <- l_ratdb[scan_lists]
        
        # convert to tibble
        l_ratbd_lists_dfs <- 1:length(l_ratdb_lists) %>% 
            purrr::map(~l_ratdb_lists[[.x]] %>% 
                           unlist() %>% 
                           as_tibble_col(column_name = scan_lists[.x]))
        
        # concatenate into tibble
        df_ratdb_lists_dfs <- l_ratbd_lists_dfs[-c(2,7)] %>% bind_cols()
        
        # combine with other scan info
        df_scans <- df_scans %>% bind_cols(df_ratdb_lists_dfs) %>% mutate(dates = as_date(dates))
    
    # identify scan type
        # add character vector as type (still indexed with type0)
        df_scans <- df_scans %>% mutate(type = strsplit(scan_type, "") %>% unlist())
        
        # key: MTs = 'm', diffusions = 'd', Ts = 't', PDs = 'p', a = actual flip angle image (calibration thing), b = b1 map, f = field map

        df_scans <- df_scans %>% 
            mutate(
                
                gsk = ifelse(optid == 3, 1, 0),
                mt = ifelse(type == "m", 1, 0), 
                dif = ifelse(proc == 1 & type == "d", 1, 0), 
                fmri = ifelse(nreps > 50, 1, 0)
                
            ) %>% 
            
            # reorder columns here
            select(c("did", "ids", "gsk", "mt", "dif", "fmri", "type", "type0", "optid", "proc", "nreps", "sind", "dages", "dates", everything()))
    

    ######
        df_scans %>% 
            select(ids, dates) %>% 
            distinct() %>% 
            filter(grepl("JWD|EDAA", ids)) %>% 
            left_join(df_scans %>% 
                          filter(fmri == 1 & trs %in% c(1832, 1312)) %>% 
                          select(ids, dates, dfiles), 
                      by = c("ids", "dates")) %>% 
            arrange(ids, dates) %>% 
            distinct() %>% 
            
            filter(is.na(dfiles)) %>% select(ids) %>% distinct()
        
df_scans %>% filter(ids == "C108088-JWD01" & fmri == 1) %>% select(ids, dates, dfiles)



# subject metadata --------------------------------------------------------

l_subjs <- l_subjs[-4]

df_subjs <- 1:length(l_subjs) %>% map(
    
    ~l_subjs[[.x]] %>% 
        as_tibble() %>% 
        filter(grepl("SUBJECT", OriginalVariableNames)) %>%
        mutate(OriginalVariableNames = str_remove(OriginalVariableNames, "SUBJECT_")) %>%
        mutate(subject = .$Var1[.$OriginalVariableNames == "id"],
               subject = str_remove(subject, "<"),
               subject = str_remove(subject, ">"),
               .before = 1) %>% 
        pivot_wider(id_cols = subject, names_from = OriginalVariableNames, values_from = Var1)
    
    
) %>% 
    bind_rows() %>% 
    mutate(subject = gsub(".*-", "", subject),
           subject = gsub(".* ", "", subject),
           subject = gsub("\\.", "", subject)) %>% 
    arrange(subject)


save(df_subjs, file = "data/20221107_all_subject_data.RDS")

df_subjs %>% 
    count(subject) %>% 
    print(n = nrow(.))


# Identify rats in GSK study and scans that exist -------------------------

df_scans %>% 
    filter(gsk == 1) %>% 
    count(did, sind, ids, dages, dates, type)

df_subjs$ages %>% unique() %>% sort()
df_scans$dages %>% unique() %>% sort()
df_scans$sind %>% unique() %>% sort()

# also....

# all the EPIs (fMRIs) in the GSK study     
df_scans %>% filter(fmri == 1, gsk == 1)

# ages
df_scans %>% filter(gsk == 1) %>% 
    ggplot(aes(x = dages)) +
    geom_histogram()


    # extract diffusion scans
    diff_scans <- df_scans %>% filter(gsk == 1, dif == 1) %>% .$dfiles %>% 
        str_replace("sjs80", "rls83") # %>% 
    # str_replace("/home/rls83", "")
    # write.table(diff_scans, "CamRat/data/gsk_diff_scans.txt")
    
    # extract MPM scans
    df_scans %>% filter(gsk == 1 & type %in% c("p", "m", "t")) %>% .$dfiles %>% 
        str_replace("sjs80", "rls83") %>% 
        str_replace("rds/q10014", "rat_data")
    
    
# export ratdb as Rdata object
    
    save(df_scans, file = "CamRat/camrat/objects/ratdb_scans.Rdata")
    
write.csv(df_scans, file = "outputs/ratdb_scans.csv")




# dcm2nii -----------------------------------------------------------------


setwd("/Users/rachelsmith/OneDrive - National Institutes of Health/CamRat/camrat/")
load("objects/ratdb_scans.Rdata")

df_scans %>% 
    filter(grepl("JWD19", ids)) %>% 
    filter(fmri == 1) %>% 
    dplyr::select(ids, fmri, nreps, dages, dates, ser, trs, dfiles)




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
        
    )) %>% 
    
    write.csv("outputs/20230124_scan_registration_qc.csv", row.names = FALSE)
    
    

# where are our scans going -----------------------------------------------



df_edaa_info <- read.csv("data/subject_info.csv") %>% 
    as_tibble() %>% 
    na.omit() %>% 
    clean_names() %>% 
    dplyr::rename("sex" = "gender") %>% 
    mutate(id = as.character(id) %>% ifelse(nchar(.) < 2, paste0("0", id), .)) %>% 
    mutate(subject = paste0("EDAA", id))

df_steve_data <- read.table("~/OneDrive - National Institutes of Health/rat_mri_preproc/afni/steve_data.txt") %>% 
    as_tibble() %>% 
    filter(grepl("EDAA|JWD", V1)) %>% 
    separate(V1, into = c("ids", "dates"), sep = "/") %>% 
    mutate(dates = str_remove(dates, "_:") %>% as.Date(format = "%Y%m%d")) %>% 
    left_join(df_scans %>%
                  filter(mt == 1) %>% 
                  distinct, by = c("ids", "dates")) %>% 
    select(ids, dates, dages) %>% 
    distinct() %>% 
    mutate(ids = str_remove(ids, ".*-"),
           ids = str_remove(ids, ".*_")) %>% 
    mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310))) %>% 
    dplyr::select(-dates, -dages) %>% 
    dplyr::rename("subject" = "ids") %>% 
    left_join(df_edaa_info) %>% 
    dplyr::select(-id) %>% 
    distinct() %>% 
    mutate(timepoint = case_when(
        
        timepoint == 1 ~ 20,
        timepoint == 2 ~ 35,
        timepoint == 3 ~ 63,
        timepoint == 4 ~ 200
        
    ))


df_raw_scans <- read.table("~/OneDrive - National Institutes of Health/rat_mri_preproc/afni/raw_scans.txt") %>% 
    as_tibble() %>% 
    filter(grepl("sub", V1)) %>% 
    separate(V1, into = c("subject", "timepoint"), sep = "/") %>% 
    mutate(subject = str_remove(subject, "sub-"),
           timepoint = str_remove(timepoint, "ses-PND"),
           timepoint = str_remove(timepoint, ":") %>% as.numeric()) %>% 
    left_join(df_edaa_info) %>% 
    dplyr::select(-id) %>% 
    distinct()

df_registered_scans <- read.table("~/OneDrive - National Institutes of Health/rat_mri_preproc/afni/registered_scans.txt") %>% 
    as_tibble() %>% 
    mutate(subject = str_remove(V1, "derived/sub-|"),
           timepoint = str_remove(V1, "ses-PND")) %>% 
    dplyr::select(-V1) %>% 
    mutate(subject = ifelse(grepl("ses-", subject), NA, subject),
           timepoint = ifelse(grepl("derived", timepoint), NA, timepoint)) %>% 
    mutate(timepoint = .[row_number() + 1, ]$timepoint %>% as.numeric) %>% 
    filter(!is.na(timepoint)) %>% 
    fill(subject) %>% 
    mutate(subject = str_remove(subject, ":")) %>% 
    left_join(df_edaa_info) %>% 
    dplyr::select(-id) %>% 
    distinct()
 
df_steve_n <- df_steve_data %>% filter(grepl("EDAA", subject) & !is.na(timepoint)) %>% distinct %>% dplyr::count(timepoint, sex, group) %>% dplyr::rename("steve_n" = "n")
df_camrat_bids_n <- df_raw_scans %>% filter(grepl("EDAA", subject)) %>% distinct %>% dplyr::count(timepoint, sex, group) %>% dplyr::rename("camrat_bids_n" = "n")
df_registered_n <- df_registered_scans %>% filter(grepl("EDAA", subject)) %>% distinct %>% dplyr::count(timepoint, sex, group) %>% dplyr::rename("registered_n" = "n")

df_steve_n %>% 
    left_join(df_camrat_bids_n) %>% 
    left_join(df_registered_n)

df_steve_data %>% 
    filter(grepl("EDAA", subject) & !is.na(timepoint)) %>% 
    left_join(df_raw_scans %>% 
                  filter(grepl("EDAA", subject) & !is.na(timepoint)),
              by = c("subject", "timepoint")) %>% 
    filter(is.na(group.y))

read.table("~/OneDrive - National Institutes of Health/rat_mri_preproc/afni/steve_data.txt") %>% 
    as_tibble() %>% 
    filter(grepl("EDAA|JWD", V1)) %>% 
    separate(V1, into = c("ids", "dates"), sep = "/") %>% 
    mutate(dates = str_remove(dates, "_:") %>% as.Date(format = "%Y%m%d")) %>% 
    left_join(df_scans %>% 
    dplyr::select(ids, dates, fmri,nreps,dages, trs, dfiles) %>% 
    filter(fmri==1 & (trs == 1832.0 | trs == 1312.0))) %>% 
    filter(is.na(fmri))

read.table("~/OneDrive - National Institutes of Health/rat_mri_preproc/afni/steve_data.txt") %>% 
    as_tibble() %>% 
    filter(grepl("EDAA|JWD", V1)) %>% 
    separate(V1, into = c("ids", "dates"), sep = "/") %>% 
    mutate(dates = str_remove(dates, "_:") %>% as.Date(format = "%Y%m%d")) %>% 
    left_join(df_scans) %>% filter(is.na(did)) %>% dplyr::select(ids, dates) 
    
read.table("~/OneDrive - National Institutes of Health/rat_mri_preproc/afni/steve_data.txt") %>% 
    as_tibble() %>% 
    filter(grepl("EDAA|JWD", V1)) %>% 
    separate(V1, into = c("ids", "dates"), sep = "/") %>% 
    mutate(dates = str_remove(dates, "_:") %>% as.Date(format = "%Y%m%d")) %>% 
    left_join(df_scans %>% 
                  dplyr::select(ids, dates, mt, fmri, trs, dages, dfiles) %>% 
                  filter(mt == 1 | (fmri == 1 & (trs == 1832.0 | trs == 1312.0)))) %>% 
    filter(is.na(mt))

# C136150,_EDAA20 2019-04-02 not in ratdb
# C126523-JWD30   2018-10-23 only b map
# C138072,_EDAA49 2019-03-14 only b map
# C138651,_EDAA53 2019-03-20 only b map

df_edaa_info %>% filter(subject == "EDAA49")
df_raw_scans %>% filter(subject == "EDAA49")


df_scans %>% 
    filter(ids == "C138072,_EDAA49" & dates == "2019-03-14")

####

df_edaa_info %>% dplyr::count(sex, group)

df_scans %>% 
    filter(grepl("EDAA", ids)) %>% 
    dplyr::select(ids, dages) %>% 
    distinct() %>% 
    mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310))) %>% 
    mutate(subject = str_remove(ids, ".*_"), .before = 1) %>% 
    dplyr::select(-ids, -dages) %>% 
    left_join(df_edaa_info) %>% 
    
    dplyr::select(subject, sex, group) %>% 
    distinct() %>% 
    count(sex, group)

df_scans %>% 
    mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310)) %>% as.factor) %>% 
    ggplot(aes(x = dates, y = dages, color = timepoint)) + 
    geom_point()


# find fMRI scans ---------------------------------------------------------

    load("objects/ratdb_scans.Rdata")
    df_scans %>% filter(fmri == 1) %>% pull(dfiles)

    
    
# registration chain .csv file --------------------------------------------
    
    

    df_scans %>% mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310)), .before = gsk)
    
    df_scans %>% mutate(ids = str_replace_all(ids, "-", "_"),
                        ids = str_remove(ids, ","))
    
    df_scans %>% filter(dages != 0 & dages < 500 & type %in% c('p', 'm', 't')) 
    
    df_scans %>% filter(grepl("JWD", ids) | grepl("EDA", ids)) %>% filter(dages != 0 & dages < 500 & type %in% c('p', 'm', 't')) 
    
    # data 
    setwd("CamRat/camrat/")
    df_mt_jwd_eda <- read.csv("data/mt.subject.files.csv") %>% as_tibble

    # add is_common column
    df_timepoint_avail <- expand_grid(subject_id = df_mt_jwd_eda$subject_id %>% unique, timepoint_avail = seq(1:4))
    
    # df_is_common <- df_mt_jwd_eda %>% 
    #     group_by(subject_id) %>% 
    #     mutate(timepoint_avail = toString(timepoint)) %>% 
    #     unique %>% 
    #     mutate(is_common = case_when(
    #         
    #         grepl("3", timepoint_avail) ~ 3,
    #         !grepl("3", timepoint_avail) & grepl("4", timepoint_avail) ~ 4,
    #         !grepl("3", timepoint_avail) & !grepl("4", timepoint_avail) & grepl("2", timepoint_avail) ~ 2,
    #         !grepl("3", timepoint_avail) & !grepl("4", timepoint_avail) & !grepl("2", timepoint_avail) & grepl("1", timepoint_avail) ~ 1,
    #         
    #     )) %>% 
    #     
    #     ungroup %>% 
    #     dplyr::select(-timepoint_avail)
        
    # using 1/0 nomenclature
    df_is_common <- 
        df_mt_jwd_eda %>%
        group_by(subject_id) %>%
        mutate(timepoint_avail = toString(timepoint)) %>%
        mutate(is_common = case_when(
            
            grepl("3", timepoint_avail) & timepoint == 3 ~ 1,
            !grepl("3", timepoint_avail) & grepl("4", timepoint_avail) & timepoint == 4 ~ 1,
            !grepl("3", timepoint_avail) & !grepl("4", timepoint_avail) & grepl("2", timepoint_avail) & timepoint == 2 ~ 1,
            !grepl("3", timepoint_avail) & !grepl("4", timepoint_avail) & !grepl("2", timepoint_avail) & grepl("1", timepoint_avail) & timepoint == 1 ~ 1,
            TRUE ~ 0
            
        )) %>%
        
        ungroup %>% 
        dplyr::select(-timepoint_avail)
    
    df_is_common %>% count(timepoint, is_common) %>% filter(is_common == 1)
    
    df_is_common %>% write.csv("data/mt_subject_files.csv", row.names = FALSE)

    # select sample size: 5 rats from each study (JWD & EDAA) and all their respective timepoints
    jwd <- df_is_common %>% filter(grepl("JWD", subject_id)) %>% pull(subject_id) %>% unique
    jwd_sample <- sample(jwd, 5, replace = FALSE)
    
    edaa <- df_is_common %>% filter(grepl("EDA", subject_id)) %>% pull(subject_id) %>% unique
    edaa_sample <- sample(edaa, 5, replace = FALSE)
    
    df_is_common_sample <- df_is_common %>% filter(subject_id %in% c(edaa_sample, jwd_sample))
    df_is_common_sample %>% write.csv("data/mt_subject_files_sample.csv", row.names = FALSE)
    
    df_is_common_sample %>% print(n = nrow(.))
    
    # select only rats with all 4 timepoints
    complete_rats <- df_is_common %>% 
        group_by(subject_id) %>% 
        mutate(tps_combine = toString(timepoint)) %>% 
        filter(tps_combine == "1, 2, 3, 4") %>% 
        pull(subject_id) %>% unique
    
    df_is_common_sample2 <- df_is_common %>% filter(subject_id %in% complete_rats[1:3])
    df_is_common_sample2 %>% write.csv("data/mt_subject_files_sample2.csv", row.names = FALSE)
    
    
    

# count number of rats per timepoint --------------------------------------


    load("objects/ratdb_scans.Rdata")
    
    df_gsk_scans <- df_scans %>% 
        mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310)), .before = gsk) %>% 
        filter(gsk == 1 & type %in% c("t", "m", "p", "d"))
    
    df_gsk_scans %>% 
        select(ids, timepoint) %>% 
        unique %>% 
        count(ids) %>% # 51 GSK rats
        
        ggplot(aes(x = n)) +
        geom_histogram()
    
    df_rat_pnd <- df_scans %>% 
        mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310)), .before = gsk) %>% 
        mutate(pnd = case_when(
            
            timepoint == 1 ~ 20,
            timepoint == 2 ~ 35,
            timepoint == 3 ~ 60,
            timepoint == 4 ~ 300
            
        )) %>% 
        
        filter(grepl("EDAA|JWD", ids)) %>% 
        
        mutate(study = ifelse(grepl("JWD", ids), "MRC", "GSK")) %>% 
        
        select(study, did, pnd) %>% 
        unique %>% 
        arrange(study, did, pnd) %>% 
        unite(col = "rat_id", c(study, did)) %>% 
        mutate(rat_id = factor(rat_id, levels = unique(.$rat_id)))
    
    df_rat_pnd %>% 
        mutate(study = substr(rat_id, start = 1, stop = 3)) %>% 
        
        ggplot(aes(x = pnd, y = rat_id, color = study)) +
        geom_point(shape = 1, size = 2) +
        geom_line(aes(group = rat_id), alpha = 0.3) +
        scale_x_continuous(breaks = c(20, 35, 60, 300)) +
        theme_bw() +
        theme(axis.text.y = element_text(size = 7))

    df_rat_pnd %>% filter(grepl("GSK", rat_id)) %>% pull(rat_id) %>% unique %>% length
       
    df_rat_pnd %>% 
        group_by(rat_id) %>% 
        mutate(pnd = toString(pnd)) %>% 
        mutate(pnd20 = ifelse(grepl("20", pnd), 1, 0),
               pnd35 = ifelse(grepl("35", pnd), 1, 0),
               pnd60 = ifelse(grepl("60", pnd), 1, 0),
               pnd300 = ifelse(grepl("300", pnd), 1, 0)) %>% 
        dplyr::select(-pnd) %>% 
        unique %>% 
        ungroup %>% 
        
        filter(pnd20 == 1 &
                   pnd35 == 1 &
                   pnd60 == 1 &
                   pnd300 == 1)
    
    
    
    
    
    
    

