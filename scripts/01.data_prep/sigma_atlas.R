
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

# get image data
    img <-  paste0(atlas_dir, "SIGMA_Anatomical_Brain_Atlas.nii")
    m_img_data <- asNifti(img)

# get labels & atlas hierarchy
    structures <- paste0(atlas_dir, "SIGMA_Anatomical_Brain_Atlas_ListOfStructures.csv")
    df_sigma_labels <- read_csv(structures) %>% 
        left_join( readxl::read_xlsx(paste0(base_dir, "data/SIGMA_analysis/SIGMA_abbreviations.xlsx"))) %>% 
        clean_names %>% 
        mutate_at(c("matter", "territories", "system", "region_of_interest"), ~tolower(.x) %>% 
                      str_replace_all(" ", "_") %>% 
                      str_replace_all("-", "_")) %>% 
        pivot_longer(2:3, names_to = "hemi", values_to = "idx") %>% 
        mutate(hemisphere = gsub("_.*", "", hemi)) %>% 
        select(-hemi) %>% 
      dplyr::select(-ends_with("_2"))
    
    # remove multiple instances of fi
    df_sigma_labels <- df_sigma_labels %>% 
      filter(!(region_of_interest == "fimbria_of_the_hippocampus" & system == "medulla") &
               !(region_of_interest == "fimbria_of_the_hippocampus" & idx %in% c(941, 942))) %>% 
      distinct()
    
    df_sigma_labels %>% filter(region_of_interest == "fimbria_of_the_hippocampus")
    
    # save for later use
    save(df_sigma_labels, file = paste0(base_dir, "objects/21June2023_sigma_labels.RDS"))
    
    
# data format -------------------------------------------------------------


    # format image data
    l_img_data <- 1:dim(m_img_data)[1] %>% map( ~ m_img_data[.x, , ] %>% 
                                                   as.data.frame %>% 
                                                   as_tibble %>% 
                                                   mutate(x = .x,
                                                          .before = 1)
                                               )
    
    l_img_data <- 1:length(l_img_data) %>% map( ~  l_img_data[[.x]] %>%     
                                                   mutate(y = row_number(), .before = 2) %>% 
                                                   pivot_longer(3:ncol(.), names_to = "z", values_to = "idx") %>% 
                                                   mutate(z = str_remove(z, "V") %>% 
                                                              as.numeric()
                                                          )
    )
    
    df_img_data <- l_img_data %>% 
        bind_rows() %>% 
        filter(idx != 0)
    
    # combine image data with labels
    doParallel::registerDoParallel()
    df_sigma <- df_sigma_labels %>% 
        left_join(df_img_data) %>% 
        distinct()
  
    # save atlas!!
    save(df_sigma, file = paste0(base_dir, "objects/23June2023_sigma_atlas_full.RDS"))
    
    
# exploratory plots -------------------------------------------------------------------

    load(paste0(base_dir, "objects/23June2023_sigma_atlas_full.RDS"))
    
    x_dim <- dim(m_img_data)[1] # 260
    y_dim <- dim(m_img_data)[2] # 342
    z_dim <- dim(m_img_data)[3] # 184
    
    # axial slice
    df_sigma %>%     
        filter(z == z_dim/2) %>%  
        ggplot(aes(x = x, y = y, color = region_of_interest)) +
        geom_point() +
        # scale_color_gradientn(colors = brewer.pal(9, "Spectral")) +
        theme_classic() +
        guides(color = guide_legend(ncol = 2))
    
    # sagittal slice
    df_sigma %>%     
        filter(x == x_dim/2) %>%  
        ggplot(aes(x = y, y = z, color = region_of_interest)) +
        geom_point() +
        # scale_color_gradientn(colors = brewer.pal(9, "Spectral")) +
        theme_classic() +
        #guides(color = guide_legend(nrow = 2)) +
        theme(legend.position = "bottom")
    
    # coronal slice
    df_sigma %>%     
        filter(y == y_dim/2) %>%  
        ggplot(aes(x = x, y = z, color = region_of_interest)) +
        geom_point() +
        # scale_color_gradientn(colors = brewer.pal(9, "Spectral")) +
        theme_classic() +
        guides(color = guide_legend(ncol = 2))
    
    # 3D
    df_sigma %>% 
        plot_ly(x = ~x, 
                y = ~y, 
                z = ~z,
                color = ~region_of_interest,
                type = "scatter3d",
                mode = "markers")

    

# find region borders -----------------------------------------------------

    load(paste0(base_dir, "objects/23June2023_sigma_atlas_full.RDS"))
    
    borders <- df_sigma %>%
        na.omit() %>% 
        st_as_sf(coords = c("x", "y", "z")) %>% 
        group_by(hemi, region_of_interest) %>%
        summarize(geometry = st_union(geometry)) %>%
        st_convex_hull()
    
    p <- borders %>% 
        filter(hemi == "left") %>% 
        ggplot() +
        geom_sf(aes(fill = region_of_interest))
    
    print(p)
    ggplotly(p)
    
    ## AXIAL
    axial <- df_sigma %>%
        filter(z == z_dim/2) %>% 
        st_as_sf(coords = c("x", "y")) %>% 
        group_by(hemisphere, region_of_interest) %>%
        summarize(geometry = st_union(geometry)) %>%
        st_convex_hull()
    p_axial <- axial %>% 
        # filter(hemisphere == "left") %>% 
        ggplot() +
        geom_sf(aes(fill = region_of_interest)) +
        theme_void() +
        theme(legend.position = "none")
    ggplotly(p_axial)
    
    ## SAGITTAL
    sagittal <- df_sigma %>%
        filter(x == x_dim/2) %>% 
        st_as_sf(coords = c("y", "z")) %>% 
        group_by(hemisphere, region_of_interest) %>%
        summarize(geometry = st_union(geometry)) %>%
        st_convex_hull()
    p_saggital <- sagittal %>% 
        filter(hemisphere == "left") %>% 
        ggplot() +
        geom_sf(aes(fill = region_of_interest)) +
        theme_void() +
        theme(legend.position = "none")
    ggplotly(p_saggital)
    
    ## CORONAL
    coronal <- df_sigma %>%
        filter(y == y_dim/2) %>% 
        st_as_sf(coords = c("x", "z")) %>% 
        group_by(hemisphere, region_of_interest) %>%
        summarize(geometry = st_union(geometry)) %>%
        st_convex_hull()
    p_coronal <- coronal %>% 
        #filter(hemisphere == "left") %>% 
        ggplot() +
        geom_sf(aes(fill = region_of_interest)) +
        theme_void() +
        theme(legend.position = "none")
    ggplotly(p_coronal)
    
    
    ##
    ## create atlas for all coronal slices
    ##
    y_slices <- df_sigma %>% filter(!is.na(y)) %>% arrange(y) %>% pull(y) %>% unique
    df_coronal_atlas <- y_slices %>% 
        map_df( ~ df_sigma %>%
                    filter(y == .x) %>% 
                    st_as_sf(coords = c("x", "z")) %>% 
                    group_by(hemisphere, region_of_interest) %>%
                    summarize(geometry = st_union(geometry)) %>%
                    st_convex_hull() %>% 
                    dplyr::mutate(side = "coronal",
                                  y_slice = .x,
                                  .before = 1)
        )
    
    df_coronal_atlas %>% 
        ggplot() +
        geom_sf(aes(fill = region_of_interest)) +
        facet_wrap(~ y_slice) +
        theme_void() +
        theme(legend.position = "none")
    
    ##
    ## create atlas for all sagittal slices
    ##
    x_slices <- df_sigma %>% filter(!is.na(x)) %>% arrange(x) %>% pull(x) %>% unique
    df_sagittal_atlas <- x_slices %>% 
        map_df( ~ df_sigma %>%
                    filter(x == .x) %>% 
                    st_as_sf(coords = c("y", "z")) %>% 
                    group_by(hemisphere, region_of_interest) %>%
                    summarize(geometry = st_union(geometry)) %>%
                    st_convex_hull() %>% 
                    dplyr::mutate(side = "sagittal",
                                  x_slice = .x,
                                  .before = 1)
        )
    
    df_sagittal_atlas %>% 
        ggplot() +
        geom_sf(aes(fill = region_of_interest)) +
        facet_wrap(~ x_slice) +
        theme_void() +
        theme(legend.position = "none")
    
    ##
    ## create atlas for all axial slices
    ##
    z_slices <- df_sigma %>% filter(!is.na(z)) %>% arrange(z) %>% pull(z) %>% unique
    df_axial_atlas <- z_slices %>% 
        map_df( ~ df_sigma %>%
                    filter(z == .x) %>% 
                    st_as_sf(coords = c("x", "y")) %>% 
                    group_by(hemisphere, region_of_interest) %>%
                    summarize(geometry = st_union(geometry)) %>%
                    st_convex_hull() %>% 
                    dplyr::mutate(side = "axial",
                                  z_slice = .x,
                                  .before = 1)
        )
    
    df_axial_atlas %>% 
        ggplot() +
        geom_sf(aes(fill = region_of_interest)) +
        facet_wrap(~ z_slice) +
        theme_void() +
        theme(legend.position = "none")
    
    ##
    ## COMBINE AND SAVE
    ##
    
    df_sigma_atlas <- df_axial_atlas %>% 
        bind_rows(df_sagittal_atlas) %>% 
        bind_rows(df_coronal_atlas) %>% 
        dplyr::select(hemisphere, side, region_of_interest, x_slice, y_slice, z_slice, geometry)
    
    save(df_sigma_atlas, file = "objects/20220909_sigma_atlas_polygons.RDS")
    
    ## 2022-11-24 add ROI abbreviations
    load("objects/20220909_sigma_atlas_polygons.RDS")
    df_sigma_atlas <- df_sigma_atlas %>%
        mutate(region_of_interest = sub("_$", "", region_of_interest)) %>% 
        left_join(df_sigma_labels) %>% 
        dplyr::select(side, matter, territories, system, region_of_interest, roi_abbreviation, hemisphere, x_slice, y_slice, z_slice, geometry)
    
    save(df_sigma_atlas, file = "objects/20230111_sigma_atlas_polygons.RDS")
    
    ## 3D plot?
    # 3D
    df_sigma %>% 
        plot_ly(x = ~x, 
                y = ~y, 
                z = ~z,
                color = ~region_of_interest,
                opacity = 0.7,
                type = "mesh3d")
  

# create maps to study
    library(ggrepel)
    library(gridExtra)
    
    df_color <- df_sigma_labels %>% 
        select(roi_abbreviation) %>% 
        distinct() %>% 
        mutate(color = hue_pal()(114) %>% sample)
    
    pdf("SIGMA_sagittal_slices.pdf", onefile = TRUE)
    x_slices <- df_sigma_atlas %>% filter(!is.na(x_slice)) %>% pull(x_slice) %>% unique
    for(i in x_slices) {
        
        print(
            
            df_sigma_atlas %>% 
                filter(side == "sagittal" & x_slice == i) %>%
                mutate(label = roi_abbreviation) %>% 
                #mutate(label = ifelse(hemisphere == "right", roi_abbreviation, "")) %>% 
                left_join(df_color, by = "roi_abbreviation") %>% 
                
                # plot
                ggplot() +
                geom_sf(aes(fill = I(color), geometry = geometry), alpha = 0.7) +
                geom_text_repel(aes(label = label, geometry = geometry),
                                stat = "sf_coordinates",
                                min.segment.length = 0,
                                size = 4) +
                theme_void() +
                ggtitle(paste0("slice ", i)) +
                theme(legend.position = "none")
            
        )

    }
    dev.off()



    