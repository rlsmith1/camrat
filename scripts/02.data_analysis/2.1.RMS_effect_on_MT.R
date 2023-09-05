
###
### Baseline MT and MT early developmental slope case-control differences & their relation to behavior
###

# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(tidymodels)
library(sf)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(gghalves)

# set plot theme ----------------------------------------------------------

theme_set(theme_light() +
            theme(plot.title = element_text(size = 10),
                  axis.title = element_text(size = 10),
                  axis.text = element_text(size = 10),
                  strip.text = element_text(size = 10, color = "black"),
                  strip.background = element_rect(fill = "white", color = "gray"),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 8)
                  #legend.key.width = unit(2, "cm")
            )
)

group_cols <- c("#FFA500", "#4682B4")
names(group_cols) <- c("control", "MS")


# helper functions --------------------------------------------------------

# RUN LME
f_run_lme <- function(df) {
  
  df %>% 
    mutate(lme_res = 
             
             map(
               
               .x = data, 
               .f = ~ nlme::lme(normalized_feature ~ group*age + sex + tbv, 
                                random = ~1 | subject, 
                                data = .x) %>% 
                 summary %>% 
                 .$tTable %>% 
                 as.data.frame %>% 
                 rownames_to_column("contrast") %>% 
                 as_tibble %>% 
                 filter(contrast != "(Intercept)") %>% 
                 clean_names
             )
           
    ) %>% 
    unnest(cols = c(lme_res)) %>% 
    ungroup %>% # ungroup to correct for multiple comparisons *across* regions (not correcting within ROI)
    arrange(p_value) %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
    dplyr::rename("term" = "contrast")
  
}

# RUN LM
f_run_lm <- function(df, timepoint = c(20, 63, 300)) {
  
  model <- ifelse(timepoint == 300,
                  "normalized_feature ~ group + sex + tbv + age",
                  "normalized_feature ~ group + sex + tbv")
  
  df_out <- df %>% 
    mutate(lm_res = 
             
             map(
               
               .x = data,
               .f = ~ lm(formula(model), data = .x) %>% 
                 tidy %>% 
                 clean_names() %>%
                 filter(term != "(Intercept)") %>%
                 dplyr::rename("t_value" = "statistic")
               
             )) %>%
    
    unnest(cols = c(lm_res)) %>%
    ungroup %>% # ungroup to correct for multiple comparisons *across* regions (not correcting within ROI)
    arrange(p_value) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  return(df_out)
  
}

# BOX PLOT
f_plot_box <- function(df, p_adj = TRUE, xlab, ylab, title, legend = TRUE) {
  
  df <- df %>% 
    unnest(cols = c(data)) %>% 
    mutate(system = str_replace_all(system, "_", " "))
  
  if (p_adj == TRUE) {
    
    df <- df %>% filter(term == "groupMS" & p_adj < 0.05) 
    
  } else if (p_adj == FALSE) {
    
    df <- df %>% filter(term == "groupMS" & p_value < 0.05)
    
  }
  
  p <- df %>%
    #mutate(system = factor(system, levels = xaxis_order)) %>% 
    
    ggplot(aes(x = fct_reorder(str_wrap(system, 10), t_value),
               y = feature_resids,
               color = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = group), size = 1.5, alpha = 0.3,
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
    scale_color_manual(values = group_cols)  +
    labs(x = paste0(xlab), y = paste0(ylab)) +
    ggtitle(paste0(title)) +
    theme(legend.position = "bottom")
  
  if (legend == TRUE) {
    
    return(p)
    
  } else if (legend == FALSE) {
    
    return(p +  theme(legend.position = "none"))
    
  }
  
}

# BRAIN MAP PLOT
f_plot_brain <- function(df, 
                         contrast = "groupMS",
                         p_adj = TRUE, 
                         labels = TRUE, 
                         column = FALSE,
                         slices = list("sagittal" = 125, "coronal" = 150, "axial" = 95),
                         legend = TRUE
) {
  
  if (p_adj == TRUE) {
    
    df <- df %>% 
      mutate(significant = ifelse(p_adj < 0.05, "yes", "no")) %>% 
      filter(term == contrast) 
    
  } else if (p_adj == FALSE) {
    
    df <- df %>% 
      mutate(significant = ifelse(p_value < 0.05, "yes", "no")) %>% 
      filter(term == contrast) 
    
  }
  
  df_slices <- map_dfr(.x = 1:length(slices),
                       .f = ~ df_sigma_atlas %>%
                         filter(region_of_interest != "striatum") %>% 
                         filter(side == names(slices[.x]) & slice == slices[.x])
                       
  ) %>% 
    mutate(side = factor(side, levels = c("axial", "sagittal", "coronal")))
  
  p <- df_slices %>% 
    left_join(df, by = "system") %>% 
    mutate(label = ifelse(significant == "yes" & hemisphere == "left", roi_abbreviation, NA)) %>% 
    arrange(!is.na(p_adj), -p_adj) %>% 
    dplyr::rename("t-value" = "t_value") %>% 
    
    ggplot() +
    geom_sf(aes(fill = `t-value`, color = significant, geometry = geometry, group = -1), 
            lwd = 0.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                         na.value = "grey", limits = c(-7.5, 7.5)) +
    scale_color_manual(values = c("darkgrey", "yellow"), na.value = "darkgrey") +
    #guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme_void() +
    theme(strip.text = element_blank(),
          legend.direction = "horizontal")
  
  if (column == TRUE) {
    
    p <- p +
      facet_wrap( ~ side + slice, ncol = 1)
    
  } else if (column == FALSE) {
    
    p <- p +
      facet_wrap( ~ side + slice, nrow = 1)
    
  }
  
  if (labels == TRUE) {
    
    p <- p +
      geom_text_repel(
        mapping = aes(label = label, geometry = geometry),
        stat = "sf_coordinates",
        min.segment.length = 0,
        size = 3
      )
    
  }
  
  if (legend == TRUE) {
    
    return(p)
    
  } else if (legend == FALSE) {
    
    return(p +  theme(legend.position = "none"))
    
  }
  
}

# TRAJECTORIES PLOT
f_plot_traj <- function(df, 
                        timepoints = c(20, 63),
                        xlab, 
                        ylab, 
                        title, 
                        dim = c(4, 3), 
                        legend = TRUE,
                        include_tails = TRUE,
                        lnwdth = 1,
                        include_individuals = TRUE) {
  
  df_plot <- df %>% 
    filter(term == "groupMS:age" & p_adj < 0.05) %>% 
    unnest(cols = c(data)) %>% 
    mutate(system = str_replace_all(system, "_", " ")) %>% 
    mutate(label = interaction(subject, region_of_interest), .before = 3) 
  
  factor_levels <- df_plot %>% arrange(t_value) %>% pull(system) %>% unique()
  
  p <- df_plot %>% 
    ggplot(aes(x = age, y = feature_resids, color = group)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = lnwdth) +
    
    geom_half_violin(data = df_plot %>% filter(age == 20) %>% mutate(age = 18),
                     aes(x = age, y = feature_resids, fill = group, color = group),
                     width = 8, alpha = 0.5, position = "identity", side = "l", 
                     trim = ifelse(include_tails == TRUE, FALSE, TRUE)
    ) +
    stat_summary(data = df_plot %>% filter(age == 20) %>% mutate(age = 18 - 2),
                 aes(x = age, y = feature_resids, fill = group),
                 geom = "crossbar", fun = "mean", width = 4, linewidth = 0.33
    ) +
    
    scale_color_manual(values = group_cols) +
    scale_fill_manual(values = group_cols) +
    
    facet_wrap(~ factor(system, levels = factor_levels), 
               scales = "free", ncol = dim[2], labeller = label_wrap_gen(width = 10)) +
    labs(x = paste0(xlab), y = paste0(ylab)) +
    ggtitle(paste0(title)) + 
    #ylim(c(-4, 4)) +
    theme_light() +
    theme(
      strip.text = element_text(color = "black", size = 14), 
      strip.background = element_rect(fill = "white", color = "gray"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14)
      
    )
  
  if (timepoints[2] == 63) {
    
    p <- p +
      geom_half_violin(data = df_plot %>% filter(age > 60) %>% mutate(age = 64),
                       aes(x = age, y = feature_resids, fill = group, color = group),
                       width = 8, alpha = 0.5, position = "identity", side = "r",  
                       trim = ifelse(include_tails == TRUE, FALSE, TRUE)
      ) +
      stat_summary(data = df_plot %>% filter(age > 60) %>% mutate(age = 64 + 2),
                   aes(x = age, y = feature_resids, fill = group),
                   geom = "crossbar", fun = "mean", width = 4, linewidth = 0.33
      ) +
      scale_x_continuous(breaks = c(20, 63))
    
  } else if (timepoints[2] == 300) {
    
    p <- p +
      geom_half_violin(data = df_plot %>% filter(age > 100) %>% mutate(age = 310),
                       aes(x = age, y = feature_resids, fill = group, color = group),
                       width = 8, alpha = 0.5, position = "identity", side = "r", 
                       trim = ifelse(include_tails == TRUE, FALSE, TRUE)
      ) +
      stat_summary(data = df_plot %>% filter(age > 100) %>% mutate(age = 312),
                   aes(x = age, y = feature_resids, fill = group),
                   geom = "crossbar", fun = "mean", width = 4, linewidth = 0.33
      ) +
      scale_x_continuous(breaks = c(20, 63, 300))
    
  }
  
  if (include_individuals == TRUE) {
    
    p <- p +
      geom_point(size = 0.5, alpha = 0.15) +
      geom_point(size = 0.5, shape = 1) +
      geom_line(aes(group = label), alpha = 0.3, lty = 2, linewidth = 0.5) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 1)
    
  }
  
  if (legend == TRUE) {
    
    return(p)
    
  } else if (legend == FALSE) {
    
    return(p +  
             theme(legend.position = "none",
                   strip.text = element_text(color = "black", size = 14), 
                   strip.background = element_rect(fill = "white", color = "gray"),
                   axis.title = element_text(size = 14),
                   axis.text = element_text(size = 14) 
             )
           
    )
    
  }
  
  
}

# MODEL LABEL
f_model_label <- function(model, type = c("linear", "quadratic", "cubic", "absolute")) {
  
  adj_r_sq <- summary(model)$adj.r.squared
  f <- summary(model)$fstatistic
  p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  coefs <- coefficients(model) %>% round(2)
  
  if (type == "linear") { formula <- paste0("y = ", coefs[1], " + ", coefs[2], "*x")
  } else if (type == "quadratic") {formula <- paste0("y = ", coefs[1], " + ", coefs[2], "*x + ", coefs[3], "*x^2")
  } else if (type == "cubic") {formula <- paste0("y = ", coefs[1], " + ", coefs[2], "*x + ", coefs[3], "*x^2", coefs[4], "*x^3")
  } else if (type == "absolute") {formula <-paste0("y = ", coefs[1], " + ", coefs[2], "*|x|")}
  
  return(paste0(formula, 
                "\nR^2 = ", round(adj_r_sq, 2), 
                "; p = ", format(p_val, scientific = TRUE, digits = 3)
  )
  )
  
}

# load data ---------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

load(paste0(base_dir, "objects/25Aug2023_df_data.RDS")) # df_data
load(paste0(base_dir, "objects/25Aug2023_sigma_labels.RDS")) # df_sigma_labels
load(paste0(base_dir, "objects/25Aug2023_sigma_atlas_for_plotting.RDS")) # df_sigma_atlas
load(paste0(base_dir, "objects/29Aug2023_df_mind_GM.RDS")) # df_mind_GM
load(paste0(base_dir, "objects/mt_degree_slopes.RDS")) # df_mt_degree_slopes

# info on if scan had gadolinium
df_gad <- read_xlsx(paste0(base_dir, "data/gadolinium_use_records.xlsx")) %>% 
  dplyr::select(subject, scan_date, gad_by_record)

# Behavior data
df_behavior_raw <- read_xlsx(paste0(base_dir, "data/GSK_E1_all_behaviour__collapsed_wide.xlsx"), 
                             sheet = 1) %>% 
  clean_names() %>% 
  pivot_longer(2:ncol(.), names_to = "task", values_to = "value") %>% 
  mutate(value = as.numeric(value)) %>% 
  
  # add task meta data
  left_join(
    
    read_xlsx(paste0(base_dir, "data/GSK_E1_all_behaviour__collapsed_wide.xlsx"), sheet = 2) %>% 
      dplyr::rename("task" = "ethan code") %>% 
      dplyr::select(-ID, -meaning) %>% 
      mutate(task = tolower(task)),
    by = join_by(task)
    
  ) %>% 
  
  # add subject info
  left_join(
    df_data %>% 
      filter(study == "GSK") %>% 
      dplyr::select(subject, sex, group) %>% 
      distinct() %>% 
      mutate(id = str_remove(subject, "EDAA") %>% as.numeric),
    by = join_by(id)
  ) %>% 
  dplyr::select(subject, sex, group, task_group, label_for_plotting, task, value) 

# VT = Video Tracking (anxiety tests)
# VS = Video Scoring (videos of animals in shock boxes)
# SPT = Sucrose Preference Tests 
#   (0.5%, 1%, 2% sucrose pre-stress, 
#   as well as an average of these three, 
#   plus 1% sucrose post-stress under "SPT_pref_s")
# FRPR = Fixed Ratio / Progressive Ratio
# "PRL_first_7" refers to the first 7 PRL sessions i.e. pre-stress sessions
# "PRL_stress" refers to the three PRL sessions during the stress period. 
# Mean is the average over sessions while "beta" is just the slope of a linear model over those sessions

# SIGMA atlas ROI centroids -----------------------------------------------

# CALCULATE CENTROID OF EACH SYSTEM FOR PLOTTING
df_centroids <- df_sigma_atlas %>% 
  filter(side == "axial" & hemisphere == "right") %>% # axial for x & y coords
  group_by(region_of_interest) %>%
  summarize(geometry = st_union(geometry)) %>% 
  st_centroid %>% 
  mutate(point = map(geometry, ~st_coordinates(.x) %>% as.data.frame %>% as_tibble)) %>% 
  unnest(cols = c(point)) %>% 
  as_tibble() %>% 
  dplyr::select(-geometry) %>% 
  clean_names() %>% 
  
  left_join(
    df_sigma_atlas %>% 
      filter(side == "coronal" & hemisphere == "right") %>% # coronal for x & z coords
      group_by(region_of_interest) %>%
      summarize(geometry = st_union(geometry)) %>% 
      st_centroid %>% 
      mutate(point = map(geometry, ~st_coordinates(.x) %>% as.data.frame %>% as_tibble)) %>% 
      unnest(cols = c(point)) %>% 
      as_tibble() %>% 
      dplyr::select(-geometry) %>% 
      clean_names() %>% 
      dplyr::rename("z" = "y"),
    by = join_by(region_of_interest)
  ) %>% 
  
  left_join(
    df_sigma_atlas %>% 
      filter(side == "sagittal" & hemisphere == "right") %>% # sagittal for y & z coords
      group_by(region_of_interest) %>%
      summarize(geometry = st_union(geometry)) %>% 
      st_centroid %>% 
      mutate(point = map(geometry, ~st_coordinates(.x) %>% as.data.frame %>% as_tibble)) %>% 
      unnest(cols = c(point)) %>% 
      as_tibble() %>% 
      dplyr::select(-geometry) %>% 
      clean_names() %>% 
      dplyr::rename("y" = "x", "z" = "y"),
    by = join_by(region_of_interest)
  ) %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  mutate(dim = substr(dim, start = 1, stop = 1)) %>% 
  group_by(region_of_interest, dim) %>%
  summarise(coord = mean(coord)) %>% 
  pivot_wider(id_cols = region_of_interest, names_from = dim, values_from = coord) %>% 
  ungroup


# format data -------------------------------------------------------------

# LIST ROIS THAT DID NOT REGISTER WELL TO REMOVE FROM ANALYSIS
bad_reg_rois <- c("spinal_cord", "brainstem", "cerebell", "commissural_stria_terminalis", "central_canal")

# CONVERT ROI TO SYSTEM LEVEL
df_sys_to_roi <- df_sigma_labels %>% 
  dplyr::select(matter, system, region_of_interest) %>% 
  distinct()

# LIST GRAY VS WHITE MATTER REGIONS
gray_matter_rois <- df_sigma_labels %>% 
  filter(matter == "grey_matter") %>% 
  pull(region_of_interest)

white_matter_rois <- df_sigma_labels %>% 
  filter(matter == "white_matter") %>% 
  pull(region_of_interest)

# ADD SYSTEM LEVEL TO df_mind_GM & REMOVE SCANS WITH GADOLINIUM
df_mind_GM <- df_mind_GM %>% 
  left_join(df_sys_to_roi %>% filter(region_of_interest %in% gray_matter_rois) %>% dplyr::rename("R1" = "region_of_interest")) %>% 
  dplyr::rename("S1" = "system") %>% 
  left_join(df_sys_to_roi %>% filter(region_of_interest %in% gray_matter_rois) %>% dplyr::rename("R2" = "region_of_interest")) %>% 
  dplyr::rename("S2" = "system") %>% 
  
  # remove scans with gadolinium stain
  anti_join(df_gad %>% filter(gad_by_record == 1), 
            by = join_by(subject, scan_date)
  ) %>% 
  filter(!is.na(sex))

# PULL ROI ABBREVIATIONS
df_roi_to_abbrev <- df_sigma_labels %>% 
  dplyr::select(region_of_interest, roi_abbreviation) %>% 
  distinct() #%>% 
# mutate(region_of_interest = factor(region_of_interest, levels = c(roi_order))) %>% 
# arrange(region_of_interest)

# FORMAT MT DATA
df_mt <- df_data %>% 
  
  # filter for ROI MT values calculated from MT scans
  filter(metric == "MTR") %>% 
  dplyr::select(-metric) %>% 
  dplyr::rename("mt" = "value") %>% 
  
  # reorder df
  dplyr::select(study, subject, sex, group, timepoint, age, tbv, scan_date,
                hemisphere, matter, territories, system, region_of_interest, roi_abbreviation, 
                mt) %>% 
  
  # average across hemispheres
  group_by(subject, timepoint, region_of_interest) %>% 
  mutate(mt = mean(mt)) %>% 
  dplyr::select(-hemisphere) %>% 
  distinct() %>% 
  
  # change grouping to correct & normalize within study, ROI, and timepoint
  group_by(study, region_of_interest, timepoint) %>% 
  nest() %>% 
  
  # take residuals (for plotting)
  mutate(
    
    model = case_when(
      
      study == "MRC" & timepoint != 300 ~ "mt ~ tbv",
      study == "MRC" & timepoint == 300 ~  "mt ~ age + tbv",
      
      study == "GSK" & timepoint == 300 ~ "mt ~ sex + age + tbv",
      TRUE ~ "mt ~ sex + tbv" 
      
    ),
    
    data = map(
      
      .x = data,
      .f = ~ .x %>% 
        mutate(feature_resids = lm(formula(model), data = .x)$residuals + mean(mt))
      
    )
    
  ) %>% 
  
  # normalize and remove outliers (defined as z_score > 4)
  mutate(data = 
           
           map(
             .x = data, 
             .f = ~ .x %>% 
               mutate(normalized_feature = (mt - mean(mt))/sd(mt))
           )
         
  ) %>% 
  unnest(cols = c(data)) %>% 
  filter(abs(normalized_feature) <= 4) %>% 
  
  # make sure timepoint is a factor ordered correctly
  mutate(timepoint = factor(timepoint, levels = c(20, 35, 63, 300)))

# GENERATE GRAPH FROM EACH CMN AND CALCULATE DEGREE & HUB SCORE
df_degree <- df_mind_GM %>% 
  group_by(study, subject, timepoint, sex, group, scan_date, age, tbv) %>% 
  nest() %>% 
  mutate(mind_graph = map(data, ~ graph_from_data_frame(.x, directed = FALSE)),
         hubs = map(
           .x = mind_graph,
           .f = ~ hub_score(.x)$vector %>% 
             enframe %>% 
             dplyr::rename("hub_score" = "value") %>% 
             mutate(weighted_degree = strength(.x)) %>% 
             arrange(-hub_score) %>% 
             dplyr::rename("region_of_interest" = "name")
         )
  ) %>% 
  unnest(cols = c(hubs)) %>% 
  
  # change grouping to correct & normalize within study
  group_by(study, region_of_interest, timepoint) %>% 
  nest() %>% 
  
  # take residuals for sex + age + tbv (for plotting)
  mutate(
    
    model = case_when(
      study == "MRC" & timepoint != 300 ~ "weighted_degree ~ tbv",
      study == "MRC" & timepoint == 300 ~ "weighted_degree ~ age + tbv",
      
      study == "GSK" & timepoint == 300 ~ "weighted_degree ~ sex + age + tbv",
      TRUE ~ "weighted_degree ~ sex + tbv" 
    ),
    
    data = map(
      .x = data,
      .y = model,
      .f = ~ .x %>% 
        mutate(feature_resids = lm(formula(.y), data = .x)$residuals + mean(weighted_degree))
    ) 
  ) %>% 
  
  # normalize and remove outliers
  mutate(data = 
           map(
             .x = data, 
             .f = ~ .x %>% 
               mutate(normalized_feature = (weighted_degree - mean(weighted_degree))/sd(weighted_degree))
           )
  ) %>% 
  unnest(cols = c(data)) %>% 
  filter(abs(normalized_feature) <= 4) %>% 
  
  # add system-level annotation
  left_join(df_sys_to_roi, by = join_by(region_of_interest)) %>% 
  
  # make sure timepoint is a factor ordered correctly
  mutate(timepoint = factor(timepoint, levels = c(20, 35, 63, 300)))

# BEHAVIOR DATA
f_impute_median <- function(x, na.rm = TRUE) (replace(x, is.na(x), median(x, na.rm = na.rm)))

df_behavior <- df_behavior_raw %>% 
  
  filter(subject != "EDAA22") %>% # no behavior data available
  group_by(task) %>% 
  
  # impute missing values
  mutate(value = ifelse(is.na(value), f_impute_median(value), value)) %>% 
  
  nest() %>% 
  mutate(
    data = map(
      .x = data,
      .f = ~ .x %>%
        
        # correct for sex
        mutate(feature_resids = lm(value ~ sex, data = .x)$residuals + mean(value)) %>% 
        
        # normalize
        mutate(normalized_feature = (value - mean(value))/sd(value))
    )
  ) %>%
  unnest(cols = c(data)) %>% 
  
  # remove outliers
  filter(abs(normalized_feature) <= 4) %>% 
  
  # remove 'beta' tasks (don't know what these are)
  filter(!str_detect(task, "beta")) %>% 
  
  # remove total_dist task
  filter(!str_detect(task, "total_dist")) %>% 
  
  # select only the most extreme FRPR
  filter(!str_detect(task, "pr4|pr8"))


# MT analysis ----------------------------------------------------------------

df_mt_gsk <- df_mt %>% filter(study == "GSK")

# BASELINE
df_mt_pnd20 <- df_mt_gsk %>%
  filter(timepoint == 20 & region_of_interest %in% gray_matter_rois) %>%
  group_by(system) %>%
  nest() %>%
  f_run_lm(20) %>% 
  filter(str_detect(term, "group"))

# DELTA
df_mt_lme <- df_mt_gsk %>% 
  filter(timepoint %in% c(20, 63) & region_of_interest %in% gray_matter_rois) %>% 
  group_by(system) %>% 
  nest() %>% 
  f_run_lme %>% 
  filter(str_detect(term, "group"))

# RELATE CASE-CONTROL DIFFERENCE IN BASELINE MT TO MRC BASELINE MT
df_baseline_baseline_comparison <- df_mt_degree_slopes %>% 
  unnest(cols = c(data)) %>% 
  ungroup %>% 
  filter(phenotype == "mt" & timepoint == 20 & study == "MRC" & region_of_interest %in% gray_matter_rois) %>% 
  dplyr::select(system, region_of_interest, subject, feature_resids) %>% 
  distinct() %>% 
  group_by(region_of_interest) %>%
  summarise(feature_resids = median(feature_resids)) %>%
  left_join(df_mt_gsk %>%
              filter(timepoint == 20 & region_of_interest %in% gray_matter_rois) %>%
              group_by(region_of_interest) %>%
              nest() %>%
              f_run_lm(20) %>% 
              filter(term == "groupMS") %>% 
              dplyr::select(region_of_interest, t_value, p_adj),
            by = join_by(region_of_interest)
  ) %>% 
  arrange(t_value) %>% 
  mutate(feature_resids2 = feature_resids^2,
         feature_resids3 = feature_resids^3)

lm(t_value ~ feature_resids, data = df_baseline_baseline_comparison) %>% summary
lm(t_value ~ feature_resids + feature_resids2, data = df_baseline_baseline_comparison) %>% summary
lm(t_value ~ feature_resids + feature_resids2 + feature_resids3, data = df_baseline_baseline_comparison) %>% summary

# RELATE CASE-CONTROL DIFFERENCE IN BASELINE MT TO MRC BASELINE MT
df_baseline_slopes_comparison <- df_mt_degree_slopes %>% 
  unnest(cols = c(data)) %>% 
  ungroup %>% 
  filter(phenotype == "mt" & period == "early" & study == "MRC" & region_of_interest %in% gray_matter_rois) %>% 
  dplyr::select(system, region_of_interest, subject, slope) %>% 
  distinct() %>% 
  group_by(region_of_interest) %>%
  summarise(slope = median(slope)) %>%
  left_join(df_mt_gsk %>%
              filter(timepoint == 20 & region_of_interest %in% gray_matter_rois) %>%
              group_by(region_of_interest) %>%
              nest() %>%
              f_run_lm(20) %>% 
              filter(term == "groupMS") %>% 
              dplyr::select(region_of_interest, t_value, p_adj),
            by = join_by(region_of_interest)
  ) %>% 
  arrange(t_value) %>% 
  mutate(t_value2 = t_value^2,
         t_value3 = t_value^3)

lm(slope ~ t_value, data = df_baseline_slopes_comparison) %>% summary
lm(slope ~ t_value + t_value2, data = df_baseline_slopes_comparison) %>% summary
lm(slope ~ t_value + t_value2 + t_value3, data = df_baseline_slopes_comparison) %>% summary


# behavior analysis -------------------------------------------------------

# CALCULATE GROUP DIFFERENCES IN EACH TASK
df_behavior_effects <- df_behavior %>% 
  group_by(task_group, label_for_plotting, task) %>% 
  nest() %>% 
  
  # run LM on normalized behavior score, correcting for sex
  mutate(
    
    lm_res = map(
      
      .x = data,
      .f = ~ lm(normalized_feature ~ group + sex, data = .x) %>% 
        tidy %>% 
        filter(term != "(Intercept)") %>% 
        clean_names
      
    )
    
  ) %>% 
  unnest(cols = c(lm_res)) %>% 
  ungroup %>% 
  
  # adjust p-values across all tasks & contrasts
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
  arrange(p_value) %>% 
  
  # look specifically at group differences (don't care about sex here)
  filter(term == "groupMS") %>% 
  mutate(label = case_when(
    
    p_value < 0.05 & p_adj > 0.05 ~ "*",
    p_adj <= 0.05 ~ "**",
    TRUE ~ ""
    
  )
  )

# POST-HOC: PRL differences pre and during stress in control vs MS
prl_control_p <- wilcox.test(
  df_behavior %>% filter(task == "prl_first_7_mean_latency_to_respond" & group == "control") %>% pull(feature_resids),
  df_behavior %>% filter(task == "prl_stress_session_3_predicted_latency_to_respond_mean" & group == "control") %>% pull(feature_resids)
) %>% .$p.value

prl_ms_p <- wilcox.test(
  df_behavior %>% filter(task == "prl_first_7_mean_latency_to_respond" & group == "MS") %>% pull(feature_resids),
  df_behavior %>% filter(task == "prl_stress_session_3_predicted_latency_to_respond_mean" & group == "MS") %>% pull(feature_resids)
) %>% .$p.value

# CORRELATE SIGNIFICANT BEHAVIOR WITH FEATURE OF EACH REGION

sig_tasks <- df_behavior_effects %>% 
  filter(p_value < 0.05 & !str_detect(task, "vs")) %>% 
  pull(task) %>% 
  unique

## correlate MT, CMN degree in each system of each group at each timepoint with each behavior task

df_metric_behavior_cor <- df_mt_gsk %>% 
  filter(timepoint == 20) %>% # focus on baseline associations
  
  # take median of residuals in each region_of_interest in each subject at each timepoint
  group_by(timepoint, region_of_interest, group, subject) %>% 
  dplyr::rename("phenotype_resids" = "feature_resids") %>% 
  mutate(phenotype = "mt", .before = 1) %>% 
  
  # add behavior data
  left_join(df_behavior, by = join_by(subject, group)) %>% 
  
  # only look at tasks that are significantly different between groups
  filter(task %in% sig_tasks) %>% 
  
  # correlate across subjects
  group_by(timepoint, group, region_of_interest, phenotype, task_group, label_for_plotting) %>% 
  nest() %>% 
  
  mutate(
    
    r = map(
      .x = data,
      .f = ~cor.test(.x$feature_resids, .x$phenotype_resids)$estimate
    ),
    
    p_value = map(
      .x = data,
      .f = ~cor.test(.x$feature_resids, .x$phenotype_resids)$p.value
    )

  ) %>% 
  
  unnest(cols = c(r, p_value)) %>% 
  ungroup %>% 
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
  arrange(p_value)


# A: Case-control baseline MT ---------------------------------------------

f_plot_brain(df = df_mt_pnd20,
             contrast = "groupMS",
             labels = FALSE,
             slices = list("sagittal" = 95, 
                           "axial" = 115, 
                           "coronal" = 150),
             legend = TRUE#,
             #column = TRUE
) +
  labs(title = "RMS-control differences at baseline (PND 20)") +
  theme(legend.position = "bottom")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6a.case_control_baseline", .x), width = 7, height = 5)
)

# B: Case-control delta MT ------------------------------------------------

f_plot_brain(df = df_mt_lme,
             contrast = "groupMS:age",
             p_adj = TRUE,
             labels = FALSE,
             slices = list("sagittal" = 95, 
                           "axial" = 115, 
                           "coronal" = 150),
             legend = FALSE#,
             #column = TRUE
) +
  labs(title = "RMS-control differences in early development MT decay") +
  theme(legend.position = "bottom")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6b.case_control_MT_decay", .x), width = 7, height = 5)
)

# B2: Anatomical patterning of baseline MT and MT slope effect size  --------

labels <- c("x (medial --> lateral)", "y (posterior --> anterior)", "z (inferior --> superior)")
names(labels) <- c("x", "y", "z")

df_centroids %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  left_join(df_mt_gsk %>%
              filter(timepoint == 20 & region_of_interest %in% gray_matter_rois) %>%
              group_by(region_of_interest) %>%
              nest() %>%
              f_run_lm(20) %>% 
              filter(str_detect(term, "group")) %>% 
              dplyr::select(region_of_interest, t_value) %>% 
              mutate(feature = "Baseline MT"), 
            by = join_by(region_of_interest)
  ) %>% 
  bind_rows(
    df_centroids %>% 
      pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
      left_join(df_mt_gsk %>% 
                  filter(timepoint %in% c(20, 63) & region_of_interest %in% gray_matter_rois) %>% 
                  group_by(region_of_interest) %>% 
                  nest() %>% 
                  f_run_lme %>% 
                  filter(term == "groupMS:age") %>%
                  dplyr::select(region_of_interest, t_value) %>% 
                  mutate(feature = "MT slope"), 
                by = join_by(region_of_interest)
      )
  ) %>% 
  mutate(sign = ifelse(t_value < 0, "negative", "positive"),
         hub = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor) %>% 
  filter(!is.na(feature)) %>% 
  
  ggplot(aes(x = coord, y = t_value)) +
  geom_point(aes(fill = t_value, color = sign, shape = hub, size = hub)) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_cor(color = "black", size = 3, label.y = -Inf, vjust = 0) +
  #geom_text_repel(aes(label = system)) +
  facet_grid(~ feature + dim, scales = "free",
             labeller = labeller(dim = labels)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(1, 3)) +
  labs(x = "Coordinate position", y = "Pearsons' r",
       title = "Baseline MT and MT slope effect size relative to anatomical position (RH)") +
  theme(legend.position = "none")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6b2.MT_effect_size_coordinates", .x), width = 12, height = 3)
)


# C: Case-control baseline MT vs MRC baseline MT --------------------------

model <- lm(t_value ~ feature_resids + feature_resids2, data = df_baseline_baseline_comparison)

df_baseline_baseline_comparison %>% 
  mutate(sign = ifelse(t_value < 0, "negative", "positive"),
         formula = f_model_label(model, type = "quadratic"),
         hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor(),
         label = ifelse(hubs == 1, str_replace_all(region_of_interest, "_", " "), NA)
  ) %>% 
  
  ggplot(aes(x = t_value, y = feature_resids)) +
  geom_point(aes(fill = t_value, color = sign, shape = hubs, size = hubs)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black") +
  geom_text(
    x = 1, y = 0.9, color = "blue", size = 2.5, hjust = 0.5,
    aes(label = formula), check_overlap = TRUE
  ) +
  geom_text_repel(aes(label = label), size = 3) +
  geom_vline(aes(xintercept = 0)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", guide = "none") +
  scale_color_manual(values = c("#2166AC", "#B2182B"), guide = "none") +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(1, 3)) +
  labs(y = "Baseline MT in developmental cohort \n(corrected for TBV)",
       x = "RMS-control t-value at baseline",
       title = "Relationship between group effect size and baseline MT (at ROI level)") +
  coord_flip()

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6c.case_control_baseline_MRC_baseline", .x), width = 5, height = 4)
)

# D: Case-control baseline MT vs MRC slopes MT --------------------------

model <- lm(slope ~ t_value + t_value2, data = df_baseline_slopes_comparison)

df_baseline_slopes_comparison %>% 
  mutate(sign = ifelse(t_value < 0, "negative", "positive"),
         formula = f_model_label(model, type = "quadratic"),
         hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor(),
         label = ifelse(hubs == 1, str_replace_all(region_of_interest, "_", " "), NA)
  ) %>% 
  
  ggplot(aes(x = t_value, y = slope)) +
  geom_point(aes(fill = t_value, color = sign, shape = hubs, size = hubs)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black") +
  geom_text(
    x = -0.5, y = -0.004, color = "blue", size = 2.5, hjust = 0.5,
    aes(label = formula), check_overlap = TRUE
  ) +
  geom_text_repel(aes(label = label), size = 3) +
  geom_vline(aes(xintercept = 0)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", guide = "none") +
  scale_color_manual(values = c("#2166AC", "#B2182B"), guide = "none") +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(1, 3)) +
  labs(y = "Early developmental MT slope in developmental cohort \n(corrected for TBV)",
       x = "RMS-control t-value at baseline",
       title = "Relationship between group effect size and MT decay slope (at ROI level)") +
  coord_flip()

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6d.case_control_baseline_MRC_slope", .x), width = 5, height = 4)
)



# E: Case-control behavior differences ------------------------------------

# case-control differences
df_behavior_effects %>% 
  unnest(cols = c(data)) %>% 
  filter(str_detect(task_group, "FRPR|PRL")) %>% 
  mutate(label_for_plotting = str_remove(label_for_plotting, "PR16")) %>% 
  filter(p_value < 0.05) %>% 
  
  ggplot(aes(x = normalized_feature, 
             y = task_group, 
             color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             alpha = 0.5) +
  labs(y = "", x = "Normalized score",
       title = "Differential response to FRPR and PRL latency to respond in RMS vs controls") +
  scale_color_manual(values = group_cols) +
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12)
  )

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6e.behavior_task_case-control_differences", .x), width = 5, height = 3)
)

# correlate tasks
df_cor_tasks <- df_behavior_effects %>% 
  unnest(cols = c(data)) %>% 
  filter(str_detect(task_group, "FRPR|PRL") & p_value < 0.05) %>% 
  pivot_wider(id_cols = c(subject, group), names_from = task_group, values_from = feature_resids)

(df_cor_tasks %>% 
    ggplot(aes(x = FRPR, y = PRL)) +
    geom_point(aes(fill = group), shape = 21, size = 2) +
    geom_smooth(method = "lm", color = "black") +
    stat_cor() +
    scale_fill_manual(values = group_cols, guide = "none")
  ) /
    
    (df_cor_tasks %>% 
       ggplot(aes(x = FRPR, y = `PRL post-stress`)) +
       geom_point(aes(fill = group), shape = 21, size = 2) +
       geom_smooth(method = "lm", color = "black") +
       stat_cor() +
       scale_fill_manual(values = group_cols, guide = "none")
    )  /
    
    (df_cor_tasks %>% 
       ggplot(aes(x = `PRL post-stress`, y = PRL)) +
       geom_point(aes(fill = group), shape = 21, size = 2) +
       geom_smooth(method = "lm", color = "black") +
       stat_cor() +
       scale_fill_manual(values = group_cols, guide = "none")
    ) + 
  plot_annotation(title = "Correlation between latency to respond across tasks")
  
# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6f.correlate_behavior_tasks", .x), width = 3, height = 8)
)


# F: task correlations with system MT at baseline -------------------------

sig_systems <- df_metric_behavior_cor %>% filter(p_adj < 0.05) %>% pull(region_of_interest) %>% unique

df_metric_behavior_cor_plot <- df_metric_behavior_cor %>% 
  filter(task_group == "FRPR") %>% 
  unnest(cols = c(data)) %>% 
  dplyr::select(group, region_of_interest, r, p_adj) %>% 
  distinct()

slices <- list("sagittal" = 95, 
               "axial" = 115, 
               "coronal" = 150)

df_slices <- map_dfr(.x = 1:length(slices),
                     .f = ~ df_sigma_atlas %>%
                       filter(region_of_interest != "striatum") %>% 
                       filter(side == names(slices[.x]) & slice == slices[.x])
                     
) %>% 
  dplyr::select(region_of_interest, side, slice, geometry) %>% 
  distinct()

df_na_regs <- df_slices %>% 
  filter(str_detect(region_of_interest, paste0(bad_reg_rois, collapse = "|"))) %>% 
  expand_grid(group = c("control", "MS"))

df_slices %>% 
  left_join(df_metric_behavior_cor_plot, 
            by = "region_of_interest") %>%
  arrange(p_adj) %>% 
  filter(!is.na(group)) %>% 
  bind_rows(df_na_regs, .) %>% 
  filter(region_of_interest %in% gray_matter_rois) %>% 
  
  ggplot() +
  geom_sf(aes(fill = r, color = r, geometry = geometry, group = -1), 
          lwd = 0.5) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  scale_color_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  facet_wrap(group ~ side, nrow = 2) +
  labs(title = "FRPR latency to respond - MT correlation map in RMS vs controls") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "bottom"
  )

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6g.behavior_task_case-control_relationships_map", .x), width = 10, height = 6)
)


# H: Relationship between case-control effect size at baseline and FRPR-MT correlation in controls vs RMS -----------------------

df_fig6h <- df_metric_behavior_cor %>% 
  filter(task_group == "FRPR") %>% 
  dplyr::select(group, region_of_interest, r) %>% 
  left_join(df_mt_gsk %>%
              filter(timepoint == 20 & region_of_interest %in% gray_matter_rois) %>%
              group_by(region_of_interest) %>%
              nest() %>%
              f_run_lm(20) %>% 
              filter(str_detect(term, "group")) %>% 
              dplyr::select(region_of_interest, t_value), 
            by = join_by(region_of_interest)
  ) %>% 
  mutate(sign = ifelse(t_value < 0, "neg", "pos"))

control_model <- lm(r ~ abs(t_value), data = df_fig6h %>% filter(group == "control"))
ms_model <- lm(r ~ abs(t_value), data = df_fig6h %>% filter(group == "MS"))

# BASELINE
df_fig6h %>% 
  
  mutate(
    label = ifelse(group == "control",
                   f_model_label(control_model, type = "absolute"),
                   f_model_label(ms_model, type = "absolute"))
  ) %>% 
  
  ggplot(aes(x = t_value, y = r)) +
  geom_point(aes(fill = t_value, color = sign), shape = 21, size = 2) +
  facet_wrap(vars(group), nrow = 1) +
  geom_smooth(method = "lm", formula = y ~ x + I(abs(x)), color = "black") +
  geom_text(aes(label = label), x = 0, y = - 0.75, hjust = 0.5, check_overlap = TRUE) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = "RMS-control baseline MT t-value", y = "Pearson's r",
       title = "Relationship between MT baseline effect size and \nFRPR latency to respond - MT correlation (ROI-level)") +
  theme(legend.position = "none")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6h.baseline_effect_size_FRPR_cor", .x), width = 6, height = 3)
)

# I: Relationship between case-control slope effect size and FRPR-MT correlation in controls vs RMS -----------------------

df_fig6i <- df_metric_behavior_cor %>% 
  filter(task_group == "FRPR") %>% 
  dplyr::select(group, region_of_interest, r) %>% 
  left_join(df_mt_gsk %>%
              filter(timepoint %in% c(20, 63) & region_of_interest %in% gray_matter_rois) %>%
              group_by(region_of_interest) %>%
              nest() %>%
              f_run_lme() %>% 
              filter(term == "groupMS:age") %>% 
              dplyr::select(region_of_interest, t_value), 
            by = join_by(region_of_interest)
  ) %>% 
  mutate(sign = ifelse(t_value < 0, "neg", "pos"))

control_model <- lm(r ~ abs(t_value), data = df_fig6i %>% filter(group == "control"))
ms_model <- lm(r ~ abs(t_value), data = df_fig6i %>% filter(group == "MS"))

# BASELINE
df_fig6i %>% 
  
  mutate(
    label = ifelse(group == "control",
                   f_model_label(control_model, type = "absolute"),
                   f_model_label(ms_model, type = "absolute"))
  ) %>% 
  
  ggplot(aes(x = t_value, y = r)) +
  geom_point(aes(fill = t_value, color = sign), shape = 21, size = 2) +
  facet_wrap(vars(group), nrow = 1) +
  geom_smooth(method = "lm", formula = y ~ x + I(abs(x)), color = "black") +
  geom_text(aes(label = label), x = -0.5, y = - 0.65, hjust = 0.5, check_overlap = TRUE) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = "RMS-control baseline MT t-value", y = "Pearson's r",
       title = "Relationship between MT slope effect size and \nFRPR latency to respond - MT correlation (ROI-level)") +
  theme(legend.position = "none")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6i.slope_effect_size_FRPR_cor", .x), width = 6, height = 3)
)


# J: Anatomical patterning of behavior-MT correlations in control  --------

labels <- c("x (medial --> lateral)", "y (posterior --> anterior)", "z (inferior --> superior)")
names(labels) <- c("x", "y", "z")

df_centroids %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  left_join(df_metric_behavior_cor %>% 
              filter(task_group == "FRPR") %>% 
              dplyr::select(group, region_of_interest, r), 
            by = join_by(region_of_interest)
  ) %>% 
  mutate(sign = ifelse(r < 0, "negative", "positive")) %>% 
  filter(!is.na(group)) %>% 
  
  ggplot(aes(x = coord, y = r)) +
  geom_point(aes(fill = r, color = sign), shape = 21, size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_cor(color = "black", size = 3, label.y = -Inf, vjust = 0) +
  #geom_text_repel(aes(label = system)) +
  facet_grid(group ~ dim, scales = "free",
             labeller = labeller(dim = labels)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = "Coordinate position", y = "Pearsons' r",
       title = "ROI-level MT correlation with FRPR latency to respond relative to anatomical position (RH)") +
  theme(legend.position = "none")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/6j.FRPR_cor_coordinates", .x), width = 6, height = 4)
)

# S5a: Baseline MT boxplot -----------------------------------------------------------------------

f_plot_box(df = df_mt_pnd20,
           xlab = "", 
           ylab = "MT \n(corrected for sex + TBV)", 
           title = "Baseline regional MT profiles", 
           legend = FALSE) +
  coord_flip()

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S5a.baseline_MT_boxplot", .x), width = 4, height = 3)
)

# S5b: Overall MT decay differences ---------------------------------------

df_mt_gsk %>% 
  group_by(timepoint, group) %>% 
  mutate(median_feaure_resids = median(feature_resids),
         median_age = median(age)) %>% 
  
  ggplot() +
  geom_point(aes(x = age, y = feature_resids), color = "black", size = 2, shape = 1) +
  geom_point(aes(x = median_age, y = median_feaure_resids, color = group), size = 2) +
  geom_line(aes(x = median_age, y = median_feaure_resids, color = group), linewidth = 0.5, lty = 2) +
  scale_color_manual(values = group_cols) +
  labs(x = "age (PND)", y = "Regional MT \n(correct for TBV + sex)",
       title = "MT across ROIs in each group \n(line shows overall median)") +
  theme(legend.position = "none")

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S5b.overall_MT_decay", .x), width = 4, height = 3)
)

# S5c: ROI-specific MT decay differences ----------------------------------

f_plot_traj(df = df_mt_lme,
                        timepoints = c(20, 63),
                        xlab = "age (PND)", 
                        ylab = "MT (corrected for sex + TBV)", 
                        title = "Differential MT trajectories in early development", 
                        dim = c(3, 4),
                        lnwdth = 0.75,
                        legend = FALSE,
                        include_individuals = TRUE
)

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S5c.ROI_MT_decay", .x), width = 8, height = 7)
)


# S5d: Group differences in behavior tasks --------------------------------

df_behavior_effects %>% 
  unnest(cols = c(data)) %>% 
  mutate(label_for_plotting = str_remove(label_for_plotting, "PR16")) %>% 
  
  ggplot(aes(x = normalized_feature, 
             y = reorder(str_wrap(label_for_plotting, 20), abs(statistic)), 
             color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             alpha = 0.5, size = 2) +
  geom_text(aes(label = label, x = -3), color = "black", size = 10, vjust = 0.75) +
  
  facet_wrap(vars(task_group), scales = "free_y") +
  labs(y = "", x = "Normalized score",
       title = "Normalized behavior task scores in RMS vs controls") +
  scale_color_manual(values = group_cols) +
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12)
  )

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S5d.group_diffs_behavior", .x), width = 12, height = 8)
)

# S5e: behavior-feature correlations --------------------------------------

sig_systems <- df_metric_behavior_cor %>% filter(p_adj < 0.05) %>% pull(system) %>% unique

df_metric_behavior_cor %>% 
  filter(system %in% sig_systems & task_group == "FRPR") %>% 
  unnest(cols = c(data)) %>% 
  
  # take feature median across ROIs within system
  group_by(system, subject, group, task_group, feature_resids) %>% 
  summarise(phenotype_resids = median(phenotype_resids)) %>% 
  mutate(task_system = paste0(task_group, sep = " ", system)) %>% 
  
  ggplot(aes(x = phenotype_resids, y = feature_resids, color = group)) +
  geom_point(aes(fill = group), color = "black", size = 3, shape = 21) +
  geom_smooth(method = "lm", lty = 2) +
  stat_cor() +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  facet_wrap(~ factor(system, levels = sig_systems), scales = "free", nrow = 2) +
  labs(x = "MT \n(corrected for sex + TBV)", 
       y = "Behavior score \n(corrected for sex)",
       title = "MT and FRPR latency to respond correlation at PND 20"
  ) + theme(legend.position = "none")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S5e.behavior_task_case-control_relationships", .x), width = 12, height = 6)
)


