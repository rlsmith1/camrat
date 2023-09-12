
###
### Post-adult stress case-control differences in MT and degree & their relation to behavior
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
library(viridis)

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

### EDGES

# in each subject-scan, take the mean of edge weights and treat it as the value for that edge (like an ROI)
df_edge <- df_mind_GM %>% 
  group_by(study, subject, group, timepoint) %>% 
  
  filter(R1 != R2 & R1 %in% gray_matter_rois & R2 %in% gray_matter_rois) %>% 
  mutate(roi_edge = paste0(pmin(R1, R2), sep = "-", pmax(R1, R2)),
         sys_edge = paste0(pmin(S1, S2), sep = "-", pmax(S1, S2))) %>%
  group_by(subject, timepoint, roi_edge) %>% 
  slice(1) %>%
  ungroup %>%
  
  # take average weight for each subject-scan
  group_by(study, subject, timepoint, group, sex, age, tbv, sys_edge, roi_edge) %>% 
  summarise(weight = mean(weight)) %>% 
  
  # change grouping to correct & normalize within study
  group_by(study, roi_edge, timepoint) %>% 
  nest() %>% 
  
  # take residuals for sex + age + tbv (for plotting)
  mutate(
    
    model = case_when(
      
      study == "MRC" & timepoint != 300 ~ "weight ~ tbv",
      study == "MRC" & timepoint == 300 ~  "weight ~ age + tbv",
      
      study == "GSK" & timepoint == 300 ~ "weight ~ sex + age + tbv",
      TRUE ~ "weight ~ sex + tbv" 
      
    ),
    
    data = map(
      
      .x = data,
      .f = ~ .x %>% 
        mutate(feature_resids = lm(formula(model), data = .x)$residuals + mean(weight))
      
    )
    
  ) %>% 
  
  # normalize and remove outliers
  mutate(data = 
           
           map(
             .x = data, 
             .f = ~ .x %>% 
               mutate(normalized_feature = (weight - mean(weight))/sd(weight))
           )
         
  ) %>% 
  unnest(cols = c(data)) %>% 
  filter(abs(normalized_feature) <= 4)


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



# MT analysis -------------------------------------------------------------

df_mt_pnd300_SYS <- df_mt_gsk %>%
  filter(timepoint == 300 & region_of_interest %in% gray_matter_rois) %>%
  group_by(system) %>%
  nest() %>%
  f_run_lm(300) %>% 
  filter(str_detect(term, "group"))

df_mt_pnd300_ROI <- df_mt_gsk %>%
  filter(timepoint == 300 & region_of_interest %in% gray_matter_rois) %>%
  group_by(region_of_interest) %>%
  nest() %>%
  f_run_lm(300) %>% 
  filter(str_detect(term, "group"))


# MIND analysis -----------------------------------------------------------

# DEGREE
df_degree_pnd300_SYS <- df_degree_gsk %>%
  filter(timepoint == 300 & region_of_interest %in% gray_matter_rois) %>%
  group_by(system) %>%
  nest() %>%
  f_run_lm(300) %>% 
  filter(str_detect(term, "group"))

df_degree_pnd300_ROI <- df_degree_gsk %>%
  filter(timepoint == 300 & region_of_interest %in% gray_matter_rois) %>%
  group_by(region_of_interest) %>%
  nest() %>%
  f_run_lm(300) %>% 
  filter(str_detect(term, "group"))

# EDGE
df_edge_pnd300_SYS <- df_edge_gsk %>%
  filter(timepoint == 300) %>%
  group_by(sys_edge) %>%
  nest() %>%
  f_run_lm(300) %>% 
  filter(str_detect(term, "group"))

df_edge_pnd300_ROI <- df_edge_gsk %>%
  filter(timepoint == 300) %>%
  group_by(roi_edge) %>%
  nest() %>%
  f_run_lm(300) %>% 
  filter(str_detect(term, "group"))


# relate results to MRC cohort --------------------------------------------

df_MRC_relationships <- df_mt_degree_slopes %>% 
  unnest(cols = c(data)) %>% 
  ungroup %>% 
  filter(timepoint == 20 & period == "early" & study == "MRC" & region_of_interest %in% gray_matter_rois) %>% 
  group_by(phenotype, system, region_of_interest, slope) %>% 
  select(subject, feature_resids) %>% 
  distinct() %>% 
  summarise(feature_resids = median(feature_resids)) %>%
  pivot_longer(cols = c(slope, feature_resids), names_to = "feature", values_to = "value") %>% 
  mutate(feature = ifelse(feature == "feature_resids", "baseline", feature),
         feature = paste0(phenotype, sep = "_", feature)) %>% 
  pivot_wider(id_cols = c(system, region_of_interest), names_from = feature, values_from = value) %>% 
  
  left_join(
    df_mt_pnd300_ROI %>% 
      dplyr::select(region_of_interest, t_value) %>% 
      dplyr::rename("mt_t_value" = "t_value"),
    by = join_by(region_of_interest)
  ) %>% 
  left_join(
    df_degree_pnd300_ROI %>% 
      dplyr::select(region_of_interest, t_value) %>% 
      dplyr::rename("degree_t_value" = "t_value"),
    by = join_by(region_of_interest)
  )


# behavior analysis -------------------------------------------------------

df_metric_behavior_cor_SYS <- df_mt_gsk %>% 
  filter(timepoint == 300) %>% # focus on adult associations
  mutate(phenotype = "mt", .before = 1) %>% 
  
  bind_rows(
    df_degree_gsk %>% 
      filter(timepoint == 300) %>% # focus on adult associations
      mutate(phenotype = "degree", .before = 1)
  ) %>% 
  
  # take median of residuals in each region_of_interest in each subject at each timepoint
  group_by(phenotype, timepoint, system, group, subject) %>% 
  dplyr::rename("phenotype_resids" = "feature_resids") %>% 

  # add behavior data
  left_join(df_behavior, by = join_by(subject, group)) %>% 
  
  # only look at tasks that are significantly different between groups
  filter(task %in% sig_tasks) %>% 
  
  # correlate across subjects
  group_by(timepoint, group, system, phenotype, task_group, label_for_plotting) %>% 
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


# A: case-control adult MT --------------------------------------------

f_plot_brain(df = df_mt_pnd300_SYS,
             contrast = "groupMS",
             labels = FALSE,
             slices = list("sagittal" = 95, 
                           "axial" = 115, 
                           "coronal" = 180),
             legend = TRUE#,
             #column = TRUE
) +
  labs(title = "RMS-control MT differences following stress (PND 300)") +
  theme(legend.position = "bottom")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/8a.case_control_adult", .x), width = 7, height = 5)
)

# B: case-control adult degree --------------------------------------------

f_plot_brain(df = df_degree_pnd300_SYS,
             contrast = "groupMS",
             labels = FALSE,
             slices = list("sagittal" = 95, 
                           "axial" = 115, 
                           "coronal" = 180),
             legend = TRUE#,
             #column = TRUE
) +
  labs(title = "RMS-control degree differences following stress (PND 300)") +
  theme(legend.position = "bottom")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/8b.case_control_adult_degree", .x), width = 7, height = 5)
)


# A2B2: Anatomical patterning of baseline degree and degree slope effect size -------------------------------


labels <- c("x (medial --> lateral)", "y (posterior --> anterior)", "z (inferior --> superior)")
names(labels) <- c("x", "y", "z")

df_centroids %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  left_join(df_mt_pnd300_ROI %>% mutate(phenotype = "mt"), 
            by = join_by(region_of_interest)
  ) %>% 
  bind_rows(
    df_centroids %>% 
      pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
      left_join(df_degree_pnd300_ROI %>% mutate(phenotype = "degree"), 
                by = join_by(region_of_interest)
      )
  ) %>% 
  mutate(sign = ifelse(t_value < 0, "negative", "positive"),
         hub = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor) %>% 
  filter(!is.na(phenotype)) %>% 
  mutate(phenotype = factor(phenotype, levels = c("mt", "degree"))) %>% 
  
  ggplot(aes(x = coord, y = t_value)) +
  geom_point(aes(fill = t_value, color = sign, shape = hub, size = hub)) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_cor(color = "black", size = 3, label.y = -Inf, vjust = 0) +
  #geom_text_repel(aes(label = system)) +
  facet_grid(~ phenotype + dim, scales = "free",
             labeller = labeller(dim = labels)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(1, 3)) +
  labs(x = "Coordinate position", y = "Effect size",
       title = "Post-stress MT and degree effect size relative to anatomical position (RH)") +
  theme(legend.position = "none")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/8a2b2.adult_degree_mt_effect_size_coordinates", .x), width = 12, height = 3)
)


# C & D analysis: relationships with MRC -----------------------------------------------

### MT

# MT effect size and MT baseline (none)
df_MRC_relationships %>% 
  ggplot(aes(x = mt_baseline, y = mt_t_value)) +
  geom_point()
lm(mt_t_value ~ mt_baseline, data = df_MRC_relationships) %>% summary
lm(mt_t_value ~ mt_baseline + mt_baseline2, 
   data = df_MRC_relationships %>% mutate(mt_baseline2 = mt_baseline^2)) %>% summary
lm(mt_t_value ~ mt_baseline + mt_baseline2 + mt_baseline3, 
   data = df_MRC_relationships %>% mutate(mt_baseline2 = mt_baseline^2, mt_baseline3 = mt_baseline^3)) %>% summary

# MT effect size and MT slope (linear)
df_MRC_relationships %>% 
  ggplot(aes(x = mt_slope, y = mt_t_value)) +
  geom_point()
lm(mt_t_value ~ mt_slope, data = df_MRC_relationships) %>% summary # significant
lm(mt_t_value ~ mt_slope + mt_slope2, 
   data = df_MRC_relationships %>% mutate(mt_slope2 = mt_slope^2)) %>% summary
lm(mt_t_value ~ mt_slope + mt_slope2 + mt_slope3, 
   data = df_MRC_relationships %>% mutate(mt_slope2 = mt_slope^2, mt_slope3 = mt_slope^3)) %>% summary

# MT effect size and degree baseline (none)
df_MRC_relationships %>% 
  ggplot(aes(x = degree_baseline, y = mt_t_value)) +
  geom_point()
lm(mt_t_value ~ degree_baseline, data = df_MRC_relationships) %>% summary
lm(mt_t_value ~ degree_baseline + degree_baseline2, 
   data = df_MRC_relationships %>% mutate(degree_baseline2 = degree_baseline^2)) %>% summary

# MT effect size and degree slope (none)
df_MRC_relationships %>% 
  ggplot(aes(x = degree_slope, y = mt_t_value)) +
  geom_point()
lm(mt_t_value ~ degree_slope, data = df_MRC_relationships) %>% summary
lm(mt_t_value ~ degree_slope + degree_slope2, 
   data = df_MRC_relationships %>% mutate(degree_slope2 = degree_slope^2)) %>% summary

### DEGREE

# degree effect size and MT baseline (none)
df_MRC_relationships %>% 
  ggplot(aes(x = mt_baseline, y = degree_t_value)) +
  geom_point()
lm(degree_t_value ~ mt_baseline, data = df_MRC_relationships) %>% summary
lm(degree_t_value ~ mt_baseline + mt_baseline2, 
   data = df_MRC_relationships %>% mutate(mt_baseline2 = mt_baseline^2)) %>% summary # significant
lm(degree_t_value ~ mt_baseline + mt_baseline2 + mt_baseline3, 
   data = df_MRC_relationships %>% mutate(mt_baseline2 = mt_baseline^2, mt_baseline3 = mt_baseline^3)) %>% summary

# degree effect size and MT slope (none)
df_MRC_relationships %>% 
  ggplot(aes(x = mt_slope, y = degree_t_value)) +
  geom_point()
lm(degree_t_value ~ mt_slope, data = df_MRC_relationships) %>% summary
lm(degree_t_value ~ mt_slope + mt_slope2, 
   data = df_MRC_relationships %>% mutate(mt_slope2 = mt_slope^2)) %>% summary
lm(degree_t_value ~ mt_slope + mt_slope2 + mt_slope3, 
   data = df_MRC_relationships %>% mutate(mt_slope2 = mt_slope^2, mt_slope3 = mt_slope^3)) %>% summary

# degree effect size and degree baseline (none)
df_MRC_relationships %>% 
  ggplot(aes(x = degree_baseline, y = degree_t_value)) +
  geom_point()
lm(degree_t_value ~ degree_baseline, data = df_MRC_relationships) %>% summary
lm(degree_t_value ~ degree_baseline + degree_baseline2, 
   data = df_MRC_relationships %>% mutate(degree_baseline2 = degree_baseline^2)) %>% summary

# degree effect size and degree slope (quadratic)
df_MRC_relationships %>% 
  ggplot(aes(x = degree_slope, y = degree_t_value)) +
  geom_point()
lm(degree_t_value ~ degree_slope, data = df_MRC_relationships) %>% summary
lm(degree_t_value ~ degree_slope + degree_slope2, 
   data = df_MRC_relationships %>% mutate(degree_slope2 = degree_slope^2)) %>% summary

# USE LINEAR RELATIONSHIPS

f_model_label(lm(mt_t_value ~ mt_slope, 
                 data = df_MRC_relationships %>% mutate(mt_slope2 = mt_slope^2)), type = "linear")


# C: MT & MRC and MT slope ------------------------------------------------

mt_model <- lm(mt_t_value ~ mt_slope, data = df_MRC_relationships)
summary(mt_model)

df_MRC_relationships %>% 
  mutate(sign = ifelse(mt_t_value < 0, "negative", "positive"),
         formula = f_model_label(mt_model, type = "linear"),
         hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor(),
         label = ifelse(hubs == 1, str_replace_all(region_of_interest, "_", " "), NA)
  ) %>% 
  
  ggplot(aes(x = mt_t_value, y = mt_slope)) +
  geom_point(aes(fill = mt_t_value, color = sign, shape = hubs, size = hubs)) +
  # geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  # geom_text(
  #   x = -1.5, y = -0.003, color = "blue", size = 2.5, hjust = 0.5,
  #   aes(label = formula), check_overlap = TRUE
  # ) +
  geom_text_repel(aes(label = label), size = 3) +
  geom_vline(aes(xintercept = 0)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", guide = "none") +
  scale_color_manual(values = c("#2166AC", "#B2182B"), guide = "none") +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(1, 3)) +
  labs(y = "MT decay slope in developmental cohort \n(corrected for TBV)",
       x = "RMS-control MT t-value post-stress",
       title = "Relationship between MT group effect size post-stress \nand MT decay slope in MRC (at ROI level)") +
  coord_flip()

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/8c.case_control_adult_MRC_baseline_MT", .x), width = 5, height = 4)
)


# D: degree - degree slope relationship -----------------------------------

degree_model <- lm(degree_t_value ~ degree_slope + degree_slope2, 
                   data = df_MRC_relationships %>% mutate(degree_slope2 = degree_slope^2))

df_MRC_relationships %>% 
  mutate(sign = ifelse(degree_t_value < 0, "negative", "positive"),
         formula = f_model_label(degree_model, type = "linear"),
         hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor(),
         label = ifelse(hubs == 1, str_replace_all(region_of_interest, "_", " "), NA)
  ) %>% 
  
  ggplot(aes(x = degree_t_value, y = degree_slope)) +
  geom_point(aes(fill = degree_t_value, color = sign, shape = hubs, size = hubs)) +
  # geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black") +
  # geom_text(
  #   x = -1.5, y = -0.003, color = "blue", size = 2.5, hjust = 0.5,
  #   aes(label = formula), check_overlap = TRUE
  # ) +
  geom_text_repel(aes(label = label), size = 3) +
  geom_vline(aes(xintercept = 0)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", guide = "none") +
  scale_color_manual(values = c("#2166AC", "#B2182B"), guide = "none") +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(1, 3)) +
  labs(y = "degree slope in developmental cohort \n(corrected for TBV)",
       x = "RMS-control degree t-value post-stress",
       title = "Relationship between degree group effect size post-stress \nand degree slope in MRC (at ROI level)") +
  coord_flip()

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/8d.case_control_adult_MRC_baseline_degree", .x), width = 5, height = 4)
)

# E: case-control edge weights at baseline ---------------------------------------

# PREP FOR BALL AND STICK PLOT

# calculate node position
df_node_coords <- df_sigma_atlas %>% 
  filter(side == "sagittal", hemisphere == "right" & roi_abbreviation != "mcp" & region_of_interest %in% gray_matter_rois) %>% 
  group_by(side, system) %>%
  summarize(geometry = st_union(geometry)) %>% 
  dplyr::select(side, system) %>% 
  st_centroid %>% 
  mutate(point = map(geometry, ~st_coordinates(.x) %>% as.data.frame %>% as_tibble)) %>% 
  unnest(cols = c(point)) %>% 
  as_tibble() %>% 
  dplyr::select(-geometry) %>% 
  clean_names() %>% 
  
  # add degree and module to nodes
  left_join(df_degree_gsk %>% 
              filter(timepoint == 20) %>% 
              group_by(system) %>% 
              summarise(degree = median(feature_resids))
  ) %>% 
  mutate(system = str_replace_all(system, "_", " ")) %>% 
  filter(system != "sensory motor")

# edge weights
df_edge_weights <- df_edge_pnd300_SYS %>% 
  separate(sys_edge, into = c("S1", "S2"), sep = "-") %>% 
  dplyr::select(S1, S2, t_value, p_adj) %>% 
  
  left_join(df_node_coords %>% dplyr::rename("S1" = "system"),
            by = join_by(S1)) %>% 
  dplyr::rename("x1" = "x", "y1" = "y") %>% 
  left_join(df_node_coords %>% dplyr::rename("S2" = "system"), 
            by = join_by(S2)) %>% 
  dplyr::rename("x2" = "x", "y2" = "y") %>% 
  filter(side.x == side.y) %>% 
  dplyr::rename("side" = "side.x")

# PLOT
ggplot() +
  geom_segment(aes(x = -x1, y = y1,
                   xend = -x2, yend = y2, 
                   color = t_value, alpha = t_value, linewidth = abs(t_value)), 
               data = df_edge_weights %>% arrange(abs(t_value))) +
  geom_point(aes(x = -x, y = y, size = degree, fill = degree), 
             data = df_node_coords, shape = 21) +
  geom_text_repel(aes(x = -x, y = y, label = system),
                  data = df_node_coords, size = 3.5) +
  
  scale_fill_viridis() +
  scale_color_gradient2("low" = "#2166AC", "mid" = "white", "high" = "#B2182B") +
  
  scale_size_continuous(range = c(3, 6), guide = "none") +
  scale_linewidth_continuous(range = c(1, 3), guide = "none") +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  
  facet_wrap(vars(side)) +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "bottom",
        legend.justification = "left")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/8e.adult_edge_case_control", .x), width = 8, height = 4) 
)


# S7a: MT box plot -------------------------------------------------------------

f_plot_box(df = df_mt_pnd300_SYS,
           xlab = "", 
           ylab = "MT \n(corrected for sex + TBV)", 
           title = "Post-stress system-level MT differences", 
           legend = FALSE) +
  coord_flip()

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S7a.adult_MT_boxplot", .x), width = 4, height = 3)
)

# S7b: degree box plot ---------------------------------------------------------

f_plot_box(df = df_degree_pnd300_SYS,
           xlab = "", 
           ylab = "degree \n(corrected for sex + TBV)", 
           title = "Post-stress system-level degree differences", 
           legend = FALSE) +
  coord_flip()

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S7b.adult_degree_boxplot", .x), width = 4, height = 3)
)


# S7c: distinguish edges from PND 63 -------------------------------

df_edge_pnd063_SYS <- df_edge_gsk %>%
  filter(timepoint == 63) %>%
  group_by(sys_edge) %>%
  nest() %>%
  f_run_lm(63) %>% 
  filter(str_detect(term, "group"))

# edge weights
df_edge_weights <- df_edge_pnd063_SYS %>% 
  separate(sys_edge, into = c("S1", "S2"), sep = "-") %>% 
  dplyr::select(S1, S2, t_value, p_adj) %>% 
  
  left_join(df_node_coords %>% dplyr::rename("S1" = "system"),
            by = join_by(S1)) %>% 
  dplyr::rename("x1" = "x", "y1" = "y") %>% 
  left_join(df_node_coords %>% dplyr::rename("S2" = "system"), 
            by = join_by(S2)) %>% 
  dplyr::rename("x2" = "x", "y2" = "y") %>% 
  filter(side.x == side.y) %>% 
  dplyr::rename("side" = "side.x")

# PLOT
ggplot() +
  geom_segment(aes(x = -x1, y = y1,
                   xend = -x2, yend = y2, 
                   color = t_value, alpha = t_value, linewidth = abs(t_value)), 
               data = df_edge_weights %>% arrange(abs(t_value))) +
  geom_point(aes(x = -x, y = y, size = degree, fill = degree), 
             data = df_node_coords, shape = 21) +
  geom_text_repel(aes(x = -x, y = y, label = system),
                  data = df_node_coords, size = 3.5) +
  
  scale_fill_viridis() +
  scale_color_gradient2("low" = "#2166AC", "mid" = "white", "high" = "#B2182B") +
  
  scale_size_continuous(range = c(3, 6), guide = "none") +
  scale_linewidth_continuous(range = c(1, 3), guide = "none") +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  
  facet_wrap(vars(side)) +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "bottom",
        legend.justification = "left")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/S7c.PND063_edge_case_control", .x), width = 8, height = 4) 
)







