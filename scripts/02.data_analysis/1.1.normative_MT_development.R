
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(viridis)
library(sf)
library(ggh4x)
library(ggpubr)
library(patchwork)
library(janitor)
library(tidytext)

# set plot theme ----------------------------------------------------------

theme_set(theme_light() +
            theme(plot.title = element_text(size = 14),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  strip.text = element_text(size = 12, color = "black"),
                  strip.background = element_rect(fill = "white", color = "gray"),
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 12)
                  #legend.key.width = unit(2, "cm")
            )
)

group_cols <- c("#FFA500", "#4682B4")
names(group_cols) <- c("control", "MS")


# helper functions --------------------------------------------------------

### T-TEST USING MEAN AND STANDARD ERROR

# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the sample sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 

t.test2 <- function(m1, m2, s1, s2, n1, n2, m0 = 0, equal.variance = FALSE) {
  
  if ( equal.variance == FALSE ) {
    
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    
  } else {
    
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1 - 1)*s1^2 + (n2 - 1)*s2^2)/(n1 + n2 - 2) ) 
    df <- n1 + n2 - 2
  }     
  
  t <- (m1 - m2 - m0)/se 
  dat <- c(m1 - m2, se, t, 2*pt(-abs(t), df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
  
}

### BRAIN MAP PLOT

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
                       
  )
  
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

# load data ---------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

load(paste0(base_dir, "objects/25Aug2023_df_data.RDS")) # df_data
load(paste0(base_dir, "objects/25Aug2023_sigma_labels.RDS")) # df_sigma_labels
load(paste0(base_dir, "objects/25Aug2023_sigma_atlas_for_plotting.RDS")) # df_sigma_atlas

# format data -------------------------------------------------------------

# LIST ROIS THAT DID NOT REGISTER WELL TO REMOVE FROM ANALYSIS
bad_reg_rois <- c("spinal_cord", "brainstem", "cerebell", "commissural_stria_terminalis", "central_canal")

# CONVERT ROI TO SYSTEM LEVEL
df_sys_to_roi <- df_sigma_labels %>% 
  dplyr::select(matter, system, region_of_interest) %>% 
  distinct()

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

# FILTER FOR JUST DEVELOPMENTAL STUDY (MRC)
df_mt_mrc <- df_mt %>% filter(study == "MRC")

# analysis ----------------------------------------------------------------

# CALCULATE MT DECAY SLOPES AT ROI LEVEL

# generate models and calculate slopes and confidence intervals for each ROI
ROIs <- df_mt_mrc %>% pull(region_of_interest) %>% unique

df_confint_roi <- tibble()
for (r in ROIs) {
  
  print(r)
  
  lme_model <- lme4::lmer(feature_resids ~ age + (1 | subject), 
                          data = df_mt_mrc %>% filter(region_of_interest == r & timepoint %in% c(20, 63))
  )
  
  df_tmp <- confint(lme_model) %>% 
    as.data.frame() %>% 
    rownames_to_column("term") %>% 
    filter(term == "age") %>% 
    mutate(region_of_interest = r,
           slope = lme_model@beta[2],
           standard_error = `97.5 %` - `2.5 %` / 3.92
    )
  
  df_confint_roi <- df_confint_roi %>% 
    bind_rows(df_tmp)
  
}
df_confint_roi <- df_confint_roi %>% arrange(slope)

# run t-test on the pairwise comparison of each ROI slope with each other slope
# df_confint_roi_pvals <- expand_grid(
#   R1 = ROIs,
#   R2 = ROIs
# ) %>% 
#   filter(R1 != R2) %>% 
#   
#   mutate(
#     p_val = map2(
#       .x = R1,
#       .y = R2,
#       .f = ~  t.test2(
#         
#         m1 = df_confint_roi %>% filter(region_of_interest == .x) %>% pull(slope), 
#         m2 = df_confint_roi %>% filter(region_of_interest == .y) %>% pull(slope),
#         
#         s1 = df_confint_roi %>% filter(region_of_interest == .x) %>% pull(standard_error),
#         s2 = df_confint_roi %>% filter(region_of_interest == .y) %>% pull(standard_error),
#         
#         n1 = df_mt_mrc %>% ungroup %>% filter(region_of_interest == .x) %>% pull(subject) %>% unique %>% length,
#         n2 = df_mt_mrc %>% ungroup %>% filter(region_of_interest == .y) %>% pull(subject) %>% unique %>% length
#         
#       ) %>% 
#         .["p-value"]
#     )
#   ) %>% 
#   unnest(cols = c(p_val)) %>% 
#   mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
#   arrange(p_val)

# Run t-test comparisons at system level
df_confint_sys <- df_confint_roi %>% 
  left_join(df_sys_to_roi) %>% 
  group_by(system) %>% 
  summarise(slope = median(slope),
            standard_error = median(standard_error)
  )

sys_order <- df_confint_sys %>% arrange(slope) %>% pull(system)

df_confint_sys_pvals <- expand_grid(
  S1 = df_confint_sys %>% pull(system) %>% unique,
  S2 = df_confint_sys %>% pull(system) %>% unique
) %>%

  # order systems by slope
  arrange(match(S1, sys_order), match(S2, sys_order)) %>% 
  
  # remove duplicates in reverse order
  group_by(grp = paste0(pmin(S1, S2), sep = "_", pmax(S1, S2))) %>% 
  slice(1) %>% 
  ungroup %>% dplyr::select(-grp) %>% 
  
  # run pairwise t-test
  mutate(
    p_val = map2(
      .x = S1,
      .y = S2,
      .f = ~  t.test2(
        
        m1 = df_confint_sys %>% filter(system == .x) %>% pull(slope),
        m2 = df_confint_sys %>% filter(system == .y) %>% pull(slope),
        
        s1 = df_confint_sys %>% filter(system == .x) %>% pull(standard_error),
        s2 = df_confint_sys %>% filter(system == .y) %>% pull(standard_error),
        
        n1 = df_mt_mrc %>% ungroup %>% filter(system == .x) %>% pull(subject) %>% unique %>% length,
        n2 = df_mt_mrc %>% ungroup %>% filter(system == .y) %>% pull(subject) %>% unique %>% length
        
      ) %>%
        .["p-value"]
    )
  ) %>%
  unnest(cols = c(p_val)) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH")) %>%
  arrange(p_val)

# CALCULATE CENTROID OF EACH SYSTEM FOR PLOTTING
df_centroids <- df_sigma_atlas %>% 
  filter(side == "axial" & hemisphere == "right") %>% # axial for x & y coords
  group_by(system) %>%
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
      group_by(system) %>%
      summarize(geometry = st_union(geometry)) %>% 
      st_centroid %>% 
      mutate(point = map(geometry, ~st_coordinates(.x) %>% as.data.frame %>% as_tibble)) %>% 
      unnest(cols = c(point)) %>% 
      as_tibble() %>% 
      dplyr::select(-geometry) %>% 
      clean_names() %>% 
      dplyr::rename("z" = "y"),
    by = join_by(system)
  ) %>% 
  
  left_join(
    df_sigma_atlas %>% 
      filter(side == "sagittal" & hemisphere == "right") %>% # sagittal for y & z coords
      group_by(system) %>%
      summarize(geometry = st_union(geometry)) %>% 
      st_centroid %>% 
      mutate(point = map(geometry, ~st_coordinates(.x) %>% as.data.frame %>% as_tibble)) %>% 
      unnest(cols = c(point)) %>% 
      as_tibble() %>% 
      dplyr::select(-geometry) %>% 
      clean_names() %>% 
      dplyr::rename("y" = "x", "z" = "y"),
    by = join_by(system)
  ) %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  mutate(dim = substr(dim, start = 1, stop = 1)) %>% 
  group_by(system, dim) %>%
  summarise(coord = mean(coord)) %>% 
  pivot_wider(id_cols = system, names_from = dim, values_from = coord) %>% 
  ungroup


# A: Brain map of normative MT development ------------------------------------

slices <- list("sagittal" = 125, "coronal" = 150, "axial" = 115)
df_slices <- map_dfr(.x = 1:length(slices),
                     .f = ~ df_sigma_atlas %>%
                       filter(region_of_interest != "striatum") %>% 
                       filter(side == names(slices[.x]) & slice == slices[.x])
                     
) %>% 
  dplyr::select(region_of_interest, side, geometry) %>% 
  distinct()

df_na_regs <- df_slices %>% 
  filter(str_detect(region_of_interest, paste0(bad_reg_rois, collapse = "|"))) %>% 
  expand_grid(timepoint = unique(df_mt$timepoint))

p_fig1a <- df_slices %>% 
  filter(!str_detect(region_of_interest, paste0(bad_reg_rois, collapse = "|"))) %>% 
  left_join(df_mt_mrc %>% 
              group_by(region_of_interest, roi_abbreviation, timepoint) %>% 
              summarise(feature_resids = median(feature_resids)), 
            by = join_by(region_of_interest)) %>% 
  arrange(roi_abbreviation != "cc") %>% 
  bind_rows(df_na_regs, .) %>% 

  ggplot() +
  geom_sf(aes(fill = feature_resids, color = feature_resids, geometry = geometry, group = -1), 
          lwd = 0.5) +
  scale_fill_viridis(na.value = "grey") +
  scale_color_viridis(na.value = "grey") +
  facet_grid(side ~ timepoint) +
  labs(title = "Median MT map through development",
       fill = "MT \n(corrected for TBV)", color = "MT \n(corrected for TBV)") +
  theme_void() +
  theme(strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.25, "cm")
  )

# B: Overall MT decay ------------------------------------------------------------

p_fig1b <- df_mt_mrc %>% 
  group_by(timepoint) %>% 
  mutate(median_feaure_resids = median(feature_resids),
         median_age = median(age)) %>% 
  filter(timepoint != 35) %>% 
  
  ggplot() +
  geom_point(aes(x = age, y = feature_resids), size = 2, shape = 1) +
  geom_point(aes(x = median_age, y = median_feaure_resids), size = 3, color = "red") +
  geom_line(aes(x = median_age, y = median_feaure_resids), linewidth = 2, color = "red", lty = 2) +
  labs(x = "age (PND)", y = "Regional MT (corrected for TBV)",
       title = "MT across ROIs and subjects \n(red shows overall median)") +
  theme(legend.position = "none")

# C: ROI-specific MT decay in early development ------------------------------------------------------------

# PLOT SLOPES AND ERRORS
p_fig1c_slopes <- df_confint_roi %>% 
  left_join(df_sys_to_roi) %>% 
  group_by(system) %>% 
  summarise(slope = median(slope),
            standard_error = median(standard_error)
  ) %>% 
    mutate(system = str_replace_all(system, "_", "\n")) %>% 

  ggplot(aes(x = reorder(system, -slope), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(alpha = abs(slope)), color = "forestgreen", size = 3) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  ylim(c(-0.015, 0)) +
  labs(x = "", y = "PND 20 --> 63 slope (+/- standard error)",
       title = "Rate of early developmental MT decay in each brain system") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )

# MATRIX OF SLOPE DIFFERENCES
p_fig1c_pvals <- df_confint_sys_pvals %>%   
  mutate(S1 = factor(str_replace_all(S1, "_", "\n"), levels = rev(str_replace_all(sys_order, "_", "\n"))),
         S2 = factor(str_replace_all(S2, "_", "\n"), levels = rev(str_replace_all(sys_order, "_", "\n"))),
         significant = p_adj < 0.05
  ) %>% 
  
  ggplot(aes(x = S2, y = S1)) +
  geom_tile(aes(fill = -log10(p_adj), color = significant), linewidth = 1) +
  scale_color_manual(values = c("white", "black"), guide = "none") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "", y = "") +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top",
        legend.justification = "right"
  )


# D: Relationship between early developmental MT decay and anatomical position ----------------------------

# PLOT SLOPES ON BRAINS
slices <- list("sagittal" = 125, "coronal" = 150)
df_slices <- map_dfr(.x = 1:length(slices),
                     .f = ~ df_sigma_atlas %>%
                       filter(region_of_interest != "striatum") %>% 
                       filter(side == names(slices[.x]) & slice == slices[.x])
                     
) %>% 
  dplyr::select(region_of_interest, side, geometry) %>% 
  distinct()

p_fig1d_brain <- df_slices %>% 
  left_join(df_confint_roi, 
            by = "region_of_interest") %>% 
  filter(region_of_interest != "central_canal") %>% 

  ggplot() +
  geom_sf(aes(color = slope, fill = slope, geometry = geometry, group = -1)) +
  scale_fill_gradient(low = "forestgreen", high = "white", na.value = "grey") +
  scale_color_gradient(low = "forestgreen", high = "white", na.value = "grey") +
  facet_wrap(vars(side), nrow = 1) +
  #labs(title = "Rate of early developmental MT decay across regions") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.position = "bottom"
  )

# SLOPE VS X, Y, Z COORDS
labels <- c("x (left --> right)", "y (posterior --> anterior)", "z (inferior --> superior)")
names(labels) <- c("x", "y", "z")

p_fig1d_xyz <- df_centroids %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  left_join(df_confint_sys, by = join_by(system)) %>%
  
  ggplot(aes(x = coord, y = slope)) +
  geom_point(aes(alpha = abs(slope)), size = 3, color = "forestgreen") +
  geom_smooth(method = "lm", se = FALSE, lty = 2) +
  stat_cor(color = "blue") +
  #geom_text_repel(aes(label = system)) +
  facet_wrap(vars(dim), scales = "free_x", labeller = as_labeller(labels)) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  ylim(c(-0.015, 0)) +
  labs(x = "Coordinate position", y = "MT decay slope",
       title = "System-level MT decay slope relative to anatomical position (RH)") +
  theme(legend.position = "none")


# layout + patch ------------------------------------------------------------

doParallel::registerDoParallel()

layout <- c(
  
  # A
  patchwork::area(t = 1, b = 100, l = 1, r = 120),
  
  # B
  patchwork::area(t = 1, b = 100, l = 121, r = 180),
  
  # C
  patchwork::area(t = 101, b = 250, l = 2, r = 75),
  patchwork::area(t = 101, b = 250, l = 76, r = 180),
  
  # D
  patchwork::area(t = 251, b = 350, l = 2, r = 60),
  patchwork::area(t = 251, b = 350, l = 61, r = 180)
  
)

p_fig1a +
  p_fig1b + 
  p_fig1c_slopes + (p_fig1c_pvals + plot_layout(tag_level = "new")) +
  p_fig1d_brain + (p_fig1d_xyz + plot_layout(tag_level = "new")) +
  
  plot_layout(design = layout) + 
  plot_annotation(title = "Regional MT development",
                  theme = theme(plot.title = element_text(size = 20, face = "bold")),
                  tag_levels =  "A")  &
  theme(plot.tag = element_text(face = 'bold', size = 16))

# Patchwork taking too long to render... save each subpanel row separately
p_fig1a + p_fig1b +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(face = 'bold', size = 16))
ggsave(paste0(base_dir, "outputs/figures/1.1a_b.pdf"), width = 10, height = 4)
ggsave(paste0(base_dir, "outputs/figures/1.1a_b.png"), width = 10, height = 4)

p_fig1c_slopes + (p_fig1c_pvals + plot_layout(tag_level = "new")) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = list(c("C"))) &
  theme(plot.tag = element_text(face = 'bold', size = 16))
ggsave(paste0(base_dir, "outputs/figures/1.1c.pdf"), width = 14, height = 12)
ggsave(paste0(base_dir, "outputs/figures/1.1c.png"), width = 14, height = 12)

p_fig1d_brain / (p_fig1d_xyz + plot_layout(tag_level = "new"))
ggsave(paste0(base_dir, "outputs/figures/1.1c_inset.pdf"), width = 8, height = 6)
ggsave(paste0(base_dir, "outputs/figures/1.1c_inset.png"), width = 8, height = 6)


# supplement analysis -----------------------------------------------------

# CALCULATE MT DECAY SLOPES AT ROI LEVEL

# generate models and calculate slopes and confidence intervals for each ROI
ROIs <- df_mt_mrc %>% pull(region_of_interest) %>% unique

df_confint_roi_sup <- tibble()
for (r in ROIs) {
  
  print(r)
  
  lme_model <- lme4::lmer(mt ~ age + (1 | subject), 
                          data = df_mt_mrc %>% filter(region_of_interest == r & timepoint %in% c(20, 63))
  )
  
  df_tmp <- confint(lme_model) %>% 
    as.data.frame() %>% 
    rownames_to_column("term") %>% 
    filter(term == "age") %>% 
    mutate(region_of_interest = r,
           slope = lme_model@beta[2],
           standard_error = `97.5 %` - `2.5 %` / 3.92
    )
  
  df_confint_roi_sup <- df_confint_roi_sup %>% 
    bind_rows(df_tmp)
  
}
df_confint_roi_sup <- df_confint_roi_sup %>% arrange(slope)

# supplement --------------------------------------------------------------

# S1 MT decay without TBV correction

# A
p_figS1a <- df_mt_mrc %>% 
  group_by(timepoint) %>% 
  mutate(median_mt = median(mt),
         median_age = median(age)) %>% 
  filter(timepoint != 35) %>% 
  
  ggplot() +
  geom_point(aes(x = age, y = mt), size = 2, shape = 1) +
  geom_point(aes(x = median_age, y = median_mt), size = 3, color = "red") +
  geom_line(aes(x = median_age, y = median_mt), linewidth = 2, color = "red", lty = 2) +
  labs(x = "age (PND)", y = "Regional MT (*not* corrected for TBV)",
       title = "MT across ROIs and subjects \n(red shows overall median)") +
  theme(legend.position = "none")

p_figS1a
ggsave(paste0(base_dir, "outputs/figures/S1a.pdf"), width = 8, height = 6)
ggsave(paste0(base_dir, "outputs/figures/S1a.png"), width = 8, height = 6)

# B slopes
p_figS1b <- df_confint_roi %>% 
  left_join(df_sys_to_roi %>% filter(!str_detect(region_of_interest, paste0(bad_reg_rois, collapse = "|")))) %>% 
  group_by(system) %>% 
  summarise(slope = median(slope),
            standard_error = median(standard_error),
            system = str_replace_all(system, "_", " ")
  ) %>% 
  
  ggplot(aes(x = reorder(str_wrap(system, 35), -slope), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(alpha = abs(slope)), color = "forestgreen", size = 3) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  ylim(c(-0.015, 0)) +
  labs(x = "", y = "PND 20 --> 63 slope (+/- standard error)",
       title = "Rate of early developmental MT decay in each brain system \n(no TBV correction)") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )

p_figS1b
ggsave(paste0(base_dir, "outputs/figures/S1b.pdf"), width = 6, height = 8)
ggsave(paste0(base_dir, "outputs/figures/S1b.png"), width = 6, height = 8)

# S2 MT decay (with TBV correction) split by white vs gray matter
n_sys <- df_sys_to_roi %>% dplyr::select(matter, system) %>% distinct %>% count(matter) %>% filter(matter != "csf") %>% pull(n)

df_confint_roi %>% 
  left_join(df_sys_to_roi %>% filter(!str_detect(region_of_interest, paste0(bad_reg_rois, collapse = "|")))) %>% 
  group_by(matter, system) %>% 
  summarise(slope = median(slope),
            standard_error = median(standard_error)
  ) %>% 
  
  ggplot(aes(x = reorder_within(str_wrap(system, 35), -slope, matter), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(alpha = abs(slope)), color = "forestgreen", size = 3) +
  facet_wrap(~ matter, scales = "free_y", ncol = 1) +
  force_panelsizes(rows = c(n_sys[1], n_sys[2])) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  scale_x_reordered() +
  ylim(c(-0.015, 0)) +
  labs(x = "", y = "PND 20 --> 63 slope (+/- standard error)",
       title = "Rate of early developmental MT decay in each brain system") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )
ggsave(paste0(base_dir, "outputs/figures/S2.pdf"), width = 6, height = 8)
ggsave(paste0(base_dir, "outputs/figures/S2.png"), width = 6, height = 8)
