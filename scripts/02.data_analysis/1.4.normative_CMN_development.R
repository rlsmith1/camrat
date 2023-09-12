
###
### Characterize development of CMN and relate that to development of MT
###

# libraries ---------------------------------------------------------------

library(tidyverse)
library(igraph)

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

# load data ---------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

load(paste0(base_dir, "objects/25Aug2023_df_data.RDS")) # df_data
load(paste0(base_dir, "objects/25Aug2023_sigma_labels.RDS")) # df_sigma_labels
load(paste0(base_dir, "objects/25Aug2023_sigma_atlas_for_plotting.RDS")) # df_sigma_atlas
load(paste0(base_dir, "objects/29Aug2023_df_mind_GM.RDS")) # df_mind_GM

# info on if scan had gadolinium
df_gad <- read_xlsx(paste0(base_dir, "data/gadolinium_use_records.xlsx")) %>% 
  dplyr::select(subject, scan_date, gad_by_record)


# format data -------------------------------------------------------------

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



# calculate median degree for each module through development -------------

df_median_module_degree <- df_degree %>% 
  ungroup %>% 
  filter(study == "MRC") %>% 
  dplyr::select(timepoint, age, region_of_interest, subject, weighted_degree, feature_resids) %>% 
  distinct() %>% 
  
  # add modules
  left_join(df_mind_mods %>% dplyr::select(-k), by = join_by(region_of_interest)) %>% 
  
  # add sample size (n)
  left_join(
    df_degree %>% 
      ungroup %>% 
      filter(study == "MRC") %>% 
      dplyr::select(timepoint, subject) %>%
      distinct %>% 
      count(timepoint),
    by = join_by(timepoint)
  ) %>% 
  
  # calculate median degree and 95% CIs for each module
  group_by(timepoint, module) %>%
  reframe(
    age = median(age),
    median_degree = median(feature_resids),
    sd_degree = sd(feature_resids),
    ci_degree = qnorm(0.975)*sd_degree/sqrt(n)
  ) %>% 
  distinct()


# relate MT and degree across development ---------------------------------


# CALCULATE MT DECAY SLOPES AT ROI LEVEL
df_mt_mrc <- df_mt %>% filter(study == "MRC")
ROIs <- df_mt_mrc %>% pull(region_of_interest) %>% unique

df_slopes_mt <- tibble()
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
           mt_slope = lme_model@beta[2],
           mt_standard_error = `97.5 %` - `2.5 %` / 3.92
    )
  
  df_slopes_mt <- df_slopes_mt %>% 
    bind_rows(df_tmp)
  
}
df_slopes_mt <- df_slopes_mt %>% arrange(mt_slope)

# CALCULATE DEGREE SLOPES AT ROI LEVEL
df_degree_mrc <- df_degree %>% filter(study == "MRC")
ROIs <- df_degree_mrc %>% pull(region_of_interest) %>% unique

df_slopes_degree <- tibble()
for (r in ROIs) {
  
  print(r)
  
  lme_model <- lme4::lmer(feature_resids ~ age + (1 | subject), 
                          data = df_degree_mrc %>% filter(region_of_interest == r & timepoint %in% c(20, 63))
  )
  
  df_tmp <- confint(lme_model) %>% 
    as.data.frame() %>% 
    rownames_to_column("term") %>% 
    filter(term == "age") %>% 
    mutate(region_of_interest = r,
           degree_slope = lme_model@beta[2],
           degree_standard_error = `97.5 %` - `2.5 %` / 3.92
    )
  
  df_slopes_degree <- df_slopes_degree %>% 
    bind_rows(df_tmp)
  
}
df_slopes_degree <- df_slopes_degree %>% arrange(-abs(degree_slope))

# COMBINE RESULTS
df_mt_degree_development <- df_mt %>% 
  filter(study == "MRC" ) %>% 
  group_by(region_of_interest, timepoint) %>% 
  summarise(mt = median(feature_resids, na.rm = TRUE)) %>% 
  
  left_join(
    df_degree %>% 
      filter(study == "MRC") %>% 
      group_by(region_of_interest, timepoint) %>% 
      summarise(degree = median(feature_resids, na.rm = TRUE)),
    by = join_by(region_of_interest, timepoint)
  ) %>% 
  
  filter(region_of_interest %in% gray_matter_rois) %>% 
  pivot_longer(c("mt", "degree"), names_to = "feature", values_to = "value") %>% 
  unite("feature", c(feature, timepoint), sep = "_") %>% 
  pivot_wider(id_cols = region_of_interest, names_from = feature, values_from = value) %>% 
  
  left_join(df_slopes_mt %>% dplyr::select(region_of_interest, mt_slope, mt_standard_error), 
            by = join_by(region_of_interest)
  ) %>% 
  left_join(df_slopes_degree %>% dplyr::select(region_of_interest, degree_slope, degree_standard_error), 
            by = join_by(region_of_interest)
  )


# A: CMN through development -----------------------------------------------

p_fig_5a <- df_mind_GM %>% 
  filter(study == "MRC") %>% 
  group_by(timepoint, R1, R2) %>% 
  summarise(weight = median(weight)) %>%  
  
  mutate(R1 = factor(R1, levels = roi_order),
         R2 = factor(R2, levels = roi_order)) %>%
  
  # plot
  ggplot(aes(x = R1, y = R2, fill = weight)) +
  geom_tile() +
  facet_wrap(vars(timepoint), nrow = 1) +
  
  geom_vline(xintercept = module_lines, color = "white", linewidth = 0.5) +
  geom_hline(yintercept = module_lines, color = "white", linewidth = 0.5) +
  
  labs(x = "", y = "", title = "Median CMN at each timepoint") +
  scale_fill_viridis(na.value = "white", limits = c(0, 1)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm")
  )

color_grob <- (df_mind_mods %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = roi_order)) %>% 
  arrange(region_of_interest) %>% 
  left_join(
    enframe(module_colors) %>% 
      dplyr::rename("color" = "value") %>% 
      mutate(module = factor(name))
  ) %>% 
  ggplot() +
  geom_tile(aes(x = 1, y = region_of_interest, fill = I(color))) +
  theme_void()) %>% ggplotGrob()

p_fig_5a

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5a.CMN_development.pdf"), width = 10, height = 4)
ggsave(paste0(base_dir, "outputs/figures/5a.CMN_development.png"), width = 10, height = 4)


# B: modular degree through development -----------------------------------

p_fig_5b <- df_median_module_degree %>% 
  mutate(label = ifelse(timepoint == 300, module, NA)) %>% 
  
  ggplot(aes(x = age, y = median_degree, color = module)) +
  geom_point() +
  geom_line(aes(group = module)) +
  geom_errorbar(aes(ymin = median_degree - ci_degree, ymax = median_degree + ci_degree)) +
  geom_label_repel(aes(label = label, fill = module), color = "black") +
  scale_color_manual(values = module_colors) +
  scale_fill_manual(values = module_colors) +
  labs(y = "Median module degree (+/- 95% CI)",
       title = "Modular degree through development") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 12))

p_fig_5b

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5b.module_degree_development.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/5b.module_degree_development.png"), width = 5, height = 4)


# C: relationship between MT and degree and changes in both ------------------

# describe model parameters
adj_r_sq <- round((parabolic_model %>% summary %>% .$adj.r.squared), 2)
coefs <- coefficients(parabolic_model) %>% round(2)
model_formula <- paste0("y = ", coefs[1], " + ", coefs[2], "*x + ", coefs[3], "*x^2")

# find MT at max(degree)
model_function <- function(x) {
  coefs[1] + coefs[2]*x + coefs[3]*x^2
}
mt_at_max_degree <- optimize(f = model_function, interval = c(min(df_mt_degree_cor$median_mt), max(df_mt_degree_cor$median_mt)), maximum = TRUE)

# plot
p_fig4c <- df_mt_degree_cor %>% 
  ggplot(aes(x = median_mt, y = median_degree)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(formula = y ~ x + I(x^2), method = "lm") +
  annotate(
    geom = "text", x = 0.35, y = 30, color = "blue", size = 3, hjust = 0.5, align = "center",
    label = paste0(model_formula, "\nR^2 = ", adj_r_sq)
  ) +
  labs(x = "Median MT", y = "Median degree",
       title = "Relationship between MT and degree in control adult subjects",
       caption = paste0("MT at max(degree) = ", round(mt_at_max_degree$maximum, 2))
  )

df_mt_degree_development
df_slopes_degree
df_slopes_mt %>% arrange(abs(mt_slope))

#### MT AND DEGREE RELATIONSHIP THROUGH DEVELOPMENT
df_mt_degree_development

parabolic_model <- lm(median_degree ~ median_mt + median_mt2, data = df_mt_degree_cor) 

### MOST DEVELOPMENTALLY ACTIVE IN MT AND DEGREE
df_mt_degree_slopes <- df_mt_degree_development %>% 
  dplyr::select(region_of_interest, contains(c("slope", "standard_error"))) %>% 
  pivot_longer(contains(c("slope", "standard_error")), names_to = "feature", values_to = "value") %>% 
  separate(feature, into = c("phenotype", "feature"), sep = "_", extra = "merge") %>% 
  
  left_join(df_sys_to_roi) %>% 
  filter(region_of_interest %in% gray_matter_rois) %>% 
  
  group_by(system, phenotype, feature) %>% 
  summarise(value = median(value))

# different from 0?
n <- df_mt %>% filter(study == "MRC" & timepoint %in% c(20, 63)) %>% pull(subject) %>% unique %>% length
df_mt_degree_slopes_pvals <- df_mt_degree_slopes %>% 
  pivot_wider(id_cols = c(system, phenotype), names_from = feature, values_from = value) %>% 
  summarise(t = (slope - 0)/(standard_error/sqrt(n))) %>% 
  group_by(phenotype) %>% 
  mutate(p_val = pt(t, df = n -1, lower.tail = FALSE),
         p_adj = p.adjust(p_val, method = "fdr")
  )

# plot
df_mt_degree_slopes %>% 
  left_join(df_mt_degree_slopes_pvals) %>% 
  mutate(system = str_replace(system, "_", " "),
         system = str_replace(system, "_", "\n")
  ) %>% 
  pivot_wider(id_cols = c(system, phenotype, p_adj), names_from = feature, values_from = value) %>% 
  mutate(sign = ifelse(slope < 0, "negative", "positive"),
         label = ifelse(p_adj < 0.05, "*", "")
  ) %>% 
  
  ggplot(aes(x = reorder_within(system, -slope, phenotype), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(fill = slope, color = sign), size = 3, shape = 21) +
  geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
  geom_text(aes(label = label), size = 5, vjust = 0) +
  facet_wrap(vars(phenotype), scales = "free") +
  scale_x_reordered() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = "", y = "PND 20 --> 63 slope (+/- standard error)",
       title = "Rate of early developmental MT decay in gray matter brain systems") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )


# actual B: regional degree over time (maps -----------------------------------------

df_degree_mrc <- df_degree %>% filter(study == "MRC" & region_of_interest %in% gray_matter_rois)

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

df_slices %>% 
  filter(!str_detect(region_of_interest, paste0(bad_reg_rois, collapse = "|"))) %>% 
  left_join(df_degree_mrc %>% 
              group_by(region_of_interest, timepoint) %>% 
              summarise(feature_resids = median(feature_resids)), 
            by = join_by(region_of_interest)) %>% 
  arrange(region_of_interest != "corpus_callosum_and_associated_subcortical_white_matter") %>% 
  bind_rows(df_na_regs, .) %>% 
  filter(timepoint != 35) %>% 
  
  ggplot() +
  geom_sf(aes(fill = feature_resids, color = feature_resids, geometry = geometry, group = -1), 
          lwd = 0.5) +
  scale_fill_viridis(na.value = "grey") +
  scale_color_viridis(na.value = "grey") +
  facet_grid(side ~ timepoint) +
  labs(title = "Median degree map through development",
       fill = "Degree \n(corrected for TBV)", color = "Degree \n(corrected for TBV)") +
  theme_void() +
  theme(strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.25, "cm")
  )

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/5b.degree_maps_development", .x), width = 5, height = 4)
)


# actual C: slopes --------------------------------------------------------

load(paste0(base_dir, "objects/mt_degree_slopes.RDS")) # df_mt_degree_slopes

df_mt_degree_slopes_SYS <- df_mt_degree_slopes %>% 
  filter(phenotype == "degree" & period == "early" & region_of_interest %in% gray_matter_rois) %>% 
  left_join(df_sys_to_roi) %>% 
  group_by(system) %>% 
  summarise(slope = median(slope),
            standard_error = median(standard_error)
  ) %>% 
  mutate(system = str_replace(system, "_", " "),
         system = str_replace(system, "_", "\n"),
         sign = ifelse(slope < 0, "neg", "pos")
  )

df_mt_degree_slopes_SYS %>% 
  ggplot(aes(x = reorder(system, -slope), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(fill = slope, color = sign), shape = 21, size = 3) +
  geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  scale_fill_gradient2(low = "#2166AC", mid ="white", high = "#B2182B") +
  #ylim(c(-0.015, 0)) +
  labs(x = "", y = "PND 20 --> 63 slope (+/- standard error)",
       title = "Rate of early developmental degree change in \ngray matter brain systems") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/5c.degree_early_development_slopes", .x), width = 5, height = 4)
)

# D: Relationship between early developmental degree slope and anatomical position ----------------------------

# PLOT SLOPES ON BRAINS
slices <- list("sagittal" = 125, 
               "axial" = 115,
               "coronal" = 150
)
df_slices <- map_dfr(.x = 1:length(slices),
                     .f = ~ df_sigma_atlas %>%
                       filter(region_of_interest != "striatum") %>% 
                       filter(side == names(slices[.x]) & slice == slices[.x])
                     
) %>% 
  dplyr::select(region_of_interest, side, geometry) %>% 
  distinct()

p_fig5d_brain <- df_slices %>% 
  left_join( df_mt_degree_slopes %>% 
               filter(phenotype == "degree" & period == "early" & region_of_interest %in% gray_matter_rois), 
            by = "region_of_interest") %>% 
  mutate(slope = ifelse(region_of_interest %in% white_matter_rois, NA, slope)) %>% 
  arrange(!is.na(slope)) %>% 
  
  ggplot() +
  geom_sf(aes(color = slope, fill = slope, geometry = geometry, group = -1)) +
  scale_fill_gradient2(low = "#2166AC", mid ="white", high = "#B2182B") +
  scale_color_gradient2(low = "#2166AC", mid ="white", high = "#B2182B") +
  facet_wrap(vars(side), nrow = 1) +
  #labs(title = "Rate of early developmental degree change across regions") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.position = "bottom"
  )

# SLOPE VS X, Y, Z COORDS
labels <- c("x (left --> right)", "y (posterior --> anterior)", "z (inferior --> superior)")
names(labels) <- c("x", "y", "z")

p_fig5d_coords <- df_centroids %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  left_join(df_mt_degree_slopes %>% 
              filter(phenotype == "degree" & period == "early" & region_of_interest %in% gray_matter_rois), by = join_by(region_of_interest)
            ) %>%
  mutate(sign = ifelse(slope < 0, "neg", "pos")) %>% 
  
  ggplot(aes(x = coord, y = slope)) +
  geom_point(aes(fill = slope, color = sign), shape = 21, size = 2) +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  scale_fill_gradient2(low = "#2166AC", mid ="white", high = "#B2182B") +
  geom_smooth(method = "lm", color = "black", se = FALSE, lty = 2) +
  stat_cor(color = "black", size = 3, label.y = -0.13) +
  #geom_text_repel(aes(label = system)) +
  facet_wrap(vars(dim), scales = "free_x", labeller = as_labeller(labels)) +
  labs(x = "Coordinate position", y = "MT decay slope",
       title = "ROI-level degree slope relative to anatomical position (RH)") +
  theme(legend.position = "none")

p_fig5d_brain / p_fig5d_coords +
  plot_annotation(title = "Anatomical patterning of degree change in gray matter regions")

# SAVE
map(
  .x = c(".pdf", ".png"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/5d.anat_degree_slope", .x), width = 6, height = 4)
)





