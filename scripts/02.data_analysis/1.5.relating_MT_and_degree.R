
###
###
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

# COMBINE MT AND DEGREE IN DEVELOPMENTAL (MRC) COHORT
df_mt_degree <- df_mt %>% 
  filter(study == "MRC") %>% 
  dplyr::rename("feature" = "mt") %>% 
  mutate(phenotype = "mt") %>% 
  bind_rows(
    df_degree %>% 
      filter(study == "MRC") %>% 
      dplyr::rename("feature" = "weighted_degree") %>% 
      dplyr::select(-c(data, mind_graph, hub_score)) %>% 
      mutate(phenotype = "degree")
  )


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


# relationship between MT and degree at each timepoint --------------------

df_mt_degree_median <- df_mt_degree %>% 
  group_by(phenotype, timepoint, region_of_interest) %>% 
  summarise(feature = median(feature)) %>% 
  pivot_wider(id_cols = c(timepoint, region_of_interest), names_from = phenotype, values_from = feature) %>% 
  filter(timepoint != 35) %>% 
  mutate(
    degree2 = degree^2,
    mt2 = mt^2
  )

df_mt_degree_models <- df_mt_degree_median %>% 
  group_by(timepoint) %>% 
  nest() %>% 
  mutate(
    
    model_res = map(
      .x = data,
      .f = ~ list(
        
        # run linear model
        linear_model <- lm(degree ~ mt, data = .x), 
        lm_adj_r_sq <- summary(linear_model)$adj.r.squared,
        lm_f <- summary(linear_model)$fstatistic,
        lm_p <- pf(lm_f[1], lm_f[2], lm_f[3], lower.tail = FALSE),
        lm_coefs <- coefficients(linear_model) %>% round(2),
        lm_formula <- paste0("y = ", lm_coefs[1], " + ", lm_coefs[2], "*x"),

        # run parabolic model
        parabolic_model <- lm(degree ~ mt + mt2, data = .x), 
        pm_adj_r_sq <- summary(parabolic_model)$adj.r.squared,
        pm_f <- summary(parabolic_model)$fstatistic,
        pm_p <- pf(pm_f[1], pm_f[2], pm_f[3], lower.tail = FALSE),
        pm_coefs <- coefficients(parabolic_model) %>% round(2),
        pm_formula <- paste0("y = ", pm_coefs[1], " + ", pm_coefs[2], "*x + ", pm_coefs[3], "*x^2"),
        
        pm_function <- function(x) {pm_coefs[1] + pm_coefs[2]*x + pm_coefs[3]*x^2},
        pm_mt_at_max_degree <- optimize(f = pm_function, interval = c(0, 1), maximum = TRUE)$maximum,
        
        # return fit
        return(
          tibble(
            model = "linear", adj_r_sq = lm_adj_r_sq, p_val = lm_p, formula = lm_formula, mt_at_max_degree = NA
          ) %>% 
            bind_rows(
              tibble(
                model = "parabolic", adj_r_sq = pm_adj_r_sq, p_val = pm_p, formula = pm_formula, mt_at_max_degree = pm_mt_at_max_degree
              )
            )
        )
        
      )
    )
    
  ) %>% unnest(cols = c(model_res))

df_mt_degree_models # parabolic model fits best at each timepoint

# plot
df_mt_degree_models %>% 
  filter(model == "parabolic") %>% 
  mutate(label = paste0(formula, "\nR^2 = ", round(adj_r_sq, 2), "; p = ", format(p_val, scientific = TRUE, digits = 3))) %>% 
  unnest(cols = c(data)) %>%
  mutate(hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor) %>% 

  ggplot(aes(x = mt, y = degree)) +
  geom_point(aes(color = hubs), shape = 1, size = 2) +
  geom_smooth(formula = y ~ x + I(x^2), method = "lm") +
  
  # model formula and fit
  geom_text(
    x = 0.65, y = 20, color = "blue", size = 3, hjust = 0.5, align = "center",
    aes(label = label), check_overlap = TRUE
  ) +
  
  # MT at max(degree)
  geom_text(
    y = 80, color = "red", hjust = 0.5, check_overlap = TRUE, size = 3,
    aes(x = mt_at_max_degree, label= paste0("MT at max(degree) = ", round(mt_at_max_degree, 2))), 
    
  ) +
  
  # aesthetics
  facet_wrap(vars(timepoint)) +
  ylim(c(15, 80)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Median MT", y = "Median degree",
       title = "Relationship between MT and degree through development"
  )

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5c.mt_degree_relationship_development.pdf"), width = 8, height = 3)
ggsave(paste0(base_dir, "outputs/figures/5c.mt_degree_relationship_development.png"), width = 8, height = 3)

# calculate slope for each combination of timepoints ----------------------

df_mt_degree_slopes <- df_mt_degree %>% 
  expand_grid(period = c("early", "adult", "lifespan")) %>% 
  group_by(region_of_interest, phenotype, period) %>% 
  nest() %>% 
  
  mutate(
    slope = map2(
      .x = data,
      .y = period,
      .f = ~list(
        
        if(.y == "early") {
          pnds <- c(20, 63)
        } else if (.y == "adult") {
          pnds <- c(63, 300)
        } else if (.y == "lifespan") {
          pnds <- c(20, 300)
        },
        
        lme_model <- lme4::lmer(feature_resids ~ age + (1 | subject), 
                                data = .x %>% filter(timepoint %in% pnds)
        ),
        
        slope <- confint(lme_model) %>% 
          as.data.frame() %>% 
          rownames_to_column("term") %>% 
          as_tibble() %>% 
          filter(term == "age") %>% 
          mutate(slope = lme_model@beta[2],
                 standard_error = `97.5 %` - `2.5 %` / 3.92
          ),
        
        return(slope)
        
      )
    )
  ) %>% 
 unnest(cols = c(slope))

save(df_mt_degree_slopes, file = paste0(base_dir, "objects/mt_degree_slopes.RDS"))

# determine which slopes are statistically different from 0
n <- df_mt %>% filter(study == "MRC" & timepoint %in% c(20, 63)) %>% pull(subject) %>% unique %>% length

df_mt_degree_slopes_sys <- df_mt_degree_slopes %>% 
  filter(region_of_interest %in% gray_matter_rois) %>% 
  left_join(df_sys_to_roi) %>% 
  group_by(system, phenotype, period) %>% 
  summarise(
    slope = median(slope),
    standard_error = median(standard_error)
  )
  
df_mt_degree_slopes_pvals <- df_mt_degree_slopes_sys %>% 
  mutate(t = (slope - 0)/(standard_error/sqrt(n))) %>% 
  group_by(phenotype, period) %>% 
  mutate(p_val = pt(t, df = n -1, lower.tail = FALSE),
         p_adj = p.adjust(p_val, method = "fdr")
  )

# plot --------------------------------------------------------------------

labels <- c("Early (PND 20 --> 63)", "Adult (PND 63 -> 300)", "Lifespan (PND 20 --> 300)")
names(labels) <- c("early", "adult", "lifespan")

# plot
p_mt_slopes <- df_mt_degree_slopes_sys %>% 
  filter(phenotype == "mt") %>% 
  mutate(period = factor(period, levels = c("early", "adult", "lifespan"))) %>% 
  mutate(system = str_replace(system, "_", " "),
         system = str_replace(system, "_", "\n")
  ) %>% 
  mutate(sign = ifelse(slope < 0, "negative", "positive")) %>% 
  
  ggplot(aes(x = reorder_within(system, -slope, period), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(fill = slope, color = sign), size = 3, shape = 21) +
  geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
  facet_wrap(vars(period), scales = "free_y", labeller = as_labeller(labels)) +
  scale_x_reordered() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = "", y = "slope (+/- standard error)",
       title = "MT") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )

p_degree_slopes <- df_mt_degree_slopes_sys %>% 
  filter(phenotype == "degree") %>% 
  mutate(period = factor(period, levels = c("early", "adult", "lifespan"))) %>% 
  mutate(system = str_replace(system, "_", " "),
         system = str_replace(system, "_", "\n")
  ) %>% 
  mutate(sign = ifelse(slope < 0, "negative", "positive")) %>% 
  
  ggplot(aes(x = reorder_within(system, slope, period), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(fill = slope, color = sign), size = 3, shape = 21) +
  geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
  facet_wrap(vars(period), scales = "free_y", labeller = as_labeller(labels)) +
  scale_x_reordered() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = "", y = "slope (+/- standard error)",
       title = "Degree") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )

p_mt_slopes / p_degree_slopes + plot_annotation(title = "Rate of MT and degree change in gray matter brain systems across the lifespan")

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S1d.mt_degree_slopes_all.pdf"), width = 10, height = 8)
ggsave(paste0(base_dir, "outputs/figures/S1d.mt_degree_slopes_all.png"), width = 10, height = 8)


# relationship between MT and degree slopes -------------------------------

df_mt_degree_slopes_models <- df_mt_degree_slopes %>% 
  pivot_wider(id_cols = c(region_of_interest, period), names_from = phenotype, values_from = slope) %>% 
  mutate(
    degree2 = degree^2,
    mt2 = mt^2
  ) %>% 
  group_by(period) %>% 
  nest() %>% 
  mutate(
    
    model_res = map(
      .x = data,
      .f = ~ list(
        
        # run linear model
        linear_model <- lm(degree ~ mt, data = .x), 
        lm_adj_r_sq <- summary(linear_model)$adj.r.squared,
        lm_f <- summary(linear_model)$fstatistic,
        lm_p <- pf(lm_f[1], lm_f[2], lm_f[3], lower.tail = FALSE),
        lm_coefs <- coefficients(linear_model) %>% round(2),
        lm_formula <- paste0("y = ", lm_coefs[1], " + ", lm_coefs[2], "*x"),
        
        # run parabolic model
        parabolic_model <- lm(degree ~ mt + mt2, data = .x), 
        pm_adj_r_sq <- summary(parabolic_model)$adj.r.squared,
        pm_f <- summary(parabolic_model)$fstatistic,
        pm_p <- pf(pm_f[1], pm_f[2], pm_f[3], lower.tail = FALSE),
        pm_coefs <- coefficients(parabolic_model) %>% round(2),
        pm_formula <- paste0("y = ", pm_coefs[1], " + ", pm_coefs[2], "*x + ", pm_coefs[3], "*x^2"),
        
        pm_function <- function(x) {pm_coefs[1] + pm_coefs[2]*x + pm_coefs[3]*x^2},
        pm_mt_at_min_degree <- optimize(f = pm_function, interval = c(-0.01, 0), maximum = FALSE)$minimum,
        
        # return fit
        return(
          tibble(
            model = "linear", adj_r_sq = lm_adj_r_sq, p_val = lm_p, formula = lm_formula, mt_at_min_degree = NA
          ) %>% 
            bind_rows(
              tibble(
                model = "parabolic", adj_r_sq = pm_adj_r_sq, p_val = pm_p, formula = pm_formula, mt_at_min_degree = pm_mt_at_min_degree
              )
            )
        )
        
      )
    )
    
  ) %>% unnest(cols = c(model_res))

df_mt_degree_slopes_models

# plot
df_mt_degree_slopes_models %>% 
  filter(model == "parabolic") %>% 
  mutate(label = paste0(formula, "\nR^2 = ", round(adj_r_sq, 2), "; p = ", format(p_val, scientific = TRUE, digits = 3))) %>% 
  unnest(cols = c(data)) %>%
  mutate(hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor,
         period = factor(period, levels = c("early", "adult", "lifespan"))
  ) %>% 
  
  # plot aspects
  ggplot(aes(x = mt, y = degree)) +
  geom_point(aes(color = hubs), shape = 1, size = 2) +
  geom_smooth(formula = y ~ x + I(x^2), method = "lm") +
  geom_vline(aes(xintercept = 0), lty = 2, color = "black") +
  geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
  
  # model formula and fit
  geom_text(
    x = -0.0004, y = -0.1, color = "blue", size = 2.5, hjust = 1, align = "center",
    aes(label = label), check_overlap = TRUE
  ) +
  
  # MT at max(degree)
  geom_text(
    y = 0.275, color = "red", hjust = 0.5, check_overlap = TRUE, size = 2.5,
    aes(x = mt_at_min_degree, label= paste0("MT slope at \nmin(degree slope) = \n", format(mt_at_min_degree, scientific = TRUE, digits = 3))), 
    
  ) +
  
  # aesthetics
  facet_wrap(vars(period), labeller = as_labeller(labels), scales = "free_x") +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "MT slope", y = "Degree slope",
       title = "Relationship between change in MT and degree during different developmental periods"
  ) +
  theme(axis.text.x = element_text(size = 5))

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5d.mt_degree_slope_relationship.pdf"), width = 8, height = 3)
ggsave(paste0(base_dir, "outputs/figures/5d.mt_degree_slope_relationship.png"), width = 8, height = 3)


# relationship between baseline feature and feature development -----------

# set up data frame
df_mt_degree_development_and_slopes <- df_mt_degree %>% 
  filter(timepoint != 35) %>% 
  group_by(phenotype, timepoint, region_of_interest) %>% 
  summarise(value = median(feature)) %>% 
  unite("feature", c(phenotype, timepoint), sep = "_") %>% 
  pivot_wider(id_cols = region_of_interest, names_from = feature, values_from = value) %>% 
  
  left_join(
    df_mt_degree_slopes %>% 
      filter(period == "early") %>% 
      mutate(phenotype = paste0(phenotype, "_slope")) %>% 
      pivot_wider(id_cols = region_of_interest, names_from = phenotype, values_from = slope),
    by = join_by(region_of_interest)
  )

# run linear and parabolic models
df_mt_degree_development_and_slopes_models <- df_mt_degree_development_and_slopes %>% 
  pivot_longer(2:ncol(.), names_to = "feature", values_to = "value") %>% 
  separate(feature, into = c("phenotype", "feature"), sep = "_") %>% 
  filter(feature != 63) %>% 
  mutate(feature = case_when(
    feature == 20 ~ "baseline",
    feature == 300 ~ "adult",
    TRUE ~ feature
  )) %>%
  pivot_wider(id_cols = c(region_of_interest, phenotype), names_from = feature, values_from = value) %>% 
  mutate(
    baseline2 = baseline^2,
    adult2 = adult^2,
    slope2 = slope^2
  ) %>% 
  group_by(phenotype) %>% 
  nest() %>% 
  mutate(
    
    model_res = map(
      .x = data,
      .f = ~ list(
        
        # BASELINE
        
        # run linear model
        baseline_linear_model <- lm(slope ~ baseline, data = .x), 
        baseline_lm_adj_r_sq <- summary(baseline_linear_model)$adj.r.squared,
        baseline_lm_f <- summary(baseline_linear_model)$fstatistic,
        baseline_lm_p <- pf(baseline_lm_f[1], baseline_lm_f[2], baseline_lm_f[3], lower.tail = FALSE),
        baseline_lm_coefs <- coefficients(baseline_linear_model) %>% round(2),
        baseline_lm_formula <- paste0("y = ", baseline_lm_coefs[1], " + ", baseline_lm_coefs[2], "*x"),
        
        # run parabolic model
        baseline_parabolic_model <- lm(slope ~ baseline + baseline2, data = .x), 
        baseline_pm_adj_r_sq <- summary(baseline_parabolic_model)$adj.r.squared,
        baseline_pm_f <- summary(baseline_parabolic_model)$fstatistic,
        baseline_pm_p <- pf(baseline_pm_f[1], baseline_pm_f[2], baseline_pm_f[3], lower.tail = FALSE),
        baseline_pm_coefs <- coefficients(baseline_parabolic_model) %>% round(2),
        baseline_pm_formula <- paste0("y = ", baseline_pm_coefs[1], " + ", baseline_pm_coefs[2], "*x + ", baseline_pm_coefs[3], "*x^2"),
        
        # ADULT
        
        # run linear model
        adult_linear_model <- lm(slope ~ adult, data = .x), 
        adult_lm_adj_r_sq <- summary(adult_linear_model)$adj.r.squared,
        adult_lm_f <- summary(adult_linear_model)$fstatistic,
        adult_lm_p <- pf(adult_lm_f[1], adult_lm_f[2], adult_lm_f[3], lower.tail = FALSE),
        adult_lm_coefs <- coefficients(adult_linear_model) %>% round(2),
        adult_lm_formula <- paste0("y = ", adult_lm_coefs[1], " + ", adult_lm_coefs[2], "*x"),
        
        # run parabolic model
        adult_parabolic_model <- lm(slope ~ adult + adult2, data = .x), 
        adult_pm_adj_r_sq <- summary(adult_parabolic_model)$adj.r.squared,
        adult_pm_f <- summary(adult_parabolic_model)$fstatistic,
        adult_pm_p <- pf(adult_pm_f[1], adult_pm_f[2], adult_pm_f[3], lower.tail = FALSE),
        adult_pm_coefs <- coefficients(adult_parabolic_model) %>% round(2),
        adult_pm_formula <- paste0("y = ", adult_pm_coefs[1], " + ", adult_pm_coefs[2], "*x + ", adult_pm_coefs[3], "*x^2"),
        
        # RETURN TIBBLE OF RESULTS
        return(
          tibble(
            time = "baseline", model = "linear", adj_r_sq = baseline_lm_adj_r_sq, p_val = baseline_lm_p, formula = baseline_lm_formula
          ) %>% 
            bind_rows(
              tibble(
                time = "baseline", model = "parabolic", adj_r_sq = baseline_pm_adj_r_sq, p_val = baseline_pm_p, formula = baseline_pm_formula
              )
            ) %>% 
            
            bind_rows(
              tibble(
                time = "adult", model = "linear", adj_r_sq = adult_lm_adj_r_sq, p_val = adult_lm_p, formula = adult_lm_formula
              )
            ) %>% 
            bind_rows(
              tibble(
                time = "adult", model = "parabolic", adj_r_sq = adult_pm_adj_r_sq, p_val = adult_pm_p, formula = adult_pm_formula
              )
            )
        )
        
      )
    )
    
  ) %>% unnest(cols = c(model_res))

# plot MT
labels <- c("Baseline (PND 20)", "Adult (PND 300)")
names(labels) <- c("baseline", "adult")

p_mt <- df_mt_degree_development_and_slopes_models %>% 
  filter(phenotype == "mt" & model == "linear") %>% 
  mutate(label = paste0(formula, "\nR^2 = ", round(adj_r_sq, 2), "; p = ", format(p_val, scientific = TRUE, digits = 3))) %>% 
  unnest(cols = c(data)) %>%
  mutate(hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor) %>% 
  pivot_longer(c(baseline, adult), names_to = "time2", values_to = "value") %>% 
  filter(time2 == time) %>% dplyr::select(-time2) %>% 
  mutate(time = factor(time, levels = c("baseline", "adult"))) %>% 
  
  # plot aspects
  ggplot(aes(x = value, y = slope)) +
  geom_point(aes(color = hubs), shape = 1, size = 2) +
  geom_smooth(formula = y ~ x, method = "lm") +

  # model formula and fit
  geom_text(
    x = Inf, y = -0.004, color = "blue", size = 2.5, hjust = 1, align = "center",
    aes(label = label), check_overlap = TRUE
  ) +

  # aesthetics
  facet_wrap(vars(time), labeller = as_labeller(labels), scales = "free_x") +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "MT value", y = "MT slope",
       title = "Relationship between MT and change in MT in early development"
  )

p_degree <- df_mt_degree_development_and_slopes_models %>% 
  filter(phenotype == "degree" & model == "linear") %>% 
  mutate(label = paste0(formula, "\nR^2 = ", round(adj_r_sq, 2), "; p = ", format(p_val, scientific = TRUE, digits = 3))) %>% 
  unnest(cols = c(data)) %>%
  mutate(hubs = ifelse(region_of_interest %in% hubs, 1, 0) %>% factor) %>% 
  pivot_longer(c(baseline, adult), names_to = "time2", values_to = "value") %>% 
  filter(time2 == time) %>% dplyr::select(-time2) %>% 
  mutate(time = factor(time, levels = c("baseline", "adult"))) %>% 
  
  # plot aspects
  ggplot(aes(x = value, y = slope)) +
  geom_point(aes(color = hubs), shape = 1, size = 2) +
  geom_smooth(formula = y ~ x, method = "lm") +
  
  # model formula and fit
  geom_text(
    x = 60, y = -0.1, color = "blue", size = 2.5, hjust = 1, align = "center",
    aes(label = label), check_overlap = TRUE
  ) +
  
  # aesthetics
  facet_wrap(vars(time), labeller = as_labeller(labels)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "degree value", y = "degree slope",
       title = "Relationship between degree and change in degree in early development"
  )

p_mt / p_degree

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5e.mt_degree_baseline_adult_slope.pdf"), width = 7, height = 5)
ggsave(paste0(base_dir, "outputs/figures/5e.mt_degree_baseline_adult_slope.png"), width = 7, height = 5)


# Top 10% plots -----------------------------------------------------------

n_ROIs <- df_mt_degree_slopes %>% pull(region_of_interest) %>% unique %>% length

df_mt_degree_slopes %>%
  filter(period == "early") %>% 
  left_join(df_sys_to_roi, by = join_by(region_of_interest)) %>% 
  filter(region_of_interest %in% gray_matter_rois) %>% 
  group_by(phenotype, region_of_interest) %>%
  summarise(
    slope = median(slope),
    standard_error = median(standard_error)
  ) %>% 
  slice_max(abs(slope), n = round(0.1*n_ROIs, 0)) %>%
  
  mutate(#region_of_interest = str_replace(region_of_interest, "_", " "),
         #region_of_interest = str_replace(region_of_interest, "_", "\n"),
         sign = ifelse(slope < 0, "negative", "positive"),
         label = ifelse(region_of_interest %in% hubs, "*", "")
  ) %>% 
  
  ggplot(aes(x = reorder_within(region_of_interest, abs(slope), phenotype), y = slope)) +
  geom_errorbar(aes(ymin = slope - standard_error, ymax = slope + standard_error)) +
  geom_point(aes(color = sign), size = 3) +
  geom_hline(aes(yintercept = 0), color = "black", lty = 2) +
  geom_text(aes(y = slope + standard_error + 0.01, label = label), size = 7, vjust = 0.75) +
  facet_wrap(vars(phenotype), scales = "free") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  scale_x_reordered() +
  labs(x = "", y = "PND 20 --> 63 slope (+/- standard error)",
       title = "Top 10% most developmentally active (highest abs(slope)) ROIs for degree & MT") +
  coord_flip() +
  theme(legend.position = "none"#,
        #axis.text.y = element_text(size = 11)
  )

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5f.mt_degree_developmentally_active.pdf"), width = 10, height = 4)
ggsave(paste0(base_dir, "outputs/figures/5f.mt_degree_developmentally_active.png"), width = 10, height = 4)


## PLOT ON BRAIN

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

p_mt_DA_brain <- df_slices %>% 
  left_join(df_mt_degree_slopes %>%
              filter(period == "early") %>% 
              filter(phenotype == "mt"), 
            by = "region_of_interest") %>% 
  mutate(slope = ifelse(region_of_interest %in% white_matter_rois, NA, slope)) %>% 
  arrange(!is.na(slope)) %>% 
  
  ggplot() +
  geom_sf(aes(color = slope, fill = slope, geometry = geometry, group = -1)) +
  scale_fill_gradient(low = "#2166AC", high = "white", na.value = "grey") +
  scale_color_gradient(low = "#2166AC", high = "white", na.value = "grey") +
  facet_wrap(vars(side), nrow = 1) +
  labs(title = "MT") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.position = "bottom"
  )

p_degree_DA_brain <- df_slices %>% 
  left_join(df_mt_degree_slopes %>%
              filter(period == "early") %>% 
              filter(phenotype == "degree"), 
            by = "region_of_interest") %>% 
  mutate(slope = ifelse(region_of_interest %in% white_matter_rois, NA, slope)) %>% 
  arrange(!is.na(slope)) %>% 
  
  ggplot() +
  geom_sf(aes(color = slope, fill = slope, geometry = geometry, group = -1)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", na.value = "grey") +
  scale_color_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", na.value = "grey") +
  facet_wrap(vars(side), nrow = 1) +
  labs(title = "Degree") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.position = "bottom"
  )


p_mt_DA_brain / p_degree_DA_brain

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5f.mt_degree_developmentally_active_brains.pdf"), width = 4, height = 4)
ggsave(paste0(base_dir, "outputs/figures/5f.mt_degree_developmentally_active_brains.png"), width = 4, height = 4)

### CENTROID COORDINATES VS SLOPE

labels <- c("x (left --> right)", "y (posterior --> anterior)", "z (inferior --> superior)")
names(labels) <- c("x", "y", "z")

df_centroids %>% 
  pivot_longer(2:ncol(.), names_to = "dim", values_to = "coord") %>% 
  left_join(df_mt_degree_slopes %>%
              filter(period == "early"), 
            by = join_by(region_of_interest)
  ) %>% 
  filter(!is.na(phenotype)) %>% 
  mutate(sign = ifelse(slope < 0, "negative", "positive")) %>% 

  ggplot(aes(x = coord, y = slope)) +
  geom_point(aes(fill = slope, color = sign), shape = 21, size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_cor(color = "black", size = 3) +
  #geom_text_repel(aes(label = system)) +
  facet_grid(phenotype ~ dim, scales = "free",
             labeller = labeller(dim = labels)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = "Coordinate position", y = "Early development slope",
       title = "ROI-level MT & degree slope relative to anatomical position (RH)") +
  theme(legend.position = "none")

# SAVE
ggsave(paste0(base_dir, "outputs/figures/5f.mt_degree_developmentally_active_coords.pdf"), width = 8, height = 5)
ggsave(paste0(base_dir, "outputs/figures/5f.mt_degree_developmentally_active_coords.png"), width = 8, height = 5)

