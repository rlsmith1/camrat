
###
###
###

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
library(cluster)
library(RColorBrewer)
library(ggside)
library(igraph)
library(ggrepel)
library(corrr)
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

# degree distribution ------------------------------------------------------------

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

# FILTER FOR ADULT ONLY, TAKE MEDIAN
df_degree_control_adult <- df_degree %>% 
  filter(study == "MRC" & timepoint == 300) %>% 
  group_by(region_of_interest) %>% 
  summarise(weighted_degree = median(weighted_degree, na.rm = TRUE))

# SW TEST FOR NORMALITY
shapiro.test(df_degree_control_adult$weighted_degree)


# hubs ---------------------------------------------------------------------------

# CALCULATE CONFIDENCE INTERVALS AROUND HUB SCORES
df_hubs <- df_degree %>% 
  filter(study == "MRC" & timepoint == 300) %>% 
  group_by(region_of_interest) %>% 
  left_join(
    df_degree %>% 
      filter(study == "MRC" & timepoint == 300) %>% 
      group_by(region_of_interest) %>% 
      count()
  ) %>% 
  summarise(
    median_hub = median(hub_score),
    sd_hub = sd(hub_score),
    ci = qnorm(0.975)*sd_hub/sqrt(n)
  ) %>% 
  distinct() %>% 
  arrange(-median_hub) %>% 
  ungroup %>% 
  left_join(df_roi_to_abbrev, by = join_by(region_of_interest))

# WHAT MODULES DO HUBS BELONG TO
hubs <- df_hubs %>% 
  top_n(n = 0.1*nrow(.), wt = median_hub) %>% 
  pull(region_of_interest)

df_mind_mods %>% 
  filter(region_of_interest %in% hubs)


# prep for ball and stick plot --------------------------------------------

# calculate node position
df_node_coords <- df_sigma_atlas %>% 
  filter(side == "sagittal", hemisphere == "right" & roi_abbreviation != "mcp" & region_of_interest %in% gray_matter_rois) %>% 
  group_by(side, region_of_interest, roi_abbreviation) %>%
  summarize(geometry = st_union(geometry)) %>% 
  dplyr::select(side, region_of_interest, roi_abbreviation) %>% 
  st_centroid %>% 
  mutate(point = map(geometry, ~st_coordinates(.x) %>% as.data.frame %>% as_tibble)) %>% 
  unnest(cols = c(point)) %>% 
  as_tibble() %>% 
  dplyr::select(-geometry) %>% 
  clean_names() %>% 
  
  # add degree and module to nodes
  left_join(df_degree %>% 
              filter(study == "MRC" & timepoint == 300) %>% 
              summarise(degree = median(feature_resids))
  ) %>% 
  left_join(df_mind_mods, by = join_by(region_of_interest)) %>% 
  mutate(label = ifelse(region_of_interest %in% hubs, roi_abbreviation, NA))

# edge weights
df_edge_weights <- df_mind_GM %>% 
  filter(study == "MRC" & timepoint == 300) %>% 
  group_by(R1, R2) %>% 
  summarise(weight = median(weight)) %>% 
  dplyr::select(R1, R2, weight) %>% 
  left_join(df_node_coords %>% dplyr::rename("R1" = "region_of_interest", "roi_abbrev1" = "roi_abbreviation"), by = join_by(R1)) %>% 
  dplyr::rename("x1" = "x", "y1" = "y") %>% 
  left_join(df_node_coords %>% dplyr::rename("R2" = "region_of_interest", "roi_abbrev2" = "roi_abbreviation"), by = join_by(R2)) %>% 
  dplyr::rename("x2" = "x", "y2" = "y") %>% 
  filter(side.x == side.y) %>% 
  dplyr::rename("side" = "side.x")


# relationship between MT and degree --------------------------------------

df_mt_degree_cor <- df_mt %>% 
  filter(study == "MRC" & timepoint == 300 & region_of_interest %in% gray_matter_rois) %>% 
  summarise(median_mt = median(feature_resids)) %>% 
  
  left_join(
    df_degree %>% 
      filter(study == "MRC" & timepoint == 300) %>% 
      group_by(region_of_interest) %>% 
      summarise(median_degree = median(feature_resids)),
    by = "region_of_interest"
  ) %>% 
  mutate(median_mt2 = median_mt^2,
         median_mt3 = median_mt^3
  )

lm(median_degree ~ median_mt, data = df_mt_degree_cor) %>% summary # simple linear model explains 5% of the variance, p = 0.03
parabolic_model <- lm(median_degree ~ median_mt + median_mt2, data = df_mt_degree_cor) 
parabolic_model %>% summary # parabolic model explains 34% of the variance, p = 7*10^-8
lm(median_degree ~ median_mt + median_mt2 + median_mt3, data = df_mt_degree_cor) %>% summary # slightly worse fit than parabolic function

# A: hubs -----------------------------------------------------------------

# HUB SCORE
p_fig4a_score <- df_hubs %>% 
  filter(region_of_interest %in% hubs) %>% 
  mutate(region_of_interest = str_replace_all(region_of_interest, "_", " ")) %>% 
  
  ggplot(aes(y = reorder(region_of_interest, median_hub))) +
  geom_errorbar(aes(xmin = median_hub - ci, xmax = median_hub + ci), width = 0.5) +
  geom_point(aes(x = median_hub, alpha = median_hub), color = module_colors["2"], size = 3) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  xlim(c(0.85, 1.0)) +
  scale_x_continuous(breaks = c(0.90, 0.95, 1.0), labels = c(0.90, 0.95, 1.0)) +
  labs(x = "hub score (+/- 95% CI)", y = "",
       title = "Median hub score of top 10% ROIs") +
  theme(legend.position = "none")

p_fig4a_score

# SAVE
ggsave(paste0(base_dir, "outputs/figures/4a.hubs.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/4a.hubs.png"), width = 5, height = 4)

# ON BRAIN
slices <- list("sagittal" = 130,
               "axial" = 105,
               "coronal" = 155
)
df_slices <- map_dfr(.x = 1:length(slices),
                     .f = ~ df_sigma_atlas %>%
                       filter(region_of_interest != "striatum") %>% 
                       filter(side == names(slices[.x]) & slice == slices[.x])
                     
) %>% 
  dplyr::select(region_of_interest, side, slice, geometry) %>% 
  distinct()

p_fig4a_brain <- df_slices %>% 
  left_join(df_degree %>% 
              filter(study == "MRC" & timepoint == 300) %>% 
              mutate(hub_score = median(hub_score)) %>% 
              mutate(hub = region_of_interest %in% hubs), 
            by = "region_of_interest") %>% 
  arrange(!is.na(hub_score), hub_score) %>% 
  
  ggplot() +
  geom_sf(aes(fill = hub_score, color = hub, geometry = geometry, group = -1), 
          lwd = 0.5) +
  scale_fill_viridis(na.value = "grey") +
  scale_color_manual(values = c("darkgrey", "red"), guide = "none") +
  facet_wrap(vars(side), nrow = 2) +
  labs(title = "Hub score map in control adult CMN") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = c(0.75, 0.25)
        )

p_fig4a_brain

# SAVE
ggsave(paste0(base_dir, "outputs/figures/4a.hubs_on_brain.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/4a.hubs_on_brain.png"), width = 5, height = 4)


# B: ball-and-stick network -----------------------------------------------

threshold <- 0.1

p_fig4b <- ggplot() +
  geom_segment(aes(x = x1, y = y1,
                   xend = x2, yend = y2, 
                   color = weight, alpha = weight, linewidth = weight), 
               data = df_edge_weights %>% arrange(-weight) %>% head(threshold*nrow(.))
  ) +
  geom_point(aes(x = x, y = y, size = degree, fill = module), 
             data = df_node_coords, shape = 21) +
  geom_text_repel(aes(x = x, y = y, label = label),
                  data = df_node_coords, size = 5) +
  scale_color_viridis(limits = c(0.63, 0.96), guide = "none") +
  scale_fill_manual(values = module_colors) +
  scale_size_continuous(range = c(3, 6)) +
  scale_linewidth_continuous(range = c(1, 3), guide = "none") +
  scale_alpha_continuous(guide = "none") +
  facet_wrap(vars(side)) +
  labs(title = "Ball-and-stick representation of top 10% of connections (hubs labeled)") +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_blank()
  )

p_fig4b

# SAVE
ggsave(paste0(base_dir, "outputs/figures/4b.hubs_ball-and-stick.pdf"), width = 10, height = 4)
ggsave(paste0(base_dir, "outputs/figures/4b.hubs_ball-and-stick.png"), width = 10, height = 4)


# C: relationship between degree and MT^2 ------------------------------------

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

p_fig4c

# SAVE
ggsave(paste0(base_dir, "outputs/figures/4c.mt-degree_relationship.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/4c.mt-degree_relationship.png"), width = 5, height = 4)

# S3a: Degree distribution of control adult CMN ---------------------------

# DISTRIBUTION
p_figS3a_dist <- df_degree_control_adult %>% 
  ggplot(aes(x = weighted_degree)) +
  geom_density(fill = "grey", trim = FALSE) +
  #facet_wrap(vars(subject)) +
  labs(x = "Median weighted degree", y = "",
       title = "Control adult CMN \ndegree distribution") +
  xlim(c(0, 105)) +
  theme(legend.position = "none")

# CDF
p_figS3a_cdf <- df_degree_control_adult %>% 
  ggplot(aes(x = weighted_degree)) +
  stat_ecdf(linewidth = 2, alpha = 0.7) +
  xlim(c(0, 105)) +
  labs(title = "ECDF", x = "Weighted degree", y = "ECDF") +
  theme(legend.position = "none")

# QQ
p_figS3a_qq <- df_degree_control_adult %>% 
  ggplot(aes(sample = weighted_degree)) +
  stat_qq(distribution = stats::qnorm, alpha = 0.7) + 
  stat_qq_line(distribution = stats::qnorm) +
  labs(title = "Normal Q-Q plot", x = "Normal theoretical quantiles", y = "Data quantiles") +
  theme(legend.position = "none")

p_figS3a_dist + p_figS3a_cdf + p_figS3a_qq

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S3a.CMN_degree_distribution.pdf"), width = 10, height = 4)
ggsave(paste0(base_dir, "outputs/figures/S3a.CMN_degree_distribution.png"), width = 10, height = 4)


# Rich club analysis ------------------------------------------------------

# GENERATE NULL NETS BY RESHUFFLING WEIGHTS

# calculate control average
df_mind_ctrl_avg <- df_mind_GM %>% 
  filter(study == "MRC" & timepoint == 300) %>% 
  group_by(R1, R2) %>% 
  summarise(weight = median(weight))

# generate 1000 null networks by reshuffling weights
set.seed(2107)
doParallel::registerDoParallel()
df_null_nets <- tibble(perm = seq(1, 1000)) %>% 
  expand_grid(df_mind_ctrl_avg) %>% 
  group_by(perm) %>% 
  nest() %>% 
  mutate(
    
    data = map(
      
      .x = data, 
      .f = ~ .x %>% 
        mutate(weight = sample(df_mind_ctrl_avg %>% pull(weight)))
      
    )
    
  ) %>% unnest(cols = c(data))

# CALCULATE RICH CLUB COEFFICIENT IN EACH NULL NETWORK
df_null_nets_richClub <- df_null_nets %>%
  group_by(perm) %>% 
  nest() %>% 
  expand_grid(threshold = c(0.05, 0.10, 0.15)) %>% 
  
  mutate(
    
    # 1. identify prominent nodes
    hubs = map2(
      .x = data,
      .y = threshold,
      .f = ~ .x %>% 
        graph_from_data_frame(directed = FALSE) %>% 
        hub_score() %>% .$vector %>% 
        enframe %>% 
        arrange(-value) %>% 
        head(ceiling(.y*nrow(.))) %>% 
        pull(name)
    ),
    
    # 2. Calculate the total sum of weights attached to ties among the prominent nodes
    hub_sum =
      
      map2(
        .x = data,
        .y = hubs,
        .f = ~ .x %>% 
          filter(R1 %in% .y & R2 %in% .y) %>% 
          ungroup %>% 
          summarise(hub_sum = sum(weight, na.rm = TRUE)) %>% 
          pull(hub_sum) 
      ),
    
    # 3. Calculate the total sum of same number of ties, but the strongest ones, in the network
    strong_sum = map2(
      .x = data,
      .y = hubs,
      .f = ~ .x %>% 
        slice_max(order_by = weight, n = length(.y)) %>% 
        ungroup %>% 
        summarise(strong_sum = sum(weight, na.rm = TRUE))
    )
    
  ) %>% unnest(cols = c(hub_sum, strong_sum)) %>% 
  mutate(perm_phi = hub_sum/strong_sum) %>% 
  dplyr::select(-data)

# CALCULATE REAL CONTROL AVG NETWORK RICH CLUB COEF
df_ctrl_avg_richClub <- df_mind_ctrl_avg %>% 
  expand_grid(threshold = c(0.05, 0.10, 0.15)) %>% 
  group_by(threshold) %>% 
  nest() %>% 
  mutate(
    
    # 1. identify prominent nodes
    hubs = map2(
      .x = data,
      .y = threshold,
      .f = ~ .x %>% 
        graph_from_data_frame(directed = FALSE) %>% 
        hub_score() %>% .$vector %>% 
        enframe %>% 
        arrange(-value) %>% 
        head(ceiling(.y*nrow(.))) %>% 
        pull(name)
    ),
    
    # 2. Calculate the total sum of weights attached to ties among the prominent nodes
    hub_sum = map2(
      .x = data,
      .y = hubs,
      .f = ~ .x %>% 
        filter(R1 %in% .y & R2 %in% .y) %>% 
        ungroup %>% 
        summarise(hub_sum = sum(weight, na.rm = TRUE)) %>% 
        pull(hub_sum) 
    ),
    
    # 3. Calculate the total sum of same number of ties, but the strongest ones, in the network
    strong_sum = map2(
      .x = data,
      .y = hubs,
      .f = ~ .x %>% 
        slice_max(order_by = weight, n = length(.y)) %>% 
        ungroup %>% 
        summarise(strong_sum = sum(weight, na.rm = TRUE))
    )
    
  ) %>% 
  
  unnest(cols = c(hub_sum, strong_sum)) %>% 
  mutate(rc_phi = hub_sum/strong_sum) %>% 
  dplyr::select(-data)

# SAVE OBJECTS
save(df_null_nets_richClub, df_ctrl_avg_richClub, 
     file = paste0(base_dir, "objects/30Aug2023_rich_club_permutations.Rdata"))


# S3b: Rich club ----------------------------------------------------------

load(paste0(base_dir, "objects/30Aug2023_rich_club_permutations.Rdata"))

# BRAIN
p_figS3b_brain <- ggplot() +
  geom_segment(aes(x = x1, y = y1,
                   xend = x2, yend = y2, 
                   color = weight, alpha = weight), 
               data = df_edge_weights %>% filter(R1 %in% hubs & R2 %in% hubs)) +
  geom_point(aes(x = x, y = y, size = degree, fill = module), 
             data = df_node_coords, shape = 21) +
  geom_text_repel(aes(x = x, y = y, label = label),
                   data = df_node_coords, size = 4) +
  scale_color_viridis(guide = "none") +
  scale_fill_manual(values = module_colors) +
  scale_size_continuous(range = c(3, 6)) +
  scale_alpha_continuous(guide = "none") +
  facet_wrap(vars(side)) +
  labs(title = "The rich club of the control adult CMN") +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_blank())

# COEFFICIENTS
p_figS3b_coefs <- df_null_nets_richClub %>% 
  left_join(df_ctrl_avg_richClub %>% dplyr::select(threshold, rc_phi)) %>% 
  
  ggplot() +
  geom_density(aes(x = perm_phi), fill = "midnightblue", alpha = 0.5) +
  geom_point(aes(x = rc_phi, y = 0), color = "midnightblue", size = 3) +
  facet_wrap(vars(threshold), nrow = 1) +
  labs(x = "Rich club coefficient", y = "", #y = "Proportion of hubs considered 'rich club'", 
       title = "Actual RC coefficient vs RC coefficients across 1000 random null networks")

p_figS3b_brain / p_figS3b_coefs + plot_layout(heights = c(1, 0.5))

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S3b.rich_club.pdf"), width = 10, height = 6)
ggsave(paste0(base_dir, "outputs/figures/S3b.rich_club.png"), width = 10, height = 6)


# generalizability analysis -----------------------------------------------

df_mind_cor_all <- df_mind_GM %>% 
  filter(study == "MRC") %>% 
  
  # Correlate each MIND network with each other network
  mutate(subject_timepoint = paste0(subject, sep = "_", timepoint)) %>% 
  pivot_wider(id_cols = c(R1, R2), names_from = subject_timepoint, values_from = weight) %>% 
  dplyr::select(-R1, -R2) %>% 
  correlate() %>% 
  
  # get metadata for each network for comparison
  dplyr::rename("net1" = "term") %>% 
  pivot_longer(2:ncol(.), names_to = "net2", values_to = "cor") %>% 
  
  # remove duplicate in reverse order
  group_by(grp = paste0(pmin(net1, net2), ".", pmax(net1, net2))) %>% 
  slice(1) %>% 
  dplyr::select(-grp) %>% 
  
  # separate to join with metadata
  separate(net1, into = c("subject1", "timepoint1"), sep = "_") %>% 
  separate(net2, into = c("subject2", "timepoint2"), sep = "_") %>% 
  
  left_join(df_mind_GM %>% 
              dplyr::select(study, subject, timepoint, sex, group, age) %>% 
              distinct() %>% 
              dplyr::rename_all(~paste0(.x, "1")),
            by = c("subject1", "timepoint1")
  ) %>% 
  left_join(df_mind_GM %>% 
              dplyr::select(study, subject, timepoint, sex, group, age) %>% 
              distinct() %>% 
              dplyr::rename_all(~paste0(.x, "2")),
            by = c("subject2", "timepoint2")
  ) %>% 
  
  mutate(label1 = paste0(study1, "_", group1),
         label2 = paste0(study2, "_", group2)
  ) %>% 
  dplyr::select(subject1, timepoint1, label1, subject2, timepoint2, label2, cor)  %>% 
  mutate(edge = paste0("PND ", timepoint1, " - ", "PND ", timepoint2))


# S3c: CMN is generalizable across subjects -------------------------------

edge_order <- df_mind_cor_all %>%
  filter(label1 == "MRC_external_control" & label2 == "MRC_external_control") %>% 
  group_by(edge) %>% 
  summarise(median_r = median(cor, na.rm = TRUE)) %>% 
  arrange(median_r) %>% 
  pull(edge)

p_figS3c <- df_mind_cor_all %>% 
  filter(label1 == "MRC_external_control" & label2 == "MRC_external_control") %>% 
  mutate(edge = factor(edge, levels = edge_order)) %>% 
  
  ggplot(aes(y = cor, x = edge, fill = edge)) +
  geom_half_point(side = 'l', alpha = 0.5) +
  geom_half_violin(side = 'r', alpha = 0.75) +
  stat_summary(fun = "median", geom = "crossbar", aes(color = edge)) +
  ylim(c(0, 1)) +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C") +
  labs(title = "Distribution of correlations of each MRC CMN with every other",
       x = "", y = "Pearon's r") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

p_figS3c

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S3c.generalizable_CMN.pdf"), width = 6, height = 4)
ggsave(paste0(base_dir, "outputs/figures/S3c.generalizable_CMN.png"), width = 6, height = 4)

