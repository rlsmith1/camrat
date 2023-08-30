
###
### Establishing modularity of adult CMN and correlating with Swanson-Sporns tract-tracing
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
library(ggalluvial)
library(ggside)

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

# load data ---------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

load(paste0(base_dir, "objects/25Aug2023_df_data.RDS")) # df_data
load(paste0(base_dir, "objects/25Aug2023_sigma_labels.RDS")) # df_sigma_labels
load(paste0(base_dir, "objects/25Aug2023_sigma_atlas_for_plotting.RDS")) # df_sigma_atlas
load(paste0(base_dir, "objects/25Aug2023_df_mind.RDS")) # df_mind

# CONVERT ROI TO SYSTEM LEVEL
df_sys_to_roi <- df_sigma_labels %>% 
  dplyr::select(matter, system, region_of_interest) %>% 
  distinct()

### ADD SYSTEM LEVEL TO DF_MIND_GM
df_mind_GM <- df_mind_GM %>% 
  left_join(df_sys_to_roi %>% dplyr::rename("R1" = "region_of_interest")) %>% 
  dplyr::rename("S1" = "system") %>% 
  left_join(df_sys_to_roi %>% dplyr::rename("R2" = "region_of_interest")) %>% 
  dplyr::rename("S2" = "system")

### DATA FOR SWANSON-SPORNS CORRELATION

# Swanson-Sporns tract-tracing data

# 2020
# df_swanson <- read_xlsx(paste0(base_dir, "data/Swanson_Sporns/swanson_sporns_2020.xlsx"), sheet = 7) %>%
#   filter_at(1, all_vars(. == 1))  %>%
#   dplyr::select(c("...4", "To side", starts_with("1")))
# 
# swanson_rois <- df_swanson %>% pull(2)

# 2022
df_swanson <- read_xlsx(paste0(base_dir, "data/swanson_sporns_2022.xlsx"), sheet = 7) %>%
  filter_at(1, all_vars(. == 1))  %>%
  dplyr::select(c("...14", "Side", starts_with("1")))

swanson_rois <- df_swanson %>% pull(2)

# SIGMA to Swanson mappings
df_sigma_to_swanson <- read_xlsx(paste0(base_dir, "data/SIGMA_swanson_mappings.xlsx"), sheet = 1)

# Modularity analysis ----------------------------------------------------------------

### DETERMINE ROI ORDER FOR PLOTTING

# convert MRC PND 300 (control adult) CMN to matrix
median_mind_MRC300 <- df_mind_GM %>% 
  filter(timepoint == 300 & study == "MRC") %>% 
  group_by(R1, R2) %>% 
  summarise(weight = median(weight, na.rm = TRUE)) %>% 
  #filter(R1 != R2) %>% 
  ungroup %>% 
  pivot_wider(id_cols = R1, names_from = R2, values_from = weight) %>% 
  column_to_rownames("R1") %>% 
  as.matrix()

# identify best ROI order based on hclust
roi_order <- median_mind_MRC300[hclust(dist(median_mind_MRC300))$order,] %>% rownames

# get ROI abbreviations
df_roi_to_abbrev <- df_sigma_labels %>% 
  dplyr::select(region_of_interest, roi_abbreviation) %>% 
  distinct() %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = c(roi_order))) %>% 
  arrange(region_of_interest)

### SILHOUETTE PERMUTATIONS (FOR CLUSTERING)

control_subjects <- df_mind_GM %>% 
  filter(study == "MRC" & timepoint == 300) %>% 
  pull(subject) %>% unique
nperm <- 1000  
n <- 10
df_sil <- tibble()   

for (i in 1:nperm) {
  
  print(paste0("perm ", i))
  
  # recalculate MIND network for random sampling of controls
  m_mind <- df_mind_GM %>% 
    filter(study == "MRC" & timepoint == 300) %>% 
    filter(subject %in% sample(control_subjects, 0.75*length(control_subjects))) %>% 
    group_by(R1, R2) %>% 
    summarise(median_weight = median(weight)) %>% 
    pivot_wider(id_cols = R1,
                names_from = R2, 
                values_from = median_weight) %>% 
    column_to_rownames("R1") 
  
  # calculate silhouette widths for these controls
  sil_width <- 2:n %>% 
    map_dbl( ~ m_mind %>% 
               dist %>% 
               pam(k = .x) %>% 
               .$silinfo %>% .$avg.width
    )
  df_tmp <- tibble(perm = i, k = 2:n, sil_width = sil_width)
  df_sil <- df_sil %>% 
    bind_rows(df_tmp)
  
}     

### HCLUST ACROSS POSSIBLE K SOLUTIONS

set.seed(0829) 
n <- 10
df_many_models_hclust <- tibble(
  k = seq(3, n)
) %>% 
  mutate(
    hclust = map(
      .x = k,
      .f = ~ cutree(hclust(dist(median_mind_MRC300)), k = .x) %>%
        enframe %>%
        dplyr::rename("region_of_interest" = "name",
                      "module" = "value") %>%
        dplyr::mutate(module = factor(module))
    )
  ) %>% 
  unnest(cols = c(hclust))

# select best HClust
best_k <- 5
df_mind_mods <- df_many_models_hclust %>% 
  filter(k == best_k)

# list regions in modules
df_mind_mods %>% 
  group_by(module) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider(id_cols = row, names_from = module, values_from = region_of_interest) %>% 
  View()

### ESTABLISH MODULE LINES FOR NETWORK PLOTS
module_lines <- df_mind_mods %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = roi_order)) %>% 
  arrange(region_of_interest) %>% 
  mutate(module = factor(module, levels = unique(.$module))) %>% 
  count(module) %>% 
  mutate(lines = cumsum(n) + 0.5) %>% 
  head(nrow(.) - 1) %>% 
  pull(lines)

### SELECT MODULE COLORS
set.seed(0830)
module_colors <- brewer.pal %>% 
  
  # filter BrewerPal for colorblind-friendly palettes
  mapply(
    brewer.pal.info %>% filter(category == "qual" & colorblind == TRUE) %>% pull(maxcolors),
    brewer.pal.info %>% filter(category == "qual" & colorblind == TRUE) %>% rownames
  ) %>% 
  unlist() %>% 
  
  # remove gray as an options
  .[!str_detect(., "#B3B3B3|#666666")] %>% 
  
  # select 1 color per module
  sample(best_k)

# name module numbers using colors
names(module_colors) <- df_mind_mods %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = roi_order)) %>% 
  arrange(region_of_interest) %>% 
  mutate(module = factor(module, levels = unique(.$module))) %>% 
  pull(module) %>% levels

# tract-tracing analysis --------------------------------------------------


### CORRELATE WITH SWANSON-SPORNS TRACT-TRACING

# pivot Swanson data long and annotate with SIGMA system-level mappings
df_swanson_long <- df_swanson %>% 
  dplyr::rename_all(~c("name", "abbr1", swanson_rois)) %>% 
  dplyr::select(-name) %>% 
  pivot_longer(2:ncol(.), names_to = "abbr2", values_to = "swanson_weight") %>% 
  left_join(df_sigma_to_swanson %>% 
              dplyr::select(abbr, system) %>% 
              dplyr::rename("abbr1" = "abbr")) %>% 
  dplyr::rename("S1" = "system") %>% 
  left_join(df_sigma_to_swanson %>% 
              dplyr::select(abbr, system) %>% 
              dplyr::rename("abbr2" = "abbr")) %>% 
  dplyr::rename("S2" = "system") %>% 
  mutate(swanson_weight = as.numeric(swanson_weight))

# bin CMN to match Swanson distribution
swanson_cdf <-  df_swanson_long %>% 
  filter(swanson_weight != 0) %>%  # remove zeros
  count(swanson_weight) %>% 
  mutate(cdf = cumsum(n)/sum(n)) %>% 
  pull(cdf)

df_median_mind_bins <- df_mind_GM %>% 
  filter(timepoint == 300 & study == "MRC" & R1 != R2) %>%
  group_by(R1, R2, S1, S2) %>% 
  summarise(median_weight = median(weight)) %>% 
  group_by(grp = paste0(pmin(R1, R2), pmax(R1, R2))) %>% 
  slice(1) %>% 
  dplyr::select(-grp) %>% 
  arrange(median_weight) %>% 
  ungroup %>% 
  mutate(perc_of_total = row_number()/nrow(.),
         bin = case_when(
           
           perc_of_total <= swanson_cdf[1] ~ 1,
           perc_of_total <= swanson_cdf[2] ~ 2,
           perc_of_total <= swanson_cdf[3] ~ 3,
           perc_of_total <= swanson_cdf[4] ~ 4,
           perc_of_total <= swanson_cdf[5] ~ 5,
           perc_of_total <= swanson_cdf[6] ~ 6,
           perc_of_total <= swanson_cdf[7] ~ 7
           
         )
         
  )

# correlate binned networks
df_swanson_mind_cor <- df_swanson_long  %>% 
  
  filter(swanson_weight != 0) %>% 
  
  # take median of each system connection
  group_by(S1, S2) %>% 
  summarise(swanson_weight = mean(swanson_weight)) %>% 
  ungroup %>% 
  
  left_join(df_median_mind_bins %>% 
              group_by(S1, S2) %>% 
              summarise(mind_weight = mean(bin)), 
            by = join_by(S1, S2)) %>% 
  
  na.omit() %>% 
  filter(S1 != S2) %>% 
  mutate(norm_swanson = (swanson_weight - mean(swanson_weight))/sd(swanson_weight),
         norm_mind = (mind_weight - mean(mind_weight))/sd(mind_weight))

# A: Median PND 300 MRC CMN -----------------------------------------------------------

p_fig3a <- df_mind_GM %>%
  filter(timepoint == 300 & study == "MRC") %>% 
  group_by(R1, R2) %>% 
  summarise(weight = median(weight)) %>%  
  
  mutate(R1 = factor(R1, levels = roi_order),
         R2 = factor(R2, levels = roi_order)) %>%
  
  # plot
  ggplot(aes(x = R1, y = R2, fill = weight)) +
  geom_tile() +
  
  # add module lines
  geom_vline(xintercept = module_lines, color = "white", linewidth = 0.5) +
  geom_hline(yintercept = module_lines, color = "white", linewidth = 0.5) +

  annotate(ymin = c(-Inf, module_lines),
           ymax = c(module_lines, Inf),
           xmin = length(roi_order) + 0.5, xmax = length(roi_order) + 2,
           geom = "rect",
           fill = module_colors) +
  annotate(xmin = c(-Inf, module_lines),
           xmax = c(module_lines, Inf),
           ymin = length(roi_order) + 0.5, ymax = length(roi_order) + 2,
           geom = "rect",
           fill = module_colors) +
  
  labs(x = "", y = "", title = "Control adult (PND 300) median CMN") +
  scale_fill_viridis(na.value = "white", limits = c(0, 1)) +
  
  # use ROI abbreviation instead of full name
  scale_x_discrete(labels = df_roi_to_abbrev %>% pull(roi_abbreviation)) +
  scale_y_discrete(labels = df_roi_to_abbrev %>% pull(roi_abbreviation)) +
  
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text = element_blank(),
    legend.justification = "top",
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(2.0, "cm")
  )

p_fig3a

# SAVE
ggsave(paste0(base_dir, "outputs/figures/3a.CMN.pdf"), width = 6, height = 5)
ggsave(paste0(base_dir, "outputs/figures/3a.CMN.png"), width = 6, height = 5)

# B: Module solution on rat brain -----------------------------------------------------------

slices <- list("sagittal" = 120,
               "axial" = 115,
               "coronal" = 170
)
df_slices <- map_dfr(.x = 1:length(slices),
                     .f = ~ df_sigma_atlas %>%
                       filter(region_of_interest != "striatum") %>% 
                       filter(side == names(slices[.x]) & slice == slices[.x])
                     
) %>% 
  dplyr::select(region_of_interest, side, slice, geometry) %>% 
  distinct()

p_fig3b <- df_slices %>% 
  left_join(df_mind_mods, by = join_by(region_of_interest)) %>% 
  arrange(!is.na(module)) %>% 
  
  ggplot() +
  geom_sf(aes(fill = module, geometry = geometry), lwd = 0.5) +
  scale_fill_manual(values = module_colors, na.value = "white") +
  facet_wrap( ~ side + slice, ncol = 2) +
  labs(title = "4 module solution") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = c(0.75, 0.25)
  )

p_fig3b

# SAVE
ggsave(paste0(base_dir, "outputs/figures/3b.brain_modules.pdf"), width = 4, height = 5)
ggsave(paste0(base_dir, "outputs/figures/3b.brain_modules.png"), width = 4, height = 5)


# C: Correlation with Swanson-Sporns  --------------------------------------------------------

p_fig3c <- df_swanson_mind_cor %>% filter(mind_weight < 6 & swanson_weight < 6) %>% 
  
  ggplot(aes(x = mind_weight, y = swanson_weight)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(lty = 2, color = "black") +
  geom_xsidedensity() +
  geom_ysidedensity() +
  stat_cor(label.y = 5.2, color = "blue") +
  theme_ggside_void() +
  labs(x = "CMN bin", y = "Swanson bin", 
       title = "Binned control adult CMN correlation with tract-tracing")

p_fig3c

# SAVE
ggsave(paste0(base_dir, "outputs/figures/3c.tract-tracing_correlation.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/3c.tract-tracing_correlation.png"), width = 5, height = 4)

# S2a: alluvial plot to assess module stability across values of k -------------------------------------

module_cols <- brewer.pal(n, "Paired")
names(module_cols) <- seq(1:n)

p_figS2a <- df_many_models_hclust %>% 
  mutate(module = factor(module, levels = seq(1, n))) %>% 
  ggplot(aes(x = k, 
             stratum = module, 
             alluvium = region_of_interest,
             fill = module, label = module)
  ) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_manual(values = module_cols) +
  scale_x_continuous(breaks = seq(3, n)) +
  labs(title ="HClust module assignments across k values", x = "K") +
  theme(legend.position = "none")

p_figS2a

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S2a.modules_alluvial.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/S2a.modules_alluvial.png"), width = 5, height = 4)


# S2b: Silhouette widths across 1000 permutations -------------------------------------

p_figS2b <- df_sil %>%
  group_by(k) %>% 
  summarise(median_width = median(sil_width),
            sd_width = sd(sil_width),
            error = qnorm(0.975)*sd_width/sqrt(nperm),
            low = median_width - error,
            high = median_width + error) %>% 
  
  ggplot(aes(x = k, y = median_width)) +
  geom_point() +
  geom_line() +
  #ylim(c(0.15, 0.35)) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.3) +
  labs(x = "K", y = "Median silhouette width (+/- 95% CI)",
       title = paste0("CMN median silhouette widths across 1000 permutations"))

p_figS2b

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S2b.modules_silhouette.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/S2b.modules_silhouette.png"), width = 5, height = 4)

# S2c: Swanson-Sporns correlation histograms ----------------------------------------

# HISTOGRAMS
p_figS2c_hist_swanson <- df_swanson_long %>% 
  filter(swanson_weight != 0) %>%  # remove zeros
  count(swanson_weight) %>% 
  ggplot(aes(x = swanson_weight, y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  scale_x_continuous(breaks = seq(1, 7)) +
  labs(x = "tract-tracting weight",
       title = "Tract-tracing network \nweight distribution")

p_figS2c_hist_mind <- df_mind_GM %>% 
  filter(timepoint == 300 & study == "MRC" & R1 != R2) %>%
  group_by(R1, R2) %>% 
  summarise(median_weight = median(weight)) %>% 
  ggplot(aes(x = median_weight), bins = 7) +
  geom_histogram() +
  #stat_bin(aes(label = ..count..), geom = "text", size = 3, angle = 45, vjust = -0.5, hjust = -0.5) +
  labs(title = "Control adult CMN weight \ndistribution", y = "", x = "CMN weight")

p_figS2c_hist_mind_binned <- df_median_mind_bins %>% 
  count(bin) %>% 
  ggplot(aes(x = bin, y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  scale_x_continuous(breaks = seq(1, 7)) +
  labs(x = "CMN bin", y = "",
       title = "Control adult CMN \nbinned weight distribution")

p_figS2c_hist_swanson + p_figS2c_hist_mind + p_figS2c_hist_mind_binned

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S2c.Swanson_MIND_histograms.pdf"), width = 10, height = 4)
ggsave(paste0(base_dir, "outputs/figures/S2c.Swanson_MIND_histograms.png"), width = 10, height = 4)


# S2d: Swanson-Sporns correlation without binning ---------------------------------------------------------------------

p_figS2d <- df_swanson_mind_cor %>% 
  filter(norm_mind < 4) %>% # remove outliers
  ggplot(aes(x = norm_mind, y = norm_swanson)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(lty = 2, color = "black") +
  stat_cor(label.y = 3.2, color = "blue") +
  labs(x = "Normalized CNM weight", 
       y = "Normalized Swanson weight", 
       title = "Control adult CNM correlation with tract-tracing (without binning CMN)")

p_figS2d

# SAVE
ggsave(paste0(base_dir, "outputs/figures/S2d.Swanson_MIND_cor_no_bins.pdf"), width = 5, height = 4)
ggsave(paste0(base_dir, "outputs/figures/S2d.Swanson_MIND_cor_no_bins.png"), width = 5, height = 4)
