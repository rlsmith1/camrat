
###
### Generate figure to demonstrate sample size in cohorts
###


# libraries ---------------------------------------------------------------

library(tidyverse)
library(tidytext)

# set plot theme ----------------------------------------------------------

theme_set(theme_light() +
            theme(plot.title = element_text(size = 14),
                  axis.title = element_text(size = 14),
                  axis.text = element_text(size = 14),
                  strip.text = element_text(size = 14, color = "black"),
                  strip.background = element_rect(fill = "white", color = "gray"),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 14)
                  #legend.key.width = unit(2, "cm")
            )
)

group_cols <- c("#FFA500", "#4682B4")
names(group_cols) <- c("control", "MS")

# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

load(paste0(base_dir, "objects/25Aug2023_df_data.RDS")) # df_data


# table of n --------------------------------------------------------------

df_data %>% 
  dplyr::select(study, timepoint, sex, group, subject) %>% 
  distinct() %>% 
  count(study, timepoint, sex, group) %>% 
  
  write.csv(paste0(base_dir, "outputs/tables/28Aug2023_table_of_n.csv"), 
            row.names = FALSE)


# figure ------------------------------------------------------------------

labels <- c("Developmental cohort", "Experimental stress cohort")
names(labels) <- c("MRC", "GSK")

df_data %>% 
  dplyr::select(study, timepoint, age, sex, group, subject) %>% 
  distinct() %>% 
  mutate(study = factor(study, levels = c("MRC", "GSK"))) %>% 
  arrange(study, group, -age) %>% 
  mutate(subject = factor(subject, levels = unique(.$subject))) %>% 
  
  ggplot(aes(x = age, y = subject, color = group)) +
  geom_point(size = 2) +
  geom_line(aes(group = subject), lty = 2, alpha = 0.5) +
  scale_color_manual(values = group_cols) +
  scale_x_continuous(breaks = c(20, 35, 63, 230, 290)) +
  labs(title = "Scan numbers and ages for each subject") +
  facet_wrap(~ study, scales = "free_y", nrow = 2, labeller = as_labeller(labels)) +
  theme(axis.text.y = element_blank(),
        legend.position = "bottom"
  )

ggsave(paste0(base_dir, "outputs/figures/0.0.sample_size.pdf"),
       height = 8, width = 6)
ggsave(paste0(base_dir, "outputs/figures/0.0.sample_size.png"),
       height = 8, width = 6)
