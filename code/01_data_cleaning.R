### wog-wog-postfire-ants

### 01_data_cleaning.R

### Matt Bitters
### matthew.bitters@colorado.edu








# ============================================================
#  0. Setup
# ============================================================

### Install and load required packages
packages <- c("here", "dplyr", "ggplot2", "tidyr")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(dplyr)
library(ggplot2)
library(tidyr)


### Create required folders

# Define folders to create
folders <- c(
  "data",
  "data/raw",
  "data/derived",
  "model_outputs",
  "figures"
)

# Create folders if they don't already exist
for (f in folders) {
  dir.create(here(f), showWarnings = FALSE, recursive = TRUE)
}



# ============================================================
#  1. Read in data and factor
# ============================================================


# Read CSV
ant_data <- read.csv(here("data", "raw", "ants_all_data_10-31-2025.csv"))

# Examine data
str(ant_data)
head(ant_data)

# Factor and set levels
ant_data$bio_year <- factor(ant_data$bio_year)                 # 1-35 (with some missing years)
ant_data$treat <- factor(ant_data$treat, levels=c(2,1,3))      # 2=controls, 1=fragments, 3=matrix
ant_data$size <- factor(ant_data$size,levels=c(4, 1, 2, 3))    # 4=controls, 1=small, 2=medium, 3=large, NA=matrix
ant_data$edge <- factor(ant_data$edge, levels=c(3, 1,2))       # 3=controls, 1=core, 2=edge, NA=matrix
ant_data$site <- factor(ant_data$site)                         # 1-72=fragments, 73-120=controls, 121-144=fragments, 145-188=matrix
ant_data$rep <- factor(ant_data$rep)                           # 1-6
ant_data$patch <- factor(ant_data$patch)                       # 1-30
ant_data$topo <- factor(ant_data$topo)                         # 1=drain, 2=slope

# Confirm factoring
str(ant_data)



# ============================================================
#  2. Visualize time series
# ============================================================


# Summarize abundance by species, year, and treatment
ant_summary <- ant_data %>%
  group_by(species, bio_year, treat) %>%
  summarize(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")

# Assign treatment labels without changing factor order
ant_summary <- ant_summary %>%
  mutate(treat_label = recode(treat,
                              "2" = "Control",
                              "1" = "Fragment",
                              "3" = "Matrix"
  ))

# Define all possible years
all_years <- 1:35

# Complete the data frame so every species × treat × bio_year combo exists
ant_summary_full <- ant_summary %>%
  mutate(bio_year = as.numeric(as.character(bio_year))) %>%  # ensure numeric
  complete(
    species,
    treat_label,
    bio_year = all_years
  ) %>%
  arrange(species, treat_label, bio_year)

# Plot
ggplot(ant_summary_full, aes(
  x = bio_year,
  y = mean_abundance,
  color = treat_label,
  group = treat_label
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ species, scales = "free_y") +
  scale_color_manual(
    values = c("Control" = "#1b9e77", "Fragment" = "#d95f02", "Matrix" = "#7570b3"),
    name = "Treatment"
  ) +
  scale_x_continuous(
    breaks = seq(1, 35, by = 2),  # tick every 2 years for readability
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Biological year",
    y = "Mean abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95", color = NA)
  )



# ============================================================
#  3. Reform dataset for analysis and plotting
# ============================================================


# Compute pre-fragmentation abundance (sum years 1-2)
pre_frag <- ant_data %>%
  filter(bio_year %in% c("1", "2")) %>%
  group_by(species, site) %>%
  summarize(yrs1_2_abund = sum(abundance, na.rm = TRUE), .groups = "drop")

# Create post-fragment dataset (years >= 3)
post_frag <- ant_data %>%
  filter(as.numeric(as.character(bio_year)) >= 3) %>%
  left_join(pre_frag, by = c("species", "site")) %>%
  mutate(
    presence = ifelse(abundance > 0, 1, 0),
    yrs1_2_abund_log = log(yrs1_2_abund + 1),
    time_period = ifelse(as.numeric(as.character(bio_year)) < 35, "post-frag", "post-fire")
  ) %>%
  arrange(species, bio_year, site)

# Save clean dataset
write.csv(post_frag, here("data", "derived", "ants_model_ready_10-31-2025.csv"),
          row.names = FALSE)
message("✅ Saved ants_model_ready_10-31-2025.csv")
