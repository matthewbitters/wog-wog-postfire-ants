### wog-wog-postfire-ants

### 03_ants_main_treats_figures.R

### Matt Bitters
### matthew.bitters@colorado.edu








# ============================================================
#  0. Setup
# ============================================================

### Install and load required packages
packages <- c("here", "dplyr", "tidyr", "brms", "emmeans", "patchwork", "forcats", "tidybayes")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(dplyr)
library(tidyr)
library(brms)
library(emmeans)
library(patchwork)
library(forcats)
library(tidybayes)


# ============================================================
#  1. Read in data
# ============================================================

# Aphan models
m4_aphan <- readRDS(here("model_outputs", "m4_aphan.rds"))
m5_years_aphan <- readRDS(here("model_outputs", "m5_years_aphan.rds"))

# Lepto models
m4_lepto <- readRDS(here("model_outputs", "m4_lepto.rds"))
m5_years_lepto <- readRDS(here("model_outputs", "m5_years_lepto.rds"))

# Summaries
summary(m4_aphan)
summary(m5_years_aphan)
summary(m4_lepto)
summary(m5_years_lepto)

# Read in main data
ant_mod_data <- read.csv(here("data", "derived", "ants_model_ready_10-31-2025.csv"))
head(ant_mod_data)
str(ant_mod_data)




# ============================================================
#  2. Main treatment effect size figure
# ============================================================


### Top panel

#  Aphaenogaster longiceps

emm_a <- emmeans(m4_aphan, ~ treat * time_period, type = "link", re_formula = NA)
emm_a_df <- as.data.frame(emm_a) %>%
  mutate(SE = (upper.HPD - lower.HPD) / (2 * 1.96))

emm_a_summ <- emm_a_df %>%
  group_by(time_period, treat) %>%
  summarize(emmean = median(emmean, na.rm = TRUE),
            SE = median(SE, na.rm = TRUE),
            .groups = "drop")

emm_a_wide <- emm_a_summ %>%
  pivot_wider(names_from = treat, values_from = c(emmean, SE))

top_a <- bind_rows(
  emm_a_wide %>%
    transmute(time_period,
              contrast = paste0("Fragments (", time_period, ")"),
              estimate = emmean_1 - emmean_2,
              se = sqrt(SE_1^2 + SE_2^2)),
  emm_a_wide %>%
    transmute(time_period,
              contrast = paste0("Matrix (", time_period, ")"),
              estimate = emmean_3 - emmean_2,
              se = sqrt(SE_3^2 + SE_2^2))
) %>%
  mutate(lower = estimate - 1.96 * se,
         upper = estimate + 1.96 * se,
         species = "A. longiceps")


#  Leptomyrmex erythrocephalus

emm_l <- emmeans(m4_lepto, ~ treat * time_period, type = "link", re_formula = NA)
emm_l_df <- as.data.frame(emm_l) %>%
  mutate(SE = (upper.HPD - lower.HPD) / (2 * 1.96))

emm_l_summ <- emm_l_df %>%
  group_by(time_period, treat) %>%
  summarize(emmean = median(emmean, na.rm = TRUE),
            SE = median(SE, na.rm = TRUE),
            .groups = "drop")

emm_l_wide <- emm_l_summ %>%
  pivot_wider(names_from = treat, values_from = c(emmean, SE))

top_l <- bind_rows(
  emm_l_wide %>%
    transmute(time_period,
              contrast = paste0("Fragments (", time_period, ")"),
              estimate = emmean_1 - emmean_2,
              se = sqrt(SE_1^2 + SE_2^2)),
  emm_l_wide %>%
    transmute(time_period,
              contrast = paste0("Matrix (", time_period, ")"),
              estimate = emmean_3 - emmean_2,
              se = sqrt(SE_3^2 + SE_2^2))
) %>%
  mutate(lower = estimate - 1.96 * se,
         upper = estimate + 1.96 * se,
         species = "L. erythrocephalus")


#  Combine both species
top_combined <- bind_rows(top_a, top_l)

# Define consistent order for y-axis
contrast_order <- c(
  "Fragments (post-frag)",
  "Matrix (post-frag)",
  "Fragments (post-fire)",
  "Matrix (post-fire)"
)

top_combined <- top_combined %>%
  mutate(contrast = factor(contrast, levels = contrast_order))


#  Faceted plot
p_top <- ggplot(top_combined,
                aes(x = estimate,
                    y = fct_rev(contrast),
                    xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray70") +
  geom_hline(yintercept = 2.5, linewidth = 1, color = "gray70") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                width = 0.15,
                linewidth = 1,
                color = "black") +
  geom_point(size = 3) +
  facet_wrap(~species, nrow = 1) +
  labs(x = "Effect size (log-odds): frag/matr compared to controls", y = NULL) +
  scale_y_discrete(labels = c(
    "Fragments (post-frag)" = "Fragments (pre-fire)",
    "Matrix (post-frag)"   = "Matrix (pre-fire)",
    "Fragments (post-fire)" = "Fragments (post-fire)",
    "Matrix (post-fire)"   = "Matrix (post-fire)"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),      # remove all grid lines
    axis.text.y = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
  
p_top


### Bottom panel

# Summarize emmeans. Ensure SEs are included
emm_a_df <- as.data.frame(emm_a) %>%
  mutate(SE = (upper.HPD - lower.HPD) / (2 * 1.96))

emm_l_df <- as.data.frame(emm_l) %>%
  mutate(SE = (upper.HPD - lower.HPD) / (2 * 1.96))

# Summarize to one row per treat x time period
emm_a_summ <- emm_a_df %>%
  group_by(time_period, treat) %>%
  summarize(emmean = median(emmean, na.rm = TRUE),
            SE = median(SE, na.rm = TRUE),
            .groups = "drop")

emm_l_summ <- emm_l_df %>%
  group_by(time_period, treat) %>%
  summarize(emmean = median(emmean, na.rm = TRUE),
            SE = median(SE, na.rm = TRUE),
            .groups = "drop")

# Compute post-pre contrasts within each treatment
make_post_pre <- function(emm_summ, species_name) {
  emm_wide <- emm_summ %>%
    pivot_wider(names_from = time_period, values_from = c(emmean, SE))
  
  # Calculate post - pre for each treatment
  df <- emm_wide %>%
    transmute(
      treat = case_when(
        treat == "1" ~ "Fragments",
        treat == "2" ~ "Control",
        treat == "3" ~ "Matrix",
        TRUE ~ treat
      ),
      estimate = `emmean_post-fire` - `emmean_post-frag`,  # adjust names if your time_period labels differ
      se = sqrt(`SE_post-fire`^2 + `SE_post-frag`^2),
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      species = species_name
    )
  
  return(df)
}

bottom_a <- make_post_pre(emm_a_summ, "A. longiceps")
bottom_l <- make_post_pre(emm_l_summ, "L. erythrocephalus")

# Combine data and order 7-axis
bottom_combined <- bind_rows(bottom_a, bottom_l)

bottom_combined <- bottom_combined %>%
  mutate(
    treat = factor(treat, levels = c("Control", "Fragments", "Matrix")),
    species = factor(species, levels = c("A. longiceps", "L. erythrocephalus"))
  )

# Plot
p_bottom <- ggplot(bottom_combined,
                   aes(x = estimate,
                       y = fct_rev(treat),
                       xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray70") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                width = 0.15,
                linewidth = 1,
                color = "black") +
  geom_point(size = 3) +
  facet_wrap(~species, nrow = 1) +
  labs(
    x = "Effect size (log-odds): post-fire compared to pre-fire",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),      # remove all grid lines
    axis.text.y = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_bottom


### Combine top and bottom panels
main_treat_es_fig <- (p_top / p_bottom)

main_treat_es_fig

# Save figure
ggsave("figures/main_treat_es_fig.png", plot = main_treat_es_fig, width = 6, height = 5, dpi = 300, bg = "white")



























###############################################################################



#### Probability (instead of log odds) - same figure as above
#### ***Not complete/final



# ============================================================
# Top panel: effect vs control
# ============================================================

# ============================================================
# 1️⃣ Prepare emmeans summaries
# ============================================================

# Aphaenogaster longiceps
emm_a <- emmeans(m4_aphan, ~ treat * time_period, type = "link", re_formula = NA)
emm_a_df <- as.data.frame(emm_a) %>%
  mutate(SE = (upper.HPD - lower.HPD) / (2 * 1.96))

# Leptomyrmex erythrocephalus
emm_l <- emmeans(m4_lepto, ~ treat * time_period, type = "link", re_formula = NA)
emm_l_df <- as.data.frame(emm_l) %>%
  mutate(SE = (upper.HPD - lower.HPD) / (2 * 1.96))

# Long-format summaries (for post-pre comparisons)
emm_a_long <- emm_a_df %>%
  group_by(time_period, treat) %>%
  summarize(emmean = median(emmean, na.rm = TRUE),
            SE = median(SE, na.rm = TRUE),
            .groups = "drop")

emm_l_long <- emm_l_df %>%
  group_by(time_period, treat) %>%
  summarize(emmean = median(emmean, na.rm = TRUE),
            SE = median(SE, na.rm = TRUE),
            .groups = "drop")

# Wide-format summaries (for top panel: vs control)
emm_a_wide <- emm_a_long %>%
  pivot_wider(names_from = treat, values_from = c(emmean, SE))

emm_l_wide <- emm_l_long %>%
  pivot_wider(names_from = treat, values_from = c(emmean, SE))

# ============================================================
# 2️⃣ Top panel: Fragments/Matrix vs Control
# ============================================================

top_a <- bind_rows(
  emm_a_wide %>%
    transmute(time_period,
              contrast = paste0("Fragments (", time_period, ")"),
              estimate = emmean_1 - emmean_2,
              se = sqrt(SE_1^2 + SE_2^2)),
  emm_a_wide %>%
    transmute(time_period,
              contrast = paste0("Matrix (", time_period, ")"),
              estimate = emmean_3 - emmean_2,
              se = sqrt(SE_3^2 + SE_2^2))
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    species = "A. longiceps",
    # convert to probability difference
    estimate_prob = 1 / (1 + exp(-estimate)) - 0.5,
    lower_prob    = 1 / (1 + exp(-lower)) - 0.5,
    upper_prob    = 1 / (1 + exp(-upper)) - 0.5
  )

top_l <- bind_rows(
  emm_l_wide %>%
    transmute(time_period,
              contrast = paste0("Fragments (", time_period, ")"),
              estimate = emmean_1 - emmean_2,
              se = sqrt(SE_1^2 + SE_2^2)),
  emm_l_wide %>%
    transmute(time_period,
              contrast = paste0("Matrix (", time_period, ")"),
              estimate = emmean_3 - emmean_2,
              se = sqrt(SE_3^2 + SE_2^2))
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    species = "L. erythrocephalus",
    estimate_prob = 1 / (1 + exp(-estimate)) - 0.5,
    lower_prob    = 1 / (1 + exp(-lower)) - 0.5,
    upper_prob    = 1 / (1 + exp(-upper)) - 0.5
  )

top_combined <- bind_rows(top_a, top_l)

# Order y-axis
contrast_order <- c(
  "Fragments (post-frag)",
  "Matrix (post-frag)",
  "Fragments (post-fire)",
  "Matrix (post-fire)"
)

top_combined <- top_combined %>%
  mutate(contrast = factor(contrast, levels = contrast_order))

# Plot top panel
p_top <- ggplot(top_combined,
                aes(x = estimate_prob,
                    y = fct_rev(contrast),
                    xmin = lower_prob, xmax = upper_prob)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray70") +
  geom_errorbar(aes(xmin = lower_prob, xmax = upper_prob),
                width = 0.15, linewidth = 1.2, color = "black") +
  geom_point(size = 4) +
  facet_wrap(~species, nrow = 1) +
  labs(x = "Probability difference relative to control",
       y = NULL) +
  scale_y_discrete(labels = c(
    "Fragments (post-frag)" = "Fragments (pre-fire)",
    "Matrix (post-frag)"   = "Matrix (pre-fire)",
    "Fragments (post-fire)" = "Fragments (post-fire)",
    "Matrix (post-fire)"   = "Matrix (post-fire)"
  )) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 11))

# ============================================================
# 3️⃣ Bottom panel: post-fire vs pre-fire
# ============================================================

make_post_pre <- function(emm_summ, species_name) {
  df <- emm_summ %>%
    pivot_wider(names_from = time_period, values_from = c(emmean, SE)) %>%
    transmute(
      treat = case_when(
        treat == "1" ~ "Fragments",
        treat == "2" ~ "Control",
        treat == "3" ~ "Matrix",
        TRUE ~ treat
      ),
      estimate = `emmean_post-fire` - `emmean_post-frag`,
      se = sqrt(`SE_post-fire`^2 + `SE_post-frag`^2),
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      species = species_name
    ) %>%
    mutate(
      estimate_prob = 1 / (1 + exp(-estimate)) - 0.5,
      lower_prob    = 1 / (1 + exp(-lower)) - 0.5,
      upper_prob    = 1 / (1 + exp(-upper)) - 0.5
    )
  
  return(df)
}

bottom_a <- make_post_pre(emm_a_long, "A. longiceps")
bottom_l <- make_post_pre(emm_l_long, "L. erythrocephalus")

bottom_combined <- bind_rows(bottom_a, bottom_l) %>%
  mutate(
    treat = factor(treat, levels = c("Control", "Fragments", "Matrix")),
    species = factor(species, levels = c("A. longiceps", "L. erythrocephalus"))
  )

p_bottom <- ggplot(bottom_combined,
                   aes(x = estimate_prob,
                       y = fct_rev(treat),
                       xmin = lower_prob, xmax = upper_prob)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray70") +
  geom_errorbar(aes(xmin = lower_prob, xmax = upper_prob),
                width = 0.15, linewidth = 1.2, color = "black") +
  geom_point(size = 4) +
  facet_wrap(~species, nrow = 1) +
  labs(x = "Probability difference: post-fire vs pre-fire",
       y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 11))

# ============================================================
# 4️⃣ Combine panels and save
# ============================================================

main_prob_fig <- (p_top / p_bottom) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

main_prob_fig

ggsave("figure_effect_sizes_prob.png",
       plot = main_prob_fig,
       width = 10, height = 8, dpi = 300, bg = "white")














# ============================================================
#  3. Main treatment effects time series figure
# ============================================================


# ============================================================
# 1️⃣ Define the range of years you want predictions for
# ============================================================
all_years <- 3:35

# ============================================================
# 2️⃣ Generate predicted probabilities for each species
# ============================================================
# Aphaenogaster longiceps
pred_a <- emmeans(
  m5_years_aphan,
  ~ treat * bio_year,
  at = list(bio_year = all_years),   # specify all years explicitly
  type = "response"
) %>%
  as.data.frame() %>%
  mutate(species = "A. longiceps")

# Leptomyrmex erythrocephalus
pred_l <- emmeans(
  m5_years_lepto,
  ~ treat * bio_year,
  at = list(bio_year = all_years),
  type = "response"
) %>%
  as.data.frame() %>%
  mutate(species = "L. erythrocephalus")

# Combine species predictions
pred_combined <- bind_rows(pred_a, pred_l)

# ============================================================
# 3️⃣ Prepare observed years per species and treatment
# ============================================================
observed_years <- data.frame(
  species = rep(c("A. longiceps", "L. erythrocephalus"), each = 3*length(unique(ant_mod_data$bio_year))),
  treat = rep(1:3, times = 2*length(unique(ant_mod_data$bio_year))),
  bio_year = rep(unique(ant_mod_data$bio_year), times = 6)
) %>%
  # Convert treat to factor with consistent labels
  mutate(
    treat = factor(treat, levels = c(1,2,3), labels = c("Control","Fragments","Matrix")),
    observed = TRUE       # helper column to identify real observations
  )

# ============================================================
# 4️⃣ Mark extrapolated years
# ============================================================
plot_data <- pred_combined %>%
  mutate(
    treat = factor(treat, levels = c(1,2,3), labels = c("Control","Fragments","Matrix"))
  ) %>%
  left_join(observed_years, by = c("species", "treat", "bio_year")) %>%
  mutate(
    extrapolated = is.na(observed)  # TRUE if year not observed
  ) %>%
  select(species, treat, bio_year, response, lower.HPD, upper.HPD, extrapolated)

# Add a grouping variable combining treat and extrapolated
plot_data <- plot_data %>%
  mutate(line_group = paste(treat, extrapolated, sep = "_"))

plot_data2 <- plot_data %>%
  arrange(species, treat, bio_year) %>%
  group_by(species, treat) %>%
  mutate(
    line_group = cumsum(extrapolated != lag(extrapolated, default = FALSE))  # separate continuous segments
  ) %>%
  ungroup()

# ============================================================
# 4️⃣ Plot faceted time series
# ============================================================


ggplot(plot_data2, aes(x = bio_year, y = response, color = treat, group = line_group)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD, fill = treat), alpha = 0.2, color = NA) +
  geom_line(aes(linetype = extrapolated), size = 1.1) +
  geom_point(data = subset(plot_data2, !extrapolated), size = 2) +
  facet_wrap(~species, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Control" = "#1b9e77", "Fragments" = "#d95f02", "Matrix" = "#7570b3")) +
  scale_fill_manual(values = c("Control" = "#1b9e77", "Fragments" = "#d95f02", "Matrix" = "#7570b3")) +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed")) +
  labs(
    x = "Year",
    y = "Probability of occurrence",
    color = "Treatment",
    fill = "Treatment",
    linetype = "Extrapolated"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())


# -------------------------------
# 5️⃣ Save figure
# -------------------------------
ggsave("figure_time_series_prob.png", plot = time_series_plot,
       width = 10, height = 8, dpi = 300, bg = "white")





