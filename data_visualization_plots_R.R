#### Visualizing our Literature Review-Identified Potential Predictor Variables
#### for hypertension
#### Author: Elizabeth Hora



library(tidyverse)
library(haven)

theme_set(theme_minimal())

bp_data <- read_xpt("data/extra/github_setup/data/extra_data/BPX_I.XPT") %>% 
  janitor::clean_names()




# Jusitifying Binary Classification ---------------------------------------
bp_data %>% 
  filter(
    bpxdi1 >= 20
  ) %>% 
  mutate(
    bp_status = as.factor(case_when(
      bpxsy1 >= 140 | bpxdi1 >= 90 ~ "High",
      bpxsy1 >= 130 | bpxdi1 >= 80 ~ "Risk of\nHigh",
      bpxsy1 < 130 | bpxdi1 < 80 ~ "Regular"
    )
  )
) %>% 
  ggplot(aes(x = bpxdi1, y = bpxsy1, color = bp_status)) +
  geom_point() +
  labs(
    x = "Diastolic (mmHg)",
    y = "Systolic (mmHg)",
    color = "Blood Pressure\nStatus"
  )
# There is no satisfactory boundary here, and adding in a 
# risk of hypertension group appears arbitrary. At this rate,
# we could make another category for those who could be 
# at the greatest potential to return to regular blood pressure.
# The arbitrary boundaries would be difficult to manage.




# Reading in Optimal Data -------------------------------------------------
best_vars <- read_csv("data/extra/github_setup/optimal_vars_revised.csv") %>% 
  select(optimal_var) %>% 
  pull()

complete <- read_rds("data/extra/github_setup/best_forward_revised.rds")

complete %>% 
  left_join(bp_data %>% select(c(seqn, bpxsy1, bpxdi1))) %>% 
  filter(
    !is.na(bpxsy1),
    !is.na(bpxdi1)
  )

# Size of sample
complete %>% 
  left_join(bp_data %>% select(c(seqn, bpxsy1, bpxdi1))) %>% 
  na.omit() %>% 
  dim()




# Blood Pressure versus Age -----------------------------------------------
#### Figure 1 of Data Visualization
regression_start <- complete %>% 
  left_join(bp_data %>% 
              select(
                c(seqn, bpxsy1, bpxdi1)
              )
  ) %>% 
  na.omit()



male_coef <- lm(
  data = regression_start %>% 
    filter(
      riagendr == 1
    ),
  formula = bpxsy1 ~ ridageyr
)$coefficients  

male_eq <- paste0("Systolic Pressure ~ ", round(male_coef[1], 2), " + ", round(male_coef[2], 2), "*Age")

female_coef <- lm(
  data = regression_start %>% 
    filter(
      riagendr == 2
    ),
  formula = bpxsy1 ~ ridageyr
)$coefficients  

female_eq <- paste0("Systolic Pressure ~ ", round(female_coef[1], 2), " + ", round(female_coef[2], 2), "*Age")


regression_eqs <- tribble(
  ~"riagendr", ~"regress_eq",
  "Male", male_eq,
  "Female", female_eq
)

regression_eqs

regression_plot <- regression_start %>% 
  mutate(
    riagendr = factor(case_when(
      riagendr == 1 ~ "Male",
      riagendr == 2 ~ "Female"
    ),
    levels = c("Male", "Female")
    ),
    regeq = ifelse(riagendr == "Male", male_eq, female_eq)
  ) %>% 
  ggplot(aes(
    x = ridageyr,
    y = bpxsy1
    ),
    size = 2
  ) + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~riagendr, ncol = 2) +
  geom_smooth(
    formula = y~x,
    se = FALSE,
    color = "red"
    ) +
  labs(
    x = "Age (Years)",
    y = "Systolic Blood Pressure (mmHg)",
    color = NULL
  ) +
  theme(
    legend.position = c(0,1),
    legend.justification = c(0,1),
    strip.text = element_text(size = 26, family = "serif"),
    axis.title = element_text(size = 30, family = "serif"),
    legend.text = element_text(family = "serif", size = 18),
    legend.title = element_text(family = "serif", size = 18),
    plot.caption = element_text(family = "serif", size = 17),
    axis.text = element_text(family = "serif", size = 21)
  ) +
  geom_text(
    data = regression_eqs,
    aes(x = 30, y = 230, label = regress_eq),
    size = 6,
    family = "serif"
  )
regression_plot


ggsave(
  regression_plot,
  filename = "saved_plots/ht_by_age_gender_eq.png",
  height = 6,
  width = 16,
  dpi = 1200
)




# Exploring Important Predictor Variables ---------------------------------
#### Figure 2 of Data Visualization
alcohol <- read_xpt("data/ALQ_I.XPT") %>% 
  janitor::clean_names()


drink_labels = c("0 or 1", seq(2, 9), "10+")


# Number of People in Total:
complete %>% 
  left_join(
    alcohol %>% 
      select(
        c(seqn, alq130)
      )
  ) %>% 
  left_join(
    bp_data %>% 
      select(
        c(seqn, bpxsy1, bpxdi1)
      )
  ) %>% 
  filter(
    !is.na(bpxsy1),
    !is.na(bpxdi1),
    !is.na(alq130) & 
      alq130 < 100
  ) %>% 
  dim()
# 3241 people


# Female-Male Breakdown:
gender_alq <- complete %>% 
  left_join(
    alcohol %>% 
      select(
        c(seqn, alq130)
      )
  ) %>% 
  left_join(
    bp_data %>% 
      select(
        c(seqn, bpxsy1, bpxdi1)
      )
  ) %>% 
  filter(
    !is.na(bpxsy1),
    !is.na(bpxdi1),
    !is.na(alq130) & 
      alq130 < 100,
  ) %>% 
  mutate(
    riagendr = as.factor(case_when(
      riagendr == 1 ~ "Male",
      riagendr == 2 ~ "Female"
    )
    )
  ) %>% 
  group_by(riagendr) %>% 
  summarize(
    n_by_gender = n()
  )
# 1511 women
# 1730 men
female_size <- gender_alq %>% 
  filter(riagendr == "Female") %>% 
  select(n_by_gender) %>% 
  pull()

male_size <- gender_alq %>% 
  filter(riagendr == "Male") %>% 
  select(n_by_gender) %>% 
  pull()


facet_label_f <- paste0("Female (n = ", female_size, ")")
facet_label_m <- paste0("Male (n = ", male_size, ")")


alq_by_gender <- complete %>% 
  left_join(alcohol %>% 
              select(c(seqn, alq130))
  ) %>% 
  left_join(bp_data %>% 
              select(c(seqn, bpxsy1, bpxdi1))
  ) %>% 
  filter(
    !is.na(bpxsy1),
    !is.na(bpxdi1),
    !is.na(alq130) & 
      alq130 < 100
  ) %>% 
  mutate(
    bp_status = case_when(
      bpxsy1 >= 140 | bpxdi1 >= 90 ~ 1,
      bpxsy1 < 140 & bpxdi1 < 90 ~ 0 
    ),
    alq130 = case_when(
      floor(alq130/10) == 1 ~ 10,
      TRUE ~ alq130
    )
  ) %>% 
  group_by(riagendr, alq130) %>% 
  summarize(
    regular = n() - sum(bp_status),
    high = sum(bp_status)
  ) %>% 
  ungroup() %>% 
  mutate(
    avg_frac_high = high / (high + regular),
    avg_frac_reg = regular / (high + regular),
    total = high + regular,
    riagendr = as.factor(case_when(
      riagendr == 1 ~ facet_label_m,
      riagendr == 2 ~ facet_label_f
    )
    )
  ) %>% 
  pivot_longer(
    cols = c(regular, high),
    values_to = "count",
    names_to = "count_type"
  ) %>% 
  mutate(
    count_type = factor(case_when(
      count_type == "high" ~ "High",
      count_type == "regular" ~ "Regular"
    ),
    levels = c("Regular", "High"))
  ) %>%
  ggplot(aes(x = alq130, y = count, fill = count_type)) +
  geom_col(position = "dodge") +
  labs(
    x = "Average Number of Alcoholic Drinks / Day",
    y = "Number of Individuals",
    fill = "Hypertension\nStatus",
    caption = "The number of individuals who drink within a certain category per gender are labeled."
  ) +
  scale_x_continuous(
    limits = c(0.5, 10.5),
    breaks = seq(1, 10),
    labels = drink_labels
  ) + 
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1,1),
    strip.text = element_text(size = 30, family = "serif"),
    strip.background = element_rect(fill = "grey90"),
    axis.title = element_text(size = 36, family = "serif"),
    legend.text = element_text(family = "serif", size = 28),
    legend.title = element_text(family = "serif", size = 30),
    plot.caption = element_text(family = "serif", size = 27),
    axis.text = element_text(family = "serif", size = 40),
    panel.grid.major.x = element_line(size = 1.5)
  ) +
  facet_wrap(
    ~riagendr, 
    ncol = 1
    ) +
  geom_text(
    aes(label = count, group = count_type),
    position = position_dodge(width = 0.9),
    vjust = -0.25,
    color = "black",
    family = "serif",
    size = 12
  ) +
  scale_fill_manual(
    values = c("#d9c3eb", "#e07070")
  ) +
  scale_y_continuous(
    limits = c(0, 600)
  )


alq_by_gender


ggsave(
  alq_by_gender,
  filename = "saved_plots/alcohol_by_gender.png",
  height = 12,
  width = 20,
  dpi = 1200
)




alq_by_gender <- complete %>% 
  left_join(alcohol %>% 
              select(c(seqn, alq130))
  ) %>% 
  left_join(bp_data %>% 
              select(c(seqn, bpxsy1, bpxdi1))
  ) %>% 
  filter(
    !is.na(bpxsy1),
    !is.na(bpxdi1),
    !is.na(alq130) & 
      alq130 < 100
  ) %>% 
  mutate(
    bp_status = case_when(
      bpxsy1 >= 140 | bpxdi1 >= 90 ~ 1,
      bpxsy1 < 140 & bpxdi1 < 90 ~ 0 
    ),
    alq130 = case_when(
      floor(alq130/10) == 1 ~ 10,
      TRUE ~ alq130
    )
  ) %>% 
  group_by(riagendr, alq130) %>% 
  summarize(
    regular = n() - sum(bp_status),
    high = sum(bp_status)
  ) %>% 
  ungroup() %>% 
  mutate(
    avg_frac_high = round(high / (high + regular), 2),
    avg_frac_reg = round(regular / (high + regular),2),
    total = high + regular,
    riagendr = as.factor(case_when(
      riagendr == 1 ~ facet_label_m,
      riagendr == 2 ~ facet_label_f
    )
    )
  ) %>% 
  pivot_longer(
    cols = c(avg_frac_high, avg_frac_reg),
    values_to = "frac",
    names_to = "frac_type"
  ) %>% 
  mutate(
    frac_type = factor(case_when(
      frac_type == "avg_frac_high" ~ "High",
      frac_type == "avg_frac_reg" ~ "Regular"
    ),
    levels = c("Regular", "High"))
  ) %>%
  ggplot(aes(x = alq130, y = frac, fill = frac_type)) +
  geom_col(position = "dodge") +
  labs(
    x = "Average Number of Alcoholic Drinks / Day",
    y = "Percentage of Individuals",
    fill = "Hypertension\nStatus",
    caption = "The number of individuals who drink within a certain category per gender are labeled."
  ) +
  scale_x_continuous(
    limits = c(0.5, 10.5),
    breaks = seq(1, 10),
    labels = drink_labels
  ) + 
  theme(
    strip.text = element_text(size = 30, family = "serif"),
    strip.background = element_rect(fill = "grey90"),
    axis.title = element_text(size = 36, family = "serif"),
    legend.text = element_text(family = "serif", size = 28),
    legend.title = element_text(family = "serif", size = 30),
    plot.caption = element_text(family = "serif", size = 27),
    axis.text = element_text(family = "serif", size = 40),
    panel.grid.major.x = element_line(size = 1.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  facet_wrap(
    ~riagendr, 
    ncol = 1
  ) +
  geom_text(
    aes(label = scales::percent(frac), group = frac_type),
    position = position_dodge(width = 0.9),
    vjust = -0.25,
    color = "black",
    family = "serif",
    size = 10
  ) +
  scale_fill_manual(
    values = c("#d9c3eb", "#e07070")
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent(seq(0, 1, 0.25))
  )


alq_by_gender


ggsave(
  alq_by_gender,
  filename = "saved_plots/alcohol_by_gender_frac.png",
  height = 12,
  width = 20,
  dpi = 1200
)
