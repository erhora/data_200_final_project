# Loading Packages --------------------------------------------------------
library(tidyverse)
library(haven)
tidymodels::tidymodels_prefer()




# Loading Downloaded Data (Predictors) ------------------------------------
downloaded_data <- dir("data/", ignore.case = TRUE, pattern = "XPT", full.names = TRUE)


# Initiating the process of including everything in one Tibble.
# This is too much to do at once, so the goal is to widdle this down.
complete <- read_xpt("data/extra_data/BPX_I.XPT") %>% 
  janitor::clean_names() %>%
  select(c(seqn, bpxsy1, bpxdi1))

for (i in seq(length(downloaded_data))){
  file <- read_xpt(downloaded_data[[i]]) %>% 
    janitor::clean_names()
  complete <- merge(complete, file)
}




# Loading in Medication- Additional Data Cleaning Required ----------------
medication <- read_xpt("data/extra_data/RXQ_RX_I.XPT") %>%
  janitor::clean_names()

medication_counts <- medication %>%
  group_by(seqn) %>%
  summarize(
    count_meds = n()
  )

medication_abridged <- medication %>%
  filter(rxdrsc1 != 99999)
medication_abridged[medication_abridged == ''] <- "None"

medication_abridged <- medication_abridged %>%
  mutate(
    section = str_sub(rxdrsc1, end = 1L)
  ) %>%
  select(c(seqn, section))

medication_abridged <- medication_abridged %>%
  dummy_cols() %>%
  select(-section)

complete <- complete %>% 
  merge(medication_counts) %>% 
  merge(medication_abridged)




# Finalizing Data ---------------------------------------------------------
# Filtering out those on high blood pressure medication
complete <- complete %>% 
  filter(section_I == 0) %>% 
  glimpse()

write_rds(
  complete,
  file = "data/complete.rds"
)
