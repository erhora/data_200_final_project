#### Running Forward Selection Iteratively
#### Author: Elizabeth Hora


# Loading data
library(tidyverse)
library(tidymodels)
library(foreign)
library(caTools)
library(ROCR) 
library(rsample)
library(Metrics)
library(gt)
library(haven)
library(olsrr)
library(fastDummies)
library(corrr)
library(corrplot)
library(ggcorrplot)
library(tictoc)
library(parallel)
library(foreach)
library(doParallel)
tidymodels::tidymodels_prefer()




# Parallel Processing ----
# Careful with this- if there's not enough memory,
# your computer can crash.
all_cores <- parallel::detectCores(logical = TRUE)
registerDoParallel(cores = all_cores - 1)


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

overall_complete <- complete

important_variables <- c()



# OLS forward regression only works on multilinear regression.
# It does not work with logistic regression.
# I will predict systolic pressure to get around this.
blood_pressure <- read_xpt("data/extra_data/BPX_I.XPT") %>% 
  janitor::clean_names()



# Defining a forward stepwise regression function.
# I will apply this to each file that I downloaded.
forward_predictors <- function(i_file){
  # Global assignment is necessary for ols_step_forward_p()
  file <<- read_xpt(i_file) %>% 
    janitor::clean_names()
  
  bp_data <<- read_xpt("data/extra_data/BPX_I.XPT") %>% 
    janitor::clean_names() %>% 
    select(c(seqn, bpxsy1))
  
  complete <<- merge(bp_data, file) %>% 
    select(-seqn) 
  
  intact_cols <<- complete %>% 
    skimr::skim_without_charts() %>% 
    as_tibble() %>% 
    filter(
      complete_rate > 0.5 &
        numeric.p0 != numeric.p100
    ) %>% 
    select(skim_variable) %>% 
    pull()
  
  
  # There should be at least one predictor variable
  if (length(intact_cols) > 0){
    complete <<- complete %>%
      select(c(intact_cols, bpxsy1)) %>%
      scale() %>%
      as_tibble()
    
    set.seed(10201)
    health_split <<- initial_split(data = complete, prop = 0.8)
    health_train <<- training(health_split)
    health_test <<- testing(health_split)
    
    
    model <<- lm(
      data = health_train,
      bpxsy1 ~ .
    )
    
    
    # Stepwise regression only works for 2 or more predictor variables
    if (ncol(health_train) > 2){
      print(ncol(health_train))
      model <<- model %>%
        ols_step_forward_p(
          progress = FALSE,
          details = FALSE,
          print_plot = FALSE
        )
      
      important_variables <<- append(important_variables, head(model$predictors, 2))
    }
    
    # Simple linear regression has only one variable
    else{
      only_variable <<- names(health_train)
      important_variables <<- append(important_variables, only_variable)
      
      model <<- lm(
        data = health_train,
        bpxsy1 ~ .
      )
      
    }
    
    return(model)
  }
  
}





# Looping Through Datasets ------------------------------------------------
# Going through my files
for(i in seq(length(downloaded_data))){
  print(forward_predictors(downloaded_data[[i]]))
}

important_variables





# Looking at Medication Alone: Counts ----
# We need to decide how we will handle people on lots of medications
# in our models. There is no shame in leaving them in, we just need to
# include that in the model to account for that. 
bp <- read_xpt("data/extra_data/BPX_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(seqn, bpxsy1, bpxdi1))


medication <- read_xpt("data/extra_data/RXQ_RX_I.XPT") %>% 
  janitor::clean_names() %>% 
  glimpse()



# What if I take the number of medications into account?
# Would that explain some variation?
medication_counts <- medication %>% 
  group_by(seqn) %>% 
  summarize(
    count_meds = n() 
  )



# Looking at Medication Alone: Types ----
medication_abridged <- medication %>% 
  filter(rxdrsc1 != 99999)
medication_abridged[medication_abridged == ''] <- "None"

# Lots of medications that perform similar functions 
# are grouped by the same starting letter in their 
# identifying code- e.g., IXXXX means something
# related to heart health.
medication_abridged <- medication_abridged %>%
  mutate(
    section = str_sub(rxdrsc1, end = 1L)
  ) %>%
  select(c(seqn, section)) 

medication_abridged <- medication_abridged %>% 
  dummy_cols() %>% 
  select(-section)


medication_abridged <- medication_abridged %>% 
  group_by(seqn) %>% 
  summarize(
    section_5 = ifelse(sum(section_5) > 1, 1, 0),
    section_7 = ifelse(sum(section_7) > 1, 1, 0),
    section_A = ifelse(sum(section_A) > 1, 1, 0),
    section_B = ifelse(sum(section_B) > 1, 1, 0),
    section_C = ifelse(sum(section_C) > 1, 1, 0),
    section_D = ifelse(sum(section_D) > 1, 1, 0),
    section_E = ifelse(sum(section_E) > 1, 1, 0),
    section_F = ifelse(sum(section_F) > 1, 1, 0),
    section_G = ifelse(sum(section_G) > 1, 1, 0),
    section_H = ifelse(sum(section_H) > 1, 1, 0),
    section_I = ifelse(sum(section_I) > 1, 1, 0),
    section_J = ifelse(sum(section_J) > 1, 1, 0),
    section_K = ifelse(sum(section_K) > 1, 1, 0),
    section_L = ifelse(sum(section_L) > 1, 1, 0),
    section_M = ifelse(sum(section_M) > 1, 1, 0),
    section_N = ifelse(sum(section_N) > 1, 1, 0),
    section_R = ifelse(sum(section_R) > 1, 1, 0),
    section_T = ifelse(sum(section_T) > 1, 1, 0),
    section_Z = ifelse(sum(section_Z) > 1, 1, 0),
  )


# Global assignment is necessary for ols_step_forward_p()
file <- medication_abridged

bp_data <- read_xpt("data/extra_data/BPX_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(seqn, bpxsy1))

complete_meds <- merge(bp_data, file) %>% 
  select(-seqn) 

intact_cols <- complete_meds %>% 
  skimr::skim_without_charts() %>% 
  as_tibble() %>% 
  filter(
    n_missing < 2500
  ) %>% 
  select(skim_variable) %>% 
  pull()

complete_meds <- complete_meds %>% 
  select(c(intact_cols, bpxsy1))



set.seed(10201)
health_split <- initial_split(data = complete_meds, prop = 0.8)
health_train <- training(health_split)
health_test <- testing(health_split)

model <- lm(
  data = health_train,
  bpxsy1 ~ .
) %>% 
  ols_step_forward_p()

meds_help <- head(model$predictors, 3)




# Looking at everything
overall_complete_subset <- overall_complete %>% 
  merge(medication_abridged) %>% 
  merge(medication_counts) %>% 
  select(c(seqn, important_variables, all_of(meds_help), count_meds))


file <- overall_complete_subset

bp_data <- read_xpt("data/extra_data/BPX_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(seqn, bpxsy1))

complete_meds <- merge(bp_data, file) %>% 
  select(-seqn) 

intact_cols <- complete_meds %>% 
  skimr::skim_without_charts() %>% 
  as_tibble() %>% 
  filter(
    n_missing < 2500
  ) %>% 
  select(skim_variable) %>% 
  pull()

complete_meds <- complete_meds %>% 
  select(c(intact_cols, bpxsy1))


set.seed(10201)
health_split <- initial_split(data = complete_meds, prop = 0.8)
health_train <- training(health_split)
health_test <- testing(health_split)

model <- lm(
  data = health_train,
  bpxsy1 ~ .
)

ols_results <- model %>% 
  ols_step_forward_p()


best_predictors <- ols_results$predictors %>% 
  head(10)



best_forward_selection <- overall_complete %>% 
  merge(medication_abridged) %>% 
  merge(medication_counts) %>% 
  select(all_of(best_predictors))

best_forward_selection %>% 
  write_rds("data/best_forward_revised.rds")

best_forward_selection %>% 
  names() %>% 
  as_tibble() %>% 
  rename(optimal_var = value) %>% 
  write_csv("data/optimal_vars_revised.csv")






# Potential correlations --------------------------------------------------
best_predictors <- read_csv("data/optimal_vars_revised.csv")
opt_predictors <- best_predictors %>%
  select(optimal_var) %>%
  pull()

opt_complete <- overall_complete %>%
  merge(medication_abridged) %>%
  merge(medication_counts) %>%
  select(c(opt_predictors, seqn))


