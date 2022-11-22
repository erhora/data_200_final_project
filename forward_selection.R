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
library(tictoc)
library(parallel)
library(foreach)
library(doParallel)
tidymodels::tidymodels_prefer()



# Potassium Intake --------------------------------------------------------
# Isaac wrote a paragraph about this, so I will include a simple predictor
diet <- read_xpt("data/DR1IFF_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(
    seqn, dr1ipota, 
    dr1isodi, dr1itfat,
    dr1icaff
  )) %>% 
  group_by(seqn) %>% 
  summarize(
    total_k = sum(dr1ipota),
    total_na = sum(dr1isodi),
    total_total_fat = sum(dr1itfat),
    total_caffeine = sum(dr1icaff)
  ) 

write_xpt(diet, "data/diet.XPT")

diet %>% 
  glimpse()




# Parallel Processing ----
# Careful with this- if there's not enough memory,
# your computer can crash.
all_cores <- parallel::detectCores(logical = TRUE)
registerDoParallel(cores = all_cores - 1)


downloaded_data <- dir("data/", ignore.case = TRUE, pattern = "XPT", full.names = TRUE)


# Initiating the process of including everything in one Tibble.
# This is too much to do at once, so the goal is to widdle this down.
complete <- read_xpt("data/BPX_I.XPT") %>% 
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
blood_pressure <- read_xpt("data/BPX_I.XPT") %>% 
  janitor::clean_names()



# Defining a forward stepwise regression function.
# I will apply this to each file that I downloaded.
forward_predictors <- function(i_file){
  # Global assignment is necessary for ols_step_forward_p()
  file <<- read_xpt(i_file) %>% 
    janitor::clean_names()
  
  bp_data <<- read_xpt("data/BPX_I.XPT") %>% 
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





# How I would go about this for logistic regression -----------------------
health_train %>% 
  mutate(
    bp_status = as.factor(case_when(
      bpxsy1 >= 140 | bpxdi1 >= 90 ~ 1,
      bpxsy1 < 140 | bpxdi1 < 90 ~ 0 
    )
    )
  ) %>% 
  select(-bpxsy1)

health_test %>% 
  mutate(
    bp_status = as.factor(case_when(
      bpxsy1 >= 140 ~ 1,
      bpxsy1 < 140 ~ 0
    )
    )
  ) %>% 
  select(-bpxsy1)





# Boutique Analysis -------------------------------------------------------
# Looking at Medication Alone: Counts ----
# We need to decide how we will handle people on lots of medications
# in our models. There is no shame in leaving them in, we just need to
# include that in the model to account for that. 
bp <- read_xpt("data/BPX_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(seqn, bpxsy1, bpxdi1))


medication <- read_xpt("data/RXQ_RX_I.XPT") %>% 
  janitor::clean_names() %>% 
  glimpse()

# There are lots of people who don't consistently take medication;
# however, there are a non-neglible amount of people who take several.
medication %>% 
  group_by(seqn) %>% 
  summarize(
    count_meds = n() 
  ) %>% 
  ggplot(
    aes(x = count_meds)
  ) +
  geom_histogram() +
  theme_minimal()


# What if I take the number of medications into account?
# Would that explain some variation?
medication_counts <- medication %>% 
  group_by(seqn) %>% 
  summarize(
    count_meds = n() 
  )


# Global assignment is necessary for ols_step_forward_p()
file <- medication_counts

bp_data <- read_xpt("data/BPX_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(seqn, bpxsy1))

complete <- merge(bp_data, file) %>% 
  select(-seqn) 

intact_cols <- complete %>% 
  skimr::skim_without_charts() %>% 
  as_tibble() %>% 
  filter(
    n_missing < 2500
  ) %>% 
  select(skim_variable) %>% 
  pull()

complete <- complete %>% 
  select(c(intact_cols, bpxsy1)) %>% 
  scale() %>% 
  as_tibble()

set.seed(10201)
health_split <- initial_split(data = complete, prop = 0.8)
health_train <- training(health_split)
health_test <- testing(health_split)

model <- lm(
  data = health_train,
  bpxsy1 ~ .
)

model %>% 
  summary()

# This doesn't necessarily account for everything,
# but people taking medication is a big factor!!



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

medication_abridged


# Global assignment is necessary for ols_step_forward_p()
file <- medication_abridged

bp_data <- read_xpt("data/BPX_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(seqn, bpxsy1))

complete <- merge(bp_data, file) %>% 
  select(-seqn) 

intact_cols <- complete %>% 
  skimr::skim_without_charts() %>% 
  as_tibble() %>% 
  filter(
    n_missing < 2500
  ) %>% 
  select(skim_variable) %>% 
  pull()

complete <- complete %>% 
  select(c(intact_cols, bpxsy1))



set.seed(10201)
health_split <- initial_split(data = complete, prop = 0.8)
health_train <- training(health_split)
health_test <- testing(health_split)

model <- lm(
  data = health_train,
  bpxsy1 ~ .
) %>% 
  ols_step_forward_p()

meds_help <- head(model$predictors, 3)


overall_complete %>% 
  merge(medication_abridged) %>% 
  select(c(seqn, important_variables, meds_help)) %>% 
  glimpse()



# Looking at everything
overall_complete_subset <- overall_complete %>% 
  merge(medication_abridged) %>% 
  merge(medication_counts) %>% 
  select(c(seqn, important_variables, meds_help, count_meds))


file <- overall_complete_subset

bp_data <- read_xpt("data/BPX_I.XPT") %>% 
  janitor::clean_names() %>% 
  select(c(seqn, bpxsy1))

complete <- merge(bp_data, file) %>% 
  select(-seqn) 

intact_cols <- complete %>% 
  skimr::skim_without_charts() %>% 
  as_tibble() %>% 
  filter(
    n_missing < 2500
  ) %>% 
  select(skim_variable) %>% 
  pull()

complete <- complete %>% 
  select(c(intact_cols, bpxsy1))

set.seed(10201)
health_split <- initial_split(data = complete, prop = 0.8)
health_train <- training(health_split)
health_test <- testing(health_split)


model <- lm(
  data = health_train,
  bpxsy1 ~ .
)

model %>% 
  summary()


model %>% 
  ols_step_forward_p()


ols_results <- model %>% 
  ols_step_forward_p()

ols_results

best_predictors <- ols_results$predictors %>% 
  head(10)

best_predictors


best_forward_selection <- overall_complete %>% 
  merge(medication_abridged) %>% 
  merge(medication_counts) %>% 
  select(c(best_predictors))

best_forward_selection %>% 
  write_rds("data/best_forward.rds")

best_forward_selection %>% 
  names() %>% 
  as_tibble() %>% 
  rename(optimal_var = value) %>% 
  write_csv("data/optimal_vars.csv")



# Figuring Out Scaling ----------------------------------------------------
problem_scale <- read_xpt("data/INQ_I.XPT") %>% 
  janitor::clean_names()

problem_scale %>% 
  distinct()

problem_scale %>% 
  skimr::skim_without_charts() %>% 
  as_tibble() %>% 
  filter(
    complete_rate > 0.5 &
      numeric.p0 != numeric.p100
  ) %>% 
  print(n = Inf)


intact_cols_problem <- problem_scale %>% 
  skimr::skim_without_charts() %>% 
  as_tibble() %>% 
  filter(
    complete_rate > 0.5 &
      numeric.p0 != numeric.p100
  ) %>% 
  select(skim_variable) %>% 
  pull()

intact_cols_problem



problem_scale <- problem_scale %>%
  select(c(intact_cols_problem)) %>% 
  merge(blood_pressure %>% select(c(seqn, bpxsy1))) %>%
  select(-seqn) %>%
  scale() %>%
  as_tibble()

problem_scale %>%
  glimpse()

model <- lm(
  data = problem_scale,
  bpxsy1 ~ .
)

model %>% 
  summary()

model %>% 
  ols_step_forward_p()

problem_scale %>%
  glimpse()


