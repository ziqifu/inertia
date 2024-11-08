
df_ch_streak <- read_csv("~/Desktop/Otto_lab/Perseveration/analysis/perseveration_exp2_all.csv")
# Step 0: Dependencies and Read-ins ------
setwd('~/Desktop/Otto_Lab/Perseveration/recency')
#rm(list = ls())

{
  library(tidyverse)
  library(R.matlab)
  library(diptest)
  library(moments)
  #library(BEST)
  library(jtools)
  library(readxl)
  library(patchwork)
  library(ggbeeswarm)
  library(cocor)
  library(ggpubr)
  library(rstatix)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(doBy)
  library(purrr)
  library(broom)
  
  
  ci = function(x, a=qnorm(.975), na.rm=T) plotrix::std.error(x, na.rm=na.rm) * a
  tab_mod = sjPlot::tab_model
}

# read-in my dataset-----
df_all = read_csv("~/Desktop/Otto_lab/Perseveration/analysis/perseveration_exp2_all.csv")

df_all = df_all[,-1]
df_all = df_all[order(df_all$id, df_all$block, df_all$trial), ]
df_all = df_all %>%
  filter(trial >= 2)

ext_rate = df_all %>%
  group_by(id) %>%
  summarise(ext_rate = mean(ext))
extreme_ext = ext_rate$id[ext_rate$ext_rate > .99 | ext_rate$ext_rate < .01]
df_all = df_all %>%
  filter(!id %in% extreme_ext)

# Alternative 1-----
# Step 1: Determine the best sample-level lambda-------

# Define possible deltas to grid search over
lambda_values = seq(0.01, 0.99, by = 0.01)

# Initialize variable to track the best overall lambda
best_overall_fit = list(deviance = Inf, lambda = NA)

# Loop through all lambda value

for (lambda in lambda_values) {
  
  # Initialize sum of devs for the current lambda across all subjects
  total_deviance = 0
  
  # Loop through all subjects
  for(id in unique(df_all$id)){
    thisID = df_all[df_all$id==id,]
    thisPreExts = thisID$pre_ext
    thisCurExts = thisID$ext
    thisAngle = thisID$angle_c
    deltas = c(0.5)
    for (pre_ext in thisPreExts) {
      deltas = c(deltas, (1 - lambda) * tail(deltas, 1) + lambda * pre_ext)
    }
    deltas = deltas[-1]  # Remove the initial value
    pre_deltas =  lag(deltas) # Use previous delta
    pre_deltas = pre_deltas[-1]
    thisCurExts = thisCurExts[-1] # start from the 3rd trial
    thisCurAngle = thisAngle[-1]
    # Fit the GLM
    log_reg = glm(thisCurExts ~ pre_deltas + thisCurAngle, family = "binomial")
    model_sum = summary(log_reg)
    print(id)
    print(model_sum)
    # Calculate the sum of squared devs
    model_deviance = deviance(log_reg)
    beta_a1 <- model_sum$coefficients['pre_deltas', "Estimate"]
    
    total_deviance = total_deviance + model_deviance
  }
  
  # Update the best overall delta
  if (total_deviance < best_overall_fit$deviance) {
    best_overall_fit$deviance = total_deviance
    best_overall_fit$lambda = lambda
    best_beta_a1 = beta_a1
  }
}

cat(sprintf("Best sample-level lambda: %f with Total Deviance: %f\n", best_overall_fit$lambda, best_overall_fit$deviance))

# Step 2: Compute delta for each trial of each subject ------
control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000), tol = 1e-4,
                        nAGQ = 10)

lambda_sample = best_overall_fit$lambda 

df_all$deltas = NA
# Loop through all subjects and compute delta for each trial of each subject
for(id in unique(df_all$id)){
  thisID = df_all[df_all$id==id,]
  thisPreExts = thisID$pre_ext
  thisCurExts = thisID$ext
  thisAngle = thisID$angle_c
  thisTrial = thisID$trial
  deltas = c(0.5)
  for (pre_ext in thisPreExts) {
    deltas = c(deltas, (1 - lambda_sample) * tail(deltas, 1) + lambda_sample * pre_ext)
  }
  deltas = deltas[-1]  
  df_all$deltas[df_all$id==id] = deltas
}



# Fit the GLMER to the whole sample
sample_glmer = df_all %>%
  group_by(id) %>%
  mutate(pre_delta = lag(deltas, default = NA),
         cur_angle = case_when(angle == 60 ~ 0,
                               angle == 120 ~ 0.5,
                               angle == 180 ~ 1,
                              )) %>%
  filter(!is.na(pre_delta))

m_subject_angle = glmer(ext ~ pre_delta*cur_angle + (cur_angle*pre_delta|id), 
                        family = "binomial", control = control, data = sample_glmer)
print(summary(m_subject_angle))
tab_mod(m_subject_angle, transform=NULL, file='./step2_glmer.html')


m_nfc_included_int = glmer(ext ~ pre_delta*nfc_c + cur_angle*nfc_c + (cur_angle*pre_delta|id), 
                        family = "binomial", control = control, data = sample_glmer)
print(summary(m_nfc_included_int))
tab_mod(m_nfc_included_int, transform=NULL, file='./step2_glmer_all.html')

m_nfc_included = glmer(ext ~ pre_delta*nfc_s + cur_angle*nfc_s + (cur_angle + pre_delta|id), 
                       family = "binomial", control = control, data = sample_glmer)
print(summary(m_nfc_included))
tab_mod(m_nfc_included, transform=NULL, file='./step2_glmer_all.html')

anova_comparison <- anova(m_subject_angle, m_nfc_included)

# Print the result
print(anova_comparison)

sample_glmer = sample_glmer[,-c()]
# Plot: ext against pre_delta per trial 
test = sample_glmer %>%
  filter(id == 101) %>%
  mutate(trial_c = 48*(block - 10) + trial)

ggplot(test, aes(x = trial_c)) + 
  geom_line(aes(y = ext, group = id, color = "ext(t)"), size = 1) +  # Line for pre_ext
  geom_line(aes(y = pre_delta, group = id, color = "Delta(t-1)"), size = 1, linetype = "dashed") +  # Dashed line for deltas
  scale_color_manual(values = c("ext(t)" = "black", "Delta(t-1)" = "blue")) +  # Color coding
  facet_wrap(~id, scales = "free_y") +  # Separate plots for each subject, with independent y scales
  labs(title = "Empirical vs Estimated for Subject id = 100",
       x = "Trial Number",
       y = "Ext") +  # Labels
  theme_minimal() +
  theme(legend.position = "bottom")  # Positioning the legend at the bottom

# Step 3: Compute beta for each participant by fitting glm -----
beta_each_subj <- sample_glmer %>%
  group_by(id) %>%
  group_split() %>%  # Splits the data into a list, each element is a dataframe for each id
  map_df(function(df) {
    # Fit the glm model for the subset
    model <- glm(ext ~ pre_delta + cur_angle, family = "binomial", data = df)
    # Use broom to tidy the model and extract coefficients
    tidy_data <- tidy(model) %>%
      filter(term == "pre_delta") %>%
      select(estimate, p.value, term) %>%  # Select relevant columns
      mutate(id = first(df$id))  # Ensure the correct id is included for each subset
    return(tidy_data)
  }, .id = NULL) 

beta_sign = beta_each_subj %>%
  filter(p.value <= .05)

beta_insign = beta_each_subj %>%
  filter(p.value > .05)

ext_rate = df_all %>%
  group_by(id) %>%
  summarise(ext_rate = mean(ext))

ext_rate
beta_ext = merge(beta_sign, ext_rate, by = 'id')

lm = lm(estimate ~ ext_rate, data = beta_ext)
summary(lm)
# Histogram for each subject's beta
# Create the histogram
png(file = "/Desktop/Otto_lab/Perseveration/recency/subj_beta_significant.png", width = 800, height = 600, res = 100)

hist(beta_sign$estimate, main = "Distribution of beta of pre_ext", xlab = "Beta")

# Calculate mean and median of lambda
mean_beta <- mean(beta_sign$estimate)
median_beta <- median(beta_sign$estimate)

# Add a red line for the mean
abline(v = mean_beta, col = "red", lwd = 2)

# Add a green line for the median
abline(v = median_beta, col = 'blue', lwd = 2)

dev.off()


# Alternative 2: Estimate best trial-level lambda for each subject -------
lambda_values = seq(0.01, 0.99, by = 0.01)

# Initialize variable to track the best overall lambda
best_lambda_per_subject = list()

# Loop through each unique subject ID
for (id in unique(df_all$id)) {
  
  # Extract data for the current subject
  thisID = df_all[df_all$id == id,]
  thisPreExts = thisID$pre_ext
  thisCurExts = thisID$ext
  thisAngle = thisID$angle_c
  # Initialize variables to track the best lambda for the current subject
  best_fit_for_subject = list(deviance = Inf, beta = NA, lambda = NA)
  
  # Loop through all lambda values
  for (lambda in lambda_values) {
    deltas = c(0.5)  # Initialize deltas with the first value
  
    for (pre_ext in thisPreExts) {
      deltas = c(deltas, (1 - lambda) * tail(deltas, 1) + lambda * pre_ext)
    }
    deltas = deltas[-1]  # Remove the initial value
    pre_deltas =  lag(deltas) # Use previous delta
    pre_deltas = pre_deltas[-1]
    CurExts = thisCurExts[-1] # start from the 3rd trial
    CurAngle = thisAngle[-1]
    # Fit the logistic regression model
    log_reg = glm(CurExts ~ pre_deltas + CurAngle, family = binomial)
    model_summary = summary(log_reg)
    # Calculate the deviance of the model
    model_deviance = deviance(log_reg)
    # Get the beta of delta of the model
    beta_delta_a2 <- model_summary$coefficients['pre_deltas', "Estimate"]
    
    # Update the best lambda if the current model has lower deviance
    if (model_deviance < best_fit_for_subject$deviance) {
      best_fit_for_subject$deviance = model_deviance
      best_fit_for_subject$beta = beta_delta_a2
      best_fit_for_subject$lambda = lambda
    }
  }
  
  # Track model fitting process
  print(paste0(id, ', done'))
  
  # Store the best fit lambda for the current subject
  best_lambda_per_subject[[id]] = best_fit_for_subject
}

best_lambda_subj <- data.frame(id = integer(), deviance = numeric(), beta = numeric(), lambda = numeric())

# Iterate through the list and extract details
for (i in seq_along(best_lambda_per_subject)) {
  if (!is.null(best_lambda_per_subject[[i]])) {
    # Append each non-null entry to the data frame
    best_lambda_subj <- rbind(best_lambda_subj, data.frame(id = i, 
                                               deviance = best_lambda_per_subject[[i]]$deviance, 
                                               beta_subj = best_lambda_per_subject[[i]]$beta,
                                               lambda_subj = best_lambda_per_subject[[i]]$lambda))
    
  }
}
best_lambda_subj = best_lambda_subj %>%
  unique()

lambda_beta_subj = lm(lambda_subj ~ beta_subj, data = best_lambda_subj)
summary(lambda_beta_subj)


# Print the resulting data frame
print(best_lambda_subj)


# Assuming 'best_lambda_subj' is your dataframe containing the lambda column
# Create the histogram
hist(best_lambda_subj$lambda, main = "Distribution of Lambda", xlab = "Lambda")

# Calculate mean and median of lambda
mean_lambda <- mean(best_lambda_subj$lambda)
median_lambda <- median(best_lambda_subj$lambda)

# Add a red line for the mean
abline(v = mean_lambda, col = "red", lwd = 2)

# Add a green line for the median
abline(v = median_lambda, col = 'blue', lwd = 2)

dev.off()

# Add lambda to the dataset
df_all = merge(best_lambda_subj, df_all, by = 'id')

# Test if lambda has anything to do with NFC or Demographic data, nothing
demo_subj = df_all %>%
  group_by(id) %>%
  summarise(ncs_tot = nfc_c, gender = gender, age = age) %>%
  unique()

var_all = merge(best_lambda_subj, demo_subj, by = 'id')

nfc_lambda = lm(lambda ~ beta*gender*age, data = var_all)
summary(nfc_lambda)

# Get the beta for delta in each participant
# Loop through all subjects and compute delta for each trial of each subject
for(id in unique(df_all$id)){
  thisID = df_all[df_all$id==id,]
  thisPreExts = thisID$pre_ext
  thisCurExts = thisID$ext
  thisAngle = thisID$angle_c
  thisTrial = thisID$trial
  deltas = c(0.5)
  lambda_subject = unique(thisID$lambda)
  for (pre_ext in thisPreExts) {
    deltas = c(deltas, (1 - lambda_subject) * tail(deltas, 1) + lambda_subject * pre_ext)
  }
  deltas = deltas[-1]  
  df_all$deltas_subject[df_all$id==id] = deltas
}

a2 = df_all %>%
  group_by(id) %>%
  mutate(pre_delta_subject = lag(deltas_subject, default = NA),
         cur_angle = case_when(angle == 60 ~ 0,
                               angle == 120 ~ 0.5,
                               angle == 180 ~ 1,
         )) %>%
  filter(!is.na(pre_delta_subject))

a2_glmer = glmer(ext ~ pre_delta_subject + cur_angle + (cur_angle|id), 
                        family = "binomial", data = subj_glmer)
print(summary(a2_glmer))
tab_mod(m_subject_angle, transform=NULL, file='./recency/step2_glmer.html')

