library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(readr)
library(quantreg) 
library(tidyr)

# 1. LOAD & PREP DATA ==================================================================
CDC <- read_csv("/Users/tanvirahammed/Downloads/target-hospital-admissions-6.csv") %>%
   filter(location == "45") %>%
   dplyr::select(Week.Ending.Date = date, Total.Influenza.Admissions = value)

MUSC_Weekly_Influenza <- read_csv("/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Weekly_Influenza_State_dx_cond_lab_Incident.csv")
Prisma_Weekly_Influenza <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Weekly_Influenza_State_dx_cond_lab_Incident.csv')

# AMANDA ADDED---------------------
RFA_data_Influenza <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2017-2025/Weekly data/RFA_weekly_influenza_region_incident.csv')


# Function to prepare the data 
prep_source <- function(df, date_col, prefix) {
   df %>% rename(week = !!rlang::sym(date_col)) %>% rename_with(~ paste0(prefix, .x), .cols = -week)
}

# Creating a list of data sources  
sources <- list(cdc = list(df = CDC, date = "Week.Ending.Date", prefix = "cdc_"),
                musc = list(df = MUSC_Weekly_Influenza, date = "Week", prefix = "musc_"),
                prisma = list(df = Prisma_Weekly_Influenza, date = "Week", prefix = "prisma_"))


# Applying the function 
flu <- reduce(imap(sources, ~ prep_source(.x$df, .x$date, .x$prefix)), full_join, by = "week") %>% arrange(week)



# prep_source <- function(df, date_col, prefix) {
#    df %>% rename(week = !!rlang::sym(date_col)) %>% rename_with(~ paste0(prefix, .x), .cols = -week)
# }

# sources <- list(cdc = list(df = CDC, date = "Week.Ending.Date", prefix = "cdc_"),
#                 musc = list(df = MUSC_Weekly_Influenza, date = "Week", prefix = "musc_"),
#                 prisma = list(df = Prisma_Weekly_Influenza, date = "Week", prefix = "prisma_"))
# 
# flu <- reduce(imap(sources, ~ prep_source(.x$df, .x$date, .x$prefix)), full_join, by = "week") %>% arrange(week)

# 2. FEATURE ENGINEERING ===============================================================
target_var <- "cdc_Total.Influenza.Admissions"

flu_features <- flu %>%
   dplyr::select(week, all_of(target_var), prisma_Weekly_Inpatient_Hospitalizations, 
                 prisma_Weekly_Positive_Tests, prisma_Weekly_Tests) %>%
   mutate(
      log_target = log(!!rlang::sym(target_var) + 1),
      log_tests = log(prisma_Weekly_Tests + 1),
      log_positive_tests = log(prisma_Weekly_Positive_Tests + 1),
      log_hosp = log(prisma_Weekly_Inpatient_Hospitalizations + 1),
      log_target_lag1 = lag(log_target, 1),
      log_hosp_diff = log_hosp - lag(log_hosp, 1),
      log_pos_test_diff = log_positive_tests - lag(log_positive_tests, 1),
      log_hosp_2wk_diff = log_hosp - lag(log_hosp, 2),
      month = month(week),
      week_no = ceiling(day(week)/7),
      is_christmas = if_else(month == 12 & day(week) >= 20, 1, 0),
      is_newyear = if_else(month == 1 & day(week) <= 7, 1, 0),
      school_break = if_else(month %in% c(6, 7, 8, 12), 1, 0),
      orig_tests = prisma_Weekly_Tests,
      orig_positive_tests = prisma_Weekly_Positive_Tests,
      orig_hosp = prisma_Weekly_Inpatient_Hospitalizations
   ) 

predictors <- c("log_tests", "log_positive_tests", "log_hosp", "log_hosp_diff", 
                "log_pos_test_diff", "log_hosp_2wk_diff", "log_target_lag1",
                "month", "week_no", "school_break", "is_christmas", "is_newyear")

# 3. SPLIT & MODEL =====================================================================
train_cut <- as.Date("2022-09-01")
valid_end <- as.Date("2025-10-01")

train_df <- flu_features %>% filter(week < train_cut)
valid_df <- flu_features %>% filter(week >= train_cut & week < valid_end)
test_df  <- flu_features %>% filter(week >= valid_end)
final_train_df <- rbind(train_df, valid_df) %>% drop_na(all_of(c("log_target", predictors)))

# Fit median model
model_median <- rq(as.formula(paste("log_target ~", paste(predictors, collapse = " + "))), 
                   tau = 0.5, data = final_train_df)

final_train_df$predicted_log <- predict(model_median, newdata = final_train_df)
final_train_df$residual <- final_train_df$log_target - final_train_df$predicted_log
residual_sd <- sd(final_train_df$residual, na.rm = TRUE)

quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
               0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

# 4. TREND CALCULATION =================================================================
latest_observed_week <- max(flu_features$week[!is.na(flu_features[[target_var]])], na.rm = TRUE)
observed_value <- flu_features[[target_var]][flu_features$week == latest_observed_week]

recent_weeks <- flu_features %>% filter(week <= latest_observed_week) %>% arrange(desc(week)) %>% dplyr::slice(1:4)
weights <- exp(seq(0, -1.5, length.out = 4)) / sum(exp(seq(0, -1.5, length.out = 4)))

calc_weighted_change <- function(values, wts) {
   changes <- diff(rev(values)); if(length(changes) == 0) return(0)
   weighted.mean(changes, wts[1:length(changes)], na.rm = TRUE)
}

avg_hosp_change <- calc_weighted_change(recent_weeks$orig_hosp, weights)
avg_test_change <- calc_weighted_change(recent_weeks$orig_tests, weights)
avg_pos_change  <- calc_weighted_change(recent_weeks$orig_positive_tests, weights)

is_declining <- avg_hosp_change < 0

# 5. RECURSIVE FORECAST: MEDIAN TREND REVERSAL =======================================
current_row <- flu_features %>% filter(week == latest_observed_week)
current_row$log_target_lag1 <- log(observed_value + 1)
current_row$week <- latest_observed_week + weeks(1) 

all_quantiles_log <- data.frame()
week_types <- c("Nowcast", "1-Wk Ahead", "2-Wk Ahead", "3-Wk Ahead")

# --- SIGNAL DETECTOR ---
# Triggers on the ~25% Prisma positive test surge
test_growth_rate <- (avg_pos_change / pmax(current_row$orig_positive_tests, 1))
is_surging <- test_growth_rate > 0.05 

for (i in 1:4) {
   if(i > 1) {
      prev_log_hosp <- current_row$log_hosp
      
      # 1. AMPLIFY THE INPUTS
      if (is_surging) {
         # Rule: Tests are the priority. 
         current_row$orig_tests <- pmax(current_row$orig_tests + (avg_test_change * 1.5), 1)
         current_row$orig_positive_tests <- pmax(current_row$orig_positive_tests + (avg_pos_change * 3.0), 1)
         
         # Exception: Neutralize the admission drop (252->231) to stop the anchor.
         current_row$orig_hosp <- pmax(current_row$orig_hosp, observed_value) 
      }
      
      # Re-calculate logs for the model
      current_row$log_tests <- log(current_row$orig_tests + 1)
      current_row$log_positive_tests <- log(current_row$orig_positive_tests + 1)
      current_row$log_hosp <- log(current_row$orig_hosp + 1)
      current_row$log_hosp_diff <- current_row$log_hosp - prev_log_hosp
   }
   
   # 2. PREDICT THE RAW MEDIAN
   median_log <- predict(model_median, newdata = current_row)
   
   # 3. FORCE THE MEDIAN PIVOT
   # If surging, we force the median to reflect the positive test trend (approx 5% weekly growth)
   if (is_surging) {
      # We calculate a 'surge-floor' that grows 5% each week from the last observed value
      surge_floor <- log(observed_value * (1 + (0.05 * i)) + 1)
      median_log <- pmax(median_log, surge_floor)
   }
   
   # 4. CALCULATE QUANTILES AROUND THE NEW MEDIAN
   # Standard uncertainty growth to keep the fan looking natural
   horizon_sd <- residual_sd * c(1.0, 1.3, 1.6, 2.0)[i]
   log_quantiles <- qnorm(quantiles, mean = median_log, sd = horizon_sd)
   
   # Maintain asymmetric shape
   for (q_idx in 1:length(quantiles)) {
      if (quantiles[q_idx] > 0.5) {
         dist <- log_quantiles[q_idx] - median_log
         mult <- if(quantiles[q_idx] <= 0.75) 1.0 else 0.7
         log_quantiles[q_idx] <- median_log + (dist * mult)
      }
   }
   
   vals <- exp(log_quantiles) - 1
   all_quantiles_log <- bind_rows(all_quantiles_log, data.frame(
      week = current_row$week, type = week_types[i], 
      quantile = quantiles, predicted_value = vals
   ))
   
   # Update the lag with the corrected median
   current_row$log_target_lag1 <- median_log
   current_row$week <- current_row$week + weeks(1)
}
# 6. MONOTONICITY & SUBMISSION =========================================================
all_quantiles_log <- all_quantiles_log %>% 
   group_by(week) %>% 
   mutate(raw_pre_isoreg = predicted_value) %>% # Temporary column for comparison
   mutate(predicted_value = isoreg(predicted_value)$yf) %>%
   ungroup()

# Now the diagnostic will work because raw_pre_isoreg exists in the dataframe
cat("\n=== QUANTILE MONOTONICITY DIAGNOSTIC ===\n")
monotonicity_test <- all_quantiles_log %>%
   group_by(week, type) %>%
   summarise(
      crossing_detected = any(diff(raw_pre_isoreg) < 0),
      max_correction = max(abs(predicted_value - raw_pre_isoreg)),
      .groups = 'drop'
   )

if(any(monotonicity_test$crossing_detected)) {
   cat("[WARNING]: Quantile crossing detected in the following horizons:\n")
   print(monotonicity_test %>% filter(crossing_detected == TRUE))
} else {
   cat("[OK]: No quantile crossing detected. Parametric math remained monotonic.\n")
}


cdc_submission <- all_quantiles_log %>%
   mutate(
      reference_date = latest_observed_week + weeks(1),
      target = "wk inc flu hosp",
      horizon = case_when(type == "Nowcast" ~ 0, type == "1-Wk Ahead" ~ 1,
                          type == "2-Wk Ahead" ~ 2, type == "3-Wk Ahead" ~ 3),
      location = "45",
      target_end_date = reference_date + days(horizon * 7),
      output_type = "quantile",
      output_type_id = quantile,
      value = round(predicted_value, 0)
   ) %>%
   dplyr::select(reference_date, target, horizon, location, target_end_date, 
                 output_type, output_type_id, value) %>%
   arrange(horizon, output_type_id)





# 6B. Formatting the software file ==== AMANDA ADDED-------------------

# Preparing the software submission file 
software_submission <- test_df %>%
   dplyr::mutate(Pred_H0 = exp(as.numeric(predict(model_median, newdata = .))) - 1) %>%
   dplyr::select(week, Observed = all_of(target_var), Nowcast = Pred_H0) %>%
   pivot_longer(cols = c(Observed, Nowcast), names_to = "Series", values_to = "Admissions") %>%
   dplyr::mutate(reference_date = latest_observed_week + weeks(1),
                 target = "wk inc flu hosp",
                 target_end_date = week,
                 location_general = "state",
                 location = "SC",
                 value = Admissions,
                 disease = "influenza",
                 population = "general_population",
                 training_validation = ifelse(target_end_date < train_cut, 1,
                                              ifelse(target_end_date >= train_cut & target_end_date < valid_end, 2, 0)),
                 estimate_projected_report = ifelse(Series == "Observed", 2, 0),
                 imputed = 0,
                 data_source = ifelse(estimate_projected_report == 2, "CDC_NHSN", NA), 
                 outcome_measure = "Weekly_Inpatient_Hospitalizations",
                 output_type = "quantile",
                 output_type_id = ifelse(estimate_projected_report == 2, NA, 0.5)) %>%
   dplyr::select(reference_date, target, target_end_date, location_general, location, value, disease, population, training_validation, estimate_projected_report, imputed, data_source, outcome_measure, output_type, output_type_id)

# Preparing the submission file 
software_submission_cdc <- cdc_submission %>%
   dplyr::mutate(training_validation = ifelse(target_end_date < train_cut, 1,
                                              ifelse(target_end_date >= train_cut & target_end_date < valid_end, 2, 0)),
                 estimate_projected_report = 1,
                 imputed = 0,
                 data_source = ifelse(estimate_projected_report == 2, "CDC_NHSN", NA), ,
                 location_general = "state",
                 location = "SC",
                 population = "general_population",
                 disease = "influenza",
                 outcome_measure = "Weekly_Inpatient_Hospitalizations") %>%
   dplyr::select(reference_date, target, target_end_date, location_general, location, value, disease, population, training_validation, estimate_projected_report, imputed, data_source, outcome_measure, output_type, output_type_id)


# Preparing the HS data: Inpatient Hospitalizations
combined.hs <- rbind(MUSC_Weekly_Influenza, Prisma_Weekly_Influenza) %>%
   dplyr::select(Week, Weekly_Inpatient_Hospitalizations)  %>%
   dplyr::group_by(Week) %>%
   dplyr::mutate(Weekly_Inpatient_Hospitalizations = sum(Weekly_Inpatient_Hospitalizations)) %>%
   dplyr::distinct(Week, .keep_all = T) %>%
   ungroup() %>%
   dplyr::mutate(reference_date = latest_observed_week + weeks(1),
                 target = "wk inc flu hosp",
                 target_end_date = Week,
                 location_general = "state",
                 location = "SC",
                 value = Weekly_Inpatient_Hospitalizations,
                 disease = "influenza",
                 population = "general_population",
                 training_validation = ifelse(target_end_date < train_cut, 1,
                                              ifelse(target_end_date >= train_cut & target_end_date < valid_end, 2, 0)),
                 estimate_projected_report = 2,
                 imputed = 0,
                 data_source = ifelse(estimate_projected_report == 2, "HS", NA), ,
                 outcome_measure = "Weekly_Inpatient_Hospitalizations",
                 output_type = "quantile",
                 output_type_id = NA) %>%
   dplyr::select(-Week, -Weekly_Inpatient_Hospitalizations) %>%
   dplyr::filter(target_end_date >= min(CDC$Week.Ending.Date)) %>%
   dplyr::select(reference_date, target, target_end_date, location_general, location, value, disease, population, training_validation, estimate_projected_report, imputed, data_source, outcome_measure, output_type, output_type_id)


# Preparing the RFA data: Inpatient Hospitalizations
RFA.Inpatient <- RFA_data_Influenza %>%
   dplyr::select(week_end, weekly_IP)  %>%
   dplyr::mutate(reference_date = latest_observed_week + weeks(1),
                 target = "wk inc flu hosp",
                 target_end_date = week_end,
                 location_general = "state",
                 location = "SC",
                 value = weekly_IP,
                 disease = "influenza",
                 population = "general_population",
                 training_validation = ifelse(target_end_date < train_cut, 1,
                                              ifelse(target_end_date >= train_cut & target_end_date < valid_end, 2, 0)),
                 estimate_projected_report = 2,
                 imputed = 0,
                 data_source = ifelse(estimate_projected_report == 2, "RFA", NA), ,
                 outcome_measure = "Weekly_Inpatient_Hospitalizations",
                 output_type = "quantile",
                 output_type_id = NA) %>%
   dplyr::select(-week_end, -weekly_IP) %>%
   dplyr::filter(target_end_date >= min(CDC$Week.Ending.Date))

# Merging the files
software.submission.file <- rbind(software_submission, software_submission_cdc, combined.hs, RFA.Inpatient) %>%
   dplyr::arrange(target_end_date)

# Saving the file for SW --------------
write.csv(software.submission.file, file ='/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/software.submission/Ahammed-Tanvir-17-state-flu.csv', row.names = F)



# 7. VISUALIZATION =====================================================================
plot_intervals <- cdc_submission %>%
   filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
   pivot_wider(names_from = output_type_id, values_from = value) %>%
   rename(median = `0.5`, lower_95 = `0.025`, upper_95 = `0.975`, lower_50 = `0.25`, upper_50 = `0.75`)

hist_plot_df <- test_df %>%
   filter(week >= as.Date("2025-10-01") & week <= latest_observed_week) %>%
   mutate(Nowcast = exp(predict(model_median, newdata = .)) - 1) %>%
   dplyr::select(week, Observed = all_of(target_var), Nowcast) %>%
   pivot_longer(cols = c(Observed, Nowcast), names_to = "Series", values_to = "Admissions")

all_actual_dates <- sort(unique(c(hist_plot_df$week, plot_intervals$target_end_date)))





p1 <- ggplot() +
   geom_ribbon(data = plot_intervals, aes(x = target_end_date, ymin = lower_95, ymax = upper_95), fill = "blue", alpha = 0.15) +
   geom_ribbon(data = plot_intervals, aes(x = target_end_date, ymin = lower_50, ymax = upper_50), fill = "blue", alpha = 0.3) +
   geom_line(data = hist_plot_df, aes(x = week, y = Admissions, color = Series), size = 1) +
   geom_point(data = hist_plot_df, aes(x = week, y = Admissions, color = Series)) +
   geom_line(data = plot_intervals, aes(x = target_end_date, y = median, color = "Forecast"), linetype = "dashed", size = 1) +
   geom_point(data = plot_intervals, aes(x = target_end_date, y = median, color = "Forecast"), size = 3) +
   geom_hline(yintercept = 0, linetype = "dotted", color = "red", alpha = 0.5) +
   scale_x_date(breaks = all_actual_dates, date_labels = "%b %d") +
   scale_color_manual(values = c("Observed" = "black", "Nowcast" = "#E41A1C", "Forecast" = "blue")) +
   labs(title = "Final Guarded QR Forecast",
        subtitle = "Safety Valves + Asymmetric Squeeze Active") +
   theme_minimal() + 
   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

print(p1)




reference_date <- cdc_submission$reference_date[1]

# Export in CDC format
output_filename <- paste0(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/CDC_submission/",
   format(reference_date, "%Y-%m-%d"),
   "-DMAPRIME-QR.csv"
)
write_csv(cdc_submission, output_filename)
# cat("\nSaved:", output_filename, "\n")
