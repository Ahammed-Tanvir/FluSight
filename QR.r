library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(readr)
library(quantreg) 
library(tidyr)

# 1. LOAD DATA =========================================================================
CDC <- read_csv("/Users/tanvirahammed/Downloads/target-hospital-admissions-2.csv") %>%
   filter(location == "45") %>%
   select(Week.Ending.Date = date, Total.Influenza.Admissions = value)

MUSC_Weekly_Influenza <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Influenza_State.csv')
Prisma_Weekly_Influenza <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Influenza_State.csv')

prep_source <- function(df, date_col, prefix) {
   df %>%
      rename(week = !!rlang::sym(date_col)) %>%
      rename_with(~ paste0(prefix, .x), .cols = -week)
}

sources <- list(
   cdc    = list(df = CDC, date = "Week.Ending.Date", prefix = "cdc_"),
   musc   = list(df = MUSC_Weekly_Influenza, date = "Week", prefix = "musc_"),
   prisma = list(df = Prisma_Weekly_Influenza, date = "Week", prefix = "prisma_")
)

prepared_list <- imap(sources, ~ prep_source(.x$df, .x$date, .x$prefix))
flu <- reduce(prepared_list, full_join, by = "week") %>% arrange(week)

# 2. FEATURE ENGINEERING WITH LOG TRANSFORMATIONS =========================================================================
target_var <- "cdc_Total.Influenza.Admissions"

flu_features <- flu %>%
   select(week, all_of(target_var), 
          prisma_Weekly_Inpatient_Hospitalizations, prisma_Weekly_Positive_Tests, prisma_Weekly_Tests) %>%
   mutate(
      # LOG TRANSFORM THE TARGET (add 1 to handle zeros)
      log_target = log(!!rlang::sym(target_var) + 1),
      
      # LOG TRANSFORM COUNT PREDICTORS (add 1 to handle zeros)
      log_tests = log(prisma_Weekly_Tests + 1),
      log_positive_tests = log(prisma_Weekly_Positive_Tests + 1),
      log_hosp = log(prisma_Weekly_Inpatient_Hospitalizations + 1),
      
      # Create 1-week lag of target (for nowcast context only)
      log_target_lag1 = lag(log_target, 1),
      
      # Differences on LOG SCALE (log differences = growth rates)
      log_hosp_diff = log_hosp - lag(log_hosp, 1),
      log_pos_test_diff = log_positive_tests - lag(log_positive_tests, 1),
      
      # Recent momentum (2-week change for stronger signal)
      log_hosp_2wk_diff = log_hosp - lag(log_hosp, 2),
      
      # Time features
      year = year(week),
      month = month(week),
      week_no = ceiling(day(week)/7),
      is_christmas = if_else(month == 12 & day(week) >= 20, 1, 0),
      is_newyear = if_else(month == 1 & day(week) <= 7, 1, 0),
      school_break = if_else(month %in% c(6, 7, 8, 12), 1, 0),
      
      # Keep original values for reference and recursive forecasting
      orig_tests = prisma_Weekly_Tests,
      orig_positive_tests = prisma_Weekly_Positive_Tests,
      orig_hosp = prisma_Weekly_Inpatient_Hospitalizations
   ) 

# *** REMOVED log_target_lag2 but kept lag1 for recent context ***
# Predictors now focus on current state, momentum, and most recent target value
predictors <- c("log_tests", "log_positive_tests", "log_hosp",
                "log_hosp_diff", "log_pos_test_diff", "log_hosp_2wk_diff",
                "log_target_lag1",
                "month", "week_no", "school_break", "is_christmas", "is_newyear")

cat("\n=== MODEL CONFIGURATION ===\n")
cat("Target variable (log scale):", "log_target", "\n")
cat("*** Using only 1-week lagged target (not 2-week) for balance ***\n")
cat("Predictors:", paste(predictors, collapse = ", "), "\n\n")

# =========================================================================
# 3. TRAIN / TEST SPLIT
# =========================================================================
train_cut   <- as.Date("2022-09-01")
valid_start <- as.Date("2022-09-01")
valid_end   <- as.Date("2025-10-01")
test_start  <- as.Date("2025-10-01")

train_df <- flu_features %>% filter(week < train_cut)
valid_df <- flu_features %>% filter(week >= valid_start & week < valid_end)
test_df  <- flu_features %>% filter(week >= test_start) 

final_train_df <- rbind(train_df, valid_df) %>% drop_na(all_of(c("log_target", predictors)))

cat("Training observations:", nrow(final_train_df), "\n")
cat("Test observations:", nrow(test_df), "\n\n")

# =========================================================================
# 4. MODEL FITTING (Multi-Quantile on LOG SCALE)
# =========================================================================
quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
               0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

# Formula now predicts LOG TARGET using LOG COUNT PREDICTORS (NO LAGS)
formula_log <- as.formula(paste("log_target ~", paste(predictors, collapse = " + ")))

cat("=== FITTING LOG-SCALE QUANTILE REGRESSION ===\n")
cat("Formula:", deparse(formula_log), "\n\n")

# Fit Quantile Regression for all tau values
model_log_multi <- rq(formula = formula_log, tau = quantiles, data = final_train_df)

cat("Model fitted successfully!\n\n")

# =========================================================================
# 5. GENERATE TEST PREDICTIONS (if test_df has data)
# =========================================================================
if(nrow(test_df) > 0) {
   quantile_preds_log <- predict(model_log_multi, newdata = test_df)
   # Back-transform from log scale
   test_df$Pred_H0 <- exp(as.numeric(quantile_preds_log[, "tau= 0.500"])) - 1
   cat("Test predictions generated\n\n")
}

# =========================================================================
# 6. TREND-BASED RECURSIVE FORECASTING WITH WEIGHTED TRENDS
# =========================================================================
latest_data <- flu_features %>% 
   arrange(desc(week)) %>% 
   dplyr::slice(1) %>%
   drop_na(all_of(predictors))

forecast_horizon <- 3

# Get recent weeks for trend calculation
recent_weeks <- flu_features %>% 
   arrange(desc(week)) %>% 
   dplyr::slice(1:4)

# *** WEIGHTED TREND CALCULATION ***
# Exponential weights: most recent week gets highest weight
weights <- exp(seq(0, -1.5, length.out = 4))
weights <- weights / sum(weights)

cat("\n=== WEIGHTED TREND PROJECTION ===\n")
cat("Week weights (newest to oldest):", paste(round(weights, 3), collapse = ", "), "\n")

# Helper function for weighted average of changes
calc_weighted_change <- function(values, wts) {
   rev_vals <- rev(values)
   changes <- diff(rev_vals)
   
   if(length(changes) == 0 || all(is.na(changes))) return(0)
   
   weighted.mean(changes, wts[1:length(changes)], na.rm = TRUE)
}

# Calculate WEIGHTED average weekly changes (original scale)
avg_test_change <- calc_weighted_change(recent_weeks$orig_tests, weights)
avg_pos_change <- calc_weighted_change(recent_weeks$orig_positive_tests, weights)
avg_hosp_change <- calc_weighted_change(recent_weeks$orig_hosp, weights)

cat("Weighted average weekly changes (original scale):\n")
cat("  Tests:", round(avg_test_change, 2), "\n")
cat("  Positive Tests:", round(avg_pos_change, 2), "\n")
cat("  Hospitalizations:", round(avg_hosp_change, 2), "\n\n")

# Initialize results dataframe for all quantiles
all_quantiles_log <- data.frame()

# Week labels
week_types <- c("Nowcast", "1-Wk Ahead", "2-Wk Ahead", "3-Wk Ahead")
week_dates <- c(latest_data$week, seq(latest_data$week + weeks(1), by = "1 week", length.out = forecast_horizon))

# Initialize recursive state
current_row <- latest_data

cat("=== STARTING RECURSIVE FORECASTING ===\n")

for (i in 1:4) {  # Nowcast + 3 forecasts
   cat("Forecasting", week_types[i], "for week", as.character(week_dates[i]), "\n")
   
   # Predict ALL quantiles on LOG SCALE using the current row state
   pred_matrix_log <- predict(model_log_multi, newdata = current_row)
   
   # BACK-TRANSFORM to original scale: exp(log_pred) - 1
   quantile_results <- data.frame(
      week = week_dates[i],
      type = week_types[i],
      quantile = quantiles,
      log_predicted = as.numeric(pred_matrix_log[1, ]),
      predicted_value = exp(as.numeric(pred_matrix_log[1, ])) - 1  # Back-transform
   )
   
   # Append to master dataframe
   all_quantiles_log <- bind_rows(all_quantiles_log, quantile_results)
   
   pred_median_log <- as.numeric(pred_matrix_log[, "tau= 0.500"])
   pred_median_original <- exp(pred_median_log) - 1
   
   if (i < 4) {  # Don't update after last prediction
      # Update Date and Time features
      next_week_date <- current_row$week + weeks(1)
      current_row$week <- next_week_date
      current_row$month <- month(next_week_date)
      current_row$year  <- year(next_week_date)
      current_row$week_no <- ceiling(day(next_week_date)/7)
      
      # Update Seasonality/Holidays
      current_row$is_christmas <- if_else(current_row$month == 12 & day(next_week_date) >= 20, 1, 0)
      current_row$is_newyear   <- if_else(current_row$month == 1 & day(next_week_date) <= 7, 1, 0)
      current_row$school_break <- if_else(current_row$month %in% c(6, 7, 8, 12), 1, 0)
      
      # Update lag1 with current prediction
      current_row$log_target_lag1 <- pred_median_log
      
      # Store previous values for diff calculation
      prev_log_tests <- current_row$log_tests
      prev_log_pos <- current_row$log_positive_tests
      prev_log_hosp <- current_row$log_hosp
      
      # PROJECT EHR DATA on ORIGINAL SCALE using weighted trends
      current_row$orig_tests <- pmax(current_row$orig_tests + avg_test_change, 1)
      current_row$orig_positive_tests <- pmax(current_row$orig_positive_tests + avg_pos_change, 1)
      current_row$orig_hosp <- pmax(current_row$orig_hosp + avg_hosp_change, 1)
      
      # Convert projected original values back to LOG SCALE for model
      current_row$log_tests <- log(current_row$orig_tests + 1)
      current_row$log_positive_tests <- log(current_row$orig_positive_tests + 1)
      current_row$log_hosp <- log(current_row$orig_hosp + 1)
      
      # Update LOG DIFFS (growth rates)
      current_row$log_hosp_diff <- current_row$log_hosp - prev_log_hosp
      current_row$log_pos_test_diff <- current_row$log_positive_tests - prev_log_pos
      
      # Update 2-week diff
      if(i == 1) {
         # For first forecast, use actual historical data
         hosp_2wk_ago <- flu_features %>% 
            filter(week == (current_row$week - weeks(2))) %>% 
            pull(log_hosp)
         current_row$log_hosp_2wk_diff <- if(length(hosp_2wk_ago) > 0) current_row$log_hosp - hosp_2wk_ago[1] else current_row$log_hosp_diff * 2
      } else {
         # For subsequent forecasts, approximate as 2x the 1-week diff
         current_row$log_hosp_2wk_diff <- current_row$log_hosp_diff * 2
      }
   }
}

cat("\n=== FORECASTING COMPLETE ===\n\n")

# =========================================================================
# 7. ENFORCE MONOTONICITY USING ISOTONIC REGRESSION
# =========================================================================

cat("=== ENFORCING QUANTILE MONOTONICITY (ISOTONIC REGRESSION) ===\n")

all_quantiles_log <- all_quantiles_log %>%
   group_by(week, type) %>%
   arrange(quantile) %>%
   mutate(
      original_value = predicted_value,
      predicted_value = isoreg(predicted_value)$yf
   ) %>%
   ungroup()

cat("Isotonic regression applied to all forecast horizons\n\n")

# Verify non-negativity
cat("=== NON-NEGATIVITY CHECK ===\n")
cat("Min predicted value:", round(min(all_quantiles_log$predicted_value), 2), "\n")
cat("All predictions >= 0:", all(all_quantiles_log$predicted_value >= 0), "\n")
cat("Number of negative predictions:", sum(all_quantiles_log$predicted_value < 0), "\n\n")

# Verify monotonicity
monotonicity_check <- all_quantiles_log %>%
   group_by(week, type) %>%
   summarise(
      is_monotonic = all(diff(predicted_value) >= -1e-10),
      .groups = "drop"
   )

cat("=== MONOTONICITY VERIFICATION ===\n")
print(monotonicity_check)

if (all(monotonicity_check$is_monotonic)) {
   cat("\n✓ All quantiles are properly ordered!\n\n")
} else {
   cat("\n✗ Warning: Some quantiles are still not monotonic\n\n")
}

# Clean up
all_quantiles_log <- all_quantiles_log %>%
   select(-original_value, -log_predicted)

# =========================================================================
# 8. DISPLAY ALL QUANTILE FORECASTS
# =========================================================================
cat("=== ALL QUANTILE FORECASTS (NO-LAG MODEL) ===\n")
print(all_quantiles_log)

# Create wide format
quantiles_wide <- all_quantiles_log %>%
   pivot_wider(names_from = quantile, values_from = predicted_value, names_prefix = "q_")

cat("\n=== QUANTILE FORECASTS (WIDE FORMAT) ===\n")
print(quantiles_wide)

# =========================================================================
# 9. SUMMARY STATISTICS BY FORECAST HORIZON
# =========================================================================
cat("\n=== SUMMARY BY FORECAST HORIZON ===\n")
summary_by_horizon <- all_quantiles_log %>%
   group_by(type) %>%
   summarise(
      week = first(week),
      q01 = predicted_value[quantile == 0.01],
      q025 = predicted_value[quantile == 0.025],
      q25 = predicted_value[quantile == 0.25],
      median = predicted_value[quantile == 0.5],
      q75 = predicted_value[quantile == 0.75],
      q975 = predicted_value[quantile == 0.975],
      q99 = predicted_value[quantile == 0.99],
      .groups = "drop"
   )

print(summary_by_horizon)

# =========================================================================
# 10. CREATE CDC SUBMISSION FORMAT
# =========================================================================
cat("\n=== CREATING CDC SUBMISSION FORMAT ===\n")

latest_data_week <- latest_data$week
reference_date <- ceiling_date(latest_data_week, "week", week_start = 6) 

cdc_submission <- all_quantiles_log %>%
   mutate(
      reference_date = reference_date,
      target = "wk inc flu hosp",
      horizon = case_when(
         type == "Nowcast" ~ 0,
         type == "1-Wk Ahead" ~ 1,
         type == "2-Wk Ahead" ~ 2,
         type == "3-Wk Ahead" ~ 3
      ),
      location = "45",
      target_end_date = reference_date + days(horizon * 7),
      output_type = "quantile",
      output_type_id = quantile,
      value = round(predicted_value, 0)  # Round to integer
   ) %>%
   select(reference_date, target, horizon, location, target_end_date, 
          output_type, output_type_id, value) %>%
   arrange(horizon, output_type_id)

# Final verification
cat("\n=== CDC SUBMISSION VERIFICATION ===\n")
cat("Min value:", min(cdc_submission$value), "\n")
cat("Max value:", max(cdc_submission$value), "\n")
cat("Any negative values:", any(cdc_submission$value < 0), "\n")
cat("Total rows:", nrow(cdc_submission), "\n\n")

# Export
output_filename <- paste0("/Users/tanvirahammed/Downloads/", 
                          format(reference_date, "%Y-%m-%d"), "-DMAPRIME-QR-NOLAG.csv")
write_csv(cdc_submission, output_filename)
cat("Submission saved to:", output_filename, "\n\n")

# =========================================================================
# 11. VISUALIZATION
# =========================================================================

# Prepare historical data
hist_plot <- test_df %>%
   select(week, Observed = all_of(target_var), Nowcast = Pred_H0) %>%
   pivot_longer(cols = c(Observed, Nowcast), names_to = "Series", values_to = "Admissions")

# Prepare forecast intervals
plot_intervals <- cdc_submission %>%
   filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
   select(target_end_date, output_type_id, value) %>%
   pivot_wider(names_from = output_type_id, values_from = value) %>%
   rename(median = `0.5`, lower_95 = `0.025`, upper_95 = `0.975`, 
          lower_50 = `0.25`, upper_50 = `0.75`) %>%
   arrange(target_end_date)

all_actual_dates <- sort(unique(c(hist_plot$week, plot_intervals$target_end_date)))

# Create plot
p1 <- ggplot() +
   # Forecast uncertainty ribbons
   geom_ribbon(data = plot_intervals, 
               aes(x = target_end_date, ymin = lower_95, ymax = upper_95), 
               fill = "blue", alpha = 0.15) +
   geom_ribbon(data = plot_intervals, 
               aes(x = target_end_date, ymin = lower_50, ymax = upper_50), 
               fill = "blue", alpha = 0.3) +
   
   # Historical & Nowcast data
   geom_line(data = hist_plot, aes(x = week, y = Admissions, color = Series), size = 1) +
   geom_point(data = hist_plot, aes(x = week, y = Admissions, color = Series)) +
   
   # Forecast median line
   geom_line(data = plot_intervals, 
             aes(x = target_end_date, y = median, color = "Forecast"), 
             linetype = "dashed", size = 1) +
   geom_point(data = plot_intervals, 
              aes(x = target_end_date, y = median, color = "Forecast"), 
              size = 3) +
   
   # Add 0 reference line
   geom_hline(yintercept = 0, linetype = "dotted", color = "red", alpha = 0.5) +
   
   scale_x_date(breaks = all_actual_dates, date_labels = "%b %d") +
   scale_color_manual(values = c("Observed" = "black", "Nowcast" = "#E41A1C", "Forecast" = "blue")) +
   labs(title = "South Carolina Flu Hospitalizations: QR Forecast (No Lags)",
        subtitle = "Model relies on EHR data and momentum with weighted recent trends",
        x = "Week Ending Date", y = "Admissions", color = "Series") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
         legend.position = "bottom")

print(p1)

# Export in CDC format
output_filename <- paste0(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/CDC_submission/",
   format(reference_date, "%Y-%m-%d"),
   "-DMAPRIME-QR.csv"
)
write_csv(cdc_submission, output_filename)
cat("\nSaved:", output_filename, "\n")
