# =========================================================================
# 1. LOAD LIBRARIES AND DATA
# =========================================================================
library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(readr)
library(quantreg) 
library(tidyr)

# Load datasets
CDC <- read_csv("/Users/tanvirahammed/Downloads/target-hospital-admissions.csv") %>%
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

# =========================================================================
# 2. FEATURE ENGINEERING
# =========================================================================
target_var <- "cdc_Total.Influenza.Admissions"

flu_features <- flu %>%
   select(week, all_of(target_var), 
          prisma_Weekly_Inpatient_Hospitalizations, prisma_Weekly_Positive_Tests, prisma_Weekly_Tests) %>%
   mutate(
      target_lag1 = lag(!!rlang::sym(target_var), 1),
      target_lag2 = lag(!!rlang::sym(target_var), 2),
      hosp_diff = prisma_Weekly_Inpatient_Hospitalizations - lag(prisma_Weekly_Inpatient_Hospitalizations, 1),
      pos_test_diff = prisma_Weekly_Positive_Tests - lag(prisma_Weekly_Positive_Tests, 1),
      
      year = year(week),
      month = month(week),
      week_no = ceiling(day(week)/7),
      is_christmas = if_else(month == 12 & day(week) >= 20, 1, 0),
      is_newyear = if_else(month == 1 & day(week) <= 7, 1, 0),
      school_break = if_else(month %in% c(6, 7, 8, 12), 1, 0)
   ) 

predictors <- c("prisma_Weekly_Tests", "prisma_Weekly_Positive_Tests", "target_lag1", "target_lag2", 
                "hosp_diff", "pos_test_diff", "year", "month", "week_no", 
                "school_break", "is_christmas", "is_newyear")

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

final_train_df <- rbind(train_df, valid_df) %>% drop_na(all_of(c(target_var, predictors)))

# =========================================================================
# 4. MODEL FITTING (Multi-Quantile)
# =========================================================================
quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
               0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

formula_h0 <- as.formula(paste(target_var, "~", paste(predictors, collapse = " + ")))

# Fit Quantile Regression for all tau values
model_h0_multi <- rq(formula = formula_h0, tau = quantiles, data = final_train_df)

# =========================================================================
# 5. GENERATE TEST PREDICTIONS
# =========================================================================
quantile_preds <- predict(model_h0_multi, newdata = test_df)
test_df$Pred_H0 <- as.numeric(quantile_preds[, "tau= 0.500"])

# =========================================================================
# 6. TREND-BASED RECURSIVE FORECASTING WITH ALL QUANTILES
# =========================================================================
latest_data <- flu_features %>% 
   arrange(desc(week)) %>% 
   dplyr::slice(1)

forecast_horizon <- 3

# Get recent trends for projection (last 4 weeks)
recent_weeks <- flu_features %>% 
   arrange(desc(week)) %>% 
   dplyr::slice(1:4)

# Calculate average weekly changes
avg_test_change <- mean(diff(rev(recent_weeks$prisma_Weekly_Tests)), na.rm = TRUE)
avg_pos_change <- mean(diff(rev(recent_weeks$prisma_Weekly_Positive_Tests)), na.rm = TRUE)
avg_hosp_change <- mean(diff(rev(recent_weeks$prisma_Weekly_Inpatient_Hospitalizations)), na.rm = TRUE)

cat("\n=== TREND-BASED PROJECTION ===\n")
cat("Average weekly change in Tests:", round(avg_test_change, 2), "\n")
cat("Average weekly change in Positive Tests:", round(avg_pos_change, 2), "\n")
cat("Average weekly change in Hospitalizations:", round(avg_hosp_change, 2), "\n\n")

# Initialize results dataframe for all quantiles
all_quantiles_trend <- data.frame()

# Week labels
week_types <- c("Nowcast", "1-Wk Ahead", "2-Wk Ahead", "3-Wk Ahead")
week_dates <- c(latest_data$week, seq(latest_data$week + weeks(1), by = "1 week", length.out = forecast_horizon))

# Initialize recursive state
current_row_trend <- latest_data

for (i in 1:4) {  # Nowcast + 3 forecasts
   # Predict ALL quantiles using the current row state
   pred_matrix <- predict(model_h0_multi, newdata = current_row_trend)
   
   # Extract all quantile predictions for this week
   quantile_results <- data.frame(
      week = week_dates[i],
      type = week_types[i],
      quantile = quantiles,
      predicted_value = as.numeric(pred_matrix[1, ])
   )
   
   # Append to master dataframe
   all_quantiles_trend <- bind_rows(all_quantiles_trend, quantile_results)
   
   # Get median for recursive update
   pred_median <- as.numeric(pred_matrix[, "tau= 0.500"])
   
   if (i < 4) {  # Don't update after last prediction
      # RECURSIVE LOGIC: Update lags
      current_row_trend$target_lag2 <- current_row_trend$target_lag1
      current_row_trend$target_lag1 <- pred_median
      
      # Update Date and Time features
      next_week_date <- current_row_trend$week + weeks(1)
      current_row_trend$week <- next_week_date
      current_row_trend$month <- month(next_week_date)
      current_row_trend$year  <- year(next_week_date)
      current_row_trend$week_no <- ceiling(day(next_week_date)/7)
      
      # Update Seasonality/Holidays
      current_row_trend$is_christmas <- if_else(current_row_trend$month == 12 & day(next_week_date) >= 20, 1, 0)
      current_row_trend$is_newyear   <- if_else(current_row_trend$month == 1 & day(next_week_date) <= 7, 1, 0)
      current_row_trend$school_break <- if_else(current_row_trend$month %in% c(6, 7, 8, 12), 1, 0)
      
      # PROJECT EHR DATA using recent trends
      current_row_trend$prisma_Weekly_Tests <- current_row_trend$prisma_Weekly_Tests + avg_test_change
      current_row_trend$prisma_Weekly_Positive_Tests <- current_row_trend$prisma_Weekly_Positive_Tests + avg_pos_change
      current_row_trend$prisma_Weekly_Inpatient_Hospitalizations <- 
         current_row_trend$prisma_Weekly_Inpatient_Hospitalizations + avg_hosp_change
      
      # Update diffs based on projected changes
      current_row_trend$hosp_diff <- avg_hosp_change
      current_row_trend$pos_test_diff <- avg_pos_change
   }
}

# Create median summary for plotting compatibility
recursive_preds_trend <- all_quantiles_trend %>%
   filter(quantile == 0.5) %>%
   select(week, type, predicted_value)

# =========================================================================
# 7. DISPLAY ALL QUANTILE FORECASTS
# =========================================================================
cat("\n=== ALL QUANTILE FORECASTS (TREND-BASED) ===\n")
print(all_quantiles_trend)

# Create wide format for easier viewing
quantiles_wide <- all_quantiles_trend %>%
   pivot_wider(names_from = quantile, values_from = predicted_value, names_prefix = "q_")

cat("\n=== QUANTILE FORECASTS (WIDE FORMAT) ===\n")
print(quantiles_wide)

# =========================================================================
# 8. VISUALIZATION WITH UNCERTAINTY INTERVALS
# =========================================================================
hist_plot <- test_df %>%
   select(week, Observed = all_of(target_var), Nowcast = Pred_H0) %>%
   pivot_longer(cols = c(Observed, Nowcast), names_to = "Series", values_to = "Admissions")

# Extract key quantiles for plotting
plot_intervals_forecast <- all_quantiles_trend %>%
   select(week, type, quantile, predicted_value) %>%
   pivot_wider(names_from = quantile, values_from = predicted_value) %>%
   rename(
      median = `0.5`,
      lower_95 = `0.025`,
      upper_95 = `0.975`,
      lower_50 = `0.25`,
      upper_50 = `0.75`
   )

p1 <- ggplot() +
   # Historical data
   geom_line(data = hist_plot, aes(x = week, y = Admissions, color = Series), size = 1) +
   geom_point(data = hist_plot, aes(x = week, y = Admissions, color = Series)) +
   
   # Forecast uncertainty ribbons
   geom_ribbon(data = plot_intervals_forecast, 
               aes(x = week, ymin = lower_95, ymax = upper_95), 
               fill = "blue", alpha = 0.2) +
   geom_ribbon(data = plot_intervals_forecast, 
               aes(x = week, ymin = lower_50, ymax = upper_50), 
               fill = "blue", alpha = 0.4) +
   
   # Forecast median line
   geom_line(data = plot_intervals_forecast, 
             aes(x = week, y = median, color = "Forecast"), 
             linetype = "dashed", size = 1) +
   geom_point(data = plot_intervals_forecast, 
              aes(x = week, y = median, color = "Forecast"), 
              size = 3) +
   
   scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
   scale_color_manual(values = c("Observed" = "black", "Nowcast" = "#E41A1C", "Forecast" = "blue")) +
   labs(title = "Trend-Based Recursive Forecast with Uncertainty Intervals",
        subtitle = "50% and 95% prediction intervals | EHR predictors projected using recent trends",
        x = "Week Ending Date", y = "Admissions", color = "Series") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

print(p1)



# =========================================================================
# 10. SUMMARY STATISTICS BY FORECAST HORIZON
# =========================================================================
cat("\n=== SUMMARY BY FORECAST HORIZON ===\n")
summary_by_horizon <- all_quantiles_trend %>%
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
# 11. CREATE CDC SUBMISSION FORMAT =========================================================================
cat("\n=== CREATING CDC SUBMISSION FORMAT ===\n")

# Get the latest data week
latest_data_week <- latest_data$week
reference_date <- ceiling_date(latest_data_week, "week", week_start = 6) 


# Create the submission format dataframe
cdc_submission <- all_quantiles_trend %>%
   mutate(
      reference_date = reference_date,
      target = "wk inc flu hosp",
      horizon = case_when(
         type == "Nowcast" ~ 0,
         type == "1-Wk Ahead" ~ 1,
         type == "2-Wk Ahead" ~ 2,
         type == "3-Wk Ahead" ~ 3
      ),
      location = "45",  # South Carolina FIPS code
      # Correct target_end_date based on horizon
      target_end_date = reference_date + days(horizon * 7),
      output_type = "quantile",
      output_type_id = quantile,
      value = round(predicted_value, 0)  # Round to nearest integer
   ) %>%
   select(reference_date, target, horizon, location, target_end_date, 
          output_type, output_type_id, value) %>%
   arrange(horizon, output_type_id)

# Export in CDC format
output_filename <- paste0(
   "/Users/tanvirahammed/Downloads/",
   format(reference_date, "%Y-%m-%d"),
   "-DMAPRIME-QR.csv"
)
write_csv(cdc_submission, output_filename)
cat("\nSaved:", output_filename, "\n")



# Summary of submission
cat("\n=== SUBMISSION SUMMARY ===\n")
cat("Reference Date:", as.character(current_date), "\n")
cat("Total Rows:", nrow(cdc_submission), "\n")
cat("Horizons:", paste(unique(cdc_submission$horizon), collapse = ", "), "\n")
cat("Quantiles:", length(unique(cdc_submission$output_type_value)), "\n")
cat("Location: South Carolina (FIPS 45)\n")
cat("Target: Weekly incident flu hospitalizations\n\n")

# Verification checks
cat("=== VERIFICATION CHECKS ===\n")
cat("✓ All horizons present:", all(0:3 %in% cdc_submission$horizon), "\n")
cat("✓ All 23 quantiles per horizon:", 
    all(table(cdc_submission$horizon) == 23), "\n")
cat("✓ Quantiles are non-decreasing within each horizon:", 
    all(cdc_submission %>% 
           group_by(horizon) %>% 
           summarise(monotonic = all(diff(value) >= 0)) %>% 
           pull(monotonic)), "\n")
cat("✓ No missing values:", sum(is.na(cdc_submission$value)) == 0, "\n")