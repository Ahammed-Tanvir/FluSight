library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(readr)
library(quantreg) 
library(tidyr)

# 1. LOAD & PREP DATA ==================================================================
CDC <- read_csv("/Users/tanvirahammed/Downloads/target-hospital-admissions-3.csv") %>%
   filter(location == "45") %>%
   select(Week.Ending.Date = date, Total.Influenza.Admissions = value)

MUSC_Weekly_Influenza <- read_csv("/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Weekly_Influenza_State_dx_cond_lab_Incident.csv")
Prisma_Weekly_Influenza <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Weekly_Influenza_State_dx_cond_lab_Incident.csv')

prep_source <- function(df, date_col, prefix) {
   df %>% rename(week = !!rlang::sym(date_col)) %>% rename_with(~ paste0(prefix, .x), .cols = -week)
}

sources <- list(cdc = list(df = CDC, date = "Week.Ending.Date", prefix = "cdc_"),
                musc = list(df = MUSC_Weekly_Influenza, date = "Week", prefix = "musc_"),
                prisma = list(df = Prisma_Weekly_Influenza, date = "Week", prefix = "prisma_"))

flu <- reduce(imap(sources, ~ prep_source(.x$df, .x$date, .x$prefix)), full_join, by = "week") %>% arrange(week)

# 2. FEATURE ENGINEERING ===============================================================
target_var <- "cdc_Total.Influenza.Admissions"

flu_features <- flu %>%
   select(week, all_of(target_var), prisma_Weekly_Inpatient_Hospitalizations, 
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

quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
               0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

model_log_multi <- rq(as.formula(paste("log_target ~", paste(predictors, collapse = " + "))), 
                      tau = quantiles, data = final_train_df)

# 4. TREND CALCULATION =================================================================
latest_observed_week <- max(flu_features$week[!is.na(flu_features[[target_var]])], na.rm = TRUE)
observed_value <- flu_features[[target_var]][flu_features$week == latest_observed_week]

recent_weeks <- flu_features %>% filter(week <= latest_observed_week) %>% arrange(desc(week)) %>% dplyr::slice(1:4)
weights <- exp(seq(0, -1.5, length.out = 4)) / sum(exp(seq(0, -1.5, length.out = 4)))

calc_weighted_change <- function(values, wts) {
   changes <- diff(rev(values)); if(length(changes) == 0) return(0)
   weighted.mean(changes, wts[1:length(changes)], na.rm = TRUE)
}

avg_test_change <- calc_weighted_change(recent_weeks$orig_tests, weights)
avg_pos_change  <- calc_weighted_change(recent_weeks$orig_positive_tests, weights)
avg_hosp_change <- calc_weighted_change(recent_weeks$orig_hosp, weights)

# 5. SMART RECURSIVE FORECAST (WITH SAFETY VALVE) ======================================
current_row <- flu_features %>% filter(week == latest_observed_week)
current_row$log_target_lag1 <- log(observed_value + 1)
current_row$week <- latest_observed_week + weeks(1) 

all_quantiles_log <- data.frame()
week_types <- c("Nowcast", "1-Wk Ahead", "2-Wk Ahead", "3-Wk Ahead")

# Stability Settings
damping <- c(1.0, 0.6, 0.3, 0.1) 
upper_bound_squeeze <- c(1.0, 0.4, 0.15, 0.05) 
max_seasonal_val <- max(flu_features[[target_var]][flu_features$week >= as.Date("2025-10-01")], na.rm = TRUE)* 1.2

for (i in 1:4) {
   if(i > 1) {
      prev_log_hosp <- current_row$log_hosp
      
      # SAFETY VALVE: If EHR data shows an actual increase, reduce damping to react
      active_damping <- if(avg_hosp_change > 0) 1.0 else damping[i]
      
      current_row$orig_hosp <- pmax(current_row$orig_hosp + (avg_hosp_change * active_damping), 1)
      current_row$log_hosp <- log(current_row$orig_hosp + 1)
      
      # THRESHOLD TRIGGER: Only allow positive momentum if weighted EHR trend is positive
      raw_diff <- current_row$log_hosp - prev_log_hosp
      if (avg_hosp_change > 0) {
         current_row$log_hosp_diff <- raw_diff 
      } else {
         current_row$log_hosp_diff <- min(raw_diff, 0) # Clip growth during decline
      }
      
      current_row$log_hosp_2wk_diff <- current_row$log_hosp_diff * 1.2
      
      # Standard EHR Feature Updates
      current_row$orig_tests <- pmax(current_row$orig_tests + (avg_test_change * active_damping), 1)
      current_row$orig_positive_tests <- pmax(current_row$orig_positive_tests + (avg_pos_change * active_damping), 1)
      current_row$log_tests <- log(current_row$orig_tests + 1)
      current_row$log_positive_tests <- log(current_row$orig_positive_tests + 1)
   }
   
   pred_matrix_log <- predict(model_log_multi, newdata = current_row)
   vals <- exp(as.numeric(pred_matrix_log[1, ])) - 1
   
   # MANUAL UNCERTAINTY SQUEEZE (Upper Quantiles Only)
   median_val <- vals[which(quantiles == 0.5)]
   # Safety: Only squeeze if the model isn't seeing a major EHR-driven rebound
   active_squeeze <- if(avg_hosp_change > 0) 1.0 else upper_bound_squeeze[i]
   
   vals <- ifelse(quantiles > 0.5, 
                  median_val + (vals - median_val) * active_squeeze, 
                  vals)
   
   # Cap and floor
   vals <- pmax(pmin(vals, max_seasonal_val), 0)
   
   quantile_results <- data.frame(
      week = current_row$week, type = week_types[i], quantile = quantiles,
      predicted_value = vals
   )
   
   all_quantiles_log <- bind_rows(all_quantiles_log, quantile_results)
   current_row$log_target_lag1 <- log(median_val + 1)
   current_row$week <- current_row$week + weeks(1)
}

# 6. MONOTONICITY & SUBMISSION =========================================================
all_quantiles_log <- all_quantiles_log %>% 
   group_by(week) %>% 
   mutate(predicted_value = isoreg(predicted_value)$yf) %>%
   ungroup()


cdc_submission <- all_quantiles_log %>%
   mutate(
      reference_date = latest_observed_week + weeks(1),
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


# 7. VISUALIZATION =====================================================================
plot_start_date <- as.Date("2025-10-01")

hist_plot_df <- test_df %>%
   filter(week >= plot_start_date & week <= latest_observed_week) %>%
   mutate(Pred_H0 = exp(as.numeric(predict(model_log_multi, newdata = .)[, which(quantiles == 0.5)])) - 1) %>%
   select(week, Observed = all_of(target_var), Nowcast = Pred_H0) %>%
   pivot_longer(cols = c(Observed, Nowcast), names_to = "Series", values_to = "Admissions")

plot_intervals <- cdc_submission %>%
   filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
   select(target_end_date, output_type_id, value) %>%
   pivot_wider(names_from = output_type_id, values_from = value) %>%
   rename(median = `0.5`, lower_95 = `0.025`, upper_95 = `0.975`, lower_50 = `0.25`, upper_50 = `0.75`)

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
   labs(title = "SC Flu Hospitalizations: QR Forecast") +
        #subtitle = "Safety Valve enabled: reacts to rebounds while stabilizing declines") +
   theme_minimal() + 
   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

print(p1)




reference_date <- cdc_submission$reference_date[1]

# Export in CDC format
output_filename <- paste0(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/CDC_submission/",
   format(reference_date, "%Y-%m-%d"),
   "-DMAPRIME-QR.csv"
)
write_csv(cdc_submission, output_filename)
# cat("\nSaved:", output_filename, "\n")
