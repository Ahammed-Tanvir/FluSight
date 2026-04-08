library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(readr)
library(quantreg) 
library(tidyr)

# ══════════════════════════════════════════════════════════════════════════════════════
# ADAPTIVE PHASE-AWARE INFLUENZA HOSPITALIZATION FORECASTING PIPELINE
# ══════════════════════════════════════════════════════════════════════════════════════

# ── CONFIGURATION ─────────────────────────────────────────────────────────────────────
USE_BOUNDS <- "adaptive"  

# 1. LOAD & PREP DATA ==================================================================
CDC <- read_csv("/Users/tanvirahammed/Downloads/target-hospital-admissions-12.csv") %>%
   filter(location == "45") %>%
   dplyr::select(Week.Ending.Date = date, Total.Influenza.Admissions = value)

MUSC_Weekly_Influenza    <- read_csv("/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Weekly_Influenza_State_dx_cond_lab_Incident.csv")
Prisma_Weekly_Influenza  <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Weekly_Influenza_State_dx_cond_lab_Incident.csv')
RFA_data_Influenza       <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2017-2025/Weekly data/RFA_weekly_influenza_region_incident.csv')

prep_source <- function(df, date_col, prefix) {
   df %>%
      rename(week = !!rlang::sym(date_col)) %>%
      rename_with(~ paste0(prefix, .x), .cols = -week)
}

sources <- list(
   cdc    = list(df = CDC,                    date = "Week.Ending.Date", prefix = "cdc_"),
   musc   = list(df = MUSC_Weekly_Influenza,  date = "Week",             prefix = "musc_"),
   prisma = list(df = Prisma_Weekly_Influenza, date = "Week",             prefix = "prisma_")
)

flu <- reduce(imap(sources, ~ prep_source(.x$df, .x$date, .x$prefix)),
              full_join, by = "week") %>%
   arrange(week)


# 2. FEATURE ENGINEERING ===============================================================
target_var <- "cdc_Total.Influenza.Admissions"

flu_features <- flu %>%
   dplyr::select(week, all_of(target_var),
                 prisma_Weekly_Inpatient_Hospitalizations,
                 prisma_Weekly_Positive_Tests,
                 prisma_Weekly_Tests) %>%
   mutate(
      log_target         = log(!!rlang::sym(target_var) + 1),
      log_tests          = log(prisma_Weekly_Tests + 1),
      log_positive_tests = log(prisma_Weekly_Positive_Tests + 1),
      log_hosp           = log(prisma_Weekly_Inpatient_Hospitalizations + 1),
      log_target_lag1    = lag(log_target, 1),
      log_target_lag2    = lag(log_target, 2),
      log_target_diff    = log_target_lag1 - log_target_lag2,
      log_hosp_diff      = log_hosp - lag(log_hosp, 1),
      log_pos_test_diff  = log_positive_tests - lag(log_positive_tests, 1),
      log_hosp_2wk_diff  = log_hosp - lag(log_hosp, 2),
      month        = month(week),
      week_no      = ceiling(day(week) / 7),
      is_christmas = if_else(month == 12 & day(week) >= 20, 1, 0),
      is_newyear   = if_else(month == 1  & day(week) <= 7,  1, 0),
      school_break = if_else(month %in% c(6, 7, 8, 12), 1, 0),
      orig_tests          = prisma_Weekly_Tests,
      orig_positive_tests = prisma_Weekly_Positive_Tests,
      orig_hosp           = prisma_Weekly_Inpatient_Hospitalizations
   )

predictors <- c(
   "log_tests", "log_positive_tests", "log_hosp",
   "log_hosp_diff", "log_pos_test_diff", "log_hosp_2wk_diff",
   "log_target_lag1", "log_target_diff"
)


# 3. SPLIT & MODEL =====================================================================
train_cut <- as.Date("2022-09-01")
valid_end <- as.Date("2025-10-01")

train_df       <- flu_features %>% filter(week < train_cut)
valid_df       <- flu_features %>% filter(week >= train_cut & week < valid_end)
test_df        <- flu_features %>% filter(week >= valid_end)
final_train_df <- rbind(train_df, valid_df) %>%
   drop_na(all_of(c("log_target", predictors)))

model_median <- rq(
   as.formula(paste("log_target ~", paste(predictors, collapse = " + "))),
   tau  = 0.5,
   data = final_train_df
)

final_train_df$predicted_log <- predict(model_median, newdata = final_train_df)
final_train_df$residual       <- final_train_df$log_target - final_train_df$predicted_log
residual_sd <- sd(final_train_df$residual, na.rm = TRUE)

quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
               0.5,
               0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)


# 4. ENHANCED PHASE-AWARE TREND SYSTEM =================================================

latest_observed_week <- max(flu_features$week[!is.na(flu_features[[target_var]])], na.rm = TRUE)
observed_value       <- flu_features[[target_var]][flu_features$week == latest_observed_week]

recent_weeks <- flu_features %>%
   filter(week <= latest_observed_week) %>%
   arrange(desc(week)) %>%
   dplyr::slice(1:4)

weights <- exp(seq(0, -1.5, length.out = 4)) / sum(exp(seq(0, -1.5, length.out = 4)))

calc_weighted_change <- function(values, wts) {
   changes <- diff(rev(values))
   if (length(changes) == 0) return(0)
   weighted.mean(changes, wts[1:length(changes)], na.rm = TRUE)
}

avg_hosp_change <- calc_weighted_change(recent_weeks$orig_hosp,          weights)
avg_test_change <- calc_weighted_change(recent_weeks$orig_tests,          weights)
avg_pos_change  <- calc_weighted_change(recent_weeks$orig_positive_tests, weights)

recent_cdc_log <- flu_features %>%
   filter(!is.na(log_target)) %>%
   arrange(desc(week)) %>%
   dplyr::slice(1:4) %>%
   arrange(week)

cdc_log_changes <- diff(recent_cdc_log$log_target)
decay_weights   <- exp(seq(-1.5, 0, length.out = length(cdc_log_changes)))
decay_weights   <- decay_weights / sum(decay_weights)
raw_weekly_log_trend <- weighted.mean(cdc_log_changes, decay_weights)

recent_4wk_values <- recent_cdc_log$log_target
recent_8wk <- flu_features %>%
   filter(!is.na(log_target)) %>%
   arrange(desc(week)) %>%
   dplyr::slice(1:min(8, n())) %>%
   pull(log_target)

trend_consistency <- if (length(cdc_log_changes) >= 2) {
   signs <- sign(cdc_log_changes)
   sum(signs == sign(raw_weekly_log_trend)) / length(signs)
} else {
   0.5
}

trend_volatility <- sd(cdc_log_changes, na.rm = TRUE)

is_at_peak <- recent_4wk_values[4] == max(recent_4wk_values, na.rm = TRUE) &&
   recent_4wk_values[4] > quantile(recent_8wk, 0.75, na.rm = TRUE)

current_magnitude <- observed_value / max(exp(recent_8wk) - 1, na.rm = TRUE)

phase <- "PLATEAU"  

if (abs(raw_weekly_log_trend) < 0.05 && trend_volatility < 0.15) {
   phase <- "PLATEAU"
} else if (raw_weekly_log_trend > 0.10 && trend_consistency > 0.6) {
   phase <- "SURGE"
} else if (is_at_peak && abs(raw_weekly_log_trend) < 0.10) {
   phase <- "PEAK"
} else if (raw_weekly_log_trend < -0.10 && trend_consistency > 0.6) {
   if (current_magnitude > 0.5) {
      phase <- "POST_PEAK_DECLINE"
   } else {
      phase <- "TAIL_DECLINE"
   }
} else if (raw_weekly_log_trend > 0.05) {
   phase <- "GRADUAL_INCREASE"
} else if (raw_weekly_log_trend < -0.05) {
   phase <- "GRADUAL_DECLINE"
}

historical_trends <- flu_features %>%
   filter(!is.na(log_target)) %>%
   arrange(desc(week)) %>%
   dplyr::slice(1:min(12, n())) %>%
   arrange(week) %>%
   mutate(log_change = log_target - lag(log_target)) %>%
   filter(!is.na(log_change))

observed_max_growth  <- max(historical_trends$log_change, na.rm = TRUE)
observed_max_decline <- min(historical_trends$log_change, na.rm = TRUE)
observed_volatility  <- sd(historical_trends$log_change, na.rm = TRUE)
safety_margin <- 2 * observed_volatility  

if (USE_BOUNDS == "none") {
   MAX_WEEKLY_SURGE   <- Inf
   MAX_WEEKLY_DECLINE <- -Inf
   bounds_applied     <- FALSE
} else if (USE_BOUNDS == "safety") {
   MAX_WEEKLY_SURGE   <- log(1.50)
   MAX_WEEKLY_DECLINE <- log(0.50)
   bounds_applied     <- TRUE
} else {  
   bounds_applied <- TRUE
   if (phase == "SURGE") {
      MAX_WEEKLY_SURGE   <- max(raw_weekly_log_trend + safety_margin, observed_max_growth, log(1.30))
      MAX_WEEKLY_DECLINE <- log(0.90)
      damping_schedule   <- c(1.0, 0.85, 0.70, 0.60)
      interval_compression <- c(0.90, 0.75)
      use_surge_logic    <- TRUE
   } else if (phase == "PEAK") {
      MAX_WEEKLY_SURGE   <- max(safety_margin, log(1.25))
      MAX_WEEKLY_DECLINE <- min(-safety_margin, log(0.75))
      damping_schedule   <- c(1.0, 0.70, 0.50, 0.35)
      interval_compression <- c(0.85, 0.65)
      use_surge_logic    <- FALSE
   } else if (phase %in% c("POST_PEAK_DECLINE", "TAIL_DECLINE", "GRADUAL_DECLINE")) {
      MAX_WEEKLY_DECLINE <- min(raw_weekly_log_trend - safety_margin, observed_max_decline - (0.5 * observed_volatility), log(0.50))
      MAX_WEEKLY_SURGE   <- log(1.15)
      if (phase == "TAIL_DECLINE") {
         damping_schedule <- c(1.0, 0.95, 0.85, 0.75)
      } else if (abs(raw_weekly_log_trend) > 0.30) {
         damping_schedule <- c(1.0, 0.90, 0.80, 0.70)
      } else {
         damping_schedule <- c(1.0, 0.85, 0.70, 0.60)
      }
      interval_compression <- c(0.92, 0.75)
      use_surge_logic    <- FALSE
   } else if (phase == "GRADUAL_INCREASE") {
      MAX_WEEKLY_SURGE   <- max(raw_weekly_log_trend + safety_margin, log(1.25))
      MAX_WEEKLY_DECLINE <- log(0.85)
      damping_schedule   <- c(1.0, 0.85, 0.70, 0.60)
      interval_compression <- c(0.88, 0.70)
      use_surge_logic    <- FALSE
   } else {  
      MAX_WEEKLY_SURGE   <- safety_margin
      MAX_WEEKLY_DECLINE <- -safety_margin
      damping_schedule   <- c(1.0, 0.70, 0.50, 0.35)
      interval_compression <- c(0.85, 0.65)
      use_surge_logic    <- FALSE
   }
}

if (!exists("damping_schedule")) damping_schedule <- c(1.0, 0.85, 0.70, 0.60)
if (!exists("interval_compression")) interval_compression <- c(0.90, 0.75)
if (!exists("use_surge_logic")) use_surge_logic <- (phase == "SURGE" && raw_weekly_log_trend > 0.10)

weekly_log_trend <- pmax(pmin(raw_weekly_log_trend, MAX_WEEKLY_SURGE), MAX_WEEKLY_DECLINE)

adaptive_bound_info <- list(
   observed_max_growth  = observed_max_growth,
   observed_max_decline = observed_max_decline,
   observed_volatility  = observed_volatility,
   safety_margin        = safety_margin,
   upper_bound          = MAX_WEEKLY_SURGE,
   lower_bound          = MAX_WEEKLY_DECLINE,
   bounds_mode          = USE_BOUNDS,
   bounds_applied       = bounds_applied && abs(raw_weekly_log_trend - weekly_log_trend) > 0.001
)

is_declining <- weekly_log_trend < -0.05
is_surging   <- use_surge_logic && weekly_log_trend > 0.10

# ── DIAGNOSTICS ───────────────────────────────────────────────────────────────────────
cat(sprintf("\n╔══════════════════════════════════════════════════════════════════╗\n"))
cat(sprintf("║          ADAPTIVE PHASE-AWARE TREND DIAGNOSTICS                  ║\n"))
cat(sprintf("║          Bounds Mode: %-44s ║\n", toupper(USE_BOUNDS)))
cat(sprintf("╚══════════════════════════════════════════════════════════════════╝\n\n"))
cat(sprintf("📊 DETECTED PHASE: %s\n", phase))
cat(sprintf("📈 OBSERVED HISTORICAL BEHAVIOR (last 12 weeks):\n"))
cat(sprintf("📉 CURRENT TREND ANALYSIS:\n"))
cat(sprintf("   ├─ Applied trend:        %+.4f  (%+.1f%% per week)\n", weekly_log_trend, (exp(weekly_log_trend) - 1) * 100))
cat(sprintf("🎯 CURRENT STATE:\n"))
cat(sprintf("   ├─ Latest observation: %s (%d admissions)\n", format(latest_observed_week, "%Y-%m-%d"), as.integer(observed_value)))
cat(sprintf("════════════════════════════════════════════════════════════════════\n\n"))


# 5. RECURSIVE FORECAST WITH PHASE-AWARE CORRECTION ===================================
current_row <- flu_features %>% filter(week == latest_observed_week)
last_two <- flu_features %>% filter(!is.na(log_target)) %>% arrange(desc(week)) %>% dplyr::slice(1:2)

current_row$log_target_lag1 <- last_two$log_target[1]
current_row$log_target_lag2 <- last_two$log_target[2]
current_row$log_target_diff  <- current_row$log_target_lag1 - current_row$log_target_lag2
current_row$week             <- latest_observed_week + weeks(1)

all_quantiles_log <- data.frame()
week_types         <- c("Nowcast", "1-Wk Ahead", "2-Wk Ahead", "3-Wk Ahead")

for (i in 1:4) {
   if (i > 1) {
      prev_log_hosp <- current_row$log_hosp
      if (is_surging) {
         current_row$orig_tests           <- pmax(current_row$orig_tests           + (avg_test_change * 1.5), 1)
         current_row$orig_positive_tests <- pmax(current_row$orig_positive_tests + (avg_pos_change  * 3.0), 1)
         current_row$orig_hosp            <- pmax(current_row$orig_hosp, observed_value)
      }
      current_row$log_tests           <- log(current_row$orig_tests           + 1)
      current_row$log_positive_tests <- log(current_row$orig_positive_tests + 1)
      current_row$log_hosp            <- log(current_row$orig_hosp           + 1)
      current_row$log_hosp_diff       <- current_row$log_hosp - prev_log_hosp
   }
   
   raw_median_log <- as.numeric(predict(model_median, newdata = current_row))
   
   if ((is_declining || phase %in% c("GRADUAL_DECLINE", "POST_PEAK_DECLINE", "TAIL_DECLINE")) && !is_surging) {
      trend_correction     <- weekly_log_trend * damping_schedule[i]
      corrected_median_log <- raw_median_log + trend_correction
   } else if (phase %in% c("SURGE", "GRADUAL_INCREASE")) {
      trend_correction     <- weekly_log_trend * damping_schedule[i]
      corrected_median_log <- raw_median_log + trend_correction
   } else {
      corrected_median_log <- raw_median_log
   }
   
   if (is_surging) {
      surge_floor           <- log(observed_value * (1 + (0.05 * i)) + 1)
      corrected_median_log <- pmax(corrected_median_log, surge_floor)
   }
   
   corrected_median_log <- pmax(corrected_median_log, 0)
   horizon_sd    <- residual_sd * c(0.7, 0.85, 1.0, 1.15)[i]
   log_quantiles <- qnorm(quantiles, mean = corrected_median_log, sd = horizon_sd)
   
   for (q_idx in seq_along(quantiles)) {
      if (quantiles[q_idx] > 0.5) {
         dist <- log_quantiles[q_idx] - corrected_median_log
         mult <- if (quantiles[q_idx] <= 0.75) interval_compression[1] else interval_compression[2]
         log_quantiles[q_idx] <- corrected_median_log + (dist * mult)
      }
   }
   
   vals <- pmax(exp(log_quantiles) - 1, 0)
   all_quantiles_log <- bind_rows(all_quantiles_log, data.frame(
      week            = current_row$week,
      type            = week_types[i],
      quantile        = quantiles,
      predicted_value = vals
   ))
   
   prev_lag1                    <- current_row$log_target_lag1
   current_row$log_target_lag2  <- prev_lag1
   current_row$log_target_lag1  <- corrected_median_log
   current_row$log_target_diff  <- current_row$log_target_lag1 - current_row$log_target_lag2
   current_row$week <- current_row$week + weeks(1)
}


# 6. MONOTONICITY & CDC SUBMISSION =====================================================
all_quantiles_log <- all_quantiles_log %>%
   group_by(week) %>%
   mutate(raw_pre_isoreg = predicted_value,
          predicted_value = isoreg(predicted_value)$yf) %>%
   ungroup()

cdc_submission <- all_quantiles_log %>%
   mutate(
      reference_date  = latest_observed_week + weeks(1),
      target          = "wk inc flu hosp",
      horizon         = case_when(
         type == "Nowcast"     ~ 0,
         type == "1-Wk Ahead"  ~ 1,
         type == "2-Wk Ahead"  ~ 2,
         type == "3-Wk Ahead"  ~ 3
      ),
      location        = "45",
      target_end_date = reference_date + days(horizon * 7),
      output_type     = "quantile",
      output_type_id  = quantile,
      value           = round(predicted_value, 0)
   ) %>%
   dplyr::select(reference_date, target, horizon, location, target_end_date,
                 output_type, output_type_id, value) %>%
   arrange(horizon, output_type_id)

cat("\n═══ FORECAST MEDIANS BY HORIZON ═══\n")
print(cdc_submission %>% filter(output_type_id == 0.5) %>% dplyr::select(horizon, target_end_date, median_forecast = value))


# 7. VISUALIZATION =====================================================================
plot_intervals <- cdc_submission %>%
   filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
   pivot_wider(names_from = output_type_id, values_from = value) %>%
   rename(median   = `0.5`,
          lower_95 = `0.025`, upper_95 = `0.975`,
          lower_50 = `0.25`,  upper_50 = `0.75`)

hist_plot_df <- test_df %>%
   filter(week >= as.Date("2025-10-01") & week <= latest_observed_week) %>%
   mutate(Nowcast = exp(predict(model_median, newdata = .)) - 1) %>%
   dplyr::select(week, Observed = all_of(target_var), Nowcast) %>%
   pivot_longer(cols = c(Observed, Nowcast), names_to = "Series", values_to = "Admissions")

all_actual_dates <- sort(unique(c(hist_plot_df$week, plot_intervals$target_end_date)))

cap_note <- if (USE_BOUNDS == "none") {
   sprintf("trend %+.1f%%/wk (NO BOUNDS - trusting raw data)", (exp(weekly_log_trend)-1)*100)
} else if (adaptive_bound_info$bounds_applied) {
   sprintf("raw %+.0f%%/wk → bounded at %+.0f%%/wk (limit: %+.0f%% to %+.0f%%)",
           (exp(raw_weekly_log_trend)-1)*100, 
           (exp(weekly_log_trend)-1)*100,
           (exp(adaptive_bound_info$lower_bound)-1)*100,
           (exp(adaptive_bound_info$upper_bound)-1)*100)
} else {
   sprintf("trend %+.1f%%/wk (within %s bounds)", 
           (exp(weekly_log_trend)-1)*100,
           USE_BOUNDS)
}

# Color scheme based on phase
phase_color <- switch(phase,
                      "SURGE"              = "#E63946",  # Red
                      "PEAK"               = "#F77F00",  # Orange
                      "POST_PEAK_DECLINE"  = "#FCBF49",  # Yellow-orange
                      "TAIL_DECLINE"       = "#06A77D",  # Green
                      "GRADUAL_INCREASE"   = "#E76F51",  # Coral
                      "GRADUAL_DECLINE"    = "#2A9D8F",  # Teal
                      "PLATEAU"            = "#457B9D"   # Blue
)

p1 <- ggplot() +
   geom_ribbon(data = plot_intervals,
               aes(x = target_end_date, ymin = lower_95, ymax = upper_95),
               fill = phase_color, alpha = 0.15) +
   geom_ribbon(data = plot_intervals,
               aes(x = target_end_date, ymin = lower_50, ymax = upper_50),
               fill = phase_color, alpha = 0.30) +
   geom_line(data  = hist_plot_df,
             aes(x = week, y = Admissions, color = Series), linewidth = 1) +
   geom_point(data = hist_plot_df,
              aes(x = week, y = Admissions, color = Series), size = 2) +
   geom_line(data  = plot_intervals,
             aes(x = target_end_date, y = median, color = "Forecast"),
             linetype = "dashed", linewidth = 1.2) +
   geom_point(data = plot_intervals,
              aes(x = target_end_date, y = median, color = "Forecast"), size = 3) +
   geom_hline(yintercept = 0, linetype = "dotted", color = "red", alpha = 0.5) +
   scale_x_date(breaks = all_actual_dates, date_labels = "%b %d") +
   scale_color_manual(values = c("Observed" = "black",
                                 "Nowcast"  = "#E41A1C",
                                 "Forecast" = phase_color)) +
   labs(
      title    = sprintf("Phase-Aware QR Forecast [Phase: %s]", phase),
      subtitle = paste0("Latest: ", format(latest_observed_week, "%Y-%m-%d"),
                        " (", observed_value, " admissions)  |  ", cap_note,
                        "  |  Damping: ", paste(sprintf("%.0f", damping_schedule*100), collapse = "→"), "%"),
      y        = "Weekly Influenza Admissions",
      x        = "Week Ending Date",
      caption  = sprintf("Model: Quantile Regression (τ=0.5) | Residual SD: %.2f | Phase Detection: Automated", residual_sd)
   ) +
   theme_minimal() +
   theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y   = element_text(size = 9),
      axis.title    = element_text(size = 11),
      legend.position = "bottom",
      legend.title    = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(hjust = 0, size = 8, color = "gray40")
   )

print(p1)


# 8. EXPORT CDC FILE ===================================================================
reference_date  <- cdc_submission$reference_date[1]
output_filename_cdc <- paste0("/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/CDC_submission/", format(reference_date, "%Y-%m-%d"), "-DMAPRIME-QR.csv")

#write_csv(cdc_submission, output_filename_cdc)










# ══════════════════════════════════════════════════════════════════════════════
# 9. GENERATE SW IMPLEMENTATION FILE
# ══════════════════════════════════════════════════════════════════════════════

sw_implementation_file <- cdc_submission %>%
   # Filter to include only the median (output_type_id == 0.5) if required by your SW,
   # or keep all quantiles if the software accepts full distributions. 
   # Based on the typical 'implementation' request, we use the median:
   mutate(
      # reference_date: Saturday of the current week 
      reference_date   = format(reference_date, "%Y-%m-%d"), 
      
      # target: Fixed character string 
      target           = "wk inc flu hosp",           
      
      # target_end_date: Future dates already in ISO format from cdc_submission [cite: 16]
      target_end_date  = format(target_end_date, "%Y-%m-%d"),    
      
      # location_general: Fixed as "state" for this specific model [cite: 19]
      location_general = "state",                    
      
      # location: "SC" for the state-level model 
      location         = "SC",      
      
      # disease: Fixed character string [cite: 35]
      disease          = "influenza",                 
      
      # population: Fixed character string [cite: 38]
      population       = "general_population",        
      
      # training_validation: Set to 0 for implementation/testing 
      training_validation = 0,                        
      
      # estimate_projected_report: Set to 1 as per your requirement [cite: 50, 52]
      estimate_projected_report = 1,                  
      
      # imputed: Numeric 0 (no imputation applied) [cite: 54, 58]
      imputed          = 0,                           
      
      # data_source: Source of observed data (CDC) 
      data_source      = NA_character_,              
      
      # outcome_measure: Fixed label [cite: 70]
      outcome_measure  = "Weekly_Inpatient_Hospitalizations", 
      
      # output_type: Fixed as quantile [cite: 74]
      output_type      = "quantile"                        
   ) %>%
   # Select and order the 15 required columns [cite: 1-81]
   dplyr::select(
      reference_date, target, target_end_date, location_general, location, 
      value, disease, population, training_validation, 
      estimate_projected_report, imputed, data_source, 
      outcome_measure, output_type, output_type_id
   )

# Software Submission
# output_filename_software <- paste0(
#    "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/software.submission/",
#    format(reference_date, "%Y-%m-%d"),
#    "-Ahammed-Tanvir-17-state-flu-implementation.csv"
# )




output_filename_software <- paste0(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Forecasting Resources/Forecast-Drop-Off/Software/Implementation/",
   format(reference_date, "%Y-%m-%d"),
   "-inpatient-hosp-DMAPRIME-QR.csv"
)





# Optional: Export files (uncomment to save)
#write_csv(sw_implementation_file, output_filename_software)


cat("\n✅ PIPELINE COMPLETE | CDC File: ", basename(output_filename_cdc), "\n")
