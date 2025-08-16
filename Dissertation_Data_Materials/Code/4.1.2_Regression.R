# =============================================================================
# 4.1.2 Panel Regression Screening – Updated R Script
# =============================================================================
# 0. Load required packages
library(tidyverse)    # data manipulation & ggplot2
library(plm)          # panel linear models
library(car)          # vif()
library(e1071)        # skewness()
library(lmtest)       # bptest(), dwtest()
library(sandwich)     # vcovHC()
library(parameters)   # model_parameters()
library(performance)  # check_model()
library(readxl) 
library(tibble)

# 1. Read in the simulated PHS panel dataset
df <- read_excel("C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/数据/Data.xlsx", sheet = 1)
# Expected columns:
# Hospital, Year, Region, TeachingHospitalFlag,
# Beds, StaffFTE,
# OperatingCost_kGBP, CapitalCost_kGBP, MedicationCost_kGBP,
# Discharges, InpatientDays, OutpatientVisits,
# Surgeries, EmergencyVisits,
# ReadmissionRate, MortalityRate, SatisfactionScore, CaseMixIndex

# 2. Compute the key dependent variable: Cost per Discharge
df <- df %>%
  mutate(
    CostPerDischarge   = OperatingCost_kGBP / Discharges,
    ln_CostPerDischarge = log(CostPerDischarge)
  )

# 3. Compute skewness for each continuous predictor as a simple vector
skewness_values <- sapply(df[cont_vars], function(x) {
  e1071::skewness(x, na.rm = TRUE)
})

# 4. Build a tibble of results so filter() works as expected
skew_tbl <- tibble(
  Variable = names(skewness_values),
  Skewness = as.numeric(skewness_values)
)
print(skew_tbl)
write_csv(skew_tbl, "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.2 回归筛选/Table4.1.2_skew_tbl.csv")

# 5. Identify and log‐transform any predictor with |skewness| > 0.5
vars_to_log <- skew_tbl %>% 
  filter(abs(Skewness) > 0.5) %>% 
  pull(Variable)

for (v in vars_to_log) {
  new_var <- paste0("ln_", v)
  df[[new_var]] <- log(df[[v]] + 1)
}

print(vars_to_log)

vars_tbl <- tibble(Variable = vars_to_log)
readr::write_csv(vars_tbl,
                 "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.2 回归筛选/Table4.1.2_vars_to_log.csv"
)
# 6. Define the predictor list for the initial “full” model
#    Use the log-version if created, otherwise the original
predictors <- map_chr(cont_vars, function(v) {
  lv <- paste0("ln_", v)
  if (lv %in% names(df)) lv else v
})

# Add factor variables
predictors <- c(predictors,
                "factor(Region)",
                "factor(TeachingHospitalFlag)",
                "factor(Year)")

# 7. Build and estimate the two‐way fixed‐effects panel model
form_full <- as.formula(
  paste("ln_CostPerDischarge ~", paste(predictors, collapse = " + "))
)
panel_full <- plm(
  form_full,
  data   = df,
  index  = c("Hospital", "Year"),
  model  = "within",
  effect = "twoways"
)

# 8. Multicollinearity check – Variance Inflation Factor
# --- after you’ve defined `form_full` and loaded `df` ---

# 8a. Build an OLS model with intercept (same regressors as the plm)
lm_for_vif <- lm(
  form_full,
  data = df
)

# 8b. Compute VIF on that lm
library(car) 
vif_vals <- vif(lm_for_vif)
print(vif_vals)


vif_df <- as.data.frame(vif_vals) %>%
  rownames_to_column(var = "Variable") %>%
  as_tibble()
out_path <- "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.2 回归筛选/Table4.1.2_vif_vals.csv"
out_dir  <- dirname(out_path)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
readr::write_csv(vif_df, out_path)

# If any vif_vals > 5 (or your chosen threshold),
# drop or merge the corresponding predictor(s), then re-run both lm_for_vif and panel_full.

# Then go back to your plm FE model:
panel_full <- plm(
  form_full,
  data   = df,
  index  = c("Hospital", "Year"),
  model  = "within",
  effect = "twoways"
)


# Define and create the output directory once
base_dir <- "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.2 回归筛选"
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

# 9. Heteroskedasticity test – Breusch–Pagan
bp_test <- bptest(panel_full, studentize = TRUE)
bp_df   <- broom::tidy(bp_test)

readr::write_csv(
  bp_df,
  file.path(base_dir, "Table4.1.2_bp_test.csv")
)

# 10. Autocorrelation test – Durbin–Watson
dw_test <- dwtest(panel_full)
dw_df   <- broom::tidy(dw_test)

readr::write_csv(
  dw_df,
  file.path(base_dir, "Table4.1.2_dw_test.csv")
)

# 11. Residual normality test – Shapiro–Wilk
resids  <- residuals(panel_full)
sw_test <- shapiro.test(resids)

# broom::tidy() also works for shapiro.test objects, but you can construct manually:
sw_df <- tibble(
  statistic = unname(sw_test$statistic),
  p.value   = sw_test$p.value
)

readr::write_csv(
  sw_df,
  file.path(base_dir, "Table4.1.2_sw_test.csv")
)

# 12. Extract coefficients with robust (HC1) standard errors
coefs <- model_parameters(
  panel_full,
  vcov          = "HC1",
  include_fixed = TRUE
)

# model_parameters() returns a tibble-like data frame already
coefs_df <- as_tibble(coefs)

readr::write_csv(
  coefs_df,
  file.path(base_dir, "Table4.1.2_coefs.csv")
)
# =============================================================================
# 13(Revised): Plot Coefficients and 95% Confidence Intervals
# =============================================================================
# 13.1:Coerce to data.frame
coef_df <- as.data.frame(coefs)

# 13.2: Inspect the names of coef_df to find the CI columns
print(names(coef_df))
# e.g. you might see something like:
#   [1] "Parameter"   "Coefficient" "CI_2.5%"     "CI_97.5%"    ...

# 13.3: Programmatically detect the lower‐ and upper‐CI column names
ci_cols <- names(coef_df)[ grepl("^CI", names(coef_df)) ]
ci_cols <- sort(ci_cols)           # ensure lower bound comes first
lower_col <- ci_cols[1]
upper_col <- ci_cols[2]

# (Optional) print to confirm
message("Lower‐CI column: ", lower_col)
message("Upper‐CI column: ", upper_col)

# 13.4:. Build the plot using aes_string to refer to dynamic column names
library(ggplot2)
coef_plot <- ggplot(coef_df, aes_string(
  x    = "Coefficient",    # the estimate column
  y    = "Parameter",      # the factor names
  xmin = lower_col,
  xmax = upper_col
)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title = "Panel Fixed‐Effects Coefficient Estimates\nwith 95% Robust Confidence Intervals",
    x     = "Estimate",
    y     = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y        = element_text(size = 10),
    plot.title         = element_text(hjust = 0.5),
    panel.grid.major.y = element_blank()
  )

# 13.5: Display the plot
print(coef_plot)

# 13.6: Save to file
ggsave(
  filename = "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.2 回归筛选/Figure4.1.2_CoefficientPlot.png",
  plot     = coef_plot,
  width    = 8,
  height   = 6,
  dpi      = 300
)
# =============================================================================
# 14. Based on VIF, tests, and p‑values, iteratively drop
#    high‑VIF and non‑significant variables, then re‑estimate
#    Repeat steps 8–13 until diagnostics are satisfied.
# =============================================================================




