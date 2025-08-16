# =============================================================================
# 4.1.1 Ratio Analysis: Complete R Script without Unnecessary Log‑Transforms
# =============================================================================
# -----------------------------------------------------------------------------
# 0. Install and load required packages
# -----------------------------------------------------------------------------
# Uncomment to install if needed:
# install.packages(c("tidyverse", "readxl", "lubridate", "moments"))

library(tidyverse)  # Data manipulation & plotting
library(readxl)     # Excel import
library(lubridate)  # Date handling
library(moments)    # Skewness calculation

# -----------------------------------------------------------------------------
# 1. Load and inspect simulated PHS panel data
# -----------------------------------------------------------------------------
# Reads "Data.xlsx" sheet 1 with columns:
# Hospital, Year, Beds, StaffFTE, OperatingCost_kGBP, Discharges, InpatientDays

df_raw <- read_excel("C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/数据/Data.xlsx", sheet = 1)
glimpse(df_raw)

# -----------------------------------------------------------------------------
# 2. Data cleaning and preprocessing (per Section 3.3.2.4)
# -----------------------------------------------------------------------------

# 2.1 Offset zeros to avoid division/log issues
df_clean <- df_raw %>%
  mutate_at(
    vars(OperatingCost_kGBP, Discharges, InpatientDays, Beds, StaffFTE),
    ~ if_else(. <= 0, 1e-6, .)
  )

# 2.2 Impute missing values: hospital mean then global median
df_clean <- df_clean %>%
  group_by(Hospital) %>%
  mutate_at(
    vars(Beds, StaffFTE, OperatingCost_kGBP, Discharges, InpatientDays),
    ~ if_else(is.na(.), mean(., na.rm = TRUE), .)
  ) %>%
  ungroup() %>%
  mutate_at(
    vars(Beds, StaffFTE, OperatingCost_kGBP, Discharges, InpatientDays),
    ~ if_else(is.na(.), median(., na.rm = TRUE), .)
  )

# 2.3 Winsorise at Tukey IQR ±1.5
winsorize_iqr <- function(x) {
  q1 <- quantile(x, .25, na.rm = TRUE)
  q3 <- quantile(x, .75, na.rm = TRUE)
  iqr <- q3 - q1
  pmin(pmax(x, q1 - 1.5 * iqr), q3 + 1.5 * iqr)
}

df_clean <- df_clean %>%
  mutate_at(
    vars(Beds, StaffFTE, OperatingCost_kGBP, Discharges, InpatientDays),
    winsorize_iqr
  )

# -----------------------------------------------------------------------------
# 3. Compute core ratio metrics (raw values)
# -----------------------------------------------------------------------------
df_ratio <- df_clean %>%
  mutate(
    CostPerDischarge  = OperatingCost_kGBP / Discharges,
    BedOccupancyRate  = InpatientDays / (Beds * 365),
    AverageLOS        = InpatientDays / Discharges,
    StaffProductivity = Discharges / StaffFTE
  )

# -----------------------------------------------------------------------------
# 4. Skewness check and distribution plots
# -----------------------------------------------------------------------------

# 4.1 Compute skewness for each ratio
skew_stats <- df_ratio %>%
  summarise(
    skew_CPD = skewness(CostPerDischarge,    na.rm = TRUE),
    skew_BOR = skewness(BedOccupancyRate,    na.rm = TRUE),
    skew_AL  = skewness(AverageLOS,          na.rm = TRUE),
    skew_SP  = skewness(StaffProductivity,   na.rm = TRUE)
  )
print(skew_stats)
# Since all |skew| < 1, no log-transform is necessary

# 4.2 Plot histograms & density for each metric
df_ratio %>%
  select(CostPerDischarge, BedOccupancyRate, AverageLOS, StaffProductivity) %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Value)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "gray80", color = "white") +
  geom_density(color = "steelblue", size = 1) +
  facet_wrap(~ Metric, scales = "free") +
  labs(
    title = "Distributions of Core Ratio Metrics",
    x     = "Value",
    y     = "Density"
  ) +
  theme_minimal() -> p_dist
print(p_dist)
ggsave("C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.1比例分析/Figure4.1.1c_RatioDistributions.png", p_dist, width = 10, height = 6)

# -----------------------------------------------------------------------------
# 5. Descriptive statistics by year
# -----------------------------------------------------------------------------
summary_table <- df_ratio %>%
  group_by(Year) %>%
  summarise(
    mean_CPD = mean(CostPerDischarge,    na.rm = TRUE),
    sd_CPD   = sd(CostPerDischarge,      na.rm = TRUE),
    med_CPD  = median(CostPerDischarge,  na.rm = TRUE),
    iqr_CPD  = IQR(CostPerDischarge,     na.rm = TRUE),
    mean_BOR = mean(BedOccupancyRate,    na.rm = TRUE),
    sd_BOR   = sd(BedOccupancyRate,      na.rm = TRUE),
    med_BOR  = median(BedOccupancyRate,  na.rm = TRUE),
    iqr_BOR  = IQR(BedOccupancyRate,     na.rm = TRUE),
    mean_AL  = mean(AverageLOS,          na.rm = TRUE),
    sd_AL    = sd(AverageLOS,            na.rm = TRUE),
    med_AL   = median(AverageLOS,        na.rm = TRUE),
    iqr_AL   = IQR(AverageLOS,           na.rm = TRUE),
    mean_SP  = mean(StaffProductivity,   na.rm = TRUE),
    sd_SP    = sd(StaffProductivity,     na.rm = TRUE),
    med_SP   = median(StaffProductivity, na.rm = TRUE),
    iqr_SP   = IQR(StaffProductivity,    na.rm = TRUE)
  )

print(summary_table)


transposed_summary <- summary_table %>%
  pivot_longer(
    cols       = -Year,
    names_to   = "Statistic",
    values_to  = "Value"
  ) %>%
  pivot_wider(
    names_from  = Year,
    values_from = Value
  )

print(transposed_summary)

write_csv(transposed_summary, "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.1比例分析/Table4.1.1_DescriptiveRatios.csv")
# -----------------------------------------------------------------------------
# 6. Boxplots to identify outliers
# -----------------------------------------------------------------------------
df_long <- df_ratio %>%
  select(Hospital, Year, CostPerDischarge, BedOccupancyRate, AverageLOS, StaffProductivity) %>%
  pivot_longer(
    cols      = CostPerDischarge:StaffProductivity,
    names_to  = "Metric",
    values_to = "Value"
  )

short_labels <- c(
  AverageLOS         = "AvgLOS",
  BedOccupancyRate   = "BedOcc",
  CostPerDischarge   = "CostPD",
  StaffProductivity  = "StaffProd"
)


p_box <- ggplot(df_long, aes(x = Metric, y = Value)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 1) +
  facet_wrap(~ Year, ncol = 3) +
  scale_x_discrete(labels = short_labels) +
  labs(
    title = "Boxplots of Core Ratio Metrics by Year",
    x     = "Metric",
    y     = "Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text  = element_text(size = 10),
    plot.title  = element_text(hjust = 0.5)
  )

print(p_box)
ggsave("C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.1比例分析/Figure4.1.1a_Boxplots_shortlabels.png", p_box, width = 12, height = 8)

# ------------------------------------------------------------
# 7. Time‐series trends for Cost per Discharge (merged approach)
# ------------------------------------------------------------

# 1) Prepare the two pieces and bind into one DF
mean_trends <- df_ratio %>%
  group_by(Year) %>%
  summarise(CostPerDischarge = mean(CostPerDischarge, na.rm = TRUE)) %>%
  mutate(Hospital = "SampleMean")

# identify top‑3 by 2024 cost
top3 <- df_ratio %>%
  filter(Year == max(Year)) %>%
  arrange(desc(CostPerDischarge)) %>%
  slice_head(n = 3) %>%
  pull(Hospital)

hosp_trends <- df_ratio %>%
  filter(Hospital %in% top3) %>%
  select(Hospital, Year, CostPerDischarge)

df_trend <- bind_rows(hosp_trends, mean_trends)

# 2) Turn Hospital into a factor with desired order
df_trend$Hospital <- factor(df_trend$Hospital,
                            levels = c(top3, "SampleMean"))

# 3) Define manual aesthetics
custom_colors    <- c(H054 = "firebrick", H065 = "royalblue", H068 = "forestgreen", SampleMean = "purple")
custom_linetypes <- c(H054 = "solid",    H065 = "solid",      H068 = "solid",         SampleMean = "dashed")
custom_shapes    <- c(H054 = 16,         H065 = 17,           H068 = 15,              SampleMean = 18)

# 4) Plot
p_trend <- ggplot(df_trend, aes(x = Year, y = CostPerDischarge,
                                color = Hospital, linetype = Hospital, shape = Hospital)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values    = custom_colors) +
  scale_linetype_manual(values = custom_linetypes) +
  scale_shape_manual(values   = custom_shapes) +
  labs(
    title = "Cost per Discharge Trends: Selected Hospitals vs. Sample Mean",
    x     = "Year",
    y     = "Cost per Discharge (£’000)",
    color = "Hospital",
    linetype = "Hospital",
    shape = "Hospital"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# 5) Render and save
print(p_trend)
ggsave("C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.1比例分析/Figure4.1.1b_CPD_Trends.png", p_trend, width = 10, height = 6)

# -----------------------------------------------------------------------------
# 8. Export outlier list for documentation
# -----------------------------------------------------------------------------
outliers <- df_long %>%
  group_by(Year, Metric) %>%
  mutate(
    Q1         = quantile(Value, 0.25, na.rm = TRUE),
    Q3         = quantile(Value, 0.75, na.rm = TRUE),
    IQR        = Q3 - Q1,
    is_outlier = Value < Q1 - 1.5 * IQR | Value > Q3 + 1.5 * IQR
  ) %>%
  filter(is_outlier) %>%
  select(Hospital, Year, Metric, Value)

write_csv(outliers, "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.1比例分析/Table4.1.1_Outliers.csv")

# =============================================================================
# End of 4.1.1 Ratio Analysis Script
# =============================================================================