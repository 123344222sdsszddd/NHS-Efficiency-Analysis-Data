# ===================================================================
# Frontier-Focused DEA Workflow in R (Years separated, Bootstrap CIs)
# ===================================================================
# ------------------------------------------------------------
# STEP 0. Libraries & Global Plot Theme
# ------------------------------------------------------------
library(dplyr)         # data manipulation
library(tidyr)         # data reshaping
library(purrr)         # iteration / functional programming
library(ggplot2)       # visualization
library(ggpubr)        # statistical annotations on ggplot
library(readxl)        # Excel import
library(RColorBrewer)  # color palettes
library(rDEA)          # DEA point estimates
library(Benchmarking)  # DEA bootstrapping

# Set a clean, minimal theme
theme_set(
  theme_minimal(base_size = 14) +
    theme(
      plot.title    = element_text(face="bold", size=16, hjust=0.5),
      plot.subtitle = element_text(size=13, hjust=0.5, color="gray40"),
      legend.position = "bottom",
      legend.title = element_text(face="bold"),
      strip.text = element_text(face="bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face="bold")
    )
)

# Single 3-class palette (Small/Medium/Large)
pal3 <- brewer.pal(3, "Set2")


# ------------------------------------------------------------
# STEP 1. Data Grouping & Preprocessing (per Year)
# ------------------------------------------------------------

# 1.1 Import raw data
data_raw <- read_excel(
  "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/数据/NHS_Scotland_Simulated_Dataset.xlsx",
  sheet = 1
) %>% mutate(Year = as.integer(Year))  # ensure integer Year

# 1.2 Winsorization helper (IQR-based)
winsorize <- function(x) {
  q <- quantile(x, c(.25, .75), na.rm = TRUE); iqr <- q[2] - q[1]
  x[x < q[1] - 1.5*iqr] <- q[1] - 1.5*iqr
  x[x > q[2] + 1.5*iqr] <- q[2] + 1.5*iqr
  x
}

# 1.3 Per-year processing
data_prep <- data_raw %>%
  group_by(Year) %>%
  # (a) Size tertiles for three proxies
  mutate(
    BedsGroup = cut(Beds,
                    breaks = quantile(Beds, probs = c(0, .33, .66, 1), na.rm = TRUE),
                    labels = c("Small","Medium","Large"), include.lowest = TRUE),
    FTEGroup = cut(StaffFTE,
                   breaks = quantile(StaffFTE, probs = c(0, .33, .66, 1), na.rm = TRUE),
                   labels = c("Small","Medium","Large"), include.lowest = TRUE),
    DiscGroup = cut(Discharges,
                    breaks = quantile(Discharges, probs = c(0, .33, .66, 1), na.rm = TRUE),
                    labels = c("Small","Medium","Large"), include.lowest = TRUE)
  ) %>%
  ungroup() %>%
  # (b) Majority vote for final SizeGroup
  mutate(
    SizeGroup = apply(select(., BedsGroup, FTEGroup, DiscGroup), 1,
                      function(r) names(which.max(table(r))))
  ) %>%
  mutate(SizeGroup = factor(SizeGroup, levels = c("Small","Medium","Large"))) %>%
  # (c) Mean imputation within Year × SizeGroup
  group_by(Year, SizeGroup) %>%
  mutate(across(
    c(Beds, StaffFTE, OperatingCost_kGBP, Discharges, OutpatientVisits, CaseMixIndex),
    ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
  )) %>%
  ungroup() %>%
  # (d) Winsorize key variables
  mutate(across(
    c(Beds, StaffFTE, OperatingCost_kGBP, Discharges, OutpatientVisits, CaseMixIndex),
    winsorize
  ))


# ------------------------------------------------------------
# STEP 2. Per-Year DEA (TE, PTE, SE, EE, AE) + Bootstrap CIs
# ------------------------------------------------------------

# 2.1 Inputs / Outputs
input_vars  <- c("Beds","StaffFTE","OperatingCost_kGBP","InpatientDays")
output_vars <- c("Discharges","OutpatientVisits")

# 2.2 Iterate over years
results_panel <- map_dfr(
  sort(unique(data_prep$Year)),
  function(yr) {
    df <- dplyr::filter(data_prep, Year == yr)
    
    # (a) Normalization for point-estimate DEA
    means_in  <- colMeans(df[input_vars])
    means_out <- colMeans(df[output_vars])
    Xn <- as.matrix(df[input_vars]  / means_in[col(df[input_vars])])
    Yn <- as.matrix(df[output_vars] / means_out[col(df[output_vars])])
    
    # Raw matrices for cost-DEA
    Xr <- as.matrix(df[input_vars])
    Yr <- as.matrix(df[output_vars])
    
    # (b) Price weights for cost-minimization
    pw <- means_in / sum(means_in)
    W  <- matrix(rep(pw, each = nrow(df)), nrow = nrow(df), byrow = FALSE)
    
    # (c) Point estimates: CRS input-oriented, VRS input-oriented, VRS cost-min
    dea_crs  <- rDEA::dea(XREF=Xn, YREF=Yn, X=Xn, Y=Yn, model="input",   RTS="constant")
    dea_bcc  <- rDEA::dea(XREF=Xn, YREF=Yn, X=Xn, Y=Yn, model="input",   RTS="variable")
    dea_cost <- rDEA::dea(XREF=Xr, YREF=Yr, X=Xr, Y=Yr, W=W, model="costmin", RTS="variable")
    
    te  <- dea_crs$thetaOpt     # Technical Efficiency (CRS)
    pte <- dea_bcc$thetaOpt     # Pure Technical Efficiency (VRS)
    se  <- te / pte             # Scale Efficiency
    ee  <- dea_cost$gammaOpt    # Economic Efficiency
    ae  <- ee / pte             # Allocative Efficiency
    
    # (d) Bootstrap (CRS TE & VRS PTE)
    boot_crs <- Benchmarking::dea.boot(
      X=Xn, Y=Yn, NREP=200, EFF=te,  RTS="crs", ORIENTATION="in",
      XREF=Xn, YREF=Yn, alpha=0.05
    )
    boot_bcc <- Benchmarking::dea.boot(
      X=Xn, Y=Yn, NREP=200, EFF=pte, RTS="vrs", ORIENTATION="in",
      XREF=Xn, YREF=Yn, alpha=0.05
    )
    
    # (e) Bias-corrected estimates & 95% CIs
    te_bc      <- boot_crs$eff.bc
    te_ci_low  <- boot_crs$conf.int[,1]
    te_ci_high <- boot_crs$conf.int[,2]
    
    pte_bc      <- boot_bcc$eff.bc
    pte_ci_low  <- boot_bcc$conf.int[,1]
    pte_ci_high <- boot_bcc$conf.int[,2]
    
    # (f) Output row-block
    tibble(
      Year       = yr,
      Hospital   = df$Hospital,
      SizeGroup  = df$SizeGroup,
      TE   = te,  PTE  = pte,  SE   = se,  EE   = ee,  AE   = ae,
      TE_bc   = te_bc,  PTE_bc  = pte_bc,
      TE_ci_low   = te_ci_low,  TE_ci_high  = te_ci_high,
      PTE_ci_low  = pte_ci_low, PTE_ci_high = pte_ci_high
    )
  }
)

# ------------------------------------------------------------
# STEP 2.6. Raincloud Plot for Cost per Discharge (CPD)
# ------------------------------------------------------------

# (Libraries already loaded above; avoid duplicates)
library(gghalves)      # half-violin/half-boxplot geoms

# 2.6.1 Compute CPD and retain grouping vars
cpd_df <- data_prep %>%
  mutate(CPD = OperatingCost_kGBP / Discharges) %>%
  select(Year, Hospital, SizeGroup, CPD)

# 2.6.2 Plot CPD by SizeGroup × Year
ggplot(cpd_df, aes(x = SizeGroup, y = CPD, fill = SizeGroup)) +
  # (i) left half-violin (density)
  geom_half_violin(side="l", trim=FALSE, alpha=0.6,
                   position=position_nudge(x=-0.2), colour=NA) +
  # (ii) right half-boxplot (summary stats)
  geom_half_boxplot(side="r", width=0.2, outlier.shape=NA,
                    position=position_nudge(x=0.2), alpha=0.7, colour="black") +
  # (iii) jittered raw points
  geom_jitter(width=0.1, size=1, alpha=0.3, colour="gray20") +
  # (iv) Kruskal–Wallis p value per facet
  stat_compare_means(method="kruskal.test", label="p.format",
                     label.x.npc="left", label.y.npc="top", size=3) +
  facet_wrap(~ Year, ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Raincloud Plot of Cost per Discharge by SizeGroup and Year",
    x     = "Hospital Size Group",
    y     = "Cost per Discharge (kGBP)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold")
  )


# ------------------------------------------------------------
# STEP 3. Grouped Summary & Rainclouds for Efficiency Metrics
# ------------------------------------------------------------

# 3.1 Summary table with CIs (wide → long → summary)
summary_stats <- results_panel %>%
  pivot_longer(
    cols = c(TE, PTE, SE, EE, AE, TE_bc, PTE_bc,
             TE_ci_low, TE_ci_high, PTE_ci_low, PTE_ci_high),
    names_to  = "Metric", values_to = "Value"
  ) %>%
  group_by(Year, SizeGroup, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD   = sd(Value,   na.rm = TRUE),
    Min  = min(Value,  na.rm = TRUE),
    Q1   = quantile(Value, .25, na.rm = TRUE),
    Med  = median(Value,    na.rm = TRUE),
    Q3   = quantile(Value, .75, na.rm = TRUE),
    Max  = max(Value,  na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# 3.2 Raincloud-style plots for TE, SE, EE, AE
metrics <- c("TE", "SE", "EE", "AE")

for (m in metrics) {
  p <- ggplot(results_panel, aes(x = SizeGroup, y = .data[[m]], fill = SizeGroup)) +
    geom_half_violin(side="l", trim=FALSE, alpha=0.6,
                     position=position_nudge(x=-0.15), colour=NA) +
    geom_boxplot(width=0.2, outlier.shape=NA,
                 position=position_nudge(x=0.15), alpha=0.7, colour="black") +
    geom_jitter(width=0.05, size=1, alpha=0.3, aes(colour = SizeGroup)) +
    stat_compare_means(method="kruskal.test", label="p.format",
                       label.x.npc="left", label.y.npc="top", size=3) +
    facet_wrap(~ Year, ncol = 3) +
    scale_fill_brewer(palette = "Set2") +
    scale_colour_brewer(palette = "Set2") +
    labs(
      title = paste0(m, " Raincloud Plot by SizeGroup and Year"),
      x     = "Hospital Size Group",
      y     = m
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position  = "none",
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold")
    )
  print(p)
}


# ------------------------------------------------------------
# STEP 4. Frontier Scatter (Cost vs Volume; mark TE=1 / EE=1)
# ------------------------------------------------------------

scatter_df <- data_prep %>%
  mutate(CostPerDischarge = OperatingCost_kGBP / Discharges) %>%
  left_join(results_panel, by = c("Year","Hospital","SizeGroup"))

ggplot(scatter_df, aes(x = CostPerDischarge, y = Discharges, color = SizeGroup)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Year) +
  # Mark TE frontiers
  geom_point(
    data = dplyr::filter(scatter_df, TE == 1),
    aes(x = CostPerDischarge, y = Discharges),
    shape = 17, color = "blue", size = 2, inherit.aes = FALSE
  ) +
  # Mark EE frontiers
  geom_point(
    data = dplyr::filter(scatter_df, EE == 1),
    aes(x = CostPerDischarge, y = Discharges),
    shape = 15, color = "red", size = 2, inherit.aes = FALSE
  ) +
  labs(
    title = "Frontier Scatter by Year\n▲ TE=1 (blue), ■ EE=1 (red)",
    x = "Cost per Discharge (kGBP)",
    y = "Annual Discharges"
  ) +
  theme_minimal() +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)


# ------------------------------------------------------------
# STEP 5. Six “All-Angles” Visualizations
# ------------------------------------------------------------

# 5.1 Parallel Coordinates of TE/PTE/SE/EE/AE (scaled within Metric)
df_pc <- results_panel %>%
  pivot_longer(cols = TE:AE, names_to = "Metric", values_to = "Value") %>%
  group_by(Year, SizeGroup, Metric) %>%
  mutate(Scaled = (Value - min(Value)) / (max(Value) - min(Value))) %>%
  ungroup()

# Overlay thin hospital lines + thick group means
df_pc_mean <- df_pc %>%
  group_by(SizeGroup, Metric) %>%
  summarise(MeanScaled = mean(Scaled, na.rm = TRUE), .groups = "drop")

ggplot(df_pc, aes(x = Metric, y = Scaled, group = Hospital, color = SizeGroup)) +
  geom_line(alpha = 0.15) +  # individual profiles
  geom_line(data = df_pc_mean, aes(x = Metric, y = MeanScaled, group = SizeGroup, color = SizeGroup),
            size = 1.2, alpha = 0.9, inherit.aes = FALSE) +
  labs(
    title = "Parallel-Coordinates of Efficiency Profiles",
    subtitle = "TE / PTE / SE / EE / AE (scaled 0–1); thin = hospitals, thick = group means"
  ) +
  theme_minimal()

# 5.2 Trend Lines with 95% CI Ribbon (Avg TE, Avg EE)
trend_stats <- results_panel %>%
  group_by(Year, SizeGroup) %>%
  summarise(
    MeanTE   = mean(TE, na.rm = TRUE),
    TE_low   = mean(TE_ci_low,  na.rm = TRUE),
    TE_high  = mean(TE_ci_high, na.rm = TRUE),
    MeanEE   = mean(EE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(MeanTE, MeanEE), names_to = "Metric", values_to = "MeanValue")

ggplot(trend_stats, aes(x = Year, y = MeanValue, color = SizeGroup, group = SizeGroup)) +
  # (i) CI ribbon for TE
  geom_ribbon(
    data = dplyr::filter(trend_stats, Metric == "MeanTE"),
    aes(x = Year, ymin = TE_low, ymax = TE_high, fill = SizeGroup),
    inherit.aes = FALSE, alpha = 0.2
  ) +
  # (ii) Mean lines & points
  geom_line(size = 1) + geom_point(size = 2) +
  # (iii) Thresholds per facet
  geom_hline(
    data = tibble(Metric = c("MeanTE", "MeanEE"), thr = c(0.90, 0.85)),
    aes(yintercept = thr), linetype = "dashed", color = "grey40", inherit.aes = FALSE
  ) +
  facet_wrap(~Metric, scales = "fixed",
             labeller = as_labeller(c(MeanTE = "Avg TE", MeanEE = "Avg EE"))) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer (palette = "Set2") +
  labs(
    title    = "Yearly Trends of Avg TE & Avg EE by SizeGroup",
    subtitle = "Dashed lines indicate TE=0.90 and EE=0.85 thresholds",
    x        = "Year",
    y        = "Average Score",
    color    = "SizeGroup",
    fill     = "SizeGroup"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

# 5.3 Cluster Scatter: TE vs EE (+ highlight TE_ci_low > 0.95)
cluster_df <- results_panel %>%
  group_by(Hospital) %>%
  summarise(AvgTE = mean(TE, na.rm = TRUE),
            AvgEE = mean(EE, na.rm = TRUE), .groups = "drop")
set.seed(42)
kmeans_res <- kmeans(cluster_df[, c("AvgTE", "AvgEE")], centers = 3)
cluster_df$Cluster <- factor(kmeans_res$cluster)

plot_cl <- results_panel %>% left_join(cluster_df %>% select(Hospital, Cluster), by = "Hospital")

ggplot(plot_cl, aes(x = TE, y = EE, color = Cluster)) +
  geom_point(shape = 16, size = 2, alpha = 0.6) +
  geom_vline(xintercept = 0.90, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 0.85, linetype = "dotted", color = "grey50") +
  geom_point(
    data = subset(plot_cl, TE_ci_low > 0.95),
    aes(x = TE, y = EE), inherit.aes = FALSE,
    shape = 8, color = "red", size = 3, show.legend = FALSE
  ) +
  facet_wrap(~ Year) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "Scatter of Technical vs Economic Efficiency by Year",
    subtitle = "Clusters colored; ★ = TE_ci_low > 0.95; dotted lines: TE=0.90, EE=0.85",
    x        = "Technical Efficiency (TE)",
    y        = "Economic Efficiency (EE)",
    color    = "Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.position="right")

# 5.4 Input–TE Scatter: Beds vs OperatingCost colored by TE, shaped by SizeGroup
input_scatter <- data_prep %>% left_join(results_panel, by = c("Year","Hospital","SizeGroup"))

ggplot(input_scatter, aes(x = Beds, y = OperatingCost_kGBP, color = TE, shape = SizeGroup)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(option = "plasma", labels = scales::percent) +
  facet_wrap(~Year) +
  labs(
    title = "Beds vs Operating Cost by Year (color = TE, shape = SizeGroup)",
    x = "Beds (Staffed)",
    y = "Operating Cost (kGBP)",
    color = "TE", shape = "SizeGroup"
  ) +
  theme_minimal()

# 5.5 Frontier Scatter: CPD vs Discharges (Year facet, SizeGroup as color/shape)
ggplot(input_scatter,
       aes(x = OperatingCost_kGBP / Discharges, y = Discharges,
           color = SizeGroup, shape = SizeGroup)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Year) +
  geom_point(
    data = dplyr::filter(input_scatter, TE == 1),
    aes(x = OperatingCost_kGBP/Discharges, y = Discharges),
    color = "blue", shape = 17, size = 2, inherit.aes = FALSE
  ) +
  geom_point(
    data = dplyr::filter(input_scatter, EE == 1),
    aes(x = OperatingCost_kGBP/Discharges, y = Discharges),
    color = "red", shape = 15, size = 2, inherit.aes = FALSE
  ) +
  labs(
    title = "Frontier Scatter: Cost per Discharge vs Discharges\n▲ TE = 1 (blue), ■ EE = 1 (red)",
    x = "Cost per Discharge",
    y = "Discharges"
  ) +
  theme_minimal()

# 5.6 Radar/Spider Chart for group means in a chosen Year
plot_year <- 2024
radar_df <- results_panel %>%
  dplyr::filter(Year == plot_year) %>%
  pivot_longer(cols = TE:AE, names_to = "Metric", values_to = "Value") %>%
  group_by(SizeGroup, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), .groups = "drop")

ggplot(radar_df, aes(x = Metric, y = Mean, group = SizeGroup, color = SizeGroup)) +
  geom_polygon(fill = NA, size = 1) +
  geom_point(size = 2) +
  coord_polar() +
  labs(
    title    = paste("Radar Chart of Mean Efficiencies in", plot_year),
    subtitle = "Metrics: TE | PTE | SE | EE | AE"
  ) +
  theme_minimal()


# ------------------------------------------------------------
# STEP 6. Slacks / Targets (VRS, input-oriented)
# ------------------------------------------------------------

## 6.0 Packages for this section
library(stringr)
library(scales)
library(ggh4x)   # per-facet axis scales
library(rlang)
library(forcats)
library(tidytext)

## 6.1 Inputs / Outputs
input_vars   <- c("Beds","StaffFTE","OperatingCost_kGBP","InpatientDays")
output_vars  <- c("Discharges","OutpatientVisits")

## 6.2 Helper: run DEA & build slack table for ONE year
run_dea_slack_year <- function(df_year, year,
                               input_vars, output_vars,
                               rts = "vrs", orientation = "in",
                               lambda_cut = 1e-4){
  X <- as.matrix(df_year[, input_vars,  drop = FALSE])
  Y <- as.matrix(df_year[, output_vars, drop = FALSE])
  
  dea_res <- Benchmarking::dea(X, Y, RTS = rts, ORIENTATION = orientation)
  theta   <- as.numeric(dea_res$eff)
  lam     <- as.matrix(dea_res$lambda)
  
  X_hat <- lam %*% X; Y_hat <- lam %*% Y
  
  sx <- sweep(X, 1, theta, "*") - X_hat;  sx[sx < 0] <- 0     # input slacks
  sy <- Y_hat - Y;                        sy[sy < 0] <- 0     # output slacks
  
  X_target <- sweep(X, 1, theta, "*") - sx
  Y_target <- Y + sy
  
  cut_pct_in  <- sx / X                    # % reductions for inputs
  add_pct_out <- sy / Y                    # % increases for outputs
  
  tibble(
    Year     = year,
    Hospital = df_year$Hospital,
    theta    = theta,
    PeerSet  = apply(lam, 1, function(v){
      idx <- which(v > lambda_cut); if(!length(idx)) return(NA_character_)
      paste(df_year$Hospital[idx], collapse = "; ")
    }),
    PeerNum  = apply(lam, 1, function(v) sum(v > lambda_cut))
  ) %>%
    bind_cols(as.data.frame(sx)          %>% set_names(paste0("CutAbs_",    input_vars))) %>%
    bind_cols(as.data.frame(cut_pct_in)  %>% set_names(paste0("CutPct_",    input_vars))) %>%
    bind_cols(as.data.frame(X_target)    %>% set_names(paste0("Target_",    input_vars))) %>%
    bind_cols(as.data.frame(sy)          %>% set_names(paste0("SlackOut_",  output_vars))) %>%
    bind_cols(as.data.frame(add_pct_out) %>% set_names(paste0("SlackPct_",  output_vars))) %>%
    bind_cols(as.data.frame(Y_target)    %>% set_names(paste0("TargetOut_", output_vars)))
}

## 6.3 Run for ALL years
years_vec <- sort(unique(data_prep$Year))
slack_all_years <- purrr::map_dfr(
  years_vec,
  ~ run_dea_slack_year(dplyr::filter(data_prep, Year == .x),
                       .x, input_vars, output_vars,
                       rts = "vrs", orientation = "in")
)

## 6.4 Helper for per-facet axis scaling
auto_y_scales <- function(df, facet_var = "Var", value_var = "Pct",
                          n_breaks = 5, accuracy = 0.1){
  stats <- df %>%
    group_by(.data[[facet_var]]) %>%
    summarise(max_v = max(.data[[value_var]], na.rm = TRUE), .groups = "drop")
  
  purrr::imap(stats$max_v, function(mx, i){
    vname <- stats[[facet_var]][i]
    brks <- pretty(c(0, mx), n = n_breaks)
    brks <- unique(c(0, brks[brks >= 0]))
    cond  <- rlang::parse_expr(paste0(facet_var, " == '", vname, "'"))
    scale <- scale_y_continuous(breaks = brks,
                                labels = percent_format(accuracy = accuracy))
    new_formula(cond, scale)
  })
}

## 6.5 Latest-year summary (table for detailed slacks)
latest_year <- max(slack_all_years$Year, na.rm = TRUE)
slack_ly    <- dplyr::filter(slack_all_years, Year == latest_year)

cols_show <- c(
  "Hospital","theta","PeerNum","PeerSet",
  grep("^CutAbs_|^CutPct_|^SlackOut_|^SlackPct_", names(slack_ly), value = TRUE)
)

table_latest <- slack_ly %>%
  select(any_of(cols_show)) %>%
  arrange(desc(CutPct_OperatingCost_kGBP))
print(table_latest)

## 6.5.1 Group-level mean slacks (latest year)
size_lookup <- results_panel %>% distinct(Hospital, SizeGroup)
slack_group_means <- slack_ly %>%
  left_join(size_lookup, by = "Hospital") %>%
  group_by(SizeGroup) %>%
  summarise(
    Mean_CutPct_Beds               = mean(CutPct_Beds,               na.rm = TRUE),
    Mean_CutPct_StaffFTE           = mean(CutPct_StaffFTE,           na.rm = TRUE),
    Mean_CutPct_OperatingCost_kGBP = mean(CutPct_OperatingCost_kGBP, na.rm = TRUE),
    Mean_CutPct_InpatientDays      = mean(CutPct_InpatientDays,      na.rm = TRUE),
    Mean_SlackPct_Discharges       = mean(SlackPct_Discharges,       na.rm = TRUE),
    Mean_SlackPct_OutpatientVisits = mean(SlackPct_OutpatientVisits, na.rm = TRUE),
    Mean_PeerNum                   = mean(PeerNum,                   na.rm = TRUE)
  )
print(slack_group_means)

## 6.5.2 Bar chart: avg operating-cost slack % by SizeGroup (latest year)
p_cost_slack <- ggplot(slack_group_means,
                       aes(x = SizeGroup, y = Mean_CutPct_OperatingCost_kGBP, fill = SizeGroup)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Average Operating-Cost Slack % by Size Group (2024)",
    x     = "Size Group",
    y     = "Mean Operating-Cost Slack %"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x= element_text(face = "bold"))
print(p_cost_slack)

## 6.5.3 Export table & bar plot
ggsave(
  filename = file.path("C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.3DEA", "4.3.2_MeanOperatingCostSlack.png"),
  plot     = p_cost_slack, width = 8, height = 5, dpi = 300
)
write.csv(
  slack_group_means,
  file = file.path("C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.3DEA", "4.3.2_SlackGroupMeans.csv"),
  row.names = FALSE
)

## 6.6 Lollipop: Top 5 per SizeGroup by Operating-Cost Slack (%) — Year 2024
topN_per_group <- 5
target_year    <- 2024

size_lookup_2024 <- data_prep %>%
  dplyr::filter(Year == target_year) %>%
  dplyr::distinct(Hospital, SizeGroup)

slack_2024 <- slack_all_years %>% dplyr::filter(Year == target_year)

lop_df <- slack_2024 %>%
  dplyr::left_join(size_lookup_2024, by = "Hospital") %>%
  dplyr::mutate(CutPct_Cost = CutPct_OperatingCost_kGBP) %>%
  dplyr::filter(!is.na(CutPct_Cost), !is.na(SizeGroup)) %>%
  dplyr::group_by(SizeGroup) %>%
  dplyr::arrange(dplyr::desc(CutPct_Cost), Hospital, .by_group = TRUE) %>%
  dplyr::slice_head(n = topN_per_group) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Hospital_lbl = tidytext::reorder_within(Hospital, CutPct_Cost, SizeGroup),
    CutPct_lab   = scales::percent(CutPct_Cost, accuracy = 0.1)
  )

y_max <- max(lop_df$CutPct_Cost, na.rm = TRUE) * 1.12

lop_plot <- ggplot(lop_df, aes(x = Hospital_lbl, y = CutPct_Cost)) +
  geom_segment(aes(xend = Hospital_lbl, y = 0, yend = CutPct_Cost),
               linewidth = 0.9, color = "grey75") +
  geom_point(size = 3.2, color = "red2") +
  geom_text(aes(label = CutPct_lab), hjust = -0.3, size = 3.2, color = "grey20") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(0, y_max), expand = expansion(mult = c(0, 0.02))) +
  tidytext::scale_x_reordered() +
  facet_wrap(~ SizeGroup, ncol = 1, scales = "free_y") +
  labs(
    title    = paste0("Top ", topN_per_group,
                      " per Size Group by Potential Operating-Cost Cut (", target_year, ")"),
    subtitle = "Input-oriented VRS DEA slack (Operating Cost)",
    x        = NULL, y = "Cut % of Operating Cost"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.major.y = element_line(color = "grey93"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
print(lop_plot)

## 6.7 Heatmap: Input slack % (latest year), faceted by SizeGroup
size_lookup <- results_panel %>% distinct(Hospital, SizeGroup)

hm_in_df <- slack_ly %>%
  left_join(size_lookup, by = "Hospital") %>%
  select(Hospital, SizeGroup, starts_with("CutPct_")) %>%
  pivot_longer(starts_with("CutPct_"), names_to = "Input", values_to = "CutPct") %>%
  mutate(Input = str_remove(Input, "^CutPct_"))

ggplot(hm_in_df, aes(Input, reorder(Hospital, CutPct), fill = CutPct)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("#f7f7f7","#fdbb84","#e34a33"),
                       labels = percent, na.value = "grey90") +
  facet_wrap(~SizeGroup, ncol = 1, scales = "free_y") +
  labs(title = paste0("Input Slack % Heatmap (", latest_year, ")"),
       x = "Input", y = "Hospital", fill = "Slack %") +
  theme(axis.text.y = element_text(size = 6))

## 6.8 Trend of mean input slack % (2019–2024)
trend_in_df <- slack_all_years %>%
  left_join(size_lookup, by = "Hospital") %>%
  select(Year, SizeGroup, starts_with("CutPct_")) %>%
  pivot_longer(starts_with("CutPct_"), names_to = "Input", values_to = "SlackPct") %>%
  mutate(Input = str_remove(Input, "^CutPct_")) %>%
  group_by(Year, SizeGroup, Input) %>%
  summarise(MeanSlackPct = mean(SlackPct, na.rm = TRUE), .groups = "drop")

p_in <- ggplot(trend_in_df, aes(Year, MeanSlackPct, color = SizeGroup)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~Input, scales = "free_y") +
  labs(title   = "Mean Input Slack% Trends (2019–2024)",
       x       = "Year",
       y       = "Mean Slack %",
       caption = "Slack % = (Observed - Target) / Observed") +
  theme_minimal(base_size = 11)

p_in +
  facetted_pos_scales(
    y = list(
      Input == "OperatingCost_kGBP" ~ scale_y_continuous(
        breaks = seq(0, 0.07, 0.01), labels = percent_format(accuracy = 1)
      ),
      Input %in% c("Beds","InpatientDays","StaffFTE") ~ scale_y_continuous(
        breaks = seq(0, 0.025, 0.005), labels = percent_format(accuracy = 0.1)
      )
    )
  )

## 6.9 Trend of mean OUTPUT slack % (2019–2024)
cols_out_pct <- grep("^SlackPct_", names(slack_all_years), value = TRUE)

trend_out_df <- slack_all_years %>%
  left_join(size_lookup, by = "Hospital") %>%
  select(Year, SizeGroup, all_of(cols_out_pct)) %>%
  pivot_longer(all_of(cols_out_pct), names_to = "Output", values_to = "SlackPct") %>%
  mutate(Output = str_remove(Output, "^SlackPct_")) %>%
  group_by(Year, SizeGroup, Output) %>%
  summarise(MeanSlackPct = mean(SlackPct, na.rm = TRUE), .groups = "drop")

p_out <- ggplot(trend_out_df, aes(Year, MeanSlackPct, color = SizeGroup)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~Output, scales = "free_y") +
  labs(title   = "Mean Output Slack% Trends (2019–2024)",
       x       = "Year",
       y       = "Mean Slack %",
       caption = "Slack % = (Target - Observed) / Observed") +
  theme_minimal(base_size = 11)

p_out_final <- p_out +
  facetted_pos_scales(
    y = list(
      Output == "Discharges" ~ scale_y_continuous(
        breaks = seq(0, 0.03, 0.005), labels = percent_format(accuracy = 0.1)
      ),
      Output == "OutpatientVisits" ~ scale_y_continuous(
        breaks = seq(0, 0.05, 0.01), labels = percent_format(accuracy = 1)
      )
    )
  )
p_out_final

## 6.10 Heatmap: Output slack % (latest year)
hm_out_df <- slack_ly %>%
  left_join(size_lookup, by = "Hospital") %>%
  select(Hospital, SizeGroup, all_of(cols_out_pct)) %>%
  pivot_longer(all_of(cols_out_pct), names_to = "Output", values_to = "SlackPct") %>%
  mutate(Output = str_remove(Output, "^SlackPct_"))

ggplot(hm_out_df, aes(Output, reorder(Hospital, SlackPct), fill = SlackPct)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("#f7f7f7","#9ecae1","#08519c"),
                       labels = percent, na.value = "grey90") +
  facet_wrap(~SizeGroup, ncol = 1, scales = "free_y") +
  labs(title = paste0("Output Slack % Heatmap (", latest_year, ")"),
       x = "Output", y = "Hospital", fill = "Slack %") +
  theme(axis.text.y = element_text(size = 6))


# ------------------------------------------------------------
# STEP 6.11. Export results to disk (tables & figures)
# ------------------------------------------------------------

# 6.11.0 Output directory
out_dir <- "C:/Users/dell/OneDrive - University of Edinburgh/桌面/毕业论文/论文草稿/图表/4.1.3DEA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 6.11.1 Tables
write.csv(results_panel,
          file = file.path(out_dir, "4.1.3_Efficiency_AllYears_ByHospital.csv"),
          row.names = FALSE)
write.csv(summary_stats,
          file = file.path(out_dir, "4.1.3_Summary_Stats.csv"),
          row.names = FALSE)
write.csv(slack_all_years,
          file = file.path(out_dir, "4.1.3_Slack_AllYears.csv"),
          row.names = FALSE)
write.csv(table_latest,
          file = file.path(out_dir, paste0("4.1.3_Slack_Table_", latest_year, ".csv")),
          row.names = FALSE)

# 6.11.2 CPD raincloud export
cpd_plot <- ggplot(cpd_df, aes(x = SizeGroup, y = CPD, fill = SizeGroup)) +
  geom_half_violin(side="l", trim=FALSE, alpha=0.6, position=position_nudge(x=-0.2), colour=NA) +
  geom_half_boxplot(side="r", width=0.2, outlier.shape=NA, position=position_nudge(x=0.2),
                    alpha=0.7, colour="black") +
  geom_jitter(width=0.1, size=1, alpha=0.3, colour="gray20") +
  stat_compare_means(method="kruskal.test", label="p.format",
                     label.x.npc="left", label.y.npc="top", size=3) +
  facet_wrap(~Year, ncol=3) +
  scale_fill_brewer(palette="Set2") +
  labs(title="Raincloud Plot of Cost per Discharge by SizeGroup and Year",
       x="Hospital Size Group", y="Cost per Discharge (kGBP)") +
  theme_minimal(base_size=13) +
  theme(legend.position="none", strip.text=element_text(face="bold"),
        panel.grid.minor=element_blank(), axis.title=element_text(face="bold"))
ggsave(cpd_plot, filename = file.path(out_dir, "4.1.3_Raincloud_CPD.png"),
       width = 10, height = 6, dpi = 300)

# 6.11.3 Metric rainclouds export (TE/SE/EE/AE)
metrics <- c("TE","SE","EE","AE")
for (m in metrics) {
  p <- ggplot(results_panel, aes(x = SizeGroup, y = .data[[m]], fill = SizeGroup)) +
    geom_half_violin(side="l", trim=FALSE, alpha=0.6,
                     position=position_nudge(x=-0.15), colour=NA) +
    geom_boxplot(width=0.2, outlier.shape=NA,
                 position=position_nudge(x=0.15), alpha=0.7, colour="black") +
    geom_jitter(width=0.05, size=1, alpha=0.3, aes(colour=SizeGroup)) +
    stat_compare_means(method="kruskal.test", label="p.format",
                       label.x.npc="left", label.y.npc="top", size=3) +
    facet_wrap(~Year, ncol=3) +
    scale_fill_brewer(palette="Set2") +
    scale_colour_brewer(palette="Set2") +
    labs(title=paste0(m, " Raincloud Plot by SizeGroup and Year"),
         x="Hospital Size Group", y=m) +
    theme_minimal(base_size=13) +
    theme(legend.position="none", strip.text=element_text(face="bold"),
          panel.grid.minor=element_blank(), axis.title=element_text(face="bold"))
  ggsave(p, filename = file.path(out_dir, paste0("4.1.3_Raincloud_", m, ".png")),
         width = 10, height = 6, dpi = 300)
}

# 6.11.4 Frontier scatter (Cost vs Volume)
fs <- ggplot(scatter_df, aes(x = CostPerDischarge, y = Discharges, color = SizeGroup)) +
  geom_point(alpha=0.5) + facet_wrap(~Year) +
  geom_point(data=dplyr::filter(scatter_df, TE==1),
             aes(x=CostPerDischarge,y=Discharges),
             shape=17, color="blue", size=2) +
  geom_point(data=dplyr::filter(scatter_df, EE==1),
             aes(x=CostPerDischarge,y=Discharges),
             shape=15, color="red", size=2) +
  labs(title="Frontier Scatter by Year (▲ TE=1, ■ EE=1)",
       x="Cost per Discharge (kGBP)", y="Annual Discharges") +
  theme_minimal() +
  scale_x_continuous(labels=scales::comma) +
  scale_y_continuous(labels=scales::comma)
ggsave(fs, filename = file.path(out_dir, "4.1.3_Frontier_Scatter_Cost_vs_Vol.png"),
       width=10, height=6, dpi=300)

# 6.11.5 Parallel-coordinates export
pc <- ggplot(df_pc, aes(x = Metric, y = Scaled, group = Hospital, color = SizeGroup)) +
  geom_line(alpha=0.15) +
  geom_line(data = df_pc_mean,
            aes(x = Metric, y = MeanScaled, group = SizeGroup, color = SizeGroup),
            size=1.2, alpha=0.9, inherit.aes=FALSE) +
  labs(title="Parallel-Coordinates of Efficiency Profiles",
       subtitle="TE / PTE / SE / EE / AE (scaled 0–1); thin=hospitals, thick=group means") +
  theme_minimal()
ggsave(pc, filename = file.path(out_dir, "4.1.3_Parallel_Coordinates.png"),
       width=10, height=6, dpi=300)

# 6.11.6 Trend (Avg TE & Avg EE) export
tr <- ggplot(trend_stats, aes(x=Year,y=MeanValue, color=SizeGroup, group=SizeGroup)) +
  geom_ribbon(data=dplyr::filter(trend_stats, Metric=="MeanTE"),
              aes(x=Year, ymin=TE_low, ymax=TE_high, fill=SizeGroup),
              inherit.aes=FALSE, alpha=0.2) +
  geom_line(size=1) + geom_point(size=2) +
  geom_hline(data=tibble(Metric=c("MeanTE","MeanEE"), thr=c(0.90,0.85)),
             aes(yintercept=thr), linetype="dashed", color="grey40", inherit.aes=FALSE) +
  facet_wrap(~Metric, scales="fixed",
             labeller=as_labeller(c(MeanTE="Avg TE",MeanEE="Avg EE"))) +
  scale_color_brewer(palette="Set2") +
  scale_fill_brewer(palette="Set2") +
  labs(title="Yearly Trends of Avg TE & Avg EE by SizeGroup",
       subtitle="Dashed: TE=0.90, EE=0.85", x="Year", y="Average Score",
       color="SizeGroup", fill="SizeGroup") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom",
        plot.title=element_text(face="bold",hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))
ggsave(tr, filename = file.path(out_dir, "4.1.3_Trend_Avg_TE_EE.png"),
       width=10, height=6, dpi=300)

# 6.11.7 Cluster scatter (TE vs EE) export
cs <- ggplot(plot_cl, aes(x=TE,y=EE,color=Cluster)) +
  geom_point(shape=16,size=2,alpha=0.6) +
  geom_vline(xintercept=0.90, linetype="dotted", color="grey50") +
  geom_hline(yintercept=0.85, linetype="dotted", color="grey50") +
  geom_point(data=subset(plot_cl,TE_ci_low>0.95),
             aes(x=TE,y=EE), shape=8, color="red", size=3) +
  facet_wrap(~Year) +
  scale_color_brewer(palette="Set2") +
  labs(title="Scatter of TE vs EE by Year",
       subtitle="Clusters; ★ if TE_ci_low>0.95; dotted: TE=0.90, EE=0.85",
       x="TE",y="EE",color="Cluster") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold",hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.position="right")
ggsave(cs, filename = file.path(out_dir, "4.1.3_Cluster_Scatter_TE_EE.png"),
       width=10, height=6, dpi=300)

# 6.11.8 Input scatter (Beds vs Cost, color=TE) export
is <- ggplot(input_scatter,
             aes(x=Beds, y=OperatingCost_kGBP, color=TE, shape=SizeGroup)) +
  geom_point(alpha=0.7,size=2) +
  scale_color_viridis_c(option="plasma", labels=scales::percent) +
  facet_wrap(~Year) +
  labs(title="Beds vs Operating Cost by Year (color=TE, shape=SizeGroup)",
       x="Beds", y="Operating Cost (kGBP)", color="TE",shape="SizeGroup") +
  theme_minimal()
ggsave(is, filename = file.path(out_dir, "4.1.3_Input_Scatter_Beds_Cost_TE.png"),
       width=10, height=6, dpi=300)

# 6.11.9 Frontier scatter (CPD vs Discharges) export
fs2 <- ggplot(input_scatter,
              aes(x=OperatingCost_kGBP/Discharges, y=Discharges,
                  color=SizeGroup, shape=SizeGroup)) +
  geom_point(alpha=0.5) +
  facet_wrap(~Year) +
  geom_point(data=dplyr::filter(input_scatter,TE==1),
             aes(x=OperatingCost_kGBP/Discharges,y=Discharges),
             color="blue",shape=17,size=2,inherit.aes=FALSE) +
  geom_point(data=dplyr::filter(input_scatter,EE==1),
             aes(x=OperatingCost_kGBP/Discharges,y=Discharges),
             color="red",shape=15,size=2,inherit.aes=FALSE) +
  labs(title="Frontier Scatter: CPD vs Discharges",
       subtitle="▲ TE=1 (blue), ■ EE=1 (red)",
       x="Cost per Discharge",y="Discharges") +
  theme_minimal()
ggsave(fs2, filename = file.path(out_dir, "4.1.3_Frontier_Scatter_CPD_vs_Discharges.png"),
       width=10, height=6, dpi=300)

# 6.11.10 Radar chart (latest year) export
rad <- ggplot(radar_df, aes(x=Metric, y=Mean, group=SizeGroup, color=SizeGroup)) +
  geom_polygon(fill=NA,size=1) + geom_point(size=2) + coord_polar() +
  labs(title=paste("Radar Chart of Mean Efficiencies in", latest_year),
       subtitle="Metrics: TE | PTE | SE | EE | AE") +
  theme_minimal()
ggsave(rad, filename = file.path(out_dir, paste0("4.1.3_Radar_", latest_year, ".png")),
       width=8, height=8, dpi=300)

# 6.11.11 Lollipop export (Top 5 per SizeGroup, 2024)
ggsave(
  lop_plot,
  filename = file.path(out_dir, paste0("4.1.3_Lollipop_Top5_PerGroup_Cost_Cut_", target_year, ".png")),
  width = 8, height = 8, dpi = 300
)

# 6.11.12 Heatmaps & Trends export (slacks)
hm_in <- ggplot(hm_in_df, aes(Input, reorder(Hospital,CutPct), fill=CutPct)) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours=c("#f7f7f7","#fdbb84","#e34a33"),
                       labels=scales::percent, na.value="grey90") +
  facet_wrap(~SizeGroup,ncol=1,scales="free_y") +
  labs(title=paste0("Input Slack % Heatmap (", latest_year, ")"),
       x="Input",y="Hospital",fill="Slack %") +
  theme(axis.text.y=element_text(size=6))
ggsave(hm_in, filename = file.path(out_dir, "4.1.3_Heatmap_Input_Slack.png"),
       width=8, height=10, dpi=300)

trend_in <- p_in + labs(title="Mean Input Slack% Trends (2019–2024)")
ggsave(trend_in, filename = file.path(out_dir, "4.1.3_Trend_Input_Slack.png"),
       width=10, height=6, dpi=300)

hm_out <- ggplot(hm_out_df, aes(Output, reorder(Hospital,SlackPct), fill=SlackPct)) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours=c("#f7f7f7","#9ecae1","#08519c"),
                       labels=scales::percent, na.value="grey90") +
  facet_wrap(~SizeGroup,ncol=1,scales="free_y") +
  labs(title=paste0("Output Slack % Heatmap (", latest_year, ")"),
       x="Output",y="Hospital",fill="Slack %") +
  theme(axis.text.y=element_text(size=6))
ggsave(hm_out, filename = file.path(out_dir, "4.1.3_Heatmap_Output_Slack.png"),
       width=8, height=10, dpi=300)

trend_out <- p_out_final + labs(title="Mean Output Slack% Trends (2019–2024)")
ggsave(trend_out, filename = file.path(out_dir, "4.1.3_Trend_Output_Slack.png"),
       width=10, height=6, dpi=300)
