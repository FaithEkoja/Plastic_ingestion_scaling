# =============================================================================
# Plastic Ingestion Window & Size-Spectrum Exposure Model
# Ekoja-Smah, Reiss, Robertson, Reuangkhuma & Perkins (in prep) — Global Change Biology
#
# Built on the TRUSTED upper edge only (longest plastic ingested ~ body size).
# The lower edge of the realised diet is the environmental floor x_min: every
# organism incidentally ingests micro- and nanoplastic (filtration, gut
# retention, trophic transfer), so the ingestion window is [x_min, w_hi(L)].
# The fitted "shortest ingested" edge is NOT used — it is detection/ascertainment
# limited (large taxa are typically scored for macroplastic only), so it cannot
# anchor the model (see Model_rebuild_plan_window_fullrange.md §11).
#
# Two elements:
#   (1) Ingestion window [x_min, w_hi(L)] vs the environmental plastic size
#       distribution -> the window lies inside the plastic range across the whole
#       body-size spectrum; large taxa ingest micro-to-macroplastic. No size class
#       escapes plastic by size-selective feeding.
#   (2) Upper-edge food exposure E(L): the factor by which the plastic-to-food
#       number ratio at the LARGEST ingested particle, w_hi(L), exceeds that at the
#       base of the size spectrum, x_min. Because plastic (lambda~2.36) is a
#       shallower spectrum than food (lambda=4), E(L) = (w_hi/x_min)^(lambda_food -
#       lambda_plastic) rises steeply with body size. Larger organisms ingest into
#       a size range where their natural food has thinned out but plastic has not.
#       E(L) is concentration-free (C_plastic, C_food cancel) and uses only the
#       upper edge + the two spectral slopes.
#
# Decisions (see plan §3, §10, §11): x_max = 724 mm; x_min = 1 um central
# {20, 0.1 um sensitivity}; lambda extrapolated full-range (envelope bounds it);
# by-number spectra.
# =============================================================================

# -- 0. Packages & paths -------------------------------------------------------

required_pkgs <- c("tidyverse", "patchwork", "readxl")
to_install    <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(readxl)
})

OUT_DIR   <- "~/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Research/PhD students/Faith/Faith meta-analysis/Modeling component"
OUT_DIR   <- normalizePath(OUT_DIR, mustWork = TRUE)
DATA_DIR  <- dirname(OUT_DIR)
DATA_XLSX <- file.path(DATA_DIR, "GCB submission/Meta-analysis data_V6_21.7.26.xlsx")

cat("Output directory:", OUT_DIR, "\nData file       :", DATA_XLSX, "\n\n")


# =============================================================================
# SECTION 1 — Build df.all and fit the upper window edge
# =============================================================================
# Mirrors Meta-analysis data_V6_21.7.26.R (canonical build). The model uses
# the longest-ingested ("upper") edge only; the shortest edge is fitted for
# reference and reporting but NOT used downstream.

df_raw <- read_excel(DATA_XLSX)

df_sp <- df_raw %>%
  group_by(Test_animal, Settings) %>%
  summarise(max_MP_size_um = max(Largest_MP_size_um, na.rm = TRUE),
            min_MP_size_um = min(Smallest_MP_um,     na.rm = TRUE),
            av_animal_size_um = mean(Animal_size_um,  na.rm = TRUE),
            .groups = "drop") %>%
  transmute(species_binomial = Test_animal,
            animal_body_length_mm        = av_animal_size_um / 1000,
            longest_ingested_plastic_mm  = max_MP_size_um   / 1000,
            shortest_ingested_plastic_mm = min_MP_size_um   / 1000,
            study = "This paper", setting = Settings)

df_jams <- read_excel(DATA_XLSX, sheet = 2) %>%
  transmute(species_binomial, animal_body_length_mm,
            longest_ingested_plastic_mm, shortest_ingested_plastic_mm,
            study = "Jams et al.", setting = "Field")

df.all <- bind_rows(df_jams, df_sp) %>%
  mutate(log10_largest_MP  = log10(longest_ingested_plastic_mm),
         log10_smallest_MP = suppressWarnings(log10(shortest_ingested_plastic_mm)),
         log10_animal_size = log10(animal_body_length_mm)) %>%
  filter(is.finite(log10_animal_size))

# Upper edge (USED): longest plastic ingested ~ body length + setting (field shown)
lm.overall <- lm(log10_largest_MP ~ log10_animal_size + setting, data = df.all)
aL <- coef(lm.overall)[["(Intercept)"]]; bL <- coef(lm.overall)[["log10_animal_size"]]

# Lower edge: shortest plastic ingested ~ body length + setting. This is the
# PREY-SIZE-MATCHED lower bound of the realised window. It is ascertainment-biased
# (high) at large body size, so used as a CONSERVATIVE lower bound on enrichment,
# bracketed against the ambient-floor case (x_min). Field intercept.
lm.overall.small <- lm(log10_smallest_MP ~ log10_animal_size + setting, data = df.all)
aS <- coef(lm.overall.small)[["(Intercept)"]]; bS <- coef(lm.overall.small)[["log10_animal_size"]]

# Bayesian major-axis upper edge (sensitivity bound; Bayesian_MA_mixed.R)
aL_MA <- -0.96; bL_MA <- 0.92   # slope 95% CrI 0.825-1.043; coincides with Jams (0.934)

cat("-- Section 1: upper window edge (USED) --\n")
cat(sprintf("  w_hi  OLS: intercept %7.4f  slope %.4f  R2 %.3f  n %d\n",
            aL, bL, summary(lm.overall)$r.squared, length(lm.overall$residuals)))
cat(sprintf("  w_hi  MA : intercept %7.4f  slope %.4f  (sensitivity)\n", aL_MA, bL_MA))
cat(sprintf("  [reference, NOT used] shortest edge slope %.4f, R2 %.3f\n",
            coef(lm.overall.small)[["log10_animal_size"]], summary(lm.overall.small)$r.squared))
cat(sprintf("  Data ranges: body %.3f - %.0f mm | longest plastic max %.1f mm\n\n",
            min(df.all$animal_body_length_mm), max(df.all$animal_body_length_mm),
            max(df.all$longest_ingested_plastic_mm)))


# =============================================================================
# SECTION 2 — Food size spectrum (Hatton et al. 2021)
# =============================================================================
# Number-density-in-length exponent: lambda_food = 3*alpha_mass + 1, where
# alpha_mass is Hatton's slope of abundance-per-log10-mass-class vs mass.
# Hatton 2021: alpha_mass = 1.04 (95% CI 1.02-1.05) -> lambda_food = 4.12
# (95% CI 4.06-4.15). Central uses Hatton's point estimate (4.12); the Sheldon
# IDEAL (slope exactly -1) corresponds to lambda_food = 4.0.
# Sensitivity envelope (cross-study, Hatton Table S3, n=325): regional mass-space
# slopes span roughly -0.88 to -1.13 -> lambda_food = 3.6 - 4.4.
LAMBDA_FOOD       <- 4.12          # Hatton 2021 point estimate (alpha_mass = 1.04)
LAMBDA_FOOD_RANGE <- c(3.6, 4.4)  # Hatton 2021 cross-study envelope (Table S3)
cat(sprintf("-- Section 2: lambda_food = %.1f (sensitivity %.1f-%.1f, Hatton 2021 Table S3)\n\n",
            LAMBDA_FOOD, LAMBDA_FOOD_RANGE[1], LAMBDA_FOOD_RANGE[2]))


# =============================================================================
# SECTION 3 — Environmental plastic size distribution (full range)
# =============================================================================
# Power law f(x) ~ x^-lambda over [x_min, x_max]. lambda_plastic from Kooi et al.
# 2021 surface-water compartments (pooled 2.36; envelope 2.07-3.25). Extrapolated
# above 5 mm (envelope bounds it; Kooi & Koelmans 2019).
lambda_MP_table <- data.frame(
  compartment = c("marine_surface","freshwater_surface","marine_sediment",
                  "freshwater_sediment","wastewater","biota_freshwater"),
  lambda_MP   = c(2.07, 2.64, 2.57, 3.25, 2.54, 2.59))

LAMBDA_MP     <- mean(lambda_MP_table$lambda_MP[
                   lambda_MP_table$compartment %in% c("marine_surface","freshwater_surface")])
LAMBDA_MP_LWR <- min(lambda_MP_table$lambda_MP)   # 2.07 flattest -> largest E
LAMBDA_MP_UPR <- max(lambda_MP_table$lambda_MP)   # 3.25 steepest -> smallest E

PLASTIC_MIN_UM <- 1        # central environmental floor (incidental nano/small-MP)
PLASTIC_MAX_UM <- 724000   # 724 mm = largest plastic ingested in the data
XMIN_SCENARIOS <- c(20, 1, 0.1, 0.005)  # 20=FTIR limit, 1=central, 0.1=nano-inclusive,
                                        # 0.005=smallest plastic ingested in our data (5 nm)

cat(sprintf("-- Section 3: lambda_plastic = %.3f (envelope %.2f-%.2f); range %.1f um - %.0f mm\n\n",
            LAMBDA_MP, LAMBDA_MP_LWR, LAMBDA_MP_UPR, PLASTIC_MIN_UM, PLASTIC_MAX_UM/1000))


# =============================================================================
# SECTION 4 — Core functions (upper edge + one-sided window)
# =============================================================================
# w_hi_um()        : largest plastic ingested (um) at body length L, clipped to x_max.
# prop_ingestible(): fraction of the plastic number-spectrum within [x_min, w_hi].
# E_of_L()         : upper-edge food exposure = plastic:food ratio at w_hi divided by
#                    the plastic:food ratio at x_min = (w_hi/x_min)^(lambda_food -
#                    lambda_plastic). Concentration-free; uses only the upper edge.

w_hi_um <- function(L_mm, a, b, x_min = PLASTIC_MIN_UM, x_max = PLASTIC_MAX_UM) {
  w <- 10^(a + b * log10(L_mm)) * 1000
  pmin(pmax(w, x_min), x_max)
}

prop_ingestible <- function(L_mm, a, b, lambda, x_min = PLASTIC_MIN_UM, x_max = PLASTIC_MAX_UM) {
  w_hi <- w_hi_um(L_mm, a, b, x_min, x_max)
  if (abs(lambda - 1) < 1e-9) {
    num <- log(w_hi) - log(x_min); den <- log(x_max) - log(x_min)
  } else {
    e <- 1 - lambda; num <- w_hi^e - x_min^e; den <- x_max^e - x_min^e
  }
  pmax(0, num / den)
}

# Upper-edge food exposure: how much more plastic-skewed (vs food) the largest
# ingested particle is, relative to the base of the spectrum. Equivalently the
# within-diet enrichment of plastic over food from x_min up to w_hi.
E_of_L <- function(L_mm, a, b, lambda_plastic = LAMBDA_MP, lambda_food = LAMBDA_FOOD,
                   x_min = PLASTIC_MIN_UM, x_max = PLASTIC_MAX_UM) {
  w_hi <- w_hi_um(L_mm, a, b, x_min, x_max)
  (w_hi / x_min)^(lambda_food - lambda_plastic)
}

# E_fixed: enrichment from the ambient floor x_min up to the ceiling w_hi
# (ambient-inclusive; UPPER bound). E_window: enrichment across the empirical
# prey-matched window, from the shortest edge w_lo(L) up to w_hi(L) (LOWER bound,
# conservative because w_lo is ascertainment-biased high at large body size).
w_lo_um <- function(L_mm) w_hi_um(L_mm, aS, bS)   # shortest-edge ceiling reuse of clip logic

E_OLS <- function(L, lam = LAMBDA_MP, lam_f = LAMBDA_FOOD, x_min = PLASTIC_MIN_UM)
  E_of_L(L, aL, bL, lam, lam_f, x_min)
E_window <- function(L, lam = LAMBDA_MP, lam_f = LAMBDA_FOOD)
  (w_hi_um(L, aL, bL) / w_lo_um(L))^(lam_f - lam)
E_MA  <- function(L, lam = LAMBDA_MP, lam_f = LAMBDA_FOOD, x_min = PLASTIC_MIN_UM)
  E_of_L(L, aL_MA, bL_MA, lam, lam_f, x_min)

cat("-- Section 4: function test --\n")
for (L in c(0.01, 0.1, 1, 10, 100, 1000)) {
  cat(sprintf("  L=%8.2f mm | w_hi=%9.1f um | prop_ingestible=%.3f | E_OLS=%.3g\n",
              L, w_hi_um(L, aL, bL), prop_ingestible(L, aL, bL, LAMBDA_MP), E_OLS(L)))
}
cat("\n")


# =============================================================================
# SECTION 5 — Run over the body-size gradient
# =============================================================================
body_sizes_mm <- 10^seq(-2, 4, length.out = 300)   # 0.01 - 10,000 mm

# For the ribbon, propagate BOTH spectra envelopes jointly:
#   upper bound = flattest plastic (LAMBDA_MP_LWR) + steepest food (LAMBDA_FOOD_RANGE[2])
#   lower bound = steepest plastic (LAMBDA_MP_UPR) + flattest food (LAMBDA_FOOD_RANGE[1])
model_df <- data.frame(
  body_size_mm    = body_sizes_mm,
  w_hi_um         = w_hi_um(body_sizes_mm, aL, bL),
  prop_ingestible = prop_ingestible(body_sizes_mm, aL, bL, LAMBDA_MP),
  prop_lwr        = prop_ingestible(body_sizes_mm, aL, bL, LAMBDA_MP_LWR),  # flattest plastic -> lowest P
  prop_upr        = prop_ingestible(body_sizes_mm, aL, bL, LAMBDA_MP_UPR),  # steepest plastic -> highest P
  E_central       = E_OLS(body_sizes_mm, LAMBDA_MP, LAMBDA_FOOD),    # ambient-floor (upper bound)
  E_hi            = E_OLS(body_sizes_mm, LAMBDA_MP_LWR, LAMBDA_FOOD_RANGE[2]),
  E_lo            = E_OLS(body_sizes_mm, LAMBDA_MP_UPR, LAMBDA_FOOD_RANGE[1]),
  E_window        = E_window(body_sizes_mm, LAMBDA_MP, LAMBDA_FOOD),  # prey-matched window (lower bound)
  E_window_lo     = E_window(body_sizes_mm, LAMBDA_MP_UPR, LAMBDA_FOOD_RANGE[1]),  # lower bound, lambda-low
  E_MA            = E_MA(body_sizes_mm, LAMBDA_MP, LAMBDA_FOOD)
)

cat("-- Section 5: upper-edge food exposure E(L): OLS vs MA --\n")
for (L in c(0.01, 0.1, 1, 10, 100, 1000, 10000)) {
  i <- which.min(abs(model_df$body_size_mm - L))
  cat(sprintf("  L = %8.2f mm | w_hi = %8.1f um | E_OLS = %10.3g | E_MA = %10.3g | prop_ingest = %.3f\n",
              L, model_df$w_hi_um[i], model_df$E_central[i], model_df$E_MA[i], model_df$prop_ingestible[i]))
}
E_fin <- model_df$E_central[is.finite(model_df$E_central)]
cat(sprintf("  Fold-increase in E across the spectrum (OLS): %.3g x (~%.0f orders)\n\n",
            tail(E_fin, 1) / E_fin[1], log10(tail(E_fin, 1) / E_fin[1])))

cat("-- Sanity checks --\n")
cat(sprintf("  E(L) monotonic increasing: %s\n", all(diff(E_fin) >= -1e-8)))
cat(sprintf("  E > 1 across spectrum (ceiling more plastic-skewed than base): %s\n", all(E_fin > 1)))
cat(sprintf("  prop_ingestible rises to %.3f at largest body size (ingest ~whole plastic pool)\n",
            tail(model_df$prop_ingestible, 1)))
# exact analytic crossing where w_hi(L) = 5 mm (not the grid approximation)
L_5mm <- 10^((log10(5) - aL) / bL)
cat(sprintf("  w_hi exceeds 5 mm above body length %.1f mm (ingests meso/macroplastic)\n", L_5mm))
cat(sprintf("    [Bayesian MA edge: 5 mm crossing at %.1f mm]\n\n", 10^((log10(5) - aL_MA) / bL_MA)))


# =============================================================================
# SECTION 6 — Sensitivity: E(L) across lambda_plastic
# =============================================================================
lambda_plastic_range <- c(2.07, 2.36, 2.66, 2.96, 3.25)
sens_lambda_df <- map_dfr(lambda_plastic_range, function(lam)
  data.frame(body_size_mm = body_sizes_mm, lambda_MP = lam, E = E_OLS(body_sizes_mm, lam)))

fig_sens_lambda <- ggplot(sens_lambda_df %>% filter(E > 0),
                          aes(body_size_mm, E, colour = factor(lambda_MP), group = factor(lambda_MP))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.9) +
  scale_x_log10(breaks = 10^(-2:4), labels = c("0.01","0.1","1","10","100","1,000","10,000")) +
  scale_y_log10() +
  scale_colour_grey(start = 0.85, end = 0.1, name = expression(lambda[plastic])) +
  labs(x = "Animal body length (mm)",
       y = "Exposure index, E(L)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")


# =============================================================================
# SECTION 6a — Sensitivity: E(L) across lambda_food
# =============================================================================
# Hatton 2021 Table S3 cross-study heterogeneity: mass-space slopes typically
# span -0.88 to -1.13, mapping to lambda_food = 3.6 - 4.4 via lambda = 3*alpha_mass + 1.
lambda_food_range <- seq(3.6, 4.4, by = 0.2)
sens_food_df <- map_dfr(lambda_food_range, function(lf)
  data.frame(body_size_mm = body_sizes_mm, lambda_food = lf,
             E = E_OLS(body_sizes_mm, LAMBDA_MP, lam_f = lf)))

fig_sens_food <- ggplot(sens_food_df %>% filter(E > 0),
                        aes(body_size_mm, E, colour = factor(lambda_food), group = factor(lambda_food))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.9) +
  scale_x_log10(breaks = 10^(-2:4), labels = c("0.01","0.1","1","10","100","1,000","10,000")) +
  scale_y_log10() +
  scale_colour_grey(start = 0.85, end = 0.1, name = expression(lambda[food])) +
  labs(x = "Animal body length (mm)",
       y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")


# =============================================================================
# SECTION 6b — Sensitivity: lower bound x_min
# =============================================================================
# x_min is the base of the diet (environmental floor). It sets the reference for
# E(L); a lower floor (more nano included) raises E. Shape (monotone rising)
# unchanged.
xmin_df <- map_dfr(XMIN_SCENARIOS, function(xm)
  data.frame(body_size_mm = body_sizes_mm, x_min_um = xm, E = E_OLS(body_sizes_mm, LAMBDA_MP, x_min = xm)))

cat("-- Section 6b: x_min sensitivity — E(L) at representative body sizes --\n")
xmin_summary <- map_dfr(c(0.1, 1, 10, 100), function(L) {
  i <- which.min(abs(body_sizes_mm - L))
  xmin_df %>% filter(abs(body_size_mm - body_sizes_mm[i]) < 1e-9) %>%
    transmute(body_size_mm = round(L, 2), x_min_um, E = signif(E, 3))
}) %>% pivot_wider(names_from = x_min_um, values_from = E, names_prefix = "xmin_")
print(xmin_summary); cat("\n")

fig_sens_xmin <- ggplot(xmin_df %>% filter(E > 0),
                        aes(body_size_mm, E, colour = factor(x_min_um), group = factor(x_min_um))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) +
  scale_x_log10(breaks = 10^(-2:4), labels = c("0.01","0.1","1","10","100","1,000","10,000")) +
  scale_y_log10() +
  scale_colour_grey(start = 0.75, end = 0.1, name = expression(x[min]~"(um)")) +
  labs(x = "Animal body length (mm)",
       y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")


# =============================================================================
# SECTION 6c — Supplementary: ingestible fraction P(L)
# =============================================================================
# Proportion of the plastic number-spectrum lying within an organism's window
# [x_min, w_hi(L)]. Because the spectrum is bottom-heavy (lambda > 1), P(L) is
# already high at the smallest body sizes (the small-particle pile-up dominates
# the integral) and rapidly saturates near 1. Supports the "no size class
# escapes plastic by particle size" Results claim.

fig_propP <- ggplot(model_df, aes(body_size_mm, prop_ingestible)) +
  geom_ribbon(aes(ymin = pmin(prop_lwr, prop_upr), ymax = pmax(prop_lwr, prop_upr)),
              fill = "firebrick", alpha = 0.15) +
  geom_line(linewidth = 1.1, colour = "firebrick") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  scale_x_log10(breaks = 10^(-2:4), labels = c("0.01","0.1","1","10","100","1,000","10,000")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(x = "Animal body length (mm)",
       y = "P(L): ingestible fraction\nof the plastic number-spectrum",
       caption = paste0("P(L) = share of the plastic number-spectrum (1 um - 724 mm) lying within ",
                        "the\nwindow [x_min, w_hi(L)]; shaded: lambda_plastic envelope ",
                        LAMBDA_MP_LWR, "-", LAMBDA_MP_UPR, ", Kooi et al. 2021.")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 7, colour = "grey50"))


# =============================================================================
# SECTION 7 — Publication figures
# =============================================================================
NANO_MM <- 0.001; DETECT_MM <- 0.02; MP_MESO <- 5; MESO_MAC <- 25

# -- Panel A: ingestion window [x_min, w_hi] vs body size ----------------------
Lseq <- 10^seq(-2, 4, length.out = 400)
win_df <- data.frame(body_mm = Lseq,
                     w_hi_mm = w_hi_um(Lseq, aL, bL) / 1000,
                     w_lo_mm = w_lo_um(Lseq) / 1000)
df_obs <- df.all %>% filter(is.finite(log10_animal_size), is.finite(log10_largest_MP))
x_lo <- min(Lseq); x_hi <- max(Lseq)

PMIN_MM <- PLASTIC_MIN_UM/1000; PMAX_MM <- PLASTIC_MAX_UM/1000   # 0.001 - 724 mm
x_arr   <- x_lo*1.10; x_lab <- x_lo*1.45                         # arrow + label x positions
band_seg <- function(ylo, yhi)
  annotate("segment", x = x_arr, xend = x_arr, y = ylo, yend = yhi,
           arrow = arrow(ends = "both", length = unit(1.1, "mm")),
           colour = "grey35", linewidth = 0.4)

fig_window <- ggplot(win_df, aes(x = body_mm)) +
  # GESAMP (2016) size-class bands (micro / meso / macro)
  annotate("rect", xmin = x_lo, xmax = x_hi, ymin = PMIN_MM, ymax = MP_MESO,  fill = "grey94") +
  annotate("rect", xmin = x_lo, xmax = x_hi, ymin = MP_MESO, ymax = MESO_MAC, fill = "grey88") +
  annotate("rect", xmin = x_lo, xmax = x_hi, ymin = MESO_MAC, ymax = PMAX_MM, fill = "grey82") +
  # realised ingestion window: floor up to the ceiling w_hi(L)
  geom_ribbon(aes(ymin = PMIN_MM, ymax = w_hi_mm), fill = "firebrick", alpha = 0.18) +
  # size-class extent arrows + labels (left margin, in empty space)
  band_seg(PMIN_MM, MP_MESO) + band_seg(MP_MESO, MESO_MAC) + band_seg(MESO_MAC, PMAX_MM) +
  annotate("text", x = x_lab, y = 1.5,  hjust = 0, size = 2.6, colour = "grey30", label = "Microplastic") +
  annotate("text", x = x_lab, y = 11.5, hjust = 0, size = 2.6, colour = "grey30", label = "Mesoplastic") +
  annotate("text", x = x_lab, y = 150,  hjust = 0, size = 2.6, colour = "grey30", label = "Macroplastic") +
  # observed data points hashed out for clarity (un-comment to show):
  # geom_point(data = df_obs, aes(x = 10^log10_animal_size, y = 10^log10_largest_MP),
  #            colour = "grey25", alpha = 0.40, size = 1.1) +
  # geom_point(data = subset(df.all, is.finite(log10_smallest_MP) & is.finite(log10_animal_size)),
  #            aes(x = 10^log10_animal_size, y = 10^log10_smallest_MP),
  #            colour = "grey60", alpha = 0.40, size = 1.1) +
  geom_line(aes(y = w_hi_mm), colour = "firebrick", linewidth = 1.1) +
  geom_line(aes(y = w_lo_mm), colour = "firebrick", linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = 9000, y = w_hi_um(9000, aL, bL)/1000*2.4, hjust = 1, size = 2.4, colour = "firebrick",
           label = "largest plastic ingested") +
  annotate("text", x = 9000, y = w_lo_um(9000)/1000*0.16, hjust = 1, size = 2.4, colour = "firebrick",
           label = "smallest plastic ingested") +
  annotate("text", x = 0.012, y = 0.0014, hjust = 0, size = 2.4, colour = "grey35",
           label = "shaded: full ingestible range") +
  scale_x_log10(breaks = 10^(-2:4), labels = c("0.01","0.1","1","10","100","1,000","10,000")) +
  scale_y_log10(breaks = 10^(-3:3), labels = c("0.001","0.01","0.1","1","10","100","1,000")) +
  labs(x = "Body length (mm)", y = "Ingestible plastic size (mm)",
       subtitle = "Ingestion window vs body size") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())

# -- Panel B: size spectra with one-sided windows ------------------------------
taxa <- tibble::tibble(
  L_mm  = c(0.05, 1, 10, 100, 1000),
  label = c("0.05 mm (ciliate)", "1 mm (zooplankton)", "10 mm (macroinvertebrate)",
            "100 mm (small fish)", "1,000 mm (mammal)")) %>%
  mutate(w_hi = w_hi_um(L_mm, aL, bL))

particle_sizes_um <- 10^seq(log10(PLASTIC_MIN_UM), log10(PLASTIC_MAX_UM), length.out = 400)
spectra_df <- data.frame(
  size_um = particle_sizes_um,
  food    = particle_sizes_um^(-LAMBDA_FOOD),
  plastic = particle_sizes_um^(-LAMBDA_MP),
  pl_lo   = particle_sizes_um^(-LAMBDA_MP_UPR),         # steepest plastic -> lower at large x
  pl_hi   = particle_sizes_um^(-LAMBDA_MP_LWR),         # flattest plastic -> higher at large x
  fd_lo   = particle_sizes_um^(-LAMBDA_FOOD_RANGE[2]),  # steepest food (lambda=4.4)
  fd_hi   = particle_sizes_um^(-LAMBDA_FOOD_RANGE[1])   # flattest food (lambda=3.6)
) %>% mutate(across(-size_um, ~ . / max(food)))

taxa <- taxa %>% mutate(y_pos = 10^seq(2, 9, length.out = n()))

  # Kooi 2021 empirical microplastic range (FTIR), 20 um - 5 mm; outside this we
  # extrapolate the power law (solid for empirical, dashed for extrapolated).
KOOI_LO_UM <- 20; KOOI_HI_UM <- 5000
spectra_emp     <- spectra_df %>% filter(size_um >= KOOI_LO_UM, size_um <= KOOI_HI_UM)
spectra_extr_lo <- spectra_df %>% filter(size_um <= KOOI_LO_UM)
spectra_extr_hi <- spectra_df %>% filter(size_um >= KOOI_HI_UM)

fig_spectra <- ggplot(spectra_df, aes(size_um)) +
  geom_ribbon(aes(ymin = pl_lo, ymax = pl_hi), fill = "firebrick", alpha = 0.12) +
  geom_ribbon(aes(ymin = fd_lo, ymax = fd_hi), fill = "steelblue", alpha = 0.12) +
  geom_line(aes(y = food),    colour = "steelblue", linewidth = 1.1) +
  geom_line(data = spectra_extr_lo, aes(y = plastic),
            colour = "firebrick", linewidth = 1.1, linetype = "dashed") +
  geom_line(data = spectra_emp,     aes(y = plastic),
            colour = "firebrick", linewidth = 1.1) +
  geom_line(data = spectra_extr_hi, aes(y = plastic),
            colour = "firebrick", linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = c(5000, 25000), linetype = "dotted", colour = "grey60", linewidth = 0.4) +
  # one-sided windows: floor (x_min) up to each taxon's ceiling w_hi
  geom_segment(data = taxa, aes(x = PLASTIC_MIN_UM, xend = w_hi, y = y_pos, yend = y_pos),
               colour = "grey25", linewidth = 1.1, lineend = "round") +
  geom_point(data = taxa, aes(x = w_hi, y = y_pos), colour = "grey25", size = 1.6) +
  geom_text(data = taxa, aes(x = PLASTIC_MIN_UM, y = y_pos, label = label),
            hjust = 0, vjust = -0.35, size = 2.2, colour = "grey25", lineheight = 0.85) +
  # in-panel slope labels: parked in empty space well clear of both curves and ribbons
  annotate("text", x = 3, y = 1e-18, label = paste0("Food (lambda~", round(LAMBDA_FOOD,2), ")"),
           colour = "steelblue", size = 3, hjust = 0) +
  annotate("text", x = 3, y = 1e0,   label = paste0("Plastic (lambda~", round(LAMBDA_MP,2), ")"),
           colour = "firebrick", size = 3, hjust = 0) +
  scale_x_log10(breaks = c(1,10,100,1000,10000,1e5,1e6),
                labels = c("1 um","10","100","1 mm","10","100","1 m"), expand = c(0.02, 0)) +
  scale_y_log10(breaks = 10^c(-20,-15,-10,-5,0,5,10),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_cartesian(ylim = c(1e-22, 1e10)) +
  labs(x = "Particle size", y = "Relative abundance\n(number-spectrum, normalised)",
       subtitle = "Size spectra & ingestion windows") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())

# -- Panel C: upper-edge food exposure E(L) ------------------------------------
# Two bounds: E_window (prey-matched window, lower) and E_central (ambient floor,
# upper). The band between them is the exposure bracket; the firebrick ribbon is
# the lambda-envelope on the upper bound. (MA edge hashed out for now.)
fig_E <- ggplot(model_df %>% filter(E_central > 0), aes(body_size_mm, E_central)) +
  # full plausible band: size-selective lower bound (E_window, with its lambda-low
  # edge) up to the ambient-floor upper bound (E_hi, with its lambda-high edge)
  geom_ribbon(aes(ymin = E_window_lo, ymax = E_hi), fill = "firebrick", alpha = 0.13) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  annotate("text", x = 0.05, y = 3, hjust = 0, size = 2.8, colour = "grey50",
           label = "E = 1 (no enrichment)") +
  # geom_line(aes(y = E_MA), colour = "firebrick", linetype = "dashed", linewidth = 0.8) +  # MA edge (hashed out)
  geom_line(aes(y = E_central), colour = "firebrick", linewidth = 1.2) +
  geom_line(aes(y = E_window),  colour = "firebrick", linewidth = 1.0, linetype = "dashed") +
  annotate("text", x = 2.5, y = 6e7, hjust = 0, size = 2.5, colour = "firebrick",
           label = "incl. ambient ingestion\n(upper bound)") +
  annotate("text", x = 40, y = 6, hjust = 0, size = 2.5, colour = "firebrick",
           label = "size-selective feeding only\n(lower bound)") +
  scale_x_log10(breaks = 10^(-2:4), labels = c("0.01","0.1","1","10","100","1,000","10,000")) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  labs(x = "Animal body length (mm)",
       y = "Exposure index, E(L)",
       subtitle = "Plastic-rich at the diet ceiling") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# -- Combined main figure: A (window) on top, B (spectra) | C (E) below --------
fig_main <- fig_window / (fig_spectra | fig_E) +
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A")

# -- Combined sensitivity panel (Appendix S1, Fig. S4): a | b | c ---------------
# Each panel keeps its own colour legend (different variables), placed at the
# bottom so the panels stay equal width; shared y-axis label on the left only.
fig_sens_panel <- (fig_sens_lambda | fig_sens_food | fig_sens_xmin) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))


# =============================================================================
# SECTION 8 — Save figures
# =============================================================================
save_fig <- function(plot, filename, width, height) {
  path <- file.path(OUT_DIR, filename)
  ggsave(plot = plot, filename = path, device = "pdf", width = width, height = height, units = "mm")
  cat(sprintf("Saved: %s\n", path))
}

cat("\n-- Section 8: saving figures --\n")
save_fig(fig_main,        "Fig_plastic_exposure_model.pdf",     width = 200, height = 200)
save_fig(fig_window,      "Fig_plastic_window.pdf",             width = 170, height = 110)
save_fig(fig_spectra,     "Fig_plastic_spectra.pdf",            width = 170, height = 110)
save_fig(fig_E,           "Fig_plastic_exposure_E.pdf",         width = 170, height = 110)
save_fig(fig_propP,       "Fig_plastic_propP.pdf",                      width = 150, height = 110)
save_fig(fig_sens_lambda, "Fig_plastic_sensitivity_lambda_plastic.pdf", width = 150, height = 110)
save_fig(fig_sens_panel,  "Fig_S4_sensitivity_panels.pdf",              width = 270, height = 100)
save_fig(fig_sens_food,   "Fig_plastic_sensitivity_lambda_food.pdf",    width = 150, height = 110)
save_fig(fig_sens_xmin,   "Fig_plastic_sensitivity_xmin.pdf",            width = 150, height = 110)

cat("\nDone.\n")
