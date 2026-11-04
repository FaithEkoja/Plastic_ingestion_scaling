# NOTE: i've updated the raw data file since there was a naming anomaly:
#lower case 'laboratory' level (instead of upper case standard) 
# one level of 'Mesocosm' which I changed to 'Field', highlighted in red -  please check this makes sense


rm(list = ls()) # clear history

########################################################
## Packages and libraries                              #
########################################################

# Function to Install and Load R Packages
Install_And_Load <- function(Required_Packages)
{
  Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];

  if(length(Remaining_Packages))
  {
    install.packages(Remaining_Packages);
  }
  for(package_name in Required_Packages)
  {
    library(package_name,character.only=TRUE,quietly=TRUE);
  }
}

# Specify the list of required packages to be installed and load
Required_Packages=c("tidyverse", "lme4", "ggplot2", "ggpubr", "plyr", "dplyr", "visreg", "lmerTest", "performance", "readxl")

# Call the Function
Install_And_Load(Required_Packages);

########################################################
## ggplot attributes                                   #
########################################################

txt.a = element_text (size = 18)
txt.b = element_text (size = 14)

pub.tweeks <- theme_bw() + theme(panel.border = element_rect(colour = "black"),
                                 strip.text.x = txt.b,
                                 strip.text.y = txt.b,
                                 axis.title.x = txt.a, axis.title.y = txt.a,
                                 axis.text.x = txt.b, axis.text.y = txt.b)

########################################################
# SET PATH TO THE DATA ON YOUR DEVICE
########################################################

aFile<-file.choose() # SET YOUR DATA PATH

data.path = dirname(aFile)
setwd(data.path)

#################################################
# Read in dataframe and summarise to species level
#################################################

df <- read_excel("Meta-analysis data_V6_21.7.26.xlsx")

# Summarise to species level (one row per species x setting)
# This avoids pseudoreplication across studies and removes
# the confounding effect of MP shape (a property of experiments, not species)

# df.sp <- ddply(df, .(Test_animal, Settings), summarise,
#                max_MP_size_um = max(Largest_MP_size_um),
#                av_largest_MP_size_um = mean(Largest_MP_size_um),
#                min_MP_size_um = min(Smallest_MP_um),
#                n_largest = length(Largest_MP_size_um),
#                n_smallest = length(Smallest_MP_um),
#                av_animal_size_um = mean(Animal_size_um),
#                organism_group = first(Organism_group_std),
#                environment = first(Environment_type_std))


# Recode MP_shape to binary sphere vs non-sphere
df$shape_binary <- ifelse(df$MP_shape == "Sphere", "Sphere", "Non-sphere")

# Then in df.sp summaries, dominant shape becomes majority vote on binary
df.sp <- ddply(df, .(Test_animal, Settings), summarise,
               max_MP_size_um        = max(Largest_MP_size_um),
               av_largest_MP_size_um = mean(Largest_MP_size_um),
               min_MP_size_um        = min(Smallest_MP_um),
               n_largest             = length(Largest_MP_size_um),
               n_smallest            = length(Smallest_MP_um),
               av_animal_size_um     = mean(Animal_size_um),
               organism_group        = first(Organism_group_std),
               environment           = first(Environment_type_std),
               # Binary shape: majority vote where mixed, else unanimous
               dominant_shape_binary = names(sort(table(shape_binary), 
                                                  decreasing = TRUE))[1],
               n_shapes              = length(unique(MP_shape)))



# Convert to mm and log-transform
df.sp$longest_ingested_plastic_mm <- df.sp$max_MP_size_um / 1000
df.sp$shortest_ingested_plastic_mm <- df.sp$min_MP_size_um / 1000
df.sp$animal_body_length_mm <- df.sp$av_animal_size_um / 1000

df.sp$log10_largest_MP <- log10(df.sp$longest_ingested_plastic_mm)
df.sp$log10_smallest_MP <- log10(df.sp$shortest_ingested_plastic_mm)
df.sp$log10_animal_size <- log10(df.sp$animal_body_length_mm)

df.sp$study <- "This paper"
df.sp$species_binomial <- df.sp$Test_animal
df.sp$setting <- df.sp$Settings

#FAith: laboratory into only one Laboratory category and mesocosm as filed
df.sp$setting <- as.character(df.sp$setting)

df.sp$setting[tolower(df.sp$setting) == "laboratory"] <- "Laboratory"
df.sp$setting[tolower(df.sp$setting) == "mesocosm"] <- "Field"

df.sp$setting <- factor(df.sp$setting)
table(df.sp$setting)
cat("Our species-level data:", nrow(df.sp), "rows\n")
cat("setting:\n"); print(table(df.sp$setting))
cat("Organism groups:\n"); print(table(df.sp$organism_group))
cat("Environments:\n"); print(table(df.sp$environment))
cat("Shapes:\n"); print(table(df.sp$dominant_shape_binary))

table(df.sp$dominant_shape_binary, df.sp$setting)

#ADD THIS TO METHODS SECTION
#MP shape was not included as a covariate as it was strongly confounded 
#with experimental setting (88% of field records used non-sphere particles; 88% of laboratory records used spheres), 
#reflecting differences in particle availability rather than independent biological effects.

#Environment type was excluded from the model as it is implicitly accounted for by the species-level random intercept, 
#which absorbs baseline differences associated with taxonomic and environmental context.

#################################################
# Plot raw species-level data
#################################################

# Largest MPs ingested
ggplot(df.sp, aes(y = log10_largest_MP, x = log10_animal_size)) +
  geom_point(size = 1.5, alpha = 1) +
  ggtitle("Body length of organism and largest ingested MP (species-level)") +
  xlab("Organism body size ([log10] mm)") +
  ylab("Longest MP ingested ([log10] mm)") +
  geom_smooth(method = "lm", colour = "black") +
  stat_regline_equation(label.x = 0.5, label.y = 2.2) +
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0.5, label.y = 2) +
  theme_minimal()

# Smallest MPs ingested
ggplot(df.sp, aes(y = log10_smallest_MP, x = log10_animal_size)) +
  geom_point(size = 1.5, alpha = 1) +
  ggtitle("Body length of organism and smallest ingested MP (species-level)") +
  xlab("Organism body size ([log10] mm)") +
  ylab("Shortest MP ingested ([log10] mm)") +
  geom_smooth(method = "lm", colour = "black") +
  stat_regline_equation(label.x = 0.5, label.y = 2.2) +
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0.5, label.y = 2) +
  theme_minimal()

###########################################################
# Model our species-level data
###########################################################

# Species occasionally appear in both field and lab, so use species as random effect
# Test interaction first
lmer.us.int <- lmer(log10_largest_MP ~ log10_animal_size * setting + (1|species_binomial), data = df.sp)
anova(lmer.us.int) # check if interaction is significant

# Final model (drop interaction if NS)
lmer.us <- lmer(log10_largest_MP ~ log10_animal_size + setting + (1|species_binomial), data = df.sp)
summary(lmer.us)
anova(lmer.us)
model_performance(lmer.us)
check_model(lmer.us)
performance::check_collinearity(lmer.us)
confint(lmer.us)

#plot patterns
visreg(lmer.us, "log10_animal_size", by = "setting")
visreg(lmer.us, "setting")

###########################################################
# Smallest MP model (same structure)
###########################################################

lmer.us.small.int <- lmer(log10_smallest_MP ~ log10_animal_size * setting + (1|species_binomial), data = df.sp)
anova(lmer.us.small.int)

lmer.us.small <- lmer(log10_smallest_MP ~ log10_animal_size + setting + (1|species_binomial), data = df.sp)
summary(lmer.us.small)
anova(lmer.us.small)
model_performance(lmer.us.small)
confint(lmer.us.small)

visreg(lmer.us.small, "log10_animal_size", by = "setting")


###########################################################
# Figure 1: Our data (field vs lab) - largest MP
###########################################################

vis_data <- visreg(lmer.us, xvar = "log10_animal_size", by = "setting", plot = FALSE)
vis_fit <- vis_data$fit %>% as.data.frame()
vis_res <- vis_data$res %>% as.data.frame()

fig1 <- ggplot() +
  geom_point(data = vis_res,
             aes(x = log10_animal_size, y = visregRes, fill = setting),
             shape = 21, size = 2.5, stroke = 0.3, colour = "black", alpha = 0.8) +
  scale_x_continuous(limits = c(-2.5, 2), expand = expansion(mult = 0.0)) +
  scale_y_continuous(limits = c(-4, 1), expand = expansion(mult = 0.0)) +
  geom_ribbon(data = subset(vis_fit, setting == "Laboratory"),
              aes(x = log10_animal_size, ymin = visregLwr, ymax = visregUpr),
              fill = "grey90", alpha = 0.3, color = NA) +
  geom_ribbon(data = subset(vis_fit, setting == "Field"),
              aes(x = log10_animal_size, ymin = visregLwr, ymax = visregUpr),
              fill = "grey70", alpha = 0.3, color = NA) +
  geom_line(data = vis_fit,
            aes(x = log10_animal_size, y = visregFit,
                color = setting, linetype = setting),
            linewidth = 1) +
  labs(
    x = expression(log[10]~"Animal body length [mm]"),
    y = expression(log[10]~"Longest microplastic ingested [mm]")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  scale_fill_manual(values = c("Field" = "black", "Laboratory" = "white")) +
  scale_color_manual(values = c("Field" = "black", "Laboratory" = "black")) +
  scale_linetype_manual(values = c("Field" = "solid", "Laboratory" = "dashed")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + theme_classic(base_size = 11) +
  scale_fill_manual(name = "Setting", values = c("Field" = "black", "Laboratory" = "white")) +
  scale_color_manual(name = "Setting", values = c("Field" = "black", "Laboratory" = "black")) +
  scale_linetype_manual(name = "Setting", values = c("Field" = "solid", "Laboratory" = "dashed")) 

fig1
#ggsave(plot = fig1, device = pdf, filename = 'Fig1_largest_MP_field_vs_lab.pdf', width = 200, height = 100, dpi = 300, units = "mm")

###########################################################
# Figure 2: Our data (field vs lab) - smallest MP
###########################################################

vis_data_sm <- visreg(lmer.us.small, xvar = "log10_animal_size", by = "setting", plot = FALSE)
vis_fit_sm <- vis_data_sm$fit %>% as.data.frame()
vis_res_sm <- vis_data_sm$res %>% as.data.frame()

fig2 <- ggplot() +
  geom_point(data = vis_res_sm,
             aes(x = log10_animal_size, y = visregRes, fill = setting),
             shape = 21, size = 2.5, stroke = 0.3, colour = "black", alpha = 0.8) +
  xlim(-2.5, 2) + ylim(-5, 1) +
  geom_ribbon(data = subset(vis_fit_sm, setting == "Laboratory"),
              aes(x = log10_animal_size, ymin = visregLwr, ymax = visregUpr),
              fill = "grey90", alpha = 0.3, color = NA) +
  geom_ribbon(data = subset(vis_fit_sm, setting == "Field"),
              aes(x = log10_animal_size, ymin = visregLwr, ymax = visregUpr),
              fill = "grey70", alpha = 0.3, color = NA) +
  geom_line(data = vis_fit_sm,
            aes(x = log10_animal_size, y = visregFit,
                color = setting, linetype = setting),
            linewidth = 1) +
  labs(
    x = expression(log[10]~"Animal body length [mm]"),
    y = expression(log[10]~"Shortest microplastic ingested [mm]")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  scale_fill_manual(values = c("Field" = "black", "Laboratory" = "white")) +
  scale_color_manual(values = c("Field" = "black", "Laboratory" = "black")) +
  scale_linetype_manual(values = c("Field" = "solid", "Laboratory" = "dashed")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + theme_classic(base_size = 11) +
  scale_fill_manual(name = "Setting", values = c("Field" = "black", "Laboratory" = "white")) +
  scale_color_manual(name = "Setting", values = c("Field" = "black", "Laboratory" = "black")) +
  scale_linetype_manual(name = "Setting", values = c("Field" = "solid", "Laboratory" = "dashed")) 

fig2
#ggsave(plot = fig2, device = pdf, filename = 'Fig2_smallest_MP_field_vs_lab.pdf', width = 150, height = 100, dpi = 300, units = "mm")

###################
# Combine with Jams data
###################

df2 <- read_excel("Meta-analysis data_V6_21.7.26.xlsx", sheet = 2)
df2 <- df2[, c("species_binomial", "animal_body_length_mm", "longest_ingested_plastic_mm", "shortest_ingested_plastic_mm")]
df2$study <- rep("Jams et al.", length(df2$animal_body_length_mm))
df2$setting <- rep("Field", length(df2$animal_body_length_mm))

# Prepare our data in same format
df.ours <- df.sp[, c("species_binomial", "animal_body_length_mm", "longest_ingested_plastic_mm", "shortest_ingested_plastic_mm", "study", "setting")]

df.all <- rbind(df2, df.ours) # combine

# log-transform variables
df.all$log10_largest_MP <- log10(df.all$longest_ingested_plastic_mm)
df.all$log10_smallest_MP <- log10(df.all$shortest_ingested_plastic_mm)
df.all$log10_animal_size <- log10(df.all$animal_body_length_mm)

cat("n :", nrow(df.all), "rows\n")
cat("Organism groups:\n"); print(table(df.all$species_binomial))
length(unique(df.all$species_binomial))

################################################
# compare our fitted values to Jams data scaling
################################################

# ── This paper (lmer.us) ────────────────────────────────────────────────────
ci.us <- confint(lmer.us, parm = "beta_")  # all fixed effects at once

slope.us       <- fixef(lmer.us)["log10_animal_size"]
slope.us.lwr   <- ci.us["log10_animal_size", 1]
slope.us.upr   <- ci.us["log10_animal_size", 2]

int.field.us     <- fixef(lmer.us)["(Intercept)"]
int.field.us.lwr <- ci.us["(Intercept)", 1]
int.field.us.upr <- ci.us["(Intercept)", 2]

int.lab.us     <- fixef(lmer.us)["(Intercept)"] + fixef(lmer.us)["settingLaboratory"]
int.lab.us.lwr <- ci.us["(Intercept)", 1] + ci.us["settingLaboratory", 1]
int.lab.us.upr <- ci.us["(Intercept)", 2] + ci.us["settingLaboratory", 2]

# ── Jams et al. (lm) ────────────────────────────────────────────────────────
df.jams  <- subset(df.all, study == "Jams et al.")
lm.jams  <- lm(log10_largest_MP ~ log10_animal_size, data = df.jams)
ci.jams  <- confint(lm.jams)

slope.jams     <- coef(lm.jams)["log10_animal_size"]
slope.jams.lwr <- ci.jams["log10_animal_size", 1]
slope.jams.upr <- ci.jams["log10_animal_size", 2]

int.jams     <- coef(lm.jams)["(Intercept)"]
int.jams.lwr <- ci.jams["(Intercept)", 1]
int.jams.upr <- ci.jams["(Intercept)", 2]

# ── Summary table ────────────────────────────────────────────────────────────
results <- data.frame(
  Model     = c("This paper (Field)", "This paper (Lab)", "Jams et al."),
  Slope     = round(c(slope.us,     slope.us,     slope.jams),     2),
  Slope_lwr = round(c(slope.us.lwr, slope.us.lwr, slope.jams.lwr), 2),
  Slope_upr = round(c(slope.us.upr, slope.us.upr, slope.jams.upr), 2),
  Intercept     = round(c(int.field.us,     int.lab.us,     int.jams),     2),
  Intercept_lwr = round(c(int.field.us.lwr, int.lab.us.lwr, int.jams.lwr), 2),
  Intercept_upr = round(c(int.field.us.upr, int.lab.us.upr, int.jams.upr), 2)
)

print(results)

# suggestive of differences...

###########################################################
# Test for differences between studies (All data)
###########################################################

# lmer approach
lmer.all <- lmer(log10_largest_MP ~ log10_animal_size + study + setting + log10_animal_size:study + (1|species_binomial), data = df.all)
summary(lmer.all)
anova(lmer.all)
confint(lmer.all)

lmer.all.small <- lmer(log10_smallest_MP ~ log10_animal_size + study + setting + log10_animal_size:study + (1|species_binomial), data = df.all)
summary(lmer.all.small)
anova(lmer.all.small)
confint(lmer.all.small)

#species are largely unique to each dataset, 
#and the random effect variance converged to zero (singular fit); 
# let's use a simplfied  fixed-effects structure with study and setting as covariates.


lm.all <- lm(log10_largest_MP ~ log10_animal_size + study + setting + log10_animal_size:study, data = df.all)
summary(lm.all)
anova(lm.all)
confint(lm.all)

lm.all.small <- lm(log10_smallest_MP ~ log10_animal_size + study + setting + log10_animal_size:study, data = df.all)
summary(lm.all.small)
anova(lm.all.small)
confint(lm.all.small)

###########################################################
# Test for differences between studies (field data only)
###########################################################

# df.field <- subset(df.all, setting == "Field")
# 
# lm.field <- lm(log10_largest_MP ~ log10_animal_size * study, data = df.field)
# summary(lm.field)
# anova(lm.field)
# model_performance(lm.field)
# check_model(lm.field)
# performance::check_collinearity(lm.field)
# 
# # Slope comparison for smallest MP too
# lm.field.small <- lm(log10_smallest_MP ~ log10_animal_size * study, data = df.field)
# summary(lm.field.small)
# anova(lm.field.small)

#############################
# Overall scaling relationship (ALL data)
##############################

# Fixed effect for setting intercept correction - uses all data
lm.overall <- lm(log10_largest_MP ~ log10_animal_size + setting, data = df.all)
summary(lm.overall)
confint(lm.overall)

lm.overall.small <- lm(log10_smallest_MP ~ log10_animal_size + setting, data = df.all)
summary(lm.overall.small)
confint(lm.overall.small)

# #############################
# # Overall scaling relationship (field data)
# ##############################
# 
# lm.all <- lm(log10_largest_MP ~ log10_animal_size, data = df.field)
# summary(lm.all)
# 
# confint(lm.all)
# model_performance(lm.all)
# check_model(lm.all)
# 
# # Smallest MP - overall field scaling
# lm.all.small <- lm(log10_smallest_MP ~ log10_animal_size, data = df.field)
# summary(lm.all.small)
# confint(lm.all.small)

######################################################
# Extract fitted values from our model for combined plot
######################################################


df.us <- subset(df.all, study == "This paper")
x_vals <- seq(
  min(df.us$log10_animal_size, na.rm = TRUE),
  max(df.us$log10_animal_size, na.rm = TRUE),
  length.out = 100
)

newdata <- data.frame(log10_animal_size = x_vals, setting = "Field")
pred <- predict(lmer.us, newdata = newdata, re.form = NA, se.fit = TRUE)

manual_fit <- data.frame(
  log10_animal_size = x_vals,
  fit = pred$fit,
  lwr = pred$fit - 1.96 * pred$se.fit,
  upr = pred$fit + 1.96 * pred$se.fit,
  study = "Our_study_field")


###########################################################
# Figure 3: Combined our data + Jams - largest MP
###########################################################

# ── Study-specific x ranges ──────────────────────────────────────────────────
x_range_jams <- seq(min(df.jams$log10_animal_size, na.rm = TRUE),
                    max(df.jams$log10_animal_size, na.rm = TRUE),
                    length.out = 200)

x_range_us <- seq(min(df.us$log10_animal_size, na.rm = TRUE),
                  max(df.us$log10_animal_size, na.rm = TRUE),
                  length.out = 200)

x_range_overall <- seq(min(df.all$log10_animal_size, na.rm = TRUE),
                       4,
                       length.out = 200)

# ── Predictions using study-specific ranges ──────────────────────────────────
pred_jams_all <- predict(lm.all,
                         newdata = data.frame(log10_animal_size = x_range_jams,
                                              study   = "Jams et al.",
                                              setting = "Field"),
                         interval = "confidence") %>%
  as.data.frame() %>%
  mutate(log10_animal_size = x_range_jams)

pred_us_all <- predict(lm.all,
                       newdata = data.frame(log10_animal_size = x_range_us,
                                            study   = "This paper",
                                            setting = "Field"),
                       interval = "confidence") %>%
  as.data.frame() %>%
  mutate(log10_animal_size = x_range_us)

pred_overall <- predict(lm.overall,
                        newdata = data.frame(log10_animal_size = x_range_overall,
                                             setting = "Field"),
                        interval = "confidence") %>%
  as.data.frame() %>%
  mutate(log10_animal_size = x_range_overall)

# ── Figure 3 ─────────────────────────────────────────────────────────────────
fig3 <- ggplot() +
  
  # Lab points - visible but subdued
  geom_point(data = subset(df.all, setting == "Laboratory"),
             aes(x = log10_animal_size, y = log10_largest_MP, shape = study),
             fill = "grey70", colour = "grey50",
             size = 2, alpha = 0.5, stroke = 0.3) +
  
  # Field points - full opacity
  geom_point(data = subset(df.all, setting == "Field"),
             aes(x = log10_animal_size, y = log10_largest_MP,
                 shape = study, fill = study),
             colour = "black", size = 2.5, alpha = 0.85, stroke = 0.5) +
  
  # Jams et al. ribbon + line - spans Jams data range only
  geom_ribbon(data = pred_jams_all,
              aes(x = log10_animal_size, ymin = lwr, ymax = upr),
              fill = "grey50", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = pred_jams_all,
            aes(x = log10_animal_size, y = fit),
            colour = "grey40", linewidth = 0.9, linetype = "dashed",
            inherit.aes = FALSE) +
  
  # This paper ribbon + line - spans this paper data range only
  geom_ribbon(data = pred_us_all,
              aes(x = log10_animal_size, ymin = lwr, ymax = upr),
              fill = "grey50", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = pred_us_all,
            aes(x = log10_animal_size, y = fit),
            colour = "grey40", linewidth = 0.9, linetype = "dotted",
            inherit.aes = FALSE) +
  
  # Overall scaling - solid black, full x range
  geom_ribbon(data = pred_overall,
              aes(x = log10_animal_size, ymin = lwr, ymax = upr),
              fill = "black", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = pred_overall,
            aes(x = log10_animal_size, y = fit),
            colour = "black", linewidth = 1.3, linetype = "solid",
            inherit.aes = FALSE) +
  
  coord_cartesian(xlim = c(-2.5, 4), ylim = c(-4, 3)) +
  labs(
    x = expression(log[10]~"Animal body length [mm]"),
    y = expression(log[10]~"Longest microplastic ingested [mm]")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    axis.title.y     = element_text(margin = margin(r = 10)),
    axis.title.x     = element_text(margin = margin(t = 10))
  ) +
  scale_fill_manual(name = "Study", values  = c("This paper" = "black", "Jams et al." = "white")) +
  scale_shape_manual(name = "Study", values = c("This paper" = 21, "Jams et al." = 23)) + theme_classic(base_size = 11)

# ── Extract equation components from lm.overall ──────────────────────────────
slope_ov     <- round(coef(lm.overall)["log10_animal_size"], 2)
intercept_ov <- round(coef(lm.overall)["(Intercept)"], 2)
r2_ov        <- round(summary(lm.overall)$r.squared, 2)

eq_label <- bquote(log[10](MP) == .(intercept_ov) + .(slope_ov)~log[10](Size)~~(R^2 == .(r2_ov)))

eq_string <- paste0("log[10](MP)==", intercept_ov,
                    "+", slope_ov, "*symbol('\\264')~log[10](Body~length)",
                    "~~~R^2==", r2_ov)

fig3 <- fig3 +
  annotate("text",
           x = -2.3, y = 2.8,
           label = eq_string,
           hjust = 0, size = 3.2, parse = TRUE)
fig3

#ggsave(plot = fig3, device = pdf,
       #filename = 'Fig3_combined_field_largest_MP.pdf',
       #width = 200, height = 100, dpi = 300, units = "mm")


###########################################################
# Figure 4: Combined our data + Jams - smallest MP
###########################################################
lm.all.small <- lm(log10_smallest_MP ~ log10_animal_size * study + study * setting, data = df.all)
lm.overall.small <- lm(log10_smallest_MP ~ log10_animal_size + setting, data = df.all)

# ── Predictions using study-specific ranges ──────────────────────────────────
pred_jams_all.small <- predict(lm.all.small,
                         newdata = data.frame(log10_animal_size = x_range_jams,
                                              study   = "Jams et al.",
                                              setting = "Field"),
                         interval = "confidence") %>%
  as.data.frame() %>%
  mutate(log10_animal_size = x_range_jams)

pred_us_all.small <- predict(lm.all.small,
                       newdata = data.frame(log10_animal_size = x_range_us,
                                            study   = "This paper",
                                            setting = "Field"),
                       interval = "confidence") %>%
  as.data.frame() %>%
  mutate(log10_animal_size = x_range_us)

pred_overall.small <- predict(lm.overall.small,
                        newdata = data.frame(log10_animal_size = x_range_overall,
                                             setting = "Field"),
                        interval = "confidence") %>%
  as.data.frame() %>%
  mutate(log10_animal_size = x_range_overall)

# ── Figure 4 ─────────────────────────────────────────────────────────────────
fig4 <- ggplot() +
  
  # Lab points - visible but subdued
  geom_point(data = subset(df.all, setting == "Laboratory"),
             aes(x = log10_animal_size, y = log10_smallest_MP, shape = study),
             fill = "grey70", colour = "grey50",
             size = 2, alpha = 0.5, stroke = 0.3) +
  
  # Field points - full opacity
  geom_point(data = subset(df.all, setting == "Field"),
             aes(x = log10_animal_size, y = log10_smallest_MP,
                 shape = study, fill = study),
             colour = "black", size = 2.5, alpha = 0.85, stroke = 0.5) +
  
  # Jams et al. ribbon + line - spans Jams data range only
  geom_ribbon(data = pred_jams_all.small,
              aes(x = log10_animal_size, ymin = lwr, ymax = upr),
              fill = "grey50", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = pred_jams_all.small,
            aes(x = log10_animal_size, y = fit),
            colour = "grey40", linewidth = 0.9, linetype = "dashed",
            inherit.aes = FALSE) +
  
  # This paper ribbon + line - spans this paper data range only
  geom_ribbon(data = pred_us_all.small,
              aes(x = log10_animal_size, ymin = lwr, ymax = upr),
              fill = "grey50", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = pred_us_all.small,
            aes(x = log10_animal_size, y = fit),
            colour = "grey40", linewidth = 0.9, linetype = "dotted",
            inherit.aes = FALSE) +
  
  # Overall scaling - solid black, full x range
  geom_ribbon(data = pred_overall.small,
              aes(x = log10_animal_size, ymin = lwr, ymax = upr),
              fill = "black", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = pred_overall.small,
            aes(x = log10_animal_size, y = fit),
            colour = "black", linewidth = 1.3, linetype = "solid",
            inherit.aes = FALSE) +
  
  coord_cartesian(xlim = c(-2.5, 4), ylim = c(-4, 3)) +
  labs(
    x = expression(log[10]~"Animal body length [mm]"),
    y = expression(log[10]~"Shortest microplastic ingested [mm]")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    axis.title.y     = element_text(margin = margin(r = 10)),
    axis.title.x     = element_text(margin = margin(t = 10))
  ) +
  scale_fill_manual(name = "Study", values  = c("This paper" = "black", "Jams et al." = "white")) +
  scale_shape_manual(name = "Study", values = c("This paper" = 21, "Jams et al." = 23)) + theme_classic(base_size = 11)

# ── Extract equation components from lm.overall ──────────────────────────────
slope_ov     <- round(coef(lm.overall.small)["log10_animal_size"], 2)
intercept_ov <- round(coef(lm.overall.small)["(Intercept)"], 2)
r2_ov        <- round(summary(lm.overall.small)$r.squared, 2)

eq_label <- bquote(log[10](MP) == .(intercept_ov) + .(slope_ov)~log[10](Size)~~(R^2 == .(r2_ov)))

eq_string <- paste0("log[10](MP)==", intercept_ov,
                    "+", slope_ov, "*symbol('\\264')~log[10](Body~length)",
                    "~~~R^2==", r2_ov)


fig4 <- fig4 +
  annotate("text",
           x = -2.3, y = 2.8,
           label = eq_string,
           hjust = 0, size = 3.2, parse = TRUE)
fig4

#ggsave(plot = fig4, device = pdf,
       #filename = 'Fig4_combined_field_smallest_MP.pdf',
       #width = 200, height = 100, dpi = 300, units = "mm")


## ============================================================================
## Bayesian model
## ============================================================================

source("Bayesian_MA_mixed.R")

