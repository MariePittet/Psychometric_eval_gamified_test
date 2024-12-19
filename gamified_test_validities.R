# Set-up ------------------------------------------------------------------

# installing packages
require(tidyverse)
require(matrixStats)
require(gridExtra)
require(ggpubr)
require(ggdist)
require(rstatix)
require(correlation)
require(Hmisc)
require(corrplot)
theme_set(theme_classic())
require(wesanderson)
require(lavaan)
require(tidySEM)
require(finalfit)
require(SimplyAgree)
require(lmerTest)
require(nlme)


# Loading and curating data -----------------------------------------------

# Load the dataset
data <- read.csv2("ACE-MEM_data.csv", header = TRUE, sep = ";")

# Convert all columns to numeric
data[] <- lapply(data, as.numeric)

# Recategorize the 'Centre' column with descriptive factor levels
data$Centre <- factor(data$Centre, levels = c("1", "2", "3"), labels = c("Lausanne", "Bern", "Nice"))

# Rename the gender column for clarity
names(data)[names(data) == "Gender..0...F...1...M."] <- "Gender"

# Convert 'Gender' and 'Education' to factors
# Gender: 0 = Female (F), 1 = Male (M)
data$Gender <- factor(data$Gender)
data$Education <- factor(data$Education)

# Keep only the first 101 rows of the dataset
data <- data[1:101, ]

# Filter rows where 'Keep' column equals 1
data <- data[data$Keep == 1, ]

#Renaming variables
data <- data %>% 
  rename(Center = Centre,
         MoCA_total = MOCA.total,
         Corsi_forward_span = corsi_forward_span,
         Corsi_backward_span = corsi_backward_span,
         ACE_forward_span = Forwardspan_countcorrect_max,
         ACE_backward_span = Backwardspan_countcorrect_max)


# Demographics (per center of data acquisition) ---------------------------

# Define demographic variables of interest
explanatory <- c("Age", "Gender", "Education", "MoCA_total", 
                 "Corsi_forward_span", "Corsi_backward_span", 
                 "ACE_forward_span", "ACE_backward_span")

# Generate a summary table with p-values, grouping by 'Center'
data %>%
  summary_factorlist("Center", explanatory, p = TRUE, na_include = FALSE)


# Inspect for extreme values (+- 3SD) -------------------------------------

identify_outliers(data = data, "Corsi_forward_span")
identify_outliers(data = data, "Corsi_backward_span")
identify_outliers(data = data, "ACE_forward_span")
identify_outliers(data = data, "ACE_backward_span")


# Convergent validity: Agreement analyses --------------------------------

# Calculate limits of agreement between ACE and Corsi spans
agreement_forward <- agreement_limit(
  data = data,
  x = "ACE_forward_span",
  y = "Corsi_forward_span"
)

agreement_backward <- agreement_limit(
  data = data,
  x = "ACE_backward_span",
  y = "Corsi_backward_span"
)

# Calculate mean and difference scores for forward and backward spans
data$mean_forward_spans <- rowMeans(data[, c("ACE_forward_span", "Corsi_forward_span")], na.rm = TRUE)
data$diff_forward_spans <- data$ACE_forward_span - data$Corsi_forward_span

data$mean_backward_spans <- rowMeans(data[, c("ACE_backward_span", "Corsi_backward_span")], na.rm = TRUE)
data$diff_backward_spans <- data$ACE_backward_span - data$Corsi_backward_span

# Plot settings for two side-by-side plots
par(mfrow = c(1, 2))

# Forward spans plot
plot(jitter(data$mean_forward_spans, 2), data$diff_forward_spans,
     xlab = expression(mu * "(gamified + standard span scores)"),
     ylab = expression(Delta * "(gamified - standard span scores)"),
     cex.lab = 0.9,
     main = substitute(paste(italic("A. Forward spans"))),
     cex.main = 1,
     col = rgb(0, 0, 0, 0.4),  # semi-transparent black
     pch = 16,
     cex = 1.2)
abline(h = c(-2.00, 2.99), col = "darkorange1", lty = 2, lwd = 2)  # Limits of Agreement
abline(h = 0.4947, col = "tomato2", lty = 1, lwd = 2)              # Bias line
legend("topright", legend = c("LoA", "Bias"),
       col = c("darkorange1", "tomato2"), 
       lty = c(2, 1),
       text.font = 4,
       cex = 0.8, 
       box.lty = 0,
       bty = "n")

# Backward spans plot
plot(jitter(data$mean_backward_spans, 2), data$diff_backward_spans,
     xlab = expression(mu * "(gamified + standard span scores)"),
     ylab = expression(Delta * "(gamified - standard span scores)"),
     cex.lab = 0.9,
     main = substitute(paste(italic("B. Backward spans"))),
     cex.main = 1,
     col = rgb(0, 0, 0, 0.4),  # semi-transparent black
     pch = 16,
     cex = 1.2)
abline(h = c(-2.00, 2.99), col = "skyblue3", lty = 2, lwd = 2)  # Limits of Agreement
abline(h = 0.5638, col = "blue4", lty = 1, lwd = 2)             # Bias line
legend("topright", legend = c("LoA", "Bias"),
       col = c("skyblue3", "blue4"), 
       lty = c(2, 1),
       text.font = 4,
       cex = 0.8, 
       box.lty = 0,
       bty = "n")

# Add a main title across both plots
mtext("Agreement between gamified and standard span scores", cex = 1.2, font = 2, side = 3, line = -1.2, outer = TRUE)


# Convergent analyses: correlations with the gold standard (Corsi test) -----------------------

# Computing scores centered around the mean of each data collection center (group-mean-centered scores) to avoid
# spurious correlations.

# Split the data by center
df_chuv <- subset(data, Center == "Lausanne")
df_bern <- subset(data, Center == "Bern")
df_nice <- subset(data, Center == "Nice")

# Select columns to adjust and subtract center-specific means
adjust_columns <- c(5, 8:ncol(data))

df_chuv_adj <- df_chuv[, adjust_columns] - colMeans(df_chuv[, adjust_columns], na.rm = TRUE)[col(df_chuv[, adjust_columns])]
df_bern_adj <- df_bern[, adjust_columns] - colMeans(df_bern[, adjust_columns], na.rm = TRUE)[col(df_bern[, adjust_columns])]
df_nice_adj <- df_nice[, adjust_columns] - colMeans(df_nice[, adjust_columns], na.rm = TRUE)[col(df_nice[, adjust_columns])]

# Combine adjusted data from all centers
df_adj <- rbind(df_chuv_adj, df_bern_adj, df_nice_adj)

# Reconstruct final data frame with additional demographic columns
df_adj <- data.frame(Patient.ID = data$Patient.ID, 
                     Center = data$Center, 
                     Age = data$Age, 
                     Gender = data$Gender, 
                     Education = data$Education, 
                     df_adj)

# Now on to the correlations using adjusted data (df_adj)

# Plot for ACE backward span vs. Corsi backward span
p1 <- ggplot(df_adj, aes(x = ACE_backward_span, y = Corsi_backward_span)) +
  geom_point(alpha = 0.4, size = 2, position = position_jitter(w = 0.1, h = 0.1), colour = "blue2") +
  geom_smooth(method = "lm", colour = "black") +
  ylim(-2.5, 4) +
  stat_cor(method = "spearman", cor.coef.name = "rho") +
  labs(
    x = "Gamified backward span",
    y = "Standard backward span",
    title = "B. Backward spans"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))

# Plot for ACE forward span vs. Corsi forward span
p2 <- ggplot(df_adj, aes(x = ACE_forward_span, y = Corsi_forward_span)) +
  geom_point(alpha = 0.5, size = 2, position = position_jitter(w = 0.1, h = 0.1), colour = "darkorange1") +
  geom_smooth(method = "lm", colour = "black") +
  ylim(-2.5, 4) +
  stat_cor(method = "spearman", cor.coef.name = "rho", digits = 2) +
  labs(
    x = "Gamified forward span",
    y = "Standard forward span",
    title = "A. Forward spans"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))

# Arrange plots side by side
convergent_val <- grid.arrange(p2, p1, ncol = 2)

# Add overall title
annotate_figure(
  convergent_val,
  top = text_grob("Convergent validity", color = "black", face = "bold", size = 14)
)


# Effects of age and modality of testing on span scores -------------------

# It can be hypothesized that older patients may struggle more with digital testing than with standard testing.
# Thus, we have to model the interaction Age*Modality (gamified vs standard) on span scores for the forward
# and the backward tests. 

# Data Wrangling: Transform data to long format for forward and backward spans
df_reg_forward <- pivot_longer(
  data = data,
  cols = c(Corsi_forward_span, ACE_forward_span),
  names_to = "Modality",
  values_to = "Scores"
)

df_reg_backward <- pivot_longer(
  data = data,
  cols = c(Corsi_backward_span, ACE_backward_span),
  names_to = "Modality",
  values_to = "Scores"
)

# Linear Mixed Models
# Forward span model: Interaction of Modality and Age with random effect of Center (to control for center differences)
model_forward <- lmer(Scores ~ Modality * Age + (1 | Center), data = df_reg_forward, na.action = na.omit, REML = FALSE)
null_forward <- lmer(Scores ~ Modality + Age + (1 | Center), data = df_reg_forward, na.action = na.omit, REML = FALSE)
anova(model_forward, null_forward)
car::Anova(model_forward, type = 3)
 
# Backward span model: Interaction of Modality and Age with random effect of Center
model_backward <- lmer(Scores ~ Modality * Age + (1 | Center), data = df_reg_backward, na.action = na.omit, REML = FALSE)
null_backward <- lmer(Scores ~ Modality + Age + (1 | Center), data = df_reg_backward, na.action = na.omit, REML = FALSE)
anova(model_backward, null_backward)
car::Anova(model_backward, type = 3)
summary(model_backward)

# Adjust Modality factor labels for plotting
df_reg_forward$Modality <- factor(df_reg_forward$Modality, levels = c("Corsi_forward_span", "ACE_forward_span"), labels = c("Standard test", "Gamified test"))
df_reg_backward$Modality <- factor(df_reg_backward$Modality, levels = c("Corsi_backward_span", "ACE_backward_span"), labels = c("Standard test", "Gamified test"))

# Plotting
# Forward span plot
p9 <- ggplot(df_reg_forward, aes(x = Age, y = Scores, colour = Modality)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("darkorange1", "tomato2"), labels = c("Standard", "Gamified")) +
  ylim(2, 9) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  labs(title = "A. Forward spans") +
  guides(colour = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -0.5, face = "italic", size = 12)) +
  annotate("text", x = 70, y = 8.5, label = expression("Age*Modality: " * italic(p) * " = 0.303"), size = 4, color = "grey60")

# Backward span plot
p10 <- ggplot(df_reg_backward, aes(x = Age, y = Scores, colour = Modality)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("skyblue3", "blue4"), labels = c("Standard", "Gamified")) +
  ylim(2, 9) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1.2) +
  labs(title = "B. Backward spans") +
  guides(colour = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -0.5, face = "italic", size = 12)) +
  annotate("text", x = 70, y = 8.5, label = expression("Age*Modality: " * italic(p) * " = 0.044*"), size = 4, color = "grey60")

# Combine plots into a single figure
p_combined <- ggarrange(p9, p10, ncol = 2, nrow = 1, common.legend = FALSE, legend = "right")

# Annotate the combined plot with a title
annotate_figure(
  p_combined,
  top = text_grob("Effect of Age and Testing Modality on Scores", color = "black", face = "bold", size = 14)
)


# Discriminant validity analyses ------------------------------------------

# Here we examinate if ACE Gem Chaser can discriminate patients with mild non clininical dementia (mNCD) and healthy
# controls of the same age range (ACE normative data). 

# Merging two age categories from the norm to have the same range as in our study
mu_forward<- (43*6.81 + 26*6.69)/69
sd_forward<- (43*0.82 + 26*0.79)/69
mu_backward<- (34*6.41 + 26*6.35)/60
sd_backward<- (34*0.86 + 26*0.89)/60

# One-sample t-test to compare mNCD patients and norms on Gem Chaser forward span
data %>%
  t_test(ACE_forward_span ~ 1, mu = mu_forward) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
data %>%
  cohens_d(ACE_forward_span ~ 1, mu = mu_forward)

# One-sample t-test to compare mNCD patients and norms on Gem Chaser backward span
data %>%
  t_test(ACE_backward_span ~ 1, mu = mu_backward) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
data %>% 
  cohens_d(ACE_backward_span ~ 1, mu = mu_backward)

# Data wrangling for side by side boxplots
patients_data <- data.frame(
  Score = c(data$ACE_forward_span, data$ACE_backward_span),  # Forward and Backward data for patients
  Modality = rep(c("Forward", "Backward"), each = 49),
  Group = rep("Patients", 98)
)

# Normative data (mean and SD for Norms)
norms_data <- data.frame(
  mean = c(mu_forward, mu_backward),  # Forward and Backward
  sd = c(sd_forward, sd_backward),  # Forward and Backward
  Group = rep("Norms", 2),
  Modality = c("Forward", "Backward")
)

# Simulate the 'Norms' data using the mean and SD for visualization (assuming normal distribution)
norms_simulated <- norms_data %>%
  rowwise() %>%
  mutate(
    Score = list(rnorm(198, mean, sd))  # Simulating 100 samples per modality
  ) %>%
  unnest(Score)

# Combine both Norms and Patients data
combined_data <- bind_rows(patients_data, norms_simulated)

# Plot: Boxplots for both Forward and Backward spans
ggplot(combined_data, aes(x = fct_rev(Modality), y = Score, fill = Group)) +
  # Boxplot for both groups (Norms and Patients)
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.75), width = 0.3) +
  scale_fill_manual(values = c("skyblue3", "tomato2")) +  # Norms (lightblue), Patients (tomato2)
  scale_color_manual(values = c("skyblue3", "tomato2")) +  
  
  labs(
    x = "Span Type",
    y = "Score",
    title = "Comparison of Forward and Backward Span Scores: Patients vs Norms"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(face = "italic", hjust = 1)
  )