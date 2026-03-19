library(mice)
library(VIM)
library(naniar)
library(ggplot2)
library(dplyr)
library(broom)

#-------------------------
# 1. Read data
#-------------------------
dat <- dget("NHANESvIDA.rda")
# Remove participant ID
dat$id <- NULL

str(dat)
summary(dat)

# Create physical activity indicator
# If average number of MVPA bouts per day > 0, classify as "yes"
dat$mvpa_any <- factor(ifelse(dat$num_mvpa_bouts > 0, "yes", "no"))

# Keep variables relevant for Task 2
dat2 <- dat %>%
  select(apoliprot, cholesterol,
         gender, age, ethnic, educ,
         waist, cpm, num_mvpa_bouts, mvpa_any)

#-------------------------
# 2. Missing-data exploration
#-------------------------

# number and proportion of missing values
miss_n <- sapply(dat2, function(x) sum(is.na(x)))
miss_prop <- sapply(dat2, function(x) mean(is.na(x)))
# summary table of missing values
data.frame(
  variable = names(miss_n),
  n_missing = miss_n,
  prop_missing = round(miss_prop, 3)
)

# Visualise missingness
aggr(dat2, prop = FALSE, numbers = TRUE)
vis_miss(dat2)
md.pattern(dat2)

# missingness indicators for informal exploration
# R_apoliprot/R_waist = 1 means apoliprot/waist is missing
dat2$R_apoliprot <- as.integer(is.na(dat2$apoliprot))
dat2$R_waist <- as.integer(is.na(dat2$waist))

# Example checks for missingness mechanism.
# If missingness differs by observed covariates, MCAR becomes less plausible.
table(dat2$gender, dat2$R_apoliprot)
prop.table(table(dat2$gender, dat2$R_apoliprot), 1)
table(dat2$ethnic, dat2$R_apoliprot)
prop.table(table(dat2$ethnic, dat2$R_apoliprot), 1)
table(dat2$educ, dat2$R_waist)
prop.table(table(dat2$educ, dat2$R_waist), 1)

# Compare distributions of observed variables across missingness groups.
tapply(dat2$age, dat2$R_apoliprot, summary, na.rm = TRUE)
tapply(dat2$cpm, dat2$R_apoliprot, summary, na.rm = TRUE)

#-------------------------
# 3. Descriptive summaries
#-------------------------

# Numerical summaries for key continuous variables
summary(dat2$apoliprot)
summary(dat2$waist)
summary(dat2$cpm)

# Frequency tables for categorical predictors
table(dat2$gender)
table(dat2$ethnic)
table(dat2$educ)
table(dat2$mvpa_any)

# Histograms
ggplot(dat2, aes(x = apoliprot)) +
  geom_histogram(bins = 30) +
  theme_minimal()

ggplot(dat2, aes(x = waist)) +
  geom_histogram(bins = 30) +
  theme_minimal()

ggplot(dat2, aes(x = cpm)) +
  geom_histogram(bins = 30) +
  theme_minimal()

# Boxplots by categorical predictors
ggplot(dat2, aes(x = gender, y = apoliprot)) +
  geom_boxplot() +
  theme_minimal()

ggplot(dat2, aes(x = ethnic, y = apoliprot)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot(dat2, aes(x = educ, y = apoliprot)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Scatterplots
ggplot(dat2, aes(x = waist, y = apoliprot)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

ggplot(dat2, aes(x = cpm, y = apoliprot)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

#-------------------------
# 4. Complete-case analysis
#-------------------------

# Response must not be imputed, so restrict to observed response
dat_cc_base <- dat2 %>% filter(!is.na(apoliprot))

# For complete-case analysis, remove any remaining rows with missing predictors.
# complete-case dataset
dat_cc <- dat_cc_base %>%
  select(apoliprot, gender, age, ethnic, educ, waist, cpm, mvpa_any) %>%
  na.omit()

# Fit the complete-case linear regression model.
fit_cc <- lm(apoliprot ~ gender + age + ethnic + educ + waist + cpm + mvpa_any,
             data = dat_cc)

summary(fit_cc)
tidy(fit_cc, conf.int = TRUE)

# Diagnostic plots
par(mfrow = c(2,2))
plot(fit_cc)

#-------------------------
# 5. Multiple imputation (predictors only)
#-------------------------

# Keep only cases with observed response
dat_mi <- dat2 %>%
  filter(!is.na(apoliprot)) %>%
  select(apoliprot, gender, age, ethnic, educ, waist, cpm, mvpa_any)

# Set up mice
ini <- mice(dat_mi, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# Do not impute response
meth["apoliprot"] <- ""

# Only waist is expected to need imputation here; PMM is good for continuous data
meth["waist"] <- "pmm"

# Fully observed variables do not need a method
meth["gender"] <- ""
meth["age"] <- ""
meth["ethnic"] <- ""
meth["educ"] <- ""
meth["cpm"] <- ""
meth["mvpa_any"] <- ""

# Run multiple imputation
set.seed(123)
imp <- mice(dat_mi, m = 20, method = meth, predictorMatrix = pred,
            maxit = 20, printFlag = FALSE)

# Check convergence
plot(imp)

# Fit the same substantive model within each imputed dataset.
fit_mi <- with(imp, lm(apoliprot ~ gender + age + ethnic + educ + waist + cpm + mvpa_any))
# Pool the estimates across all imputed datasets
pool_fit <- pool(fit_mi)

summary(pool_fit, conf.int = TRUE)

#-------------------------
# 6. Compare CC and MI
#-------------------------

cc_tab <- tidy(fit_cc, conf.int = TRUE) %>%
  mutate(method = "Complete case")

mi_tab <- summary(pool_fit, conf.int = TRUE) %>%
  as.data.frame() %>%
  select(term, estimate, std.error, statistic, df, p.value, `2.5 %`, `97.5 %`) %>%
  rename(conf.low = `2.5 %`,
         conf.high = `97.5 %`) %>%
  mutate(method = "Multiple imputation")

compare_tab <- bind_rows(cc_tab, mi_tab)
compare_tab2 <- compare_tab %>%
  select(method, term, estimate, std.error, conf.low, conf.high, p.value) %>%
  arrange(term, method)
print(compare_tab2, n = 28, width = Inf)
# MI yielded results that were very similar to CC.

# Missingness table
miss_tab <- data.frame(
  variable = names(miss_n),
  n_missing = miss_n,
  prop_missing = round(miss_prop, 3)
)
miss_tab

# Sample sizes
n_total <- nrow(dat2)
n_response_obs <- nrow(dat_cc_base)
n_cc <- nrow(dat_cc)

cat("Total sample size:", n_total, "\n")
cat("Observed response sample size:", n_response_obs, "\n")
cat("Complete-case sample size:", n_cc, "\n")

# Check missingness in analysis dataset
colSums(is.na(dat_mi))

# Key comparison table
key_compare_tab <- compare_tab2 %>%
  filter(term %in% c("genderfemale", "age", "educCollege/AA", "waist", "cpm", "mvpa_anyyes"))

key_compare_tab


