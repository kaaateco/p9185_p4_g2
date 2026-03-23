
# ──────────────────────────────────────────────────────────────────────────
# Author: Yunan Zhao
# SPRIS Project 4
# Date: 03/23/2026
# ──────────────────────────────────────────────────────────────────────────

## b.	After the company had your statistical analysis plan in place, the trial
## began, and data were collected. In addition to the primary efficacy outcome,
## the company also monitored the serious adverse events (SAE) closely by
## contacting study participants on a monthly basis for 3 months. SAE was
## defined as any untoward medical occurrence that was life-threatening or
## required inpatient hospitalization. Please use the provided data (Q2b_BL.xlsx
## and Q2b.xlsx) to determine whether the vaccine group has greater odds of
## having SAE at any of the three assessment time points.  

# -----------------------------
# Load Packages
# -----------------------------
library(readxl)
library(lme4)
library(gtsummary)
library(ggplot2)

## ---- Load the data ----------------------------------------------------------
Q2b <- readxl::read_excel("./Q2b.xlsx")
baseline.dat <- readxl::read_excel("./Q2b_BL.xlsx")

long.dat <- left_join(Q2b, baseline.dat, by = "ID") %>% 
  mutate(OBS = if_else(is.na(SAE), 0, 1)) %>% 
  mutate(TIME = factor(TIME),
         GROUP = factor(GROUP, labels = c("Control", "Vaccine"),
                        levels = c(0,1)),
         SITE = factor(SITE),
         SEX = factor(SEX, labels = c("Female", "Male"), levels = c(0,1)),
         SAE = factor(SAE, labels = c("No", "Yes"), levels = c(0,1)))

## ---- eda plot ---------------------------------------------------------------

## Number of subjects experienced serious adverse events (SAE) at each time
## point, stratified by treatment group
long.dat %>% 
  group_by(TIME, GROUP) %>% 
  mutate(SAE_TOTAL = sum(SAE == "Yes", na.rm = T),
         MISSING_TOTAL = sum(OBS == 0)) %>% 
  select(TIME, GROUP, SAE_TOTAL, MISSING_TOTAL) %>% 
  unique() %>% 
  ggplot(aes(x=TIME, y=SAE_TOTAL, fill=GROUP)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Number of Subjects Experienced SAE by Group") +
  ylab("Number of SAE") +
  xlab("Time (Months)") + 
  labs(fill = "Treatment Group") + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = c("Control" = "lightgrey",
                                "Vaccine" = "skyblue")) +
  geom_text(
    aes(label = SAE_TOTAL),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 4) +
  ylim(0, 40) 

long.dat %>% 
  group_by(TIME, GROUP) %>% 
  summarise(
    n_total = sum(!is.na(SAE)),
    n_sae = sum(SAE == "Yes", na.rm = TRUE),
    percent = round(100 * n_sae / n_total, 2),
    .groups = "drop"
  )

## Number of missing observations at each time point, stratified by treatment
## group
long.dat %>% 
  group_by(TIME, GROUP) %>% 
  mutate(SAE_TOTAL = sum(SAE == "Yes", na.rm = T),
         MISSING_TOTAL = sum(OBS == 0)) %>% 
  select(TIME, GROUP, SAE_TOTAL, MISSING_TOTAL) %>% 
  unique() %>% 
  ggplot(aes(x=TIME, y=MISSING_TOTAL, fill=GROUP)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Number of Missing Observations by Treatment Group") +
  ylab("Number of Missing Observations") +
  xlab("Time (Months)") + 
  labs(fill = "Treatment Group") + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = c("Control" = "lightgrey",
                               "Vaccine" = "skyblue")) +
  geom_text(
    aes(label = MISSING_TOTAL),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 4) +
  ylim(0, 4200) 

long.dat %>% 
  group_by(TIME, GROUP) %>% 
  summarise(
    n_total = sum(OBS),
    n_missing = sum(OBS == 0, na.rm = TRUE),
    percent = round(100 * n_missing / n_total, 2),
    .groups = "drop"
  )


## ---- glmm model -------------------------------------------------------------
glmm.fit0 <- glmer(SAE ~ TIME * GROUP + (1|ID) , data = long.dat,
                   family = binomial, nAGQ = 9)
summary(glmm.fit0)

glmm.fit1 <- glmer(SAE ~ TIME * GROUP + SEX + AGE + (1|ID) , data = long.dat,
                   family = binomial, nAGQ = 9)
# using nAGQ=0 and nAGQ=1 provides stable outcome
summary(glmm.fit.1)
anova(glmm.fit0, glmm.fit1)
isSingular(glmm.fit1, tol=1e-4) #TRUE
VarCorr(glmm.fit1)

glmm.tbl <- 
  tbl_regression(glmm.fit1, exponentiate = T, 
                 label = list(
                   TIME ~ "Time point",
                   GROUP ~ "Treatment group",
                   SEX ~ "Gender",
                   AGE ~ "Age"
                   ))

## Not using glmm because we have very little within-subject correlation
glm.fit1 <- glm(SAE ~ TIME * GROUP + SEX + AGE, data = long.dat,
                family = binomial)
tbl_regression(glm.fit1, exponentiate = T, 
               label = list(
                 TIME ~ "Time point",
                 GROUP ~ "Treatment group",
                 SEX ~ "Gender",
                 AGE ~ "Age"
                 ))



