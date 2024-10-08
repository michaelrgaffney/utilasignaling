# loading libraries -------------------------------------------------------

# note: list is cleaned up but some may no longer used
library(utiladata2023) # data package
library(tidyverse)
library(labelled)
# library(lme4)
library(effects)
library(visreg)
library(emmeans)
library(car)
library(haven)
library(MASS)
library(glmmTMB)
library(sandwich)
library(lmtest)
library(ggcorrplot)
library(marginaleffects)
library(ordinal) # clmm # polr
library(brms)
library(hagenutils)
library(localgrowth)
library(knitr)
library(kableExtra)
library(modelsummary)
library(ggalluvial)
library(patchwork)
library(mgcv)
library(gratia)
library(pvclust)

source("recode.R")
source("dictionaries.R")
source("dataprep.R")
source("regularize.R") # typo in name

# raw data plots  --------------------------------------------------

## main paper plots --------------------------------------------------------

#
signalfreqdf <- d2 |>
  dplyr::select(
    householdID,
    childHHid,
    uniqueID,
    SadFreqOF,
    CryFreqOF,
    TantrumFreqOF,
    SignalCost
  ) |>
  rename(
    Sadness = SadFreqOF,
    Crying = CryFreqOF,
    Tantrums = TantrumFreqOF
  )

freq_short <- c(
  "Never" = "0",
  "Once a month or less" = "≤1",
  "More than once a month but less than once a week" = "2-3",
  "More than once a week but not daily" = "4-29",
  "Daily" = "30",
  "Multiple times per day" = ">30"
)

freq_short2 <- c(
  "0" = "Low",
  "≤1" = "Low",
  "2-3" = "Medium",
  "4-29" = "Medium",
  "30" = "High",
  ">30" = "High"
)

sfdf <-
  signalfreqdf |>
  mutate(
    across(Sadness:Tantrums, \(x) factor(c(freq_short[x]), levels = c(freq_short)))
  ) |>
  na.omit()

signalheatmap <- function(d, signal1, signal2){
  out <- table(d[[signal1]], d[[signal2]])
  print(out)
  out2 <- matrix(out, ncol=ncol(out), dimnames = dimnames(out))
  hagenheat(out2, seriation_method = 'Identity', viridis_option = 'C') +
    ggtitle("Signal frequency (times per month)") +
    ylab(signal1) +
    xlab(signal2) +
    theme_minimal(20) +
    theme(axis.title.y = element_text(angle = 0))
}

signalheatmap(sfdf, 'Crying', 'Sadness')
signalheatmap(sfdf, 'Tantrums', 'Sadness')
signalheatmap(sfdf, 'Tantrums', 'Crying')

sfdfsum <-
  sfdf |>
  summarise(
    Freq = n(),
    .by = c(Sadness, Crying, Tantrums)
  ) |>
  arrange(desc(Freq))

sfdfsum2 <-
  sfdf |>
  mutate(
    across(Sadness:Tantrums, \(x) factor(c(freq_short2[x]), levels = rev(unique(freq_short2))))
  ) |>
  summarise(
    Freq = n(),
    MeanSignalCost = mean(SignalCost, na.rm = T),
    .by = c(Sadness, Crying, Tantrums)
  ) |>
  mutate(
    CostCategory = cut(
      MeanSignalCost,
      quantile(MeanSignalCost)[-3], # Bottom 1/4, Middle 1/2, Top 1/4
      labels = c("Low", "Medium", "High"),
      include.lowest = T
    )
  ) |>
  arrange(desc(Freq)) |>
  dplyr::filter(Freq > 2)

# sfdfsum2 |> arrange(desc(CostCategory), desc(Freq))

signal_alluvial_plot <- ggplot(sfdfsum2, aes(axis1 = Sadness, axis2 = Crying, axis3 = Tantrums, y = Freq)) +
  geom_alluvium(aes(fill = Sadness), color = "black", show.legend = F) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Sadness", "Crying", "Tantrums"), expand = c(.2, .05)) +
  ylab("Number of\nchildren") +
  theme_minimal(20) +
  theme(axis.title.y = element_text(angle = 0, hjust = 1)) +
  scale_fill_viridis_d()

signal_alluvial_plot

sfdfsum3 <- sfdfsum2 |>
  unite(SignalCombo, Sadness, Crying, Tantrums, remove = FALSE)

signal_alluvial_plot2 <- ggplot(sfdfsum3, aes(axis1 = Sadness, axis2 = Crying, axis3 = Tantrums, y = Freq)) +
  geom_alluvium(aes(fill = as.factor(SignalCombo)), color = "black", show.legend = F) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Sadness", "Crying", "Tantrums"), expand = c(.08, .05)) +
  ylab("Number of\nchildren") +
  theme_minimal(20) +
  theme(axis.title.y = element_text(angle = 0, hjust = 1)) +
  scale_fill_viridis_d(option = "B")

signal_alluvial_plot2

signaldf_long <-
  sfdf |>
  pivot_longer(cols = Sadness:Tantrums, names_to = "Signal", values_to = "Freq")  |>
  mutate(
    Signal = factor(Signal, levels = c("Sadness", "Crying", "Tantrums"))
  ) |>
  na.omit()

# barplot_SignalFreq <- ggplot(signaldf_long, aes(y=Freq, fill=Signal)) + geom_bar(position = "dodge")
# barplot_SignalFreq

ggplot(signaldf_long, aes(Signal, Freq, group = uniqueID)) + geom_line(show.legend = FALSE, alpha = .25 , position = position_jitter(h = .1, w = .1))

signaldf_long_sum <- signaldf_long |>
  group_by(Signal, Freq) |>
  summarise(N = n())

barplot_SignalFreq <- ggplot(signaldf_long_sum, aes(N, Signal, fill = fct_rev(Freq))) +
  geom_col(position = "stack", color = "grey") +
  scale_fill_viridis_d(option = "B", direction = -1) +
  labs(x = "Number of children", y = "") +
  guides(fill = guide_legend("Times per month:", nrow = 1, reverse = TRUE, position = "top")) +
  theme_minimal(15)

barplot_SignalFreq

ggplot(signaldf_long_sum, aes(Freq, N, color = Signal, group = Signal)) + geom_line()

# alluvial plot
x <- table(signalfreqdf$Sadness, signalfreqdf$Crying) |>
  as.data.frame() |>
  dplyr::filter(Freq > 3)

d2_conflict_filter <- modeldf |>
  mutate(
    ConflictFreq2 = factor(freq_short[ConflictFreqOF], levels = c(freq_short))
  )|>
  filter_at(vars(ConflictFreqOF),all_vars(!is.na(.)))

barplot_conflict <- ggplot(d2_conflict_filter, aes(x=ConflictFreq2, fill= Sex)) +
  geom_bar(position = "stack") +
  labs(title = "Conflict") +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 0.7) +
  labs(x = "Frequency of conflict (per month)", y = "Number of children") +
  theme_minimal(15)
barplot_conflict

d2_alloparent_filter <- modeldf |>
  mutate(
    AlloparentingFreq02 = factor(freq_short[AlloparentingFreq0], levels = c(freq_short))
  )|>
  filter_at(vars(AlloparentingFreq0),all_vars(!is.na(.)))

barplot_alloparenting <- ggplot(d2_alloparent_filter, aes(x= AlloparentingFreq02, fill = Sex)) +
  geom_bar(position = "stack") +
  labs(title = "Alloparenting") +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 0.7) +
  labs(x = "Frequency of child alloparenting (per month)", y = "Number of children") +
  theme_minimal(15)

barplot_alloparenting

## possible SI plots ----------------------------------------------------

#bar plot of household income
ggplot(caregivers, aes(y=IncomeCategory)) +
  geom_bar() + labs(y = "IncomeCaregory (household level)")

# scatterplot of CaregiverAge and ChildAge
# note 20 year old caregiver reporting on a 20 year old has a CaregiverRelation of "Other"
ggplot(data = modeldf, aes(CaregiverAge, ChildAge)) + geom_count() + geom_smooth() + labs(title = "Scatterplot of CaregiverAge and ChildAge")

# bar plot of ChildAge with Sex as color, but with the sample used for our main models
ggplot(modeldf, aes(y=ChildAge, fill = Sex)) +
  geom_bar()

# # bar plot of Neighborhood with ImmigrateUtila as color
caregivers2 <- caregivers |>  filter(!is.na(Neighborhood))
ggplot(caregivers2, aes(y=Neighborhood, fill=ImmigrateUtila)) +
  geom_bar() + labs(title = "Immigration status by neighborhood")

# Correlation matrix of predictor and outcome variables

maincormat <- d2[c("IllnessSusceptibilityMean", "EducationLevelYears", "OtherChildrenHH", "LogIncome", "ConflictFreqN", "number_adults", "PartnerStatus", "AlloparentingFreqN", "Sex", "ChildAge", "RelativeNeed3", "RelativeMaternalInvestment2", "SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalFreqMax", "SignalCost", "Family2")] |>
  mutate(
    PartnerStatus = as.numeric(PartnerStatus),
    Sex = as.numeric(Sex),
    RelativeNeed3 = as.numeric(RelativeNeed3),
    RelativeMaternalInvestment2 = as.numeric(RelativeMaternalInvestment2),
    Family2 = as.numeric(Family2),
  ) |>
  cor( use = "pairwise.complete.obs")

ggcorrplot(maincormat, hc.order = TRUE, lab = TRUE)


# model setup -------------------------------------------------------------

signals <- c("SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalFreqMax", "SignalCost")
std_predictors <- "ChildAge + Sex + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + NeighborhoodQuality + HouseQuality"

std_formulas <- str_glue("{signals} ~ {std_predictors} + OtherChildrenHH + IllnessSusceptibilityMean + (1|householdID)")
bodyfat_formulas <- str_glue("{signals} ~ {std_predictors} + OtherChildrenHH + BodyFat * Sex + (1|householdID)")
conflict_formula <- "ConflictFreqN ~ ChildAge + Sex + OlderKids + YoungerKids + LogIncome + AdultsNoChildcare + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + AdultsChildcare + OtherChildAlloparentingFreqN + (1|householdID)"

glmmTMBformulas <- c(std_formulas, bodyfat_formulas, conflict_formula)
names(glmmTMBformulas) <- c(paste0(signals, "illness"), paste0(signals, "bodyfat"), "ConflictFreqN")

glmmTMBmodels <- tibble(
    Outcome = c(signals, signals, "ConflictFreqN"),
    Model = map(glmmTMBformulas, \(f) glmmTMB(as.formula(f), family = nbinom2, data = d2)),
    AIC = map_dbl(Model, AIC)
)

polrNeedFormulas <- str_glue("RelativeNeed3 ~ {signals} + {std_predictors} + OtherChildrenHH + IllnessSusceptibilityMean")
polrInvestFormula <- paste0("RelativeMaternalInvestment2 ~ ", std_predictors, "+ OtherChildrenHH + IllnessSusceptibilityMean + RelativeNeed3")
polrResponseFormula <- paste0("CaregiverResponse ~ ", std_predictors, "+ OtherChildrenHH + IllnessSusceptibilityMean + UnwantedTask + Punishment2 + Family2 + DiscomfortPainInjuryIllness")

polrFormulas <- c(polrNeedFormulas, polrInvestFormula, polrResponseFormula)
names(polrFormulas) <- c(paste0("Need", signals), "Invest", "Response")

polrModels <- tibble(
    Model = map(polrFormulas, \(f) call("polr", f, data = quote(d2)) |> eval()),
    Test = map(Model, \(m) coeftest(m, vcov = vcovCL, type = "HC0", cluster = ~householdID))
)

summary(glmmTMBmodels$Model$SadFreqNillness)
summary(glmmTMBmodels$Model$CryFreqNillness)
summary(glmmTMBmodels$Model$TantrumFreqNillness)
summary(glmmTMBmodels$Model$SignalFreqillness)
summary(glmmTMBmodels$Model$SignalFreqMaxillness)
summary(glmmTMBmodels$Model$SignalCostillness)

summary(glmmTMBmodels$Model$ConflictFreqN)

summary(glmmTMBmodels$Model$SadFreqNbodyfat)
summary(glmmTMBmodels$Model$CryFreqNbodyfat)
summary(glmmTMBmodels$Model$TantrumFreqNbodyfat)
summary(glmmTMBmodels$Model$SignalFreqbodyfat)
summary(glmmTMBmodels$Model$SignalFreqMaxbodyfat)
summary(glmmTMBmodels$Model$SignalCostbodyfat)

m <- glmmTMB(Neighborhood2 ~ CaregiverAge + EducationLevelYears + NumberOfChildren + NeighborhoodQuality + HouseQuality + (1|householdID), family = binomial, data = modeldf, na.action = na.exclude)
m <- glm(Neighborhood2 ~ CaregiverAge + EducationLevelYears + NumberOfChildren + NeighborhoodQuality + HouseQuality, family = binomial, data = modeldf, na.action = na.exclude)
summary(m)
plot(allEffects(m))

x <- predict(m, type = "response")
table(x > 0.5, modeldf$Neighborhood2)

# add more plots

polrModels$Test$NeedSadFreqN
polrModels$Test$NeedCryFreqN
polrModels$Test$NeedTantrumFreqN
polrModels$Test$NeedSignalFreq
polrModels$Test$NeedSignalFreqMax
polrModels$Test$NeedSignalCost

polrModels$Test$Invest

# note: added predictors
polrModels$Test$Response


# Consideration of which controls to use
# - conflict - illness susceptibility ?+ CaregiverAge ?+ Neighborhood2
# less signaling with more adults (aligns with lower conflict freq in families with more adults)
# mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears +  Neighborhood2 + (1|householdID),data = d2, family = nbinom2)
# summary(mscH2)
# plot(allEffects(mscH2))
#
# mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + AlloparentingFreqN*Sex + Neighborhood2 + (1|householdID),data = d2, family = nbinom2)
# summary(mscH2)

# coeftest(mN1, vcov = vcovCL, type = "HC0", cluster = ~householdID)
# Shared variables
# ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears


# relative child need models
# polr(RelativeNeed3 ~ SIGNALOUTCOME(cost and freq) + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean, d2)
# Not universal: OtherChildrenHH + number_adults + IllnessSusceptibilityMean

# relative investment model
# mpr <- polr(RelativeMaternalInvestment2 ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + RelativeNeed3, d2)
# Not universal: OtherChildrenHH + number_adults + IllnessSusceptibilityMean
# Unique: RelativeNeed3

# response to last signal
# mpr <- polr(CaregiverResponse ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + UnwantedTask + IllnessSusceptibilityMean + Punishment2 + Family2 + DiscomfortPainInjuryIllness, d2)
# Not universal: OtherChildrenHH + number_adults + IllnessSusceptibilityMean
# Unique: UnwantedTask + IllnessSusceptibilityMean + Punishment2 + Family2 + DiscomfortPainInjuryIllness

# Models of signaling (6 measures) report effects of age, conflict, caregiver education (depends upon conflict being in the model), alloparenting*sex

## Disease resistance
# SIGNALMODELS_health <- glmmTMB(SIGNALOUTCOME ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
# Not universal: OtherChildrenHH + number_adults + IllnessSusceptibilityMean

## Body fat index
# SIGNALMODELS_bodyfat <- glmmTMB(SIGNALOUTCOME ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
# Not universal: OtherChildrenHH + number_adults
# Unique: BodyFat*Sex

# conflict model
# glmmTMB(ConflictFreqN ~ ChildAge + Sex + OlderKids + YoungerKids + LogIncome + AdultsNoChildcare + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + AdultsChildcare + OtherChildAlloparentingFreqN + (1|householdID),data = d2, family = nbinom2)
# Unique: OlderKids + YoungerKids + AdultsNoChildcare + AdultsChildcare + OtherChildAlloparentingFreqN

# main models ------------------------------------------------

# conflict is a strong predictor. However, signal cost could be influencing frequency of conflict
# include models without conflict

# signal cost
mscH <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
summary(mscH)
plot(allEffects(mscH))
# OtherChildAlloparents*Sex

mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus  + AlloparentingFreqN*Sex + EducationLevelYears + Neighborhood2 + (1|householdID),data = d2, family = nbinom2)
summary(mscH2)
plot(allEffects(mscH2))
# signal frequency

msfH <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
summary(msfH)
plot(allEffects(msfH))

# frequency of most common signal (running away excluded)

mmsfH <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
summary(mmsfH)
plot(allEffects(mmsfH))

# Signal specific models

msadH <- glmmTMB(SadFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
summary(msadH)
plot(allEffects(msadH))

mcryH <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
summary(mcryH)
plot(allEffects(mcryH))

mtantrumH <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
summary(mtantrumH)
plot(allEffects(mtantrumH))

# anthropometric models ---------------------------------------------------

# HeightZ, WeightZ, BMIZ, and BodyFatPercentage are not significant predictors
# (not close either) of signaling behavior in our standard models

# GripStrengthMean and GripR are not either (as main effects or interacted with Sex)

# Tricep residuals data = d2[d2$BodyFat < 6,]
mtriC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + TricepR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mtriC)
plot(allEffects(mtriC))

  # Tricep residuals = not a signficant predictor of signal frequency.
  # Effect on cost may be due to its association with tantrum frequency.
  # (Tantrum model below)

mtriCry <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + TricepR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mtriCry)
plot(allEffects(mtriCry))

mtriTantrum <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + TricepR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mtriTantrum)
plot(allEffects(mtriTantrum))

# Subscap residuals
msubC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + SubscapR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msubC)
plot(allEffects(msubC))

  # Subscap residuals = not a signficant predictor of signal frequency.
  # Effect on cost may be due to its association with tantrum frequency.
  # (Tantrum model below)

msubCry <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + SubscapR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msubCry)
plot(allEffects(msubCry))

msubTantrum <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + SubscapR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msubTantrum)
plot(allEffects(msubTantrum))

# Flexed residuals

mFlexC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2, family = nbinom2)
#mFlexC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2[d2$FlexedR < 10,], family = nbinom2)
summary(mFlexC)
plot(allEffects(mFlexC))
plot_predictions(mFlexC, condition = "FlexedR", residuals = T, vcov = T, points = 1)

# removing outliers from d2 based on FlexedR = body fat measures no longer being significant predictors
mFlexCb <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedRb + (1|householdID),data = d2, family = nbinom2)
summary(mFlexCb)
plot(allEffects(mFlexCb))

# remember this filter. Filters out children with largest bodyfat. Has impact on results that is worth noting.
# hist(d2[d2$FlexedR < 10,]$BodyFat)
# filter out those higher than 5.87100254 (8.33906184 & 8.53234405)

mFlexFb <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2[d2$FlexedR < 10,], family = nbinom2)
summary(mFlexFb)
plot(allEffects(mFlexFb))

mcryF <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR*Sex + (1|householdID),data = d2, family = nbinom2)
#mcryF <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR*Sex + (1|householdID),data = d2[d2$FlexedR < 10,], family = nbinom2)
summary(mcryF)
plot(allEffects(mcryF))

mtantrumF <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2, family = nbinom2)
summary(mtantrumF)
plot(allEffects(mtantrumF))

# do we want IllnessSusceptibilityMean in these models?

mBFf <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFf)
plot(allEffects(mBFf))

mBFfm <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFfm)
plot(allEffects(mBFfm))
# hist(d2[d2$FlexedR < 10,]$BodyFat)
# filter out those higher than 5.87100254 (8.33906184 & 8.53234405) data = d2[d2$BodyFat < 6,]
mBFc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFc)
plot(allEffects(mBFc))

mBFsad <- glmmTMB(SadFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFsad)
plot(allEffects(mBFsad))

mBFcry <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFcry)
plot(allEffects(mBFcry))

# BodyFat Estimate: -0.418049  Std. Error 0.227308 z value: -1.839  p: 0.06590
mBFtantrum <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFtantrum)
plot(allEffects(mBFtantrum))

mtest <- glmmTMB(HeightZ ~ FoodSecurityMean + Sex + (1|householdID), data = d2, family = gaussian())
summary(mtest)
plot(allEffects(mtest))

mtest <- glmmTMB(WeightZ ~ FoodSecurityMean + Sex + (1|householdID), data = d2, family = gaussian())
summary(mtest)
plot(allEffects(mtest))

mtest <- glmmTMB(BodyFat ~ FoodSecurityMean*Sex + (1|householdID), data = d2, family = gaussian())
summary(mtest)
plot(allEffects(mtest))

mtest <- glmmTMB(FlexedR ~ FoodSecurityMean*Sex + (1|householdID), data = d2, family = gaussian())
summary(mtest)
plot(allEffects(mtest))

anthropometricMeansWide <-
  anthropometricMeans0 |>
  dplyr::filter(
    measurements == 2
  ) |>
  pivot_wider(names_from = "Year", values_from = BicepMean:BMI, id_cols = ChildID) |>
  left_join(d2, by = "ChildID") |>
  dplyr::select(
    householdID,
    ChildAge,
    FoodSecurityMean,
    Sex,
    SadFreqN,
    CryFreqN,
    TantrumFreqN,
    contains("2023"),
    contains("2024"),
    -contains("Wrist"),
    -contains("Tricep"),
    -contains("Subscap"),
    -contains("Seated")
  ) |>
  dplyr::filter(if_any(SadFreqN:TantrumFreqN, ~!is.na(.x)))

explore_long <- function(signal, outcome){
  frmla <- str_glue("{outcome}_2024 ~ {outcome}_2023 + {signal} + Sex + ChildAge + FoodSecurityMean + (1|householdID)")
  glmmTMB(as.formula(frmla), data = anthropometricMeansWide, family = gaussian)
}

outcomes <-
  names(anthropometricMeansWide[8:29]) |>
  str_remove("_2023|_2024") |>
  unique()

# Commented out for now
# e <-
#   expand_grid(Signal = c("SadFreqN", "CryFreqN", "TantrumFreqN"), Outcome = outcomes) |>
#   mutate(
#     Model = map2(Signal, Outcome, explore_long),
#     Summary = map(Model, summary),
#     mainZ = map_dbl(Summary, ~.x$coefficients$cond[3,3]),
#     interactionZ = map_dbl(Summary, ~.x$coefficients$cond[6,3]),
#     possible = abs(mainZ) >= 2 | abs(interactionZ) >= 2
#   )

m <- glmmTMB(WeightMean_KG_2024 ~ WeightMean_KG_2023 + CryFreqN*Sex + ChildAge + (1|householdID), data = anthropometricMeansWide, family = gaussian)
summary(m)
plot(allEffects(m))

# age confound?
m <- glmmTMB(BMI_2024 ~ BMI_2023 + CryFreqN + Sex + (1|householdID), data = anthropometricMeansWide, family = gaussian)
summary(m)
plot(allEffects(m))

m <- glmmTMB(GripStrengthMean_2024 ~ GripStrengthMean_2023 + CryFreqN + Sex + (1|householdID), data = anthropometricMeansWide, family = gaussian)
summary(m)
plot(allEffects(m))

m <- glmmTMB(WeightMean_KG_2024 ~ WeightMean_KG_2023 + CryFreqN*ChildAge + Sex  + (1|householdID), data = anthropometricMeansWide, family = gaussian)
summary(m)
plot(allEffects(m))

weight_long_cry <- glmmTMB(WeightMean_KG_2024 ~ WeightMean_KG_2023 + CryFreqN*ChildAge + (1|householdID), data = anthropometricMeansWide, family = gaussian, na.action = na.exclude)
summary(weight_long_cry)
plot(allEffects(weight_long_cry))

anthropometricMeansWide$weightR <- residuals(weight_long_cry)
ggplot(anthropometricMeansWide, aes(ChildAge, CryFreqN, colour = weightR)) + geom_jitter(size = 3) + scico::scale_color_scico(palette = 'vik', midpoint = 0)
ggplot(anthropometricMeansWide, aes(ChildAge, weightR, colour = CryFreqN)) + geom_jitter(size = 3) + scico::scale_color_scico(palette = 'vik', midpoint = 0)

p <- plot_predictions(weight_long_cry, condition = c("CryFreqN", "ChildAge"), vcov = TRUE) +
  geom_rug(data = weight_long_cry$frame, aes(x = CryFreqN , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b") +
  scale_color_viridis_d(option = "A", end = .8) +
  guides(color = guide_legend(reverse = TRUE, title = "Child age (years)"), fill = guide_legend(reverse = TRUE, title = "Child age (years)")) +
  xlab("Crying freq. (per month)") +
  ylab("Weight (2024)") +
  ylim(40, 55) +
  theme_minimal(15)
p$layers[[2]]$aes_params$linewidth <- 2
p$layers[[1]]$aes_params$alpha <- 0
p

height_long_cry <- glmmTMB(HeightMean_2024 ~ HeightMean_2023 + CryFreqN*ChildAge + Sex + (1|householdID), data = anthropometricMeansWide, family = gaussian, na.action = na.exclude)
summary(height_long_cry)
plot(allEffects(height_long_cry))

anthropometricMeansWide$heightR <- residuals(height_long_cry)
anthropometricMeansWide$Sex2 <- as.factor(anthropometricMeansWide$Sex)
ggplot(anthropometricMeansWide, aes(ChildAge, CryFreqN, colour = heightR)) + geom_jitter(size = 3) + scico::scale_color_scico(palette = 'vik', midpoint = 0)
ggplot(anthropometricMeansWide, aes(ChildAge, heightR, colour = CryFreqN)) + geom_jitter(size = 3) + scico::scale_color_scico(palette = 'vik', midpoint = 0)

m <- glmmTMB(HeightMean_2024 ~ HeightMean_2023 + Sex2 + (1|householdID), data = anthropometricMeansWide, family = gaussian, na.action = na.exclude)
summary(m)
anthropometricMeansWide$heightR <- residuals(m)

m <- gam(HeightMean_2024 ~ Sex2 + s(HeightMean_2023) + s(ChildAge, by = Sex2, k = 3), data = anthropometricMeansWide, na.action = na.exclude)
summary(m)
draw(m)
anthropometricMeansWide$heightR <- residuals(m)

ggplot(anthropometricMeansWide, aes(SadFreqN, heightR, colour = Sex2)) +
  geom_jitter(size = 3) +
  geom_smooth(method = 'lm') +
  scale_color_binary()


m <- glmmTMB(heightR ~ SadFreqN * Sex2 + (1|householdID), data = anthropometricMeansWide, family = gaussian, na.action = na.exclude)
summary(m)

m <- gam(HeightMean_2024 ~ Sex2 * SadFreqN + s(HeightMean_2023) + s(ChildAge, by = Sex2, k = 3), data = anthropometricMeansWide, na.action = na.exclude)
summary(m)
plot(m, all.terms = T)
plot_predictions(m, condition = c("SadFreqN", "Sex2"), points = 1) +
  theme_minimal(15)


p2 <- plot_predictions(height_long_cry, condition = c("CryFreqN", "ChildAge"), vcov = TRUE) +
  geom_rug(data = height_long_cry$frame, aes(x = CryFreqN , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b") +
  scale_color_viridis_d(option = "A", end = .8) +
  guides(color = guide_legend(reverse = TRUE, title = "Child age (years)"), fill = guide_legend(reverse = TRUE, title = "Child age (years)")) +
  xlab("Crying freq. (per month)") +
  ylab("Height (2024)") +
  ylim(130, 145) +
  theme_minimal(15)
p2$layers[[2]]$aes_params$linewidth <- 2
p2$layers[[1]]$aes_params$alpha <- 0
p2

longitudinal_plot <- p + p2 +
  plot_layout(axes = "collect_x", ncol = 1, byrow = FALSE, guides = "collect") &
  theme(legend.position = "top")

# causality going the other direction (weight causes more crying?)

out <- avg_comparisons(height_long_cry, variables = list("CryFreqN" = c(0, 30), "ChildAge" = c(5,14)), cross = TRUE)
out <- avg_comparisons(height_long_cry, variables = list("CryFreqN" = c(0, 8), "ChildAge" = c(5,14)), cross = TRUE)

# ConflictFreqN model -----------------------------------------------------------

# d2$AdultsNoChildcare <- d2$number_adults - d2$AdultsChildcare

mconflict <- glmmTMB(ConflictFreqN ~ ChildAge + Sex + OlderKids + YoungerKids + LogIncome + AdultsNoChildcare + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + AdultsChildcare + OtherChildAlloparentingFreqN + (1|householdID),data = d2, family = nbinom2)
summary(mconflict)
plot(allEffects(mconflict))

conflict_scatterplot <- ggplot(modeldf, aes(ConflictFreqN, SignalCost)) + geom_count() + geom_smooth(method='lm')

# RunawayFreqN model ------------------------------------------------------

# fat measures marginally significant despite severe loss of nobs 225 -> 83

mrunawayH <- glmmTMB(RunawayFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(mrunawayH)
plot(allEffects(mrunawayH))


# childcare within the family models ----------------------------------

# keep adults and children on their original scales

# measures = or number of caretakers

# note childcare measure often not significant wihtout BodyFat in the model


#Signal Cost
table(d2$AdultsNoChildcare)

mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherChildAlloparentingFreqN + AdultsChildcare + (1|householdID),data = d2, family = nbinom2)
summary(mscH2)
plot(allEffects(mscH2))

mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + AdultsNoChildcare + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherChildAlloparentingFreqN + AdultsChildcare + (1|householdID),data = d2, family = nbinom2)
summary(mscH2)
plot(allEffects(mscH2))

modeldf2 <- d2 |>
  dplyr::select(
    householdID,
    childHHid,
    AdultsChildcare,
    AdultsNoChildcare
  ) |>
  na.omit()

#Signal Frequency

msf2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareB + (1|householdID),data = d2, family = nbinom2)
summary(msf2)
plot(allEffects(msf2))

msf2i <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareB*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(msf2i)
plot(allEffects(msf2i))

# Maximum Signal Frequency

mmsf2 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareB + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2)
plot(allEffects(mmsf2))

mmsf2i <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareB*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2i)
plot(allEffects(mmsf2i))

# Adult care lightly weighted

#Signal Cost

mscH3 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears+ BodyFat + HHChildCareA + (1|householdID),data = d2, family = nbinom2)
summary(mscH3)
plot(allEffects(mscH3))

mscH3i <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareA*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mscH3i)
plot(allEffects(mscH3i))

#Signal Frequency

msf2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareA + (1|householdID),data = d2, family = nbinom2)
summary(msf2)
plot(allEffects(msf2))

# Maximum Signal Frequency

mmsf2 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareA + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2)
plot(allEffects(mmsf2))

mmsf2i <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + HHChildCareA*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2i)
plot(allEffects(mmsf2i))

# Number of alloparents

hist(d2$OtherCare) # Few beyond 4, most 3 and below

#Signal Cost

mscH3 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherCare + (1|householdID),data = d2, family = nbinom2)
summary(mscH3)
plot(allEffects(mscH3))

mscH3i <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherCare*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mscH3i)
plot(allEffects(mscH3i))

#Signal Frequency

msf2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherCare + (1|householdID),data = d2, family = nbinom2)
summary(msf2)
plot(allEffects(msf2))

msf2i <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherCare*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(msf2i)
plot(allEffects(msf2i))

# Maximum Signal Frequency

mmsf2 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherCare + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2)
plot(allEffects(mmsf2))

mmsf2i <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat + OtherCare*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2i)
plot(allEffects(mmsf2i))

# child alloparenting effort model ----------------------------------------

# OtherChildAlloparentingFreqN < 121) included due to unrealistic predictions from the model for alloparenting frequency
d3 <- d2 |>
  dplyr::filter(YoungerKids > 1) # , OtherChildAlloparentingFreqN < 121

# model cannot handle high alloparenting family

# care (adult's heavily weighted)
malloparenting <- glmmTMB(AlloparentingFreqN ~ ChildAge + Sex + YoungerKids + OlderKids + LogIncome + AdultsChildcare + OtherChildAlloparents + (1|householdID),data = d3, family = nbinom2)
summary(malloparenting)
plot(allEffects(malloparenting))

# d3$Sex2 <- as.factor(d3$Sex)
#
# m <- mgcv::gam(AlloparentingFreqN ~ ChildAge + Sex2 + s(YoungerKids, k = 3) + s(OlderKids, by = Sex2, k = 3) + LogIncome + AdultsChildcare + s(OtherChildAlloparentingFreqN) + s(householdID, bs = "re"),data = d3, family = quasipoisson)
# summary(m)
# visreg(m, scale = "response")

d3$AlloparentingProp <- d3$AlloparentingFreqN/60 # 60 = maximum value

malloparenting2 <- glmmTMB(AlloparentingProp ~ ChildAge + Sex + YoungerKids + OlderKids + LogIncome + AdultsChildcare + OtherChildAlloparents + (1|householdID),data = d3, family = binomial)
summary(malloparenting2)
plot(allEffects(malloparenting2))
#visreg(malloparenting2, scale = "response")

ggplot(d3, aes(OtherChildAlloparents, AlloparentingFreqN)) + geom_count() + geom_smooth()
ggplot(d3, aes(OtherChildAlloparents, AlloparentingFreqN, size = CryFreqN, color = CryFreqN)) + geom_jitter() + scale_color_viridis_c(option = "A", begin = .2, end = .9)
ggplot(drop_na(d3, Sex), aes(size = OtherChildAlloparents, x = AlloparentingFreqN, y = CryFreqN, color = OtherChildAlloparents)) + geom_jitter() + scale_color_viridis_c(option = "A", begin = .2, end = .9) + geom_smooth(method = "lm") + facet_wrap(~Sex)

# Leaning against including these models of alloparenting freq
# consider putting into our models of signaling vs alloparenting the number of other child alloparents
# higher numbers of other child alloparents seems to be a protective factor at high alloparenting frequency
ggplot(drop_na(d3, Sex), aes(size = OtherChildAlloparents, x = AlloparentingFreqN, y = SignalCost, color = OtherChildAlloparents)) + geom_jitter() + scale_color_viridis_c(option = "A", begin = .2, end = .9) + geom_smooth(method = "lm") + facet_wrap(~Sex)

ggplot(d2, aes(SignalCost, color = factor(Sex))) + stat_ecdf()
ggplot(drop_na(d2, Sex), aes(AlloparentingFreqN, color = factor(Sex))) + stat_ecdf()
# relatedness within the family -------------------------------------------

# Child relatedness (relatedness relative to every kid in household) = not a significant
# predictor

# Sibling relatedness

# SignalCost and SignalFreqMax = marginally significant

# SignalFreq model
mSRF <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + MeanSiblingRelatedness + (1|householdID),data = d2, family = nbinom2)
summary(mSRF)
plot(allEffects(mSRF))

# with YoungerKids*ChildAge
mSRF3 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OlderKids + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + MeanSiblingRelatedness*number_adults + YoungerKids*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mSRF3)
plot(allEffects(mSRF3))

# ordinal logistic regression models for relative need  --------

# Results with NeighborhoodF as a random effect and vc0 seem like outliers compared to all other methods

# Using preferred model for signal frequency (without partner status)
# signal Freq and Child age still significant after controlling for education years
# IllnessSusceptibilityMean has no predictive ability when added


# mN1 <- polr(RelativeNeed3 ~ SignalFreq + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults+ ConflictFreqN + AlloparentingFreqN, d2)
# summary(mN1)
# coeftest(mN1, vcov = vcovCL, type = "HC0", cluster = ~householdID)

mN1 <- polr(RelativeNeed3 ~ SignalFreq + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean, d2)
summary(mN1)
coeftest(mN1, vcov = vcovCL, type = "HC0", cluster = ~householdID)

p_need_age_sf <- plot_predictions(mN1, condition = c("ChildAge", "group"), type = "probs", vcov = ~householdID) +
  ylab("Relative need") +
  xlab("Child age (years)") +
  scale_color_viridis_d(option = "B", end = .8) +
  scale_fill_viridis_d(option = "B", end = .8) +
  scale_linewidth_ordinal(range = c(5,10)) +
  guides(color = guide_legend(override.aes = list(linewidth=2.5))) +
  labs(color = "", fill = "") +
  theme_minimal() +
  theme(legend.position = "top", plot.tag.position = "right") +
  ylim(0,1)
p_need_age_sf

p_need_signalfreq <-
  plot_predictions(mN1, condition = c("SignalFreq", "group"), type = "probs", vcov = ~householdID) +
  ylab("Relative need") +
  xlab("Signal frequency (times per month)") +
  scale_color_viridis_d(option = "B", end = .8) +
  scale_fill_viridis_d(option = "B", end = .8) +
  labs(color = "", fill = "") +
  theme_minimal() +
  theme(plot.tag.position = "right")
p_need_signalfreq

signaling_plot_need_freq <- (p_need_age_sf) +
  (p_need_signalfreq + theme(legend.position = "none", legend.title = element_blank())) +
  plot_layout(ncol = 1, byrow = FALSE, axes = "collect_y") +
  plot_annotation(title = "", tag_levels = "A") & theme(plot.tag.position = c(.98, .92))
signaling_plot_need_freq

mN2 <- polr(RelativeNeed3 ~ SignalCost + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean, d2)
summary(mN2)
coeftest(mN2, vcov = vcovCL, type = "HC0", cluster = ~householdID)

p_need_age_cost <- plot_predictions(mN2, condition = c("ChildAge", "group"), type = "probs", vcov = ~householdID) +
  ylab("Relative need") +
  xlab("Child age (years)") +
  scale_color_viridis_d(option = "B", end = .8) +
  scale_fill_viridis_d(option = "B", end = .8) +
  scale_linewidth_ordinal(range = c(5,10)) +
  guides(color = guide_legend(override.aes = list(linewidth=2.5))) +
  labs(color = "", fill = "") +
  theme_minimal() +
  theme(legend.position = "top", plot.tag.position = "right") +
  ylim(-0.05,1)
p_need_age_cost

p_need_signal_cost <-
  plot_predictions(mN2, condition = c("SignalCost", "group"), type = "probs", vcov = ~householdID) +
  ylab("Relative need") +
  xlab("Signal cost") +
  scale_color_viridis_d(option = "B", end = .8) +
  scale_fill_viridis_d(option = "B", end = .8) +
  labs(color = "", fill = "") +
  theme_minimal() +
  theme(plot.tag.position = "right")
p_need_signal_cost

p_need_age_cost <- (p_need_age_cost) +
  (p_need_signal_cost + theme(legend.position = "none", legend.title = element_blank())) +
  plot_layout(ncol = 1, byrow = FALSE, axes = "collect_y") +
  plot_annotation(title = "", tag_levels = "A") & theme(plot.tag.position = c(.98, .92))
p_need_age_cost

p_need_age_sc <- plot_predictions(mN2, condition = c("ChildAge", "group"), type = "probs", vcov = ~householdID)
p_need_age_sc

p_need_signalcost <- plot_predictions(mN2, condition = c("SignalCost", "group"), type = "probs", vcov = ~householdID)
p_need_signalcost

(ci <- confint(mN1))
exp(cbind(coef(mN1),(ci)))

# brm
#mN1b <- brm(RelativeNeed3 ~ SignalFreq + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN + IllnessSusceptibilityMean + (1|householdID), d2, family = cumulative(link = "logit", threshold = "flexible"))
#save(mN1b, file = "mN1b.rda")
load(file = "mN1b.rda")
summary(mN1b)
plot_predictions(mN1b, condition = c("ChildAge", "group"), type = "response")

# Using preferred model for signal frequency (with partner status)
mI <- polr(RelativeMaternalInvestment2 ~ SignalFreq + ChildAge*Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + IllnessSusceptibilityMean, d2)
summary(mI)
coeftest(mI, vcov = vcovCL, type = "HC0", cluster = ~householdID)

(ci <- confint(mI))
exp(cbind(coef(mI),(ci)))

plot_predictions(mI, condition = c("OtherChildrenHH", "group"), type = "probs")
plot_predictions(mI, condition = c("ChildAge", "group"), type = "probs")

# Want/need models
summary(table(d2$WantNeed2, d2$CaregiverResponse))
plot(table(d2$WantNeed2, d2$Sex))
summary(table(d2$WantNeed2, d2$Sex))
plot(table(d2$WantNeed2, d2$SadFreqN))

m<- glmmTMB(SadFreqN ~ ChildAge*WantNeed2 + (1|householdID), family = nbinom2, d2)
summary(m)
plot(allEffects(m))

m2 <- glmmTMB(CryFreqN ~ ChildAge + WantNeed2 + (1|householdID), family = nbinom2, d2)
summary(m2)
plot(allEffects(m2))

m3 <- glmmTMB(TantrumFreqN ~ ChildAge + WantNeed2 + (1|householdID), family = nbinom2, d2)
summary(m3)
plot(allEffects(m3))

mI2 <- polr(CaregiverResponse ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + WantNeed2, d2)
summary(mI2)
coeftest(mI2, vcov = vcovCL, type = "HC0", cluster = ~householdID)

mWantNeed <- glm(WantNeedBinary ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults, data = d2, family = "binomial")
summary(mWantNeed)
coeftest(mWantNeed, vcov = vcovCL, type = "HC0", cluster = ~householdID)

# ordinal logistic regression of parental response ------------------------

# note 6 instances of both positive and negative coded as NA
table(modeldf$PositiveResponse, modeldf$NegativeResponse)

# mpr <- polr(CaregiverResponse ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + UnwantedTask + Punishment2 + Family2 + DiscomfortPainInjuryIllness, d2)
# summary(mpr)
# coeftest(mpr, vcov = vcovCL, type = "HC0", cluster = ~householdID)

polrModels$Test$Response
mpr <- polr(CaregiverResponse ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + UnwantedTask + Punishment2 + Family2 + DiscomfortPainInjuryIllness, d2)
summary(mpr)
coeftest(mpr, vcov = vcovCL, type = "HC0", cluster = ~householdID)

p_response_punishment <- plot_predictions(polrModels$Model$Response, condition = c("Punishment2", "group"), type = "probs", vcov = ~householdID) +
  xlab("Cause involved punishment") +
  ylab("Caregiver response") +
  scale_color_viridis_d(option = "B", end = .8) +
  guides(color = guide_legend(override.aes = list(linewidth=1))) +
  theme_minimal() + # theme_minimal needs to come before plot.tag.position
  theme(legend.position = "top", plot.tag.position = "right") +
  ylim(0,1)

p_response_pain <- plot_predictions(polrModels$Model$Response, condition = c("DiscomfortPainInjuryIllness", "group"), type = "probs", vcov = ~householdID) +
  xlab("Cause involved discomfort, pain, injury, illness") +
  ylab("Caregiver response") +
  scale_fill_discrete(guide="none") +
  scale_color_viridis_d(option = "B", end = .8) +
  theme_minimal() +
  theme(plot.tag.position = "right") +
  labs(color = "", fill = "")

signaling_plot_last_signal_response <- (p_response_punishment + theme(legend.position = "top", legend.title = element_blank())) +
  (p_response_pain + theme(legend.position = "none", legend.title = element_blank())) +
  plot_layout(ncol = 1, byrow = FALSE, axes = "collect_y") +
  plot_annotation(title = "", tag_levels = "A") & theme(plot.tag.position = c(.98, .92))

signaling_plot_last_signal_response

# extracting estimates from model

# converted to characters for the purpose of avg_predictions function
# d2$UnwantedTask <- as.character(d2$UnwantedTask)
# d2$DiscomfortPainInjuryIllness <- as.character(d2$DiscomfortPainInjuryIllness)
# d2$Family2 <- as.numeric(d2$Family2)

# out_mpr <- avg_predictions(mpr, newdata = datagrid(Punishment2 = "1", Sex = c("Male", "Female")))
# out_mpr$estimate[3] # 0.320612

# Note: Output is identical with LogIncome = mean(mpr$model$LogIncome) in
# the call and with it not. This suggests that it correctly handles
# numeric variables on its own?
# out_mpr_full <- avg_predictions(mpr, newdata = datagrid(Punishment2 = "1", Sex = c("Male", "Female"), Family2 = c("0", "1"), DiscomfortPainInjuryIllness = c("0", "1"), UnwantedTask = c("0", "1"), LogIncome = mean(mpr$model$LogIncome)))
# out_mpr_full
# out_mpr_full$estimate[3]

updated_out <- avg_predictions(mpr, newdata = datagrid(Punishment2 = "1", grid_type = "balanced"))
updated_out$estimate[3]
# d2$Family2 <- as.character(d2$Family2)
# out_mpr_full2 <- avg_predictions(mpr, newdata = datagrid(Punishment2 = "1", Sex = c("Male", "Female"), Family2 = c("0", "1"), DiscomfortPainInjuryIllness = c("0", "1"), UnwantedTask = c("0", "1"), LogIncome = mean(mpr$model$LogIncome)))
# out_mpr_full2
# out_mpr_full2$estimate[3]
# mean(as.numeric(mpr$model$Family2))

# The factors are almost working correctly based on this:
# Default estimate for group = 1 when family2 is removed is 0.434,
# same as family = 0. Same estimate is .310 when family 2 = 1.
# Mean of this = .372, same as when both are included in the data grid
# as follows: Family2 = c("0", "1").

# However, family2 = mean(mpr$model$Family2) and family2 = .5, both give
# slightly different (but very similar answers)

test_grid <- datagrid(model = mpr, Punishment2 = "1", grid_type = "balanced")
avg_predictions(mpr, newdata = datagrid(Punishment2 = "1", grid_type = "balanced"))
predictions(mpr, newdata = datagrid(Punishment2 = unique))
# out_mpr_full
# out_mpr_full$estimate[3]

# relative need predicts relative investment ------------------------------

# mpr <- polr(RelativeMaternalInvestment2 ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + RelativeNeed3, d2)
# summary(mpr)
# coeftest(mpr, vcov = vcovCL, type = "HC0", cluster = ~householdID)

mpI <- polr(RelativeMaternalInvestment2 ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + RelativeNeed3, d2)
summary(mpI)
coeftest(mpI, vcov = vcovCL, type = "HC0", cluster = ~householdID)

mpI_plot_data <- plot_predictions(polrModels$Model$Invest, condition = c("RelativeNeed3", "group"), type = "probs", vcov = ~householdID, draw = FALSE)

plot_relative_investment <- ggplot(mpI_plot_data, aes(RelativeNeed3, estimate, color = group, group = group, ymin = conf.low, ymax = conf.high)) +
  geom_line(position = position_dodge(width = .2), linewidth = 1.5) +
  geom_pointrange(position = position_dodge(width = .2)) +
  theme_minimal(15) +
  xlab("\nRelative need") +
  ylab("Proportion of children") +
  #scico::scale_color_scico_d(palette = "lipari", begin = .3, end = .8) +
  scale_color_viridis_d(option = "B", begin = .25, end = .85) +
  guides(color = guide_legend(title = "Relative investment:")) &
  theme(legend.position = "top")

plot_predictions(mpI, condition = c("OtherChildrenHH", "group"), type = "probs", vcov = ~householdID)

mpI2 <- polr(RelativeMaternalInvestment2 ~ ChildAge + Sex + YoungerKids + OlderKids + LogIncome + number_adults + RelativeNeed3, d2)
summary(mpI2)
coeftest(mpI2, vcov = vcovCL, type = "HC0", cluster = ~householdID)

plot_predictions(mpI2, condition = c("YoungerKids", "group"), type = "probs", vcov = ~householdID)

plot_relative_investment_kids <- plot_predictions(mpI2, condition = c("YoungerKids", "group"), type = "probs", vcov = ~householdID) +
  xlab("\nOtherChildrenHH") +
  ylab("Proportion of children") +
  scale_color_viridis_d(option = "B", end = .8) +
  scale_fill_viridis_d(option = "B", end = .8) +
  scale_linewidth_ordinal(range = c(5,10)) +
  guides(color = guide_legend(override.aes = list(linewidth=2.5))) +
  labs(color = "", fill = "") +
  theme_minimal(15) +
  theme(legend.position = "top", plot.tag.position = "right") +
  guides(color = guide_legend(title = "Relative investment:")) +
  ylim(-.03,1)

plot_relative_investment_kids

mpr_plot_data2 <- plot_predictions(mpI, condition = c("RelativeNeed3", "group"), type = "probs", vcov = ~householdID, draw = FALSE)

# NOTE error ribbons dip below 0
plot_predictions(mpI2, condition = c("YoungerKids", "group"), type = "probs", vcov = ~householdID)

plot_predictions(mpI2, condition = c("YoungerKids", "group"), type = "probs", vcov = ~householdID) +
  ylab("Proportion") +
  xlab("\nYounger children") +
  scale_color_viridis_d(option = "B", end = .8) +
  scale_fill_viridis_d(option = "B", end = .8) +
  scale_linewidth_ordinal(range = c(5,10)) +
  guides(color = guide_legend(override.aes = list(linewidth=2.5)), title = "Relative investment:") +
  labs(color = "", fill = "") +
  ylim(-.03 ,1) +
  theme_minimal() &
  theme(legend.position = "top")

#  2 kids with values for relative need but no other kids in house.
# Could be due to children elsewhere.
modeldf2 <- d2 |>
  dplyr::select(
    householdID,
    childHHid,
    Sex,
    ChildAge,
    OtherChildrenHH,
    LogIncome,
    number_adults,
    RelativeNeed3
  ) |>
  na.omit()

# plots -------------------------------------------------------


## correlation matrices  ---------------------------------------------------
maincormat <- d2[c("IllnessSusceptibilityMean", "EducationLevelYears", "OtherChildrenHH", "LogIncome", "ConflictFreqN", "number_adults", "PartnerStatus", "AlloparentingFreqN", "Sex", "ChildAge", "RelativeMaternalInvestment2", "SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalFreqMax", "SignalCost")] |>
  mutate(
    PartnerStatus = as.numeric(PartnerStatus),
    Sex = as.numeric(Sex),
    RelativeMaternalInvestment2 = as.numeric(RelativeMaternalInvestment2)
  ) |>
  cor( use = "pairwise.complete.obs")

ggcorrplot(maincormat, hc.order = TRUE, hc.method = "ward.D",lab = TRUE, lab_col = "black", lab_size = 4.5) +
  scico::scale_fill_scico(palette = "vik", midpoint = 0, begin = .1, end = .9)

signal_subset <- d2[c("ConflictFreqN", "Sex", "ChildAge", "SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalFreqMax", "SignalCost")]
names(signal_subset) <- shortform_dict[names(signal_subset)]

smallcormat <- signal_subset |>
  mutate(
    Sex = as.numeric(Sex)
  ) |>
  cor( use = "pairwise.complete.obs")

signal_corrplot <- ggcorrplot(smallcormat, hc.order = TRUE, hc.method = "ward.D",lab = TRUE, lab_col = "black", lab_size = 4.5) +
  scico::scale_fill_scico(palette = "vik", midpoint = 0, begin = .1, end = .9, limits = c(-1, 1)) +
  guides(fill = guide_colorbar(title = "Correlation coefficients"))
signal_corrplot


## age --------------------------------------------------------------------

plot_age <- function(m, ylabel){
  p <- suppressWarnings(plot_predictions(m, condition = c("ChildAge", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp")) +
    geom_rug(data = m$frame, aes(x = ChildAge , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123)) +
    ylab(ylabel) +
    xlab("Child age (years)") +
    xlim(5,20) +
    ylim(0,NA) +
    scale_color_binary() +
    scale_fill_binary() +
    theme_bw(15)

  p$layers[[2]]$aes_params$linewidth <- 2
  return(p)
}

p_age_sad <- plot_age(glmmTMBmodels$Model$SadFreqNillness, "Sadness freq.")
p_age_cry <- plot_age(glmmTMBmodels$Model$CryFreqNillness, "Crying freq.")
p_age_tantrum <- plot_age(glmmTMBmodels$Model$TantrumFreqNillness, "Tantrum freq.")
p_age_freq <- plot_age(glmmTMBmodels$Model$SignalFreqillness, "Summed freq.")
p_age_freq_max <- plot_age(glmmTMBmodels$Model$SignalCostillness, "Max freq.")
p_age_cost <- plot_age(glmmTMBmodels$Model$SignalFreqMaxillness, "Signal cost")


signaling_data_plot_age <- p_age_sad + p_age_cry + p_age_tantrum + p_age_freq + p_age_freq_max + p_age_cost +
  plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE, guides = "collect") &
  # plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.98, .92), legend.position = "top")
  # plot_annotation(title = "Signaling measures") & theme(legend.position = "top")
  theme(legend.position = "top")

ggsave("Figures/signaling_data_plot_age.pdf", signaling_data_plot_age, width = 9, height = 9)

signaling_data_plot_age

# model for rugs for plots without functions
# geom_rug(data = msadH$frame, aes(x = ChildAge , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))

## conflict (independent var)----------------------------------------------------------

plot_conflict <- function(m, ylabel){
  p <- suppressWarnings(plot_predictions(m, condition = c("ConflictFreqN", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp")) +
    geom_rug(data = m$frame, aes(x = ConflictFreqN , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123)) +
    ylab(ylabel) +
    xlab("Conflict freq. (per month)") +
    scale_color_binary() +
    scale_fill_binary() +
    theme_bw(15)

  p$layers[[2]]$aes_params$linewidth <- 2
  return(p)
}

p_conflict_sad <- plot_conflict(msadH, "Sadness freq.")
p_conflict_cry <- plot_conflict(mcryH, "Crying freq.")
p_conflict_tantrum <- plot_conflict(mtantrumH, "Tantrum freq.")
p_conflict_freq <- plot_conflict(msfH, "Summed freq.")
p_conflict_cost <- plot_conflict(mscH, "Signal cost")
p_conflict_freq_max <- plot_conflict(mmsfH, "Max freq.")

signaling_data_plot_conflict <- p_conflict_sad + p_conflict_cry + p_conflict_tantrum + p_conflict_freq + p_conflict_freq_max + p_conflict_cost +
  plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE, guides = "collect") &
  # plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.98, .92), legend.position = "none")
  # plot_annotation(title = "Signaling measures") & theme(legend.position = "none")
  theme(legend.position = "top")

ggsave("Figures/signaling_data_plot_conflict.pdf", signaling_data_plot_conflict, width = 9, height = 9)

# p_conflict_sad <- plot_predictions(msadH, condition = "ConflictFreqN", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Sadness freq.") +
#   xlab("Conflict freq. (per month)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ConflictFreqN , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_conflict_cry <- plot_predictions(mcryH, condition = "ConflictFreqN", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Crying freq.") +
#   xlab("Conflict freq. (per month)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ConflictFreqN , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_conflict_tantrum <- plot_predictions(mtantrumH, condition = "ConflictFreqN", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Tantrum freq.") +
#   xlab("Conflict freq. (per month)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ConflictFreqN , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_conflict_freq <- plot_predictions(msfH, condition = "ConflictFreqN", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Summed freq.") +
#   xlab("Conflict freq. (per month)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ConflictFreqN , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_conflict_cost <- plot_predictions(mscH, condition = "ConflictFreqN", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Signal cost") +
#   xlab("Conflict freq. (per month)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ConflictFreqN , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_conflict_freq_max <- plot_predictions(mmsfH, condition = "ConflictFreqN", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Max freq.") +
#   xlab("Conflict freq. (per month)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ConflictFreqN , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# signaling_data_plot_conflict <- p_conflict_sad + p_conflict_cry + p_conflict_tantrum + p_conflict_freq + p_conflict_freq_max + p_conflict_cost +
#   plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE) +
#   # plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.98, .92), legend.position = "none")
#   plot_annotation(title = "Signaling measures") & theme(legend.position = "none")
#
# signaling_data_plot_conflict

## conflict (dependent var) ------------------------------------------------
p_conflict_dep_var_age <-
  plot_predictions(glmmTMBmodels$Model$ConflictFreqN, condition = c("ChildAge", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_bw(15) +
  ylab("Conflict freq. (per month)") +
  xlab("Child age (years)") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = glmmTMBmodels$Model$ConflictFreqN$frame, aes(x = ChildAge , y = 0, color = Sex), position = position_jitter(width = .3, seed = 123), sides = "b") +
  theme(legend.position = "none")
  #theme(plot.tag.position = "topright")
p_conflict_dep_var_age$layers[[2]]$aes_params$linewidth <- 2

p_conflict_dep_var_age

p_conflict_dep_var_aCare <-
  plot_predictions(glmmTMBmodels$Model$ConflictFreqN, condition = c("AdultsChildcare", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_bw(15) +
  ylab("Conflict freq. (per month)") +
  xlab("Number of adults that contribute childcare") +
  scale_color_binary() +
  scale_fill_binary() +
  #theme(plot.tag.position = "topright")
  geom_rug(data = glmmTMBmodels$Model$ConflictFreqN$frame, aes(x = AdultsChildcare , y = 0, color = Sex), position = position_jitter(width = .3, seed = 123), sides = "b") +
 theme(legend.position = "none")

p_conflict_dep_var_aCare$layers[[2]]$aes_params$linewidth <- 2
p_conflict_dep_var_aCare

signaling_data_plot_conflict_dep_var <-
  p_conflict_dep_var_age + p_conflict_dep_var_aCare +
  plot_layout(axes = "collect_y", ncol = 1, byrow = FALSE, guides = "collect") +
  plot_annotation(title = "", tag_levels = "A") & theme(plot.tag.position = c(.95, .95)) &
  theme(legend.position = "top")

signaling_data_plot_conflict_dep_var


## education --------------------------------------------------------------
plot_education <- function(m, ylabel){
  p <- suppressWarnings(plot_predictions(m, condition = c("EducationLevelYears", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp")) +
    geom_rug(data = m$frame, aes(x = EducationLevelYears, y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123)) +
    ylab(ylabel) +
    xlab("Caregiver education (years)") +
    scale_color_binary() +
    scale_fill_binary() +
    theme_bw(15)

  p$layers[[2]]$aes_params$linewidth <- 2
  return(p)
}

p_education_sad <- plot_education(msadH, "Sadness freq.")
p_education_cry <- plot_education(mcryH, "Crying freq.")
p_education_tantrum <- plot_education(mtantrumH, "Tantrum freq.")
p_education_freq <- plot_education(mscH, "Summed freq.")
p_education_cost <- plot_education(msfH, "Signal cost")
p_education_freq_max <- plot_education(mmsfH, "Max freq.")

signaling_data_plot_education <- p_education_sad + p_education_cry + p_education_tantrum + p_education_freq + p_education_freq_max + p_education_cost +
  plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE, guides = "collect") &
  # plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.98, .92))
  #plot_annotation(title = "Signaling measures")
  theme(legend.position = "top")

ggsave("Figures/signaling_data_plot_education.pdf", signaling_data_plot_education, width = 9, height = 9)

# p_edu_sad <- plot_predictions(msadH, condition = "EducationLevelYears", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Sadness freq.") +
#   xlab("Caregiver education (years)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = EducationLevelYears , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_edu_cry <- plot_predictions(mcryH, condition = "EducationLevelYears", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Crying freq.") +
#   xlab("Caregiver education (years)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = EducationLevelYears , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_edu_tantrum <- plot_predictions(mtantrumH, condition = "EducationLevelYears", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Tantrum freq.") +
#   xlab("Caregiver education (years)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = EducationLevelYears , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_edu_freq <- plot_predictions(msfH, condition = "EducationLevelYears", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Summed freq.") +
#   xlab("Caregiver education (years)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = EducationLevelYears , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_edu_cost <- plot_predictions(mscH, condition = "EducationLevelYears", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Signal cost") +
#   xlab("Caregiver education (years)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = EducationLevelYears , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_edu_freq_max <- plot_predictions(mmsfH, condition = "EducationLevelYears", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Max freq.") +
#   xlab("Caregiver education (years)") +
#   geom_rug(data = drop_na(d2, Sex), aes(x = EducationLevelYears , y = 0), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# signaling_data_plot_edu <- p_edu_sad + p_edu_cry + p_edu_tantrum + p_edu_freq + p_edu_freq_max + p_edu_cost +
#   plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE) +
#   # plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.98, .92))
#   plot_annotation(title = "Signaling measures")
#
# signaling_data_plot_education

## alloparenting:sex  --------------------------------------------------------------------
plot_allo <- function(d, ylab){
  ggplot(d, aes(AlloparentingFreqN, estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Sex), alpha = .25) +
    geom_line(aes(color = Sex)) +
    geom_rug(data = d$frame, aes(x = AlloparentingFreqN, y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b") +
    scale_color_binary() +
    scale_fill_binary() +
    xlab("Alloparenting frequency (times per month)") +
    ylab(ylab) +
    theme_minimal(15)
}

pointsize <- .25

summary(glmmTMBmodels$Model$SadFreqNillness)
summary(glmmTMBmodels$Model$CryFreqNillness)
summary(glmmTMBmodels$Model$TantrumFreqNillness)
summary(glmmTMBmodels$Model$SignalFreqillness)
summary(glmmTMBmodels$Model$SignalFreqMaxillness)
summary(glmmTMBmodels$Model$SignalCostillness)

d_allo_by_sex_sad <- plot_predictions(glmmTMBmodels$Model$SadFreqNillness, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_sad <- plot_allo(d_allo_by_sex_sad, "Sadness freq.")

d_allo_by_sex_cry <- plot_predictions(glmmTMBmodels$Model$CryFreqNillness, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_cry <- plot_allo(d_allo_by_sex_cry, "Crying freq.")

d_allo_by_sex_tantrum <- plot_predictions(glmmTMBmodels$Model$TantrumFreqNillness, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_tantrum <- plot_allo(d_allo_by_sex_tantrum, "Tantrum freq.")

d_allo_by_sex_freq <- plot_predictions(glmmTMBmodels$Model$SignalFreqillness, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_freq <- plot_allo(d_allo_by_sex_freq, "Summed freq.")

d_allo_by_sex_freq_max <- plot_predictions(glmmTMBmodels$Model$SignalFreqMaxillness, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_freq_max <- plot_allo(d_allo_by_sex_freq_max, "Signal cost")

d_allo_by_sex_cost <- plot_predictions(glmmTMBmodels$Model$SignalCostillness, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_cost <- plot_allo(d_allo_by_sex_cost, "Max freq.")

signaling_data_plot_allo_by_sex <- p_allo_by_sex_sad + p_allo_by_sex_cry + p_allo_by_sex_tantrum + p_allo_by_sex_freq + p_allo_by_sex_freq_max + p_allo_by_sex_cost +
  plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE, guides = "collect") +
  # plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.98, .92), legend.position = "top")
  plot_annotation(title = "Signaling measures") & theme(legend.position = "top")

signaling_data_plot_allo_by_sex

# p_allo_by_sex_sad <- plot_predictions(msadH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp") +
#   scale_color_binary() +
#   theme_bw(15) +
#   ylab("Sadness freq.") +
#   xlab("Alloparenting freq. (per month)")
#
# p_allo_by_sex_sad
#
# p_allo_by_sex_cry <- plot_predictions(mcryH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Crying freq.") +
#   xlab("Alloparenting freq. (per month)")
#
# p_allo_by_sex_tantrum <- plot_predictions(mtantrumH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Tantrum freq.") +
#   xlab("Alloparenting freq. (per month)")
#
# p_allo_by_sex_freq <- plot_predictions(msfH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Summed freq.") +
#   xlab("Alloparenting freq. (per month)")
#
# p_allo_by_sex_cost <- plot_predictions(mscH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Signal cost") +
#   xlab("Alloparenting freq. (per month)")
#
# p_allo_by_sex_freq_max <- plot_predictions(mmsfH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Max freq.") +
#   xlab("Alloparenting freq. (per month)")
#
# signaling_data_plot_allo_by_sex <- p_allo_by_sex_sad + p_allo_by_sex_cry + p_allo_by_sex_tantrum + p_allo_by_sex_freq + p_allo_by_sex_freq_max + p_allo_by_sex_cost +
#   plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE, guides = "collect") +
#   plot_annotation(title = "Signaling measures", tag_levels = "A")
#
# signaling_data_plot_allo_by_sex


## alloparenting (dependent var) -------------------------------------------

alloparenting_dep_var_plot <- plot_predictions(malloparenting, condition = c("OlderKids"), vcov = TRUE, points = .25, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Frequency of alloparenting effort") +
  xlab("Number of older children in household")
alloparenting_dep_var_plot

## body fat ---------------------------------------------------------------

p_bodyfat_cost <- plot_predictions(glmmTMBmodels$Model$SignalCostbodyfat, condition = c("BodyFat"), vcov = TRUE, points = .25, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Signal cost") +
  xlab("Body fat composite")

p_bodyfat_cost2 <- plot_predictions(glmmTMBmodels$Model$SignalCostbodyfat, condition = c("BodyFat", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp", conf_level = .95) +
  theme_minimal(15) +
  ylab("Signal cost") +
  xlab("Body fat composite") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = glmmTMBmodels$Model$SignalCostbodyfat$frame, aes(x = BodyFat , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))
p_bodyfat_cost2$layers[[2]]$aes_params$linewidth <- 2

p_bodyfat_crying <- plot_predictions(glmmTMBmodels$Model$CryFreqNbodyfat, condition = "BodyFat", vcov = TRUE, points = .25, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Crying freq.") +
  xlab("Body fat composite")

p_bodyfat_crying2 <- plot_predictions(glmmTMBmodels$Model$CryFreqNbodyfat, condition = c("BodyFat", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Crying freq.") +
  xlab("Body fat composite") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = glmmTMBmodels$Model$CryFreqNbodyfat$frame, aes(x = BodyFat , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))
p_bodyfat_crying2$layers[[2]]$aes_params$linewidth <- 2

p_bodyfat_tantrum <- plot_predictions(glmmTMBmodels$Model$TantrumFreqNbodyfat, condition = "BodyFat", vcov = TRUE, points = .25, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Tantrum freq.") +
  xlab("Body fat composite")

p_bodyfat_tantrum2 <- plot_predictions(glmmTMBmodels$Model$TantrumFreqNbodyfat, condition = c("BodyFat", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Tantrum freq.") +
  xlab("Body fat composite") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = glmmTMBmodels$Model$TantrumFreqNbodyfat$frame, aes(x = BodyFat , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))
p_bodyfat_tantrum2$layers[[2]]$aes_params$linewidth <- 2

signaling_data_plot_bodyfat <- p_bodyfat_cost + p_bodyfat_crying + p_bodyfat_tantrum +
  plot_layout(axes = "collect_x", ncol = 1, byrow = FALSE) +
  plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.95, .95))

signaling_data_plot_bodyfat

signaling_data_plot_bodyfat2 <- (p_bodyfat_cost2 + coord_cartesian(clip = "on", ylim = c(0, 100))) + (p_bodyfat_crying2 + coord_cartesian(clip = "on", ylim = c(0, 25))) + (p_bodyfat_tantrum2 + coord_cartesian(clip = "on", ylim = c(0, 25))) +
  plot_layout(axes = "collect_x", ncol = 1, byrow = FALSE, guides = "collect") +
  plot_annotation(title = "Signaling measures") &
  theme(legend.position = "top")

signaling_data_plot_bodyfat2
# x <- summary(mBFtantrum)
# x$coefficients$cond["BodyFat", "Estimate"]
# x$coefficients$cond["BodyFat", "Pr(>|z|)"]


# tables and dataframes to extract values from ---------------------------------------
sex_table <- table(modeldf$Sex)

rel_need_table <- table(modeldf$RelativeNeed3)

response_overlap <- table(modeldf$PositiveResponse, modeldf$NegativeResponse)

modeldf_FULLSIG <- modeldf |>
  dplyr::filter(!is.na(CryFreqN) & !is.na(SadFreqN) & !is.na(TantrumFreqN))

modeldf_FULLANTHROPOMETRICS <- modeldf |>
  dplyr::filter(!is.na(SubscapR))

sex_table2 <- table(modeldf_FULLANTHROPOMETRICS$Sex)

causematrix <-
  causes |>
  mutate(
    across(Family:StatusConcerns, as.numeric)
  ) |>
  dplyr::select(Family:StatusConcerns) |>
  na.omit()

hagenheat(t(causematrix), hc_method = "ward.D2")

out <- pvclust(causematrix, method.hclust = "ward.D2")
x <- cutree(out$hclust, k = 3)

causes2 <- causes |>
  rowwise() |>
  mutate(
    FamilyConflict = if_any(all_of(names(x)[x == 1]), \(x) x == 1),
    Adversity = if_any(all_of(names(x)[x == 2]), \(x) x == 1),
    Transgression = if_any(all_of(names(x)[x == 3]), \(x) x == 1),
    CauseType = case_when(
      Transgression ~ "Transgression",
      FamilyConflict ~ "Family Conflict",
      Adversity ~ "Adversity",
      .default = "Unknown"
    )
  )

# If the cause categories were not mutually exclusive, transgression received top priority, followed by family conflict.

d2b <- left_join(d2, causes2[c("uniqueID", "CauseType")])

response_last_signal <- d2b |>
  dplyr::select(
    householdID,
    childHHid,
    ChildAge,
    Sex,
    uniqueID,
    CaregiverResponse,
    PositiveResponseBinary,
    NeutralResponseBinary,
    NegativeResponseBinary,
    DiscomfortPainInjuryIllness,
    Punishment2,
    Family2,
    CauseType
  ) |>
  mutate(
    NeutralResponseBinary2 = ifelse((NeutralResponseBinary == 1) & ((PositiveResponseBinary == 1) | (NegativeResponseBinary == 1)), 0, NeutralResponseBinary),
    MixedResponse = ifelse(PositiveResponseBinary == 1 & NegativeResponseBinary == 1, 1, 0)
  ) |>
  dplyr::filter(
    !((PositiveResponseBinary == 0) & (NeutralResponseBinary == 0) & (NegativeResponseBinary == 0)),
    MixedResponse != 1
    )

response_last_signal_long <-
  response_last_signal |>
  pivot_longer(cols =  c("PositiveResponseBinary", "NeutralResponseBinary2", "NegativeResponseBinary"), names_to = "Response", values_to = "Present") |>
  mutate(
    Response = factor(Response, levels = c("PositiveResponseBinary", "NeutralResponseBinary2", "NegativeResponseBinary"), labels = c("Positive", "Neutral", "Negative")),
    Cause = case_when(
      Punishment2 == 1 ~ "Punishment",
      DiscomfortPainInjuryIllness == 1 ~ "Pain",
      Family2 == 1 ~ "Familial conflict",
      .default = "Other"
    )
  ) |>
  dplyr::filter(Present != 0)


response_last_signal_long_sum <- response_last_signal_long |>
  mutate(
    AgeCategory = ifelse(ChildAge < 9, "Younger children", "Older children"),
    AgeCategory = factor(AgeCategory, levels = c("Younger children", "Older children"))
  ) |>
  group_by(AgeCategory, CauseType, Response) |>
  summarise(N = n())

barplot_CaregiverResponse <- ggplot(response_last_signal_long_sum, aes(N, CauseType, fill = Response)) +
  geom_col(position = "stack", color = "grey") +
  scale_fill_viridis_d(option = "B", direction = -1) +
  scale_y_discrete(limits = rev) +
  labs(x = "Number of children", y = "") +
  guides(fill = guide_legend("Caregiver response:", nrow = 1, reverse = TRUE, position = "top")) +
  theme_minimal(15) +
  facet_wrap(~AgeCategory)

# note lack of major sex difference

barplot_CaregiverResponse

last_signal_response_tbl_hurt <- table(d2$CaregiverResponse, d2$DiscomfortPainInjuryIllness)
last_signal_response_tbl_punish <- table(d2$CaregiverResponse, d2$Punishment2)

# model statistics --------------------------------------------------------

library (glue)
library (tidyverse)
library (broom)
tdy <- function (m) {
  tidy (m, conf.int = T) |>
mutate(
  # Adjust string to taste
  str_old = glue ("$\\beta={signif(estimate, 2)}$ ({signif(conf.low, 2)}, {signif(conf.high, 2)})"),
  str2 = glue("p = {signif(p.value, 2)}"),
  str = glue ("$\\beta={signif(estimate, 2)}$, 95% CI:{signif(conf.low, 2)} - {signif(conf.high, 2)}, p = {signif(p.value, 2)}"),
) %>%
split (.$term)
}

models <- list (msadH = glmmTMBmodels$Model$SadFreqNillness, mcryH = glmmTMBmodels$Model$CryFreqNillness, mtantrumH = glmmTMBmodels$Model$TantrumFreqNillness, msfH = glmmTMBmodels$Model$SignalFreqillness, mscH = glmmTMBmodels$Model$SignalCostillness, mmsfH = glmmTMBmodels$Model$SignalFreqMaxillness, mBFc = glmmTMBmodels$Model$SignalCostbodyfat, mBFcry = glmmTMBmodels$Model$CryFreqNbodyfat, mBFtantrum = glmmTMBmodels$Model$TantrumFreqNbodyfat, mconflict = glmmTMBmodels$Model$ConflictFreqN, malloparenting = malloparenting)
stats <- map (models, tdy)
stats$mBFc$BodyFat$str
stats$mBFc$BodyFat$str2


anthropometricMeansWide$Sex2 <- as.factor(anthropometricMeansWide$Sex)
# need to control for current height and age (puberty effects)
m <- mgcv::gam(HeightMean_2024 ~ Sex2 + s(HeightMean_2023, by = Sex2) + s(ChildAge, by = Sex2, k = 3), data = anthropometricMeansWide, na.action = na.exclude)
summary(m)
draw(m)
anthropometricMeansWide$heightR <- residuals(m)

# filtered out males who were sad daily or more often
# m2 <- glmmTMB(heightR ~ SadFreqN * Sex2 + FoodSecurityMean + ChildAge + (1|householdID), data = anthropometricMeansWide[anthropometricMeansWide$SadFreqN <30,], family = gaussian, na.action = na.exclude)

# no filter
m2 <- glmmTMB(heightR ~ SadFreqN * Sex2 + FoodSecurityMean + ChildAge + (1|householdID), data = anthropometricMeansWide, family = gaussian, na.action = na.exclude)
summary(m2)
plot(allEffects(m2))

# p_long_sad <- plot_predictions(m2, condition = c("SadFreqN", "Sex2"), points = 1, vcov = T, draw = T)
# p_long_sad$layers[[1]] <- geom_point(data = anthropometricMeansWide[anthropometricMeansWide$SadFreqN <30,], aes(x = SadFreqN, y = heightR, colour = Sex2), position = position_dodge(width = 0.5))

out <- plot_predictions(m2, condition = c("SadFreqN", "Sex2"), vcov = T, draw = F) |>
  dplyr::filter(SadFreqN < 9)

# caption needs to not we are not displaying the male outliders on sadfreq
p_long_sad <- ggplot(out, aes(x = SadFreqN)) +
  # geom_boxplot(data = anthropometricMeansWide[anthropometricMeansWide$SadFreqN <30,], aes(x = SadFreqN, y = heightR, color = Sex2), position = position_dodge(0.5)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Sex2), alpha = .2) +
  geom_line(aes(y = estimate, color = Sex2), linewidth = 1) +
  geom_point(data = anthropometricMeansWide[anthropometricMeansWide$SadFreqN <30,], aes(y = heightR, color = Sex2), size = 3, position = position_dodge(0.2)) +
  labs(x = "Sadness freq. (per month)", y = "Height residuals") +
  guides(color = guide_legend(title = ""), fill = guide_none()) +
  scale_x_continuous(breaks = c(0, 1, 3, 8)) +
  scale_color_binary() +
  scale_fill_binary() +
  theme_minimal(15) +
  theme(panel.grid.minor = element_blank())
p_long_sad

# how old are 0 girls and how old are 8
ggplot(anthropometricMeansWide[anthropometricMeansWide$SadFreqN <30,], aes(SadFreqN, ChildAge, color = heightR)) + geom_jitter(size = 3) + facet_wrap(~Sex) + scale_color_viridis_c(option = "A")



ggplot(drop_na(anthropometricMeansWide, Sex2), aes(SadFreqN, heightR, colour = Sex2)) +
  geom_jitter(size = 3) +
  geom_smooth(method = 'lm') +
  scale_color_binary()

ggplot(anthropometricMeansWide[anthropometricMeansWide$Sex2 == "F",], aes(SadFreqN, heightR)) +
  geom_jitter(size = 3) +
  geom_smooth(method = 'lm')


m <- mgcv::gam(WeightMean_KG_2024 ~ Sex2 + s(WeightMean_KG_2023) + s(ChildAge, by = Sex2, k = 3), data = anthropometricMeansWide, na.action = na.exclude)
summary(m)
draw(m)
anthropometricMeansWide$WeightR <- residuals(m)

m <- glmmTMB(WeightR ~ SadFreqN * Sex2 + (1|householdID), data = anthropometricMeansWide, family = gaussian, na.action = na.exclude)
summary(m)
plot(allEffects(m))

ggplot(anthropometricMeansWide, aes(SadFreqN, WeightR, colour = Sex2)) +
  geom_jitter(size = 3) +
  geom_smooth(method = 'lm') +
  scale_color_binary()

  # New comment to test ZED
