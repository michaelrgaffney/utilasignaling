# loading libraries -------------------------------------------------------

#+ message=F, warning=F, fig.height=9, fig.width=10

# note some are no longer used
library(utiladata2023) # data package
library(tidyverse)
# library(lme4)
library(effects)
library(visreg)
library(emmeans)
library(car)
library(haven)
library(MASS)
library(glmmTMB) # try this negative binomial (random effect for neighborhood)
#library(gee)
#library(geepack)
library(sandwich)
library(lmtest)
library(ggcorrplot)
#library(geeM)
library(marginaleffects)
library(ordinal) # clmm # polr
library(brms)
library(hagenutils)
library(localgrowth)
library(knitr)
library(kableExtra)
library(modelsummary)
library(ggalluvial)

source("recode.R")


# functions used ----------------------------------------------------------

# function to invert the scoring on a Likert scale
invert <- function(x){
  ((x - max(x, na.rm = TRUE)) * -1) + min(x, na.rm = TRUE)
  }

# functions to calculate values from multiple columns in ways that 1) NAs are not counted
# 2) yet do not result in NAs if values exist for other columns.
sum2 <- function(..., na.rm = F){
  v <- as.numeric(list(...))
  if(all(is.na(v))) return(NA)
  sum(v, na.rm = na.rm)
}

max2 <- function(..., na.rm = F){
  v <- as.numeric(list(...))
  if(all(is.na(v))) return(NA)
  max(v, na.rm = na.rm)
}

mean2 <- function(..., na.rm = F){
  v <- as.numeric(list(...))
  if(all(is.na(v))) return(NA)
  mean(v, na.rm = na.rm)
}


# dataframe prep ----------------------------------------------------------

anthropometricMeans <- anthropometrics |>
  dplyr::select(
    ChildID,
    AgeAtMeasurement = Age,
    Sex,
    Date.of.measurement,
    contains("Mean")
  ) |>
  mutate(
    measurements = n(),
    .by = ChildID
  ) |>
  mutate(
    BMI = WeightMean_KG / (HeightMean / 100 )^2
  ) |>
  dplyr::filter(
    measurements == 1 | year(Date.of.measurement) == 2023
  )

anthropometricMeans$WeightZ <- growthRef(AgeAtMeasurement, WeightMean_KG, Sex, anthropometricMeans, type = "Weight", pop = "CDC")
anthropometricMeans$HeightZ <- growthRef(AgeAtMeasurement, HeightMean, Sex, anthropometricMeans, type = "Height", pop = "CDC")
anthropometricMeans$BMIZ <- growthRef(AgeAtMeasurement, BMI, Sex, anthropometricMeans, type = "BMI", pop = "CDC")

m <- mgcv::gam(GripStrengthMean ~ Sex+s(AgeAtMeasurement, by = as.factor(Sex), k = 4), data = anthropometricMeans)
anthropometricMeans$GripR <- residuals(m)

m2 <- mgcv::gam(TricepMean ~ Sex+s(AgeAtMeasurement, by = as.factor(Sex)), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$TricepR <- residuals(m2)

m3 <- mgcv::gam(SubscapMean ~ Sex+s(AgeAtMeasurement, by = as.factor(Sex), k = 2), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$SubscapR <- residuals(m3)

m4 <- mgcv::gam(BodyFatPercentageMean ~ Sex+s(AgeAtMeasurement, by = as.factor(Sex), k = 3), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$BodyFatPercentageR <- residuals(m4)

anthropometricMeans$Sex <- NULL

d <- children |>
  #dplyr::filter(CompleteSurvey) %>%
  mutate(
    YoungerKids = map_int(ChildAge, \(x) sum(x > ChildAge)),
    OlderKids = map_int(ChildAge, \(x) sum(x < ChildAge)),
    .by = householdID
  ) |>
  mutate(
    SadFreqOF = factor(
      SadFreq,
      ordered = TRUE,
      levels = c("Never", "Once a month or less", "More than once a month but less than once a week", "More than once a week but not daily", "Daily", "Multiple times per day")),
    SadFreqN = case_when(
      SadFreq == "Never" ~ 0,
      SadFreq == "Once a month or less" ~ 1,
      SadFreq == "More than once a month but less than once a week" ~ 3,
      SadFreq == "More than once a week but not daily" ~ 8,
      SadFreq == "Daily" ~ 30,
      SadFreq == "Multiple times per day" ~ 60,
      #.default = SadFreq
      ),
    CryFreqOF = factor(
      CryFreq,
      ordered = TRUE,
      levels = c("Never", "Once a month or less", "More than once a month but less than once a week", "More than once a week but not daily", "Daily", "Multiple times per day")),
    CryFreqN = case_when(
      CryFreq == "Never" ~ 0,
      CryFreq == "Once a month or less" ~ 1,
      CryFreq == "More than once a month but less than once a week" ~ 3,
      CryFreq == "More than once a week but not daily" ~ 8,
      CryFreq == "Daily" ~ 30,
      CryFreq == "Multiple times per day" ~ 60,
      #.default = CryFreq
      ),
    TantrumFreqOF = factor(
      TantrumFreq,
      ordered = TRUE,
      levels = c("Never", "Once a month or less", "More than once a month but less than once a week", "More than once a week but not daily", "Daily", "Multiple times per day")),
    TantrumFreqN = case_when(
      TantrumFreq == "Never" ~ 0,
      TantrumFreq == "Once a month or less" ~ 1,
      TantrumFreq == "More than once a month but less than once a week" ~ 3,
      TantrumFreq == "More than once a week but not daily" ~ 8,
      TantrumFreq == "Daily" ~ 30,
      TantrumFreq == "Multiple times per day" ~ 60,
      #.default = TantrumFreq
      ),
    RunawayFreqOF = factor(
      RunawayFreq,
      ordered = TRUE,
      levels = c("Never", "Once a month or less", "More than once a month but less than once a week", "More than once a week but not daily", "Daily", "Multiple times per day")),
    RunawayFreqN = case_when(
      RunawayFreq == "Never" ~ 0,
      RunawayFreq == "Once a month or less" ~ 1,
      RunawayFreq == "More than once a month but less than once a week" ~ 3,
      RunawayFreq == "More than once a week but not daily" ~ 8,
      RunawayFreq == "Daily" ~ 30,
      RunawayFreq == "Multiple times per day" ~ 60,
      #.default = RunawayFreq
      ),
    ConflictFreqOF = factor(
      ConflictFreq,
      ordered = TRUE,
      levels = c("Never", "Once a month or less", "More than once a month but less than once a week", "More than once a week but not daily", "Daily", "Multiple times per day")),
    ConflictFreqN = case_when(
      ConflictFreq == "Never" ~ 0,
      ConflictFreq == "Once a month or less" ~ 1,
      ConflictFreq == "More than once a month but less than once a week" ~ 3,
      ConflictFreq == "More than once a week but not daily" ~ 8,
      ConflictFreq == "Daily" ~ 30,
      ConflictFreq == "Multiple times per day" ~ 60,
      #.default = ConflictFreq
      ),
    AlloparentingFreq0 = case_when(
      is.na(AlloparentingFreq) & ChildAge < 5 ~ "Never",
      is.na(AlloparentingFreq) & NumberOfChildren == 1 ~ "Never",
      is.na(AlloparentingFreq) & Siblings == 0 ~ "Never", # check later
      is.na(AlloparentingFreq) & YoungerKids == 0 ~ "Never",
      .default = AlloparentingFreq
    ),
    AlloparentingFreq0 = factor(
      AlloparentingFreq0,
      ordered = TRUE,
      levels = c("Never", "Once a month or less", "More than once a month but less than once a week", "More than once a week but not daily", "Daily", "Multiple times per day")),
    AlloparentingFreqN = case_when(
      AlloparentingFreq0 == "Never" ~ 0,
      AlloparentingFreq0 == "Once a month or less" ~ 1,
      AlloparentingFreq0 == "More than once a month but less than once a week" ~ 3,
      AlloparentingFreq0 == "More than once a week but not daily" ~ 8,
      AlloparentingFreq0 == "Daily" ~ 30,
      AlloparentingFreq0 == "Multiple times per day" ~ 60,
      #.default = AlloparentingFreq
    ),
    ChildAlloparent = ifelse(AlloparentingFreqN >3, 1, 0),
    RelativeMaternalInvestment = factor( #NOTE this is not always the mom
      RelativeMaternalInvestment,
      ordered = TRUE,
      levels = c("Much less", "Less", "The same amount", "More", "Much more")),
    RelativeMaternalInvestment2 = case_when(
      RelativeMaternalInvestment == "Much more" ~ "More",
      .default = RelativeMaternalInvestment),
    RelativeNeed = factor(
      RelativeNeed,
      ordered = TRUE,
      levels = c("Much less", "Less", "The same amount", "More", "Much more")),
    RelativeNeedRecode1 = case_when(
      RelativeNeed == "Much less" ~ "Less",
      .default = RelativeNeed),
    RelativeNeedRecode2 = case_when(
      RelativeNeed == "Much more" ~ "More",
      RelativeNeed == "Much less" ~ "Less",
      .default = RelativeNeed),
    RelativeMaternalInvestment2 = factor( # NOTE this is not always the mom
      RelativeMaternalInvestment2,
      ordered = TRUE,
      levels = c("Less", "The same amount", "More")),
    RelativeNeed2 = factor(
      RelativeNeedRecode1,
      ordered = TRUE,
      levels = c("Less", "The same amount", "More", "Much more")),
    RelativeNeed3 = factor(
      RelativeNeedRecode2,
      ordered = TRUE,
      levels = c("Less", "The same amount", "More")),
    Currenthealth1 = invert(Currenthealth1),
    Currenthealth2 = invert(Currenthealth2),
    Currenthealth3 = Currenthealth3 + 0, # to remove labels
    Sex = as.factor(Sex),
    ChildAge = as.numeric(ChildAge),
    Age0_2.5 = ifelse(ChildAge < 2.5, 1, 0), #Helfrecht and Meehan 2015 style risk periods
    Age2.5_5 = ifelse(ChildAge >=2.5 & ChildAge <5, 1, 0),
    Age5_10 = ifelse(ChildAge >=5 & ChildAge <10, 1, 0),
    Age10Plus = ifelse(ChildAge >= 10, 1, 0),
    Age5Minus = ifelse(ChildAge < 5, 1, 0),
    PositiveResponseBinary = ifelse(PositiveResponse >= 1, 1, 0),
    PositiveResponseF = as.factor(PositiveResponseBinary),
    NeutralResponseBinary = ifelse(NeutralResponse >= 1, 1, 0),
    NeutralResponseF = as.factor(NeutralResponseBinary),
    NegativeResponseBinary = ifelse(NegativeResponse >= 1, 1, 0),
    NegativeResponseF = as.factor(NegativeResponseBinary),
    ParentalResponse = case_when(
      PositiveResponse > 0 & NegativeResponse == 0 ~ "1",
      PositiveResponse == 0 & NegativeResponse == 0 ~ "0",
      PositiveResponse == 0 & NegativeResponse > 0 ~ "-1"
    ),
    ParentalResponse = factor(ParentalResponse, ordered = TRUE, levels = c("-1", "0", "1"))
  ) |>
  rowwise() |>
  mutate(
    SignalFreq = sum2(SadFreqN, CryFreqN, TantrumFreqN, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    SignalFreqMax = max2(SadFreqN, CryFreqN, TantrumFreqN, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    SignalCost = sum2(SadFreqN, CryFreqN * 2, TantrumFreqN * 3, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    CurrentHealthSum = sum2(Currenthealth1, Currenthealth2, Currenthealth3, na.rm = TRUE),
    CurrentHealthMean = mean2(Currenthealth1, Currenthealth2, Currenthealth3, na.rm = TRUE)
  ) |>
  ungroup() |>
  left_join(anthropometricMeans, by = "ChildID")

# add household data to d2
d <- left_join(d, dplyr::select(caregivers, householdID, CPRatio, Neighborhood, ImmigrateUtila, IncomeCategory,
                                EducationLevel, AdultsMoney, AdultsHousework, AdultsChildcare, CaregiverAge,
                                number_children2, number_adults, UserLanguage, contains("CaregiverMarital"),
                                contains("SocialSupport"), contains("FoodSecurity")), by = "householdID") |>
  left_join(dplyr::select(causes, -householdID, -childHHid, -ChildID, -Sad1), by = "uniqueID")

# prepare household variables for analyses
d2 <- d |>
  mutate(
    EducationLevel = case_when(
      EducationLevel == "Fourth grade" ~ "4th grade or less",
      EducationLevel == "Tenth grade" ~ "Some high school",
      EducationLevel == "Eleventh grade" ~ "Some high school",
      .default = EducationLevel),
    IncomeCategory = ifelse(IncomeCategory == "Prefer not to answer", NA, IncomeCategory),
    IncomeCategoryN = case_when(
      IncomeCategory == "Less than 5,000 Lempira per month" ~ 2500,
      IncomeCategory == "5,001 to 10,000 Lempira per month" ~ 7500,
      IncomeCategory == "10,001 to 20,000 Lempira per month" ~ 15000,
      IncomeCategory == "20,001 to 40,000 Lempira per month" ~ 30000,
      IncomeCategory == "40,001 to 80,000 Lempira per month" ~ 60000,
      IncomeCategory == "Over 80,000 Lempira per month" ~ 120000
      ),
    LogIncome = log10(IncomeCategoryN),
    RecentSignalSeverity = com1[as.character(Sad1)],
    RecentSignalCause = com2[as.character(Sad1)],
    AdultsHousework = as.numeric(AdultsHousework),
    AdultsChildcare = as.numeric(AdultsChildcare),
    Neighborhood2 = as.numeric(Neighborhood == "Camponado/Campolancho"),
    NeighborhoodF = as.factor(Neighborhood),
    OtherChildrenHH = number_children2 - 1,
    IncomeCategory = factor(
      IncomeCategory,
      ordered = TRUE,
      levels = c("Less than 5,000 Lempira per month", "5,001 to 10,000 Lempira per month", "10,001 to 20,000 Lempira per month", "20,001 to 40,000 Lempira per month", "40,001 to 80,000 Lempira per month", "Over 80,000 Lempira per month")),
    IncomeCategory2 = as.numeric(IncomeCategory),
    EducationLevel = factor(
      EducationLevel,
      ordered = TRUE,
      levels = c("No education", "First grade", "Second grade", "4th grade or less", "Fifth grade", "Completed 6th grade", "Seventh grade", "Completed 8th grade", "Ninth grade", "Some high school", "Completed high school", "Some university", "Professional certificate", "Completed university", "Post-graduate degree")),
    EducationLevelYears = case_when( #
      EducationLevel == "No education" ~ "0",
      EducationLevel == "First grade" ~ "1",
      EducationLevel == "Second grade" ~ "2",
      EducationLevel == "4th grade or less" ~ "3.5", # what number for this?
      EducationLevel == "Fifth grade" ~ "5",
      EducationLevel == "Completed 6th grade" ~ "6",
      EducationLevel == "Seventh grade" ~ "7",
      EducationLevel == "Completed 8th grade" ~ "8",
      EducationLevel == "Ninth grade" ~ "9",
      EducationLevel == "Some high school" ~ "10.5", # what number for this?
      EducationLevel == "Completed high school" ~ "12",
      EducationLevel == "Some university" ~ "13.5",
      EducationLevel == "Professional certificate" ~ "14",
      EducationLevel == "Completed university" ~ "16",
      EducationLevel == "Post-graduate degree" ~ "17", # what number here? 24 doctorate; 19 masters and professional degree
      .default = EducationLevel),
    EducationLevelYears = as.numeric(EducationLevelYears),
    Punishment2 = ifelse(Punishment == "Maybe", 0, Punishment),
    Family2 = ifelse(Family == "Maybe", 0, Family),
    OutsideFamily2 = ifelse(OutsideFamily == "Maybe", 0, OutsideFamily),
    DeathORIllnessInjuryHarmInOthers2 = ifelse(DeathORIllnessInjuryHarmInOthers == "Maybe", 0, DeathORIllnessInjuryHarmInOthers),
    HouseholdAdversity2 = ifelse(HouseholdAdversity == "Maybe", 0, HouseholdAdversity),
    ExplicitInvestmentDesired2 = ifelse(ExplicitInvestmentDesired == "Maybe", 0, ExplicitInvestmentDesired),
    LossOfPrivlegesOrItem2 = ifelse(LossOfPrivlegesOrItem == "Maybe", 0, LossOfPrivlegesOrItem),
    SeparationAttentionSeeking2 = ifelse(SeparationAttentionSeeking == "Maybe", 0, SeparationAttentionSeeking),
    PartnerStatus = as.factor(if_else(as.factor(CaregiverMarital_1) == "1" | as.factor(CaregiverMarital_2) == "1", "Partnered", "Unpartnered", missing = "Unpartnered")) # CHECK
  ) |>
  rowwise() |>
  mutate( # average could be better (does not let missing data bias the result as much)
    FoodSecurity = sum2(FoodSecurity_1, FoodSecurity_2, FoodSecurity_3, FoodSecurity_4, FoodSecurity_5, FoodSecurity_6, na.rm = TRUE),
    FoodSecurityMean = mean2(FoodSecurity_1, FoodSecurity_2, FoodSecurity_3, FoodSecurity_4, FoodSecurity_5, FoodSecurity_6, na.rm = TRUE),
    SocialSupportSum = sum2(SocialSupport_1, SocialSupport_2, SocialSupport_3, SocialSupport_4, SocialSupport_5, SocialSupport_6, SocialSupport_7, SocialSupport_8, SocialSupport_9, SocialSupport_10, SocialSupport_11, SocialSupport_12, na.rm = TRUE),
    SocialSupportMean = mean2(SocialSupport_1, SocialSupport_2, SocialSupport_3, SocialSupport_4, SocialSupport_5, SocialSupport_6, SocialSupport_7, SocialSupport_8, SocialSupport_9, SocialSupport_10, SocialSupport_11, SocialSupport_12, na.rm = TRUE),
  ) |>
  ungroup() |>
  rename(CSSVisits = SocialSupport_1, # C prefix = current support
        CSSHouseHelp = SocialSupport_2,
        PSSEmergencyMoney = SocialSupport_3, # P prefix = potential support
        CSSPeopleCare = SocialSupport_4,
        CSSLoveAffection = SocialSupport_5,
        PSSChildCare = SocialSupport_6,
        CSSConfidantsWork = SocialSupport_7,
        CSSConfidantsFamilySelf = SocialSupport_8,
        CSSConfidantsMoney = SocialSupport_9,
        CSSSocialInvitations = SocialSupport_10,
        CSSAdvice = SocialSupport_11,
        CSSSickHelp = SocialSupport_12
  ) |>
  mutate( # CHECK
    OtherChildAlloparentingFreqN = sum(AlloparentingFreqN, na.rm = TRUE) - AlloparentingFreqN,
    HHChildCareA = as.numeric(AdultsChildcare[1]) * 30 + OtherChildAlloparentingFreqN,
    HHChildCareB = as.numeric(AdultsChildcare[1]) * 60 + OtherChildAlloparentingFreqN,
    OtherChildAlloparents = sum(ChildAlloparent, na.rm = TRUE) - ChildAlloparent,
    HHOtherKidsAge0_2.5 = sum(Age0_2.5, na.rm = TRUE) - Age0_2.5, # Helfrecht and Meehan 2015 style risk periods
    HHOtherKidsAge2.5_5 = sum(Age2.5_5, na.rm = TRUE) - Age2.5_5,
    HHOtherKidsAge5_10 = sum(Age5_10, na.rm = TRUE) - Age5_10,
    HHOtherKidsAge10Plus = sum(Age10Plus, na.rm = TRUE) - Age10Plus,
    HHOtherKidsAge5Minus = sum(Age5Minus, na.rm = TRUE) - Age5Minus,
    .by = householdID
  ) |>
  mutate(
    OtherCare = as.numeric(AdultsChildcare) + OtherChildAlloparents #Check Change to TotalChildCare
  )

# sample characteristics --------------------------------------------------

# scatterplot of CaregiverAge and ChildAge
# note 20 year old caregiver reporting on a 20 year old has a CaregiverRelation of "Other"
ggplot(data = d2, aes(CaregiverAge, ChildAge)) + geom_count() + geom_smooth() + labs(title = "Scatterplot of CaregiverAge and ChildAge")

#bar plot of household income
ggplot(caregivers, aes(y=IncomeCategory)) +
  geom_bar() + labs(y = "IncomeCaregory (household level)")

# bar plot of ChildAge with Sex as color
ggplot(d2, aes(y=ChildAge, fill = Sex)) +
  geom_bar()

# dataframe used to plot sample characteristics using the sample from our main models
modeldf <- d2 |>
  dplyr::select(
    householdID,
    childHHid,
    Sex,
    SignalFreq,
    SignalCost,
    SignalFreqMax,
    ChildAge,
    OtherChildrenHH,
    LogIncome,
    number_adults,
    ConflictFreqN,
    AlloparentingFreqN,
    EducationLevelYears,
    CurrentHealthMean,
    ImmigrateUtila,
    UserLanguage,
    SadFreqOF,
  ) |>
  na.omit()

# bar plot of ChildAge with Sex as color, but with the sample used for our main models
ggplot(modeldf, aes(y=ChildAge, fill = Sex)) +
  geom_bar()

# # bar plot of Neighborhood with ImmigrateUtila as color
caregivers2 <- caregivers |>  filter(!is.na(Neighborhood))
ggplot(caregivers2, aes(y=Neighborhood, fill=ImmigrateUtila)) +
  geom_bar() + labs(title = "Immigration status by neighborhood")



signalfreqdf <- d2 |>
  dplyr::select(
    householdID,
    childHHid,
    uniqueID,
    # Sex,
    SadFreqOF,
    CryFreqOF,
    TantrumFreqOF
    # ConflictFreqOF,
    # AlloparentingFreq0
  ) |>
  rename(
    Sadness = SadFreqOF,
    Crying = CryFreqOF,
    Tantrums = TantrumFreqOF
    #Conflict = ConflictFreqOF,
    #Alloparenting = AlloparentingFreq0
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
    .by = c(Sadness, Crying, Tantrums)
  ) |>
  arrange(desc(Freq)) |>
  dplyr::filter(Freq > 2)


signal_alluvial_plot <- ggplot(sfdfsum2, aes(axis1 = Sadness, axis2 = Crying, axis3 = Tantrums, y = Freq)) +
  geom_alluvium(aes(fill = Sadness), color = "black", show.legend = F) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Sadness", "Crying", "Tantrums"), expand = c(.2, .05)) +
  ylab("Number of\nchildren") +
  theme_minimal(20) +
  theme(axis.title.y = element_text(angle = 0, hjust = 1))
signal_alluvial_plot

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

barplot_SignalFreq <- ggplot(signaldf_long_sum, aes(N, Signal, fill = Freq)) + geom_col(position = "stack")
barplot_SignalFreq

ggplot(signaldf_long_sum, aes(Freq, N, color = Signal, group = Signal)) + geom_line()

# alluvial plot
x <- table(signalfreqdf$Sadness, signalfreqdf$Crying) |>
  as.data.frame() |>
  dplyr::filter(Freq > 3)

ggplot(x, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1), color = "black") +
  geom_stratum() +
  scale_x_discrete(limits = c("Sadness", "Crying"), expand = c(.2, .05))

d2_conflict_filter <- d2 |>
  filter_at(vars(ConflictFreqOF,Sex),all_vars(!is.na(.)))

barplot_conflict <- ggplot(d2_conflict_filter, aes(y=ConflictFreqOF, fill= Sex)) + geom_bar(position = "dodge") + labs(title = "Frequency of conflict") + theme(axis.title.y = element_blank())
# barplot_conflict <- ggplot(d2_conflict_filter, aes(y=ConflictFreqOF, fill= Sex)) + geom_bar(position = "dodge") + labs(title = "Frequency of conflict") + theme(axis.title.y = element_blank(), legend.position = "none")
barplot_conflict

d2_alloparent_filter <- d2 |>
  filter_at(vars(AlloparentingFreq0,Sex),all_vars(!is.na(.)))

barplot_alloparenting <- ggplot(d2_alloparent_filter, aes(y=AlloparentingFreq0, fill= Sex)) + geom_bar(position = "dodge") + labs(title = "Frequency of child alloparenting") + theme(axis.title.y = element_blank())
# barplot_alloparenting <- ggplot(d2_alloparent_filter, aes(y=AlloparentingFreq0, fill= Sex)) + geom_bar(position = "dodge") + labs(title = "Frequency of child alloparenting") + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
barplot_alloparenting

# library(gridExtra)
# grid.arrange(barplot_conflict, barplot_alloparenting, ncol = 2)
# signaling frequency and cost histograms ---------------------------------

ggplot(d2, aes(x=SadFreqN, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=CryFreqN, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=TantrumFreqN, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=SignalFreq, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=SignalCost, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=SignalFreqMax, fill=Sex)) +
  geom_histogram()

barplot_sad <- ggplot(modeldf, aes(y=SadFreqOF, fill=Sex, na.rm = TRUE)) + geom_bar()
barplot_sad

#' # Correlation matrix of predictor and outcome variables

#+ message=F, warning=F, fig.height=9, fig.width=10

maincormat <- d2[c("CurrentHealthMean", "EducationLevelYears", "OtherChildrenHH", "LogIncome", "ConflictFreqN", "number_adults", "PartnerStatus", "AlloparentingFreqN", "Sex", "ChildAge", "RelativeNeed3", "RelativeMaternalInvestment2", "SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalFreqMax", "SignalCost")] |>
  mutate(
    PartnerStatus = as.numeric(PartnerStatus),
    Sex = as.numeric(Sex),
    RelativeNeed3 = as.numeric(RelativeNeed3),
    RelativeMaternalInvestment2 = as.numeric(RelativeMaternalInvestment2)
  ) |>
  cor( use = "pairwise.complete.obs")

ggcorrplot(maincormat, hc.order = TRUE, lab = TRUE)


# CurrentHealthMean models ------------------------------------------------

# signal cost
mscH <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mscH)
plot(allEffects(mscH))

# signal frequency

msfH <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(msfH)
plot(allEffects(msfH))

# frequency of most common signal (running away excluded)

mmsfH <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mmsfH)
plot(allEffects(mmsfH))

# Signal specific models

# Current health only a significant predictor of sad frequency
# little in the way of conflicts on interest with illness as long as it does not involve low survival probability at the youngest of ages?
# in line with this reasoning, CurrentHealthMean does not predict conflict

msadH <- glmmTMB(SadFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(msadH)
plot(allEffects(msadH))

mcryH <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mcryH)
plot(allEffects(mcryH))

mtantrumH <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mtantrumH)
plot(allEffects(mtantrumH))

# anthropometric models ---------------------------------------------------

# HeightZ, WeightZ, BMIZ, and BodyFatPercentage are not significant predictors
# (not close either) of signaling behavior in our standard models

# GripStrengthMean and GripR are not either (as main effects or interacted with Sex)

# Tricep residuals
mtriC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + TricepR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mtriC)
plot(allEffects(mtriC))

  # Tricep residuals = not a signficant predictor of signal frequency.
  # Effect on cost may be due to its association with tantrum frequency.
  # (Tantrum model below)

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

msubTantrum <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + SubscapR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msubTantrum)
plot(allEffects(msubTantrum))

# BodyFat
d2$BodyFat <- c(scale(d2$TricepR) + scale(d2$SubscapR))
mBFc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFc)
plot(allEffects(mBFc))

  # BodyFat scaled composite = not a signficant predictor of signal frequency.
  # Effect on cost may be due to its association with tantrum frequency.
  # (Tantrum model below)

mBFtantrum <- glmmTMB(TantrumFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFtantrum)
plot(allEffects(mBFtantrum))

mtest <- glmmTMB(HeightZ ~ FoodSecurityMean + Sex + (1|householdID), data = d2, family = gaussian())
summary(mtest)
plot(allEffects(mtest))

mtest <- glmmTMB(WeightZ ~ FoodSecurityMean + Sex + (1|householdID), data = d2, family = gaussian())
summary(mtest)
plot(allEffects(mtest))

# ConflictFreqN model -----------------------------------------------------------

# original
mconflict <- glmmTMB(ConflictFreqN ~ ChildAge+ Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(mconflict)
plot(allEffects(mconflict))

# with HHChildCareB*ChildAge
mconflict <- glmmTMB(ConflictFreqN ~ ChildAge+ Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + HHChildCareB*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mconflict)
plot(allEffects(mconflict))

mconflict <- glmmTMB(ConflictFreqN ~ ChildAge+ Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + HHChildCareA*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mconflict)
plot(allEffects(mconflict))

# number of alloparents not a significant predictor
# (as main effect or in an HHChildCareB*ChildAge interaction)

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

d2$AdultsNoChildcare <- d2$number_adults - d2$AdultsChildcare
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

# care (adult's heavily weighted)
malloparenting <- glmmTMB(AlloparentingFreqN ~ ChildAge + Sex + YoungerKids + OlderKids + LogIncome + AdultsNoChildcare + PartnerStatus + EducationLevelYears + OtherChildAlloparentingFreqN + AdultsChildcare + (1|householdID),data = d2, family = nbinom2)
summary(malloparenting)
plot(allEffects(malloparenting))

# number of others providing care
malloparenting2 <- glmmTMB(AlloparentingFreqN ~ ChildAge + Sex + YoungerKids + OlderKids + LogIncome + number_adults + PartnerStatus + EducationLevelYears + OtherCare + (1|householdID),data = d2, family = nbinom2)
summary(malloparenting2)
plot(allEffects(malloparenting2))

# relatedness within the family -------------------------------------------

# Child relatedness (relatedness relative to every kid in household) = not a significant
# predictor

# Sibling relatedness

# SignalCost and SignalFreqMax = marginally significant

# SignalFreq model
# sibling relatedness no longer significant with BodyFat
mSRF <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + (1|householdID),data = d2, family = nbinom2)
summary(mSRF)
plot(allEffects(mSRF))

# with YoungerKids*ChildAge
mSRF3 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OlderKids + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + YoungerKids*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mSRF3)
plot(allEffects(mSRF3))

# ordinal logistic regression of parental response ------------------------

# note 6 instances of both positive and negative coded as NA
table(d2$PositiveResponse, d2$NegativeResponse)

mp7 <- polr(ParentalResponse ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + UnwantedTask + Punishment2 + Family2 + DiscomfortPainInjuryIllness, d2)
summary(mp7)

coeftest(mp7, vcov = vcovCL, type = "HC0", cluster = ~householdID)
plot_predictions(mp7, condition = c("Punishment2", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mp7, condition = c("DiscomfortPainInjuryIllness", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mp7, condition = c("Family2", "group"), type = "probs", vcov = ~householdID)


# ordinal logistic regression models for relative need and materna --------

# Results with NeighborhoodF as a random effect and vc0 seem like outliers compared to all other methods

# Using preferred model for signal frequency (without partner status)
# signal Freq and Child age still significant after controlling for education years
# CurrentHealthMean has no predictive ability when added
mN1 <- polr(RelativeNeed3 ~ SignalFreq + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults+ ConflictFreqN + AlloparentingFreqN, d2)
summary(mN1)
coeftest(mN1, vcov = vcovCL, type = "HC0", cluster = ~householdID)
plot_predictions(mN1, condition = c("ChildAge", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mN1, condition = c("SignalFreq", "group"), type = "probs", vcov = ~householdID)

(ci <- confint(mN1))
exp(cbind(coef(mN1),(ci)))

mN2 <- polr(RelativeNeed3 ~ SignalCost + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults+ ConflictFreqN + AlloparentingFreqN, d2)
summary(mN2)
coeftest(mN2, vcov = vcovCL, type = "HC0", cluster = ~householdID)
plot_predictions(mN2, condition = c("ChildAge", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mN2, condition = c("SignalCost", "group"), type = "probs", vcov = ~householdID)

(ci <- confint(mN1))
exp(cbind(coef(mN1),(ci)))

# brm
# mN1b <- brm(RelativeNeed3 ~ SignalFreq + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN + CurrentHealthMean + (1|householdID), d2, family = cumulative(link = "logit", threshold = "flexible"))
# save(mN1b, file = "mN1b.rda")
load(file = "mN1b.rda")
summary(mN1b)
plot_predictions(mN1b, condition = c("ChildAge", "group"), type = "response")

# Using preferred model for signal frequency (with partner status)
mI <- polr(RelativeMaternalInvestment2 ~ SignalFreq + ChildAge*Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + CurrentHealthMean, d2)
summary(mI)
coeftest(mI, vcov = vcovCL, type = "HC0", cluster = ~householdID)

(ci <- confint(mI))
exp(cbind(coef(mI),(ci)))

plot_predictions(mI, condition = c("OtherChildrenHH", "group"), type = "probs")
plot_predictions(mI, condition = c("ChildAge", "group"), type = "probs")



# Notes:

# Survey language as a way to get at differences in culture
# - Spanish: 322
# - English: 53
# - Nothing close to significant in our existing models.

