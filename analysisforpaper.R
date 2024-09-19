# loading libraries -------------------------------------------------------

# note: list is cleaned up but some may no longer used
library(utiladata2023) # data package
library(tidyverse)
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

children2 <- children |>
  dplyr::filter(ChildID != "") |>
  dplyr::select(ChildID, householdID)

anthropometricMeans0 <- anthropometrics |>
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
    BMI = WeightMean_KG / (HeightMean / 100 )^2,
    Year = year(Date.of.measurement)
  ) |>
  relocate(measurements, .after = Year)

anthropometricMeans <- anthropometricMeans0 |>
  dplyr::filter(
    measurements == 1 | year(Date.of.measurement) == 2023
  ) |>
  left_join(children2) # note that we filter out householdID based on it's column number

# moved below
# anthropometricMeansWide <- anthropometricMeans0 |>
#   dplyr::filter(
#     measurements == 2
#   ) |>
#   pivot_wider(names_from = "Year", values_from = BicepMean:BMI, id_cols = ChildID) |>
#   left_join(d2, by = "ChildID")

anthropometricMeans$WeightZ <- growthRef(AgeAtMeasurement, WeightMean_KG, Sex, anthropometricMeans, type = "Weight", pop = "CDC")
anthropometricMeans$HeightZ <- growthRef(AgeAtMeasurement, HeightMean, Sex, anthropometricMeans, type = "Height", pop = "CDC")
anthropometricMeans$BMIZ <- growthRef(AgeAtMeasurement, BMI, Sex, anthropometricMeans, type = "BMI", pop = "CDC")
anthropometricMeans$Sex2 <- factor(anthropometricMeans$Sex)

m <- mgcv::gam(GripStrengthMean ~ Sex+s(AgeAtMeasurement, by = Sex2, k = 4), data = anthropometricMeans)
anthropometricMeans$GripR <- residuals(m)

m2 <- mgcv::gam(TricepMean ~ Sex+s(AgeAtMeasurement, by = Sex2, k = 3), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$TricepR <- residuals(m2)

m3 <- mgcv::gam(SubscapMean ~ Sex+s(AgeAtMeasurement, by = Sex2, k = 3), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$SubscapR <- residuals(m3)

m4 <- mgcv::gam(BodyFatPercentageMean ~ Sex+s(AgeAtMeasurement, by = Sex2, k = 3), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$BodyFatPercentageR <- residuals(m4)

m5 <- mgcv::gam(FlexedMean ~ Sex+s(AgeAtMeasurement, by = Sex2, k = 3), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$FlexedR <- residuals(m5)

m5b <- mgcv::gam(FlexedMean ~ Sex+s(AgeAtMeasurement, by = Sex2, k = 3) + s(householdID, bs = "re"), data = anthropometricMeans, na.action = na.exclude)
anthropometricMeans$FlexedRb <- residuals(m5b)

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
      is.na(AlloparentingFreq) & ChildAge < 5 ~ "Never", #does nothing
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
      #RelativeNeed == "The same amount" ~ "The same",
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
      levels = c("Less", "The same amount", "More")), # change back to "The same" here
    IllnessSusceptibility1 = invert(Currenthealth1),
    IllnessSusceptibility2 = invert(Currenthealth2),
    IllnessSusceptibility3 = Currenthealth3 + 0, # to remove labels
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
    CaregiverResponse = case_when(
      PositiveResponse > 0 & NegativeResponse == 0 ~ "1",
      PositiveResponse == 0 & NegativeResponse == 0 ~ "0",
      PositiveResponse == 0 & NegativeResponse > 0 ~ "-1"
    ),
    CaregiverResponse = factor(CaregiverResponse, ordered = TRUE, levels = c("-1", "0", "1"))
  ) |>
  rowwise() |>
  mutate(
    SignalFreq = sum2(SadFreqN, CryFreqN, TantrumFreqN, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    SignalFreqMax = max2(SadFreqN, CryFreqN, TantrumFreqN, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    SignalCost = sum2(SadFreqN, CryFreqN * 2, TantrumFreqN * 3, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    IllnessSusceptibilitySum = sum2(IllnessSusceptibility1, IllnessSusceptibility2, IllnessSusceptibility3, na.rm = TRUE),
    IllnessSusceptibilityMean = mean2(IllnessSusceptibility1, IllnessSusceptibility2, IllnessSusceptibility3, na.rm = TRUE)
  ) |>
  ungroup() |>
  left_join(anthropometricMeans[-21], by = "ChildID") |>  # removing householdID column
  left_join(dplyr::select(caregivers, householdID, CPRatio, Neighborhood, ImmigrateUtila, IncomeCategory, IncomeCategoryN,
                             EducationLevel, EducationLevelYears, AdultsMoney, AdultsHousework, AdultsChildcare, AdultsNoChildcare, CaregiverAge,
                             number_children2, number_adults, UserLanguage, contains("CaregiverMarital"),
                             contains("SocialSupport"), contains("FoodSecurity")), by = "householdID") |>
  left_join(dplyr::select(causes, -householdID, -childHHid, -ChildID, -Sad1), by = "uniqueID")


# # add household data to d2
# d <- left_join(d, dplyr::select(caregivers, householdID, CPRatio, Neighborhood, ImmigrateUtila, IncomeCategory, IncomeCategoryN,
#                                 EducationLevel, EducationLevelYears, AdultsMoney, AdultsHousework, AdultsChildcare, AdultsNoChildcare, CaregiverAge,
#                                 number_children2, number_adults, UserLanguage, contains("CaregiverMarital"),
#                                 contains("SocialSupport"), contains("FoodSecurity")), by = "householdID") |>
#   left_join(dplyr::select(causes, -householdID, -childHHid, -ChildID, -Sad1), by = "uniqueID")

# prepare household variables for analyses
d2 <- d |>
  mutate(
    # EducationLevel = case_when(
    #   EducationLevel == "Fourth grade" ~ "4th grade or less",
    #   EducationLevel == "Tenth grade" ~ "Some high school",
    #   EducationLevel == "Eleventh grade" ~ "Some high school",
    #   .default = EducationLevel),
    # IncomeCategory = ifelse(IncomeCategory == "Prefer not to answer", NA, IncomeCategory),
    # IncomeCategoryN = case_when(
    #   IncomeCategory == "Less than 5,000 Lempira per month" ~ 2500,
    #   IncomeCategory == "5,001 to 10,000 Lempira per month" ~ 7500,
    #   IncomeCategory == "10,001 to 20,000 Lempira per month" ~ 15000,
    #   IncomeCategory == "20,001 to 40,000 Lempira per month" ~ 30000,
    #   IncomeCategory == "40,001 to 80,000 Lempira per month" ~ 60000,
    #   IncomeCategory == "Over 80,000 Lempira per month" ~ 120000
    #   ),
    LogIncome = log10(IncomeCategoryN),
    BodyFat = c(scale(TricepR) + scale(SubscapR)),
    RecentSignalSeverity = com1[as.character(Sad1)],
    RecentSignalCause = com2[as.character(Sad1)],
    AdultsHousework = as.numeric(AdultsHousework),
    #AdultsChildcare = as.numeric(AdultsChildcare),
    Neighborhood2 = as.numeric(Neighborhood == "Camponado/Campolancho"),
    NeighborhoodF = as.factor(Neighborhood),
    OtherChildrenHH = number_children2 - 1,
    # IncomeCategory = factor(
    #   IncomeCategory,
    #   ordered = TRUE,
    #   levels = c("Less than 5,000 Lempira per month", "5,001 to 10,000 Lempira per month", "10,001 to 20,000 Lempira per month", "20,001 to 40,000 Lempira per month", "40,001 to 80,000 Lempira per month", "Over 80,000 Lempira per month")),
    # IncomeCategory2 = as.numeric(IncomeCategory),
    # EducationLevel = factor(
    #   EducationLevel,
    #   ordered = TRUE,
    #   levels = c("No education", "First grade", "Second grade", "4th grade or less", "Fifth grade", "Completed 6th grade", "Seventh grade", "Completed 8th grade", "Ninth grade", "Some high school", "Completed high school", "Some university", "Professional certificate", "Completed university", "Post-graduate degree")),
    # EducationLevelYears = case_when( #
    #   EducationLevel == "No education" ~ "0",
    #   EducationLevel == "First grade" ~ "1",
    #   EducationLevel == "Second grade" ~ "2",
    #   EducationLevel == "4th grade or less" ~ "3.5", # what number for this?
    #   EducationLevel == "Fifth grade" ~ "5",
    #   EducationLevel == "Completed 6th grade" ~ "6",
    #   EducationLevel == "Seventh grade" ~ "7",
    #   EducationLevel == "Completed 8th grade" ~ "8",
    #   EducationLevel == "Ninth grade" ~ "9",
    #   EducationLevel == "Some high school" ~ "10.5", # what number for this?
    #   EducationLevel == "Completed high school" ~ "12",
    #   EducationLevel == "Some university" ~ "13.5",
    #   EducationLevel == "Professional certificate" ~ "14",
    #   EducationLevel == "Completed university" ~ "16",
    #   EducationLevel == "Post-graduate degree" ~ "17", # what number here? 24 doctorate; 19 masters and professional degree
    #   .default = EducationLevel),
    # EducationLevelYears = as.numeric(EducationLevelYears),
    Punishment2 = ifelse(Punishment == "Maybe", 0, Punishment),
    Family2 = ifelse(Family == "Maybe", 0, Family),
    OutsideFamily2 = ifelse(OutsideFamily == "Maybe", 0, OutsideFamily),
    WantNeed2 = ifelse(WantNeed == "Want" | WantNeed == "Need", WantNeed, NA),
    WantNeedBinary = case_when(
      WantNeed2 == "Want" ~ 0,
      WantNeed2 == "Need" ~ 1),
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

# dataframe for plots -----------------------------------------------------

# filters based on having at least one value of SadFreqN, CryFreqN, or TantrumFreqN
modeldf <- d2 |>
  dplyr::select(
    householdID,
    childHHid,
    Sex,
    CryFreqN,
    SadFreqN,
    TantrumFreqN,
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
    IllnessSusceptibility1,
    IllnessSusceptibility2,
    IllnessSusceptibility3,
    IllnessSusceptibilityMean,
    SadFreqOF,
    CryFreqOF,
    TantrumFreqOF,
    ConflictFreqOF,
    AlloparentingFreq0,
    NumberOfChildren,
    RelativeNeed,
    RelativeNeed3,
    PositiveResponse,
    NegativeResponse,
    CaregiverAge,
    IncomeCategoryN,
    SubscapR,
    HeightZ,
    WeightZ,
    BMIZ,
    SubscapMean,
    SubscapR,
    TricepMean,
    TricepR,
    BodyFat,
    FlexedR

  ) |>
  dplyr::filter(!(is.na(CryFreqN) & is.na(SadFreqN) & is.na(TantrumFreqN)))

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
    TantrumFreqOF
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
  guides(fill = guide_legend("Frequency")) +
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

maincormat <- d2[c("IllnessSusceptibilityMean", "EducationLevelYears", "OtherChildrenHH", "LogIncome", "ConflictFreqN", "number_adults", "PartnerStatus", "AlloparentingFreqN", "Sex", "ChildAge", "RelativeNeed3", "RelativeMaternalInvestment2", "SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalFreqMax", "SignalCost")] |>
  mutate(
    PartnerStatus = as.numeric(PartnerStatus),
    Sex = as.numeric(Sex),
    RelativeNeed3 = as.numeric(RelativeNeed3),
    RelativeMaternalInvestment2 = as.numeric(RelativeMaternalInvestment2)
  ) |>
  cor( use = "pairwise.complete.obs")

ggcorrplot(maincormat, hc.order = TRUE, lab = TRUE)

# main models ------------------------------------------------

# signal cost
mscH <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + IllnessSusceptibilityMean + (1|householdID),data = d2, family = nbinom2)
summary(mscH)
plot(allEffects(mscH))

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

# Flexed residuals

mFlexC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2[d2$FlexedR < 10,], family = nbinom2)
#mFlexC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2[d2$FlexedR < 10,], family = nbinom2)
summary(mFlexC)
plot(allEffects(mFlexC))

mFlexCb <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedRb + (1|householdID),data = d2[d2$FlexedR < 10,], family = nbinom2)
summary(mFlexCb)
plot(allEffects(mFlexCb))

mFlexF <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2[d2$FlexedR < 10,], family = nbinom2)
#mFlexF <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + FlexedR + (1|householdID),data = d2, family = nbinom2)
summary(mFlexF)
plot(allEffects(mFlexF))

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

mBFc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFc)
plot(allEffects(mBFc))

mBFsad <- glmmTMB(SadFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFsad)
plot(allEffects(mBFsad))

mBFcry <- glmmTMB(CryFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mBFcry)
plot(allEffects(mBFcry))

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

anthropometricMeansWide <- anthropometricMeans0 |>
  dplyr::filter(
    measurements == 2
  ) |>
  pivot_wider(names_from = "Year", values_from = BicepMean:BMI, id_cols = ChildID) |>
  left_join(d2, by = "ChildID")

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

weight_long_cry <- glmmTMB(WeightMean_KG_2024 ~ WeightMean_KG_2023 + CryFreqN*ChildAge + (1|householdID), data = anthropometricMeansWide, family = gaussian)
summary(weight_long_cry)
plot(allEffects(weight_long_cry))

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

height_long_cry <- glmmTMB(HeightMean_2024 ~ HeightMean_2023 + CryFreqN*ChildAge + (1|householdID), data = anthropometricMeansWide, family = gaussian)
summary(height_long_cry)
plot(allEffects(height_long_cry))

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

d2$AdultsNoChildcare <- d2$number_adults - d2$AdultsChildcare

mconflict <- glmmTMB(ConflictFreqN ~ ChildAge + Sex + OlderKids + YoungerKids + LogIncome + AdultsNoChildcare + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + AdultsChildcare + OtherChildAlloparentingFreqN + (1|householdID),data = d2, family = nbinom2)
summary(mconflict)
plot(allEffects(mconflict))

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

d3 <- d2 |>
  dplyr::filter(YoungerKids > 1, OtherChildAlloparentingFreqN < 121)

# model cannot handle high alloparenting family

# care (adult's heavily weighted)
malloparenting <- glmmTMB(AlloparentingFreqN ~ ChildAge + Sex + YoungerKids + OlderKids + LogIncome + AdultsChildcare + OtherChildAlloparentingFreqN + (1|householdID),data = d3, family = nbinom2)
summary(malloparenting)
plot(allEffects(malloparenting))

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
mN1 <- polr(RelativeNeed3 ~ SignalFreq + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults+ ConflictFreqN + AlloparentingFreqN, d2)
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

signaling_plot_need <- (p_need_age_sf) +
  (p_need_signalfreq + theme(legend.position = "none", legend.title = element_blank())) +
  plot_layout(ncol = 1, byrow = FALSE, axes = "collect_y") +
  plot_annotation(title = "", tag_levels = "A") & theme(plot.tag.position = c(.98, .92))

signaling_plot_need

mN2 <- polr(RelativeNeed3 ~ SignalCost + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults+ ConflictFreqN + AlloparentingFreqN, d2)
summary(mN2)
coeftest(mN2, vcov = vcovCL, type = "HC0", cluster = ~householdID)

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

# converted to characters for the purpose of avg_predictions function
d2$UnwantedTask <- as.character(d2$UnwantedTask)
d2$DiscomfortPainInjuryIllness <- as.character(d2$DiscomfortPainInjuryIllness)
d2$Family2 <- as.numeric(d2$Family2)

mpr <- polr(CaregiverResponse ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + UnwantedTask + Punishment2 + Family2 + DiscomfortPainInjuryIllness, d2)
summary(mpr)
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


coeftest(mpr, vcov = vcovCL, type = "HC0", cluster = ~householdID)
p_response_punishment <- plot_predictions(mpr, condition = c("Punishment2", "group"), type = "probs", vcov = ~householdID) +
  xlab("Cause involved punishment") +
  ylab("Caregiver response") +
  scale_color_viridis_d(option = "B", end = .8) +
  guides(color = guide_legend(override.aes = list(linewidth=1))) +
  theme_minimal() + # theme_minimal needs to come before plot.tag.position
  theme(legend.position = "top", plot.tag.position = "right") +
  ylim(0,1)

p_response_pain <- plot_predictions(mpr, condition = c("DiscomfortPainInjuryIllness", "group"), type = "probs", vcov = ~householdID) +
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

# relative need predicts relative investment ------------------------------

mpr <- polr(RelativeMaternalInvestment2 ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + RelativeNeed3, d2)
summary(mpr)
coeftest(mpr, vcov = vcovCL, type = "HC0", cluster = ~householdID)

mpr_plot_data <- plot_predictions(mpr, condition = c("RelativeNeed3", "group"), type = "probs", vcov = ~householdID, draw = FALSE)

ggplot(mpr_plot_data, aes(RelativeNeed3, estimate, color = group, group = group, ymin = conf.low, ymax = conf.high)) +
  geom_line(position = position_dodge(width = .2), linewidth = 1.5) +
  geom_pointrange(position = position_dodge(width = .2)) +
  theme_minimal(15) +
  xlab("\nRelative need") +
  ylab("Proportion of children") +
  #scico::scale_color_scico_d(palette = "lipari", begin = .3, end = .8) +
  scale_color_viridis_d(option = "B", begin = .25, end = .85) +
  guides(color = guide_legend(title = "Relative investment:")) &
  theme(legend.position = "top")

plot_predictions(mpr, condition = c("OtherChildrenHH", "group"), type = "probs", vcov = ~householdID)

mpr2 <- polr(RelativeMaternalInvestment2 ~ ChildAge + Sex + YoungerKids + OlderKids + LogIncome + number_adults + RelativeNeed3, d2)
summary(mpr2)
coeftest(mpr2, vcov = vcovCL, type = "HC0", cluster = ~householdID)

mpr_plot_data2 <- plot_predictions(mpr2, condition = c("RelativeNeed3", "group"), type = "probs", vcov = ~householdID, draw = FALSE)

# NOTE error ribbons dip below 0
plot_predictions(mpr2, condition = c("YoungerKids", "group"), type = "probs", vcov = ~householdID)

plot_predictions(mpr2, condition = c("YoungerKids", "group"), type = "probs", vcov = ~householdID) +
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

p_age_sad <- plot_age(msadH, "Sadness freq.")
p_age_cry <- plot_age(mcryH, "Crying freq.")
p_age_tantrum <- plot_age(msadH, "Tantrum freq.")
p_age_freq <- plot_age(msadH, "Summed freq.")
p_age_cost <- plot_age(msadH, "Signal cost")
p_age_freq_max <- plot_age(msadH, "Max freq.")

signaling_data_plot_age <- p_age_sad + p_age_cry + p_age_tantrum + p_age_freq + p_age_freq_max + p_age_cost +
  plot_layout(axes = "collect_x", ncol = 2, byrow = FALSE, guides = "collect") &
  # plot_annotation(title = "Signaling measures", tag_levels = "A") & theme(plot.tag.position = c(.98, .92), legend.position = "top")
  # plot_annotation(title = "Signaling measures") & theme(legend.position = "top")
  theme(legend.position = "top")

ggsave("Figures/signaling_data_plot_age.pdf", signaling_data_plot_age, width = 9, height = 9)

signaling_data_plot_age

# model for rugs for plots without functions
# geom_rug(data = msadH$frame, aes(x = ChildAge , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))



# p_age_sad <- plot_predictions(msadH, condition = c("ChildAge", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Sadness freq.") +
#   xlab("Child age (years)") +
#   xlim(5,20) +
#   ylim(0,NA) +
#   scale_color_binary() +
#   scale_fill_binary() +
#   geom_rug(data = msadH$frame, aes(x = ChildAge , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))
# p_age_sad$layers[[2]]$aes_params$linewidth <- 2
#
# p_age_cry <- plot_predictions(mcryH, condition = "ChildAge", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Crying freq.") +
#   xlab("Child age (years)") +
#   xlim(5,20) +
#   ylim(0,NA) +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ChildAge , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_age_tantrum <- plot_predictions(mtantrumH, condition = "ChildAge", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Tantrum freq.") +
#   xlab("Child age (years)") +
#   xlim(5,20) +
#   ylim(0,NA) +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ChildAge , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_age_freq <- plot_predictions(msfH, condition = "ChildAge", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Summed freq.") +
#   xlab("Child age (years)") +
#   xlim(5,20) +
#   ylim(0,NA) +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ChildAge , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_age_cost <- plot_predictions(mscH, condition = "ChildAge", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Signal cost") +
#   xlab("Child age (years)") +
#   xlim(5,20) +
#   ylim(0,NA) +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ChildAge , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")
#
# p_age_freq_max <- plot_predictions(mmsfH, condition = "ChildAge", vcov = TRUE, points = 0, type = "link", transform = "exp") +
#   theme_bw(15) +
#   ylab("Max freq.") +
#   xlab("Child age (years)") +
#   xlim(5,20) +
#   ylim(0,NA) +
#   geom_rug(data = drop_na(d2, Sex), aes(x = ChildAge , y = 0, color = Sex), position = position_jitter(width = 1, seed = 123), sides = "b")



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
p_conflict_tantrum <- plot_conflict(msadH, "Tantrum freq.")
p_conflict_freq <- plot_conflict(msadH, "Summed freq.")
p_conflict_cost <- plot_conflict(msadH, "Signal cost")
p_conflict_freq_max <- plot_conflict(msadH, "Max freq.")

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
  plot_predictions(mconflict, condition = c("ChildAge", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_bw(15) +
  ylab("Conflict freq. (per month)") +
  xlab("Child age (years)") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = mconflict$frame, aes(x = ChildAge , y = 0, color = Sex), position = position_jitter(width = .3, seed = 123), sides = "b") +
  theme(legend.position = "none")
  #theme(plot.tag.position = "topright")
p_conflict_dep_var_age$layers[[2]]$aes_params$linewidth <- 2

p_conflict_dep_var_age

p_conflict_dep_var_aCare <-
  plot_predictions(mconflict, condition = c("AdultsChildcare", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_bw(15) +
  ylab("Conflict freq. (per month)") +
  xlab("Number of adults that contribute childcare") +
  scale_color_binary() +
  scale_fill_binary() +
  #theme(plot.tag.position = "topright")
  geom_rug(data = mconflict$frame, aes(x = AdultsChildcare , y = 0, color = Sex), position = position_jitter(width = .3, seed = 123), sides = "b") +
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
p_education_tantrum <- plot_education(msadH, "Tantrum freq.")
p_education_freq <- plot_education(msadH, "Summed freq.")
p_education_cost <- plot_education(msadH, "Signal cost")
p_education_freq_max <- plot_education(msadH, "Max freq.")

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

d_allo_by_sex_sad <- plot_predictions(msadH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_sad <- plot_allo(d_allo_by_sex_sad, "Sadness freq.")

d_allo_by_sex_cry <- plot_predictions(mcryH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_cry <- plot_allo(d_allo_by_sex_cry, "Crying freq.")

d_allo_by_sex_tantrum <- plot_predictions(mtantrumH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_tantrum <- plot_allo(d_allo_by_sex_tantrum, "Tantrum freq.")

d_allo_by_sex_freq <- plot_predictions(msfH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_freq <- plot_allo(d_allo_by_sex_freq, "Summed freq.")

d_allo_by_sex_cost <- plot_predictions(mscH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_cost <- plot_allo(d_allo_by_sex_cost, "Max freq.")

d_allo_by_sex_freq_max <- plot_predictions(mmsfH, condition = c("AlloparentingFreqN", "Sex"), vcov = TRUE, points = pointsize, type = "link", transform = "exp", draw = FALSE)
p_allo_by_sex_freq_max <- plot_allo(d_allo_by_sex_freq_max, "Signal cost")

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

p_bodyfat_cost <- plot_predictions(mBFc, condition = c("BodyFat"), vcov = TRUE, points = .25, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Signal cost") +
  xlab("Body fat composite")

p_bodyfat_cost2 <- plot_predictions(mBFc, condition = c("BodyFat", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp", conf_level = .95) +
  theme_minimal(15) +
  ylab("Signal cost") +
  xlab("Body fat composite") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = mBFc$frame, aes(x = BodyFat , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))
p_bodyfat_cost2$layers[[2]]$aes_params$linewidth <- 2

p_bodyfat_crying <- plot_predictions(mBFcry, condition = "BodyFat", vcov = TRUE, points = .25, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Crying freq.") +
  xlab("Body fat composite")

p_bodyfat_crying2 <- plot_predictions(mBFcry, condition = c("BodyFat", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Crying freq.") +
  xlab("Body fat composite") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = mBFcry$frame, aes(x = BodyFat , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))
p_bodyfat_crying2$layers[[2]]$aes_params$linewidth <- 2

p_bodyfat_tantrum <- plot_predictions(mBFtantrum, condition = "BodyFat", vcov = TRUE, points = .25, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Tantrum freq.") +
  xlab("Body fat composite")

p_bodyfat_tantrum2 <- plot_predictions(mBFtantrum, condition = c("BodyFat", "Sex"), vcov = TRUE, points = 0, type = "link", transform = "exp") +
  theme_minimal(15) +
  ylab("Tantrum freq.") +
  xlab("Body fat composite") +
  scale_color_binary() +
  scale_fill_binary() +
  geom_rug(data = mBFtantrum$frame, aes(x = BodyFat , y = 0, color = Sex), sides = "b", position = position_jitter(width = .3, seed = 123))
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
models <- list (msadH = msadH, mcryH = mcryH, mtantrumH = mtantrumH, msfH = msfH, mscH = mscH, mmsfH = mmsfH, mBFc = mBFc, mBFcry = mBFcry, mBFtantrum = mBFtantrum, mconflict = mconflict, malloparenting = malloparenting)
stats <- map (models, tdy)
stats$mBFc$BodyFat$str
stats$mBFc$BodyFat$str2



