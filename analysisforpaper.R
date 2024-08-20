#' ---
#' title: Child signaling
#' ---
#' # Child signaling
#' **Data**: interviews with caretakers (usually mothers) in 157 families, with data
#' on 375 children (more complete data on ~242)
#'
#' **Principle outcome variable**: child signaling frequency (Likert scale converted to times per month)
#'
#' **Secondary outcome variable**: child signaling cost (sad = low, crying = moderate, tantrums = high)
#' (SadFreq + CryFreq\*2 + TantrumFreq\*3)
#'
#' **Tertiary outcome variable**: parental response to the most recent instance of child signaling
#'
#' **Predictor variables**:
#'
#' * ChildAge
#' * Sex
#' * Other children in the household (OtherChildrenHH)
#' * LogIncome (median value of range for each option on the Likert Scale)
#' * Number of adults, including the caretaker taking the survey (number_adults)
#' * Partner status of caretaker (PartnerStatus: Partnered/Unpartnered)
#' * Frequency of conflict between child and caretaker (Likert scale converted to times per month)
#' * Frequency of the child's alloparenting contributions as a Likert scale converted to times per month (AlloparentingFreqN)
#' * Caretaker's education level in years (EducationLevelYears)
#' * Mean of caretaker's responses to three child health measures (CurrentHealthMean):
#'
#'   1) If an illness is 'going around', my child will get it.
#'
#'   2) My child’s immune system protects him/her from most illnesses that other kids get.
#'
#'   3) In general, my child is very susceptible to colds, flu, and other infectious diseases.
#' * Random effects grouping variable (householdID)
#'
#' **Results**:
#'
#' * Negative association between age and child signaling (all signaling measures)
#' * Positive association between conflict with parent and child signaling (all signaling measures)
#' * Sex and focal child alloparenting effort interact (postive for females, negative/flat for males)
#' * In most models, there is a negative association between education and child signaling?
#' * Negative association between caretaker estimates of child health and child signaling for some measures. This does not extend to the costlier signals.
#'
#' **Null Results (often not presented here)**:
#'
#' * Caretaker's reported social support not associated with child signaling
#' * Very few age by sex interactions for signaling
#' * Older models showed a LogIncome by household interaction. LogIncome is not a significant predictor of child signaling in our current models.
#'
#' **Note**:
#'
#' * Code for analyzing the effects of relatedness within households still needs work.
#'
#' # Load libraries


# loading libraries -------------------------------------------------------

#+ message=F, warning=F, fig.height=9, fig.width=10
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

source("recode.R")

#' # Functions used

# function to invert the scoring on a Likert scale
invert <- function(x){
  ((x - max(x, na.rm = TRUE)) * -1) + min(x, na.rm = TRUE)
  }

# functions to calculate values from multiple columns in ways that 1) NAs are not counted 2)
# yet do not result in NAs if values exist for other columns.
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

#' # Data frame prep

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
      RelativeMaternalInvestment, # do we need Much less if it was never picked?
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
  mutate( # should we not do this for when we convert to numeric?
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
    # AdultsChildcareN1 = AdultsChildcare*30,
    # AdultsChildcareN2 = AdultsChildcare*60,
    # OtherAlloparentingN1 = OtherChildAlloparenting + AdultsChildcareN1, # how to handle NAs?
    # OtherAlloparentingN2 = OtherChildAlloparenting + AdultsChildcareN2,
    #CaregiverMarried = as.factor(CaregiverMarital_1),
    #CaregiverPartneredLongTerm = as.factor(CaregiverMarital_2),
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
    OtherAlloparentingA = as.numeric(AdultsChildcare[1]) * 30 + OtherChildAlloparentingFreqN,
    OtherAlloparentingB = as.numeric(AdultsChildcare[1]) * 60 + OtherChildAlloparentingFreqN,
    OtherChildAlloparents = sum(ChildAlloparent, na.rm = TRUE) - ChildAlloparent,
    HHOtherKidsAge0_2.5 = sum(Age0_2.5, na.rm = TRUE) - Age0_2.5, # Helfrecht and Meehan 2015 style risk periods
    HHOtherKidsAge2.5_5 = sum(Age2.5_5, na.rm = TRUE) - Age2.5_5,
    HHOtherKidsAge5_10 = sum(Age5_10, na.rm = TRUE) - Age5_10,
    HHOtherKidsAge10Plus = sum(Age10Plus, na.rm = TRUE) - Age10Plus,
    HHOtherKidsAge5Minus = sum(Age5Minus, na.rm = TRUE) - Age5Minus,
    .by = householdID
  ) |>
  mutate(
    OtherAlloparents = as.numeric(AdultsChildcare) + OtherChildAlloparents #Check Change to TotalChildCare
  )

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
    UserLanguage, # only 19 children with english speaking
  ) |>
  na.omit()

#' # Sample characteristics

#+ message=F, warning=F, fig.height=9, fig.width=10

# Note 20 year old caregiver reporting on a 20 year old has a CaregiverRelation of "Other"
ggplot(data = d2, aes(CaregiverAge, ChildAge)) + geom_count() + geom_smooth() + labs(title = "Scatterplot of CaregiverAge and ChildAge")

ggplot(caregivers, aes(y=IncomeCategory)) +
  geom_bar() + labs(y = "IncomeCaregory (household level)")

table(d2$Sex)

ggplot(d2, aes(y=ChildAge, fill = Sex)) +
  geom_bar()

# Data for out main models without NAs
ggplot(modeldf, aes(y=ChildAge, fill = Sex)) +
  geom_bar()

caregivers3 <- caregivers |>  filter(!is.na(Neighborhood))
ggplot(caregivers3, aes(y=Neighborhood, fill=ImmigrateUtila)) +
  geom_bar() + labs(title = "Immigration status by neighborhood")

#' # Signaling frequency and cost histograms

#+ message=F, warning=F, fig.height=9, fig.width=10

ggplot(d2, aes(x=SadFreqN, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=CryFreqN, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=TantrumFreqN, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=RunawayFreqN, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=SignalFreq, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=SignalCost, fill=Sex)) +
  geom_histogram()
ggplot(d2, aes(x=SignalFreqMax, fill=Sex)) +
  geom_histogram()

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

#' # Analyses

#+ message=F, warning=F, fig.height=9, fig.width=10

# CONFIRM WE WANT ALLOPARENTINGFreqN

ggplot(data = d2, aes(ChildAge, AlloparentingFreqN, colour = Sex)) + geom_count() + geom_smooth() + labs(title = "Scatterplot of Child Age and AlloparentingFreqN")
ggplot(data = d2, aes(ChildAge, AlloparentingFreqN, colour = Sex)) + geom_jitter(width = 0.25, height = 0.25) + geom_smooth() + labs(title = "Scatterplot of Child Age and AlloparentingFreqN")

mscH <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mscH)
plot(allEffects(mscH))

# Signal frequency (sum of frequencies of sadness, crying, and temper tantrums [running away excluded])

msfH <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(msfH)
plot(allEffects(msfH))

# Frequency of most common signal (running away excluded)

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

#' # Appendix 1: Visualization of signaling by alloparenting contributions and sex

#+ message=F, warning=F, fig.height=9, fig.width=10

# visualize signaling by sex based on alloparenting contributions
ggplot(d2, aes(ChildAge, SignalFreq)) + geom_jitter(aes(colour = Sex), width = .25) + facet_wrap(vars(AlloparentingFreq))
ggplot(d2, aes(ChildAge, SignalFreqMax)) + geom_jitter(aes(colour = Sex), width = .25) + facet_wrap(vars(AlloparentingFreq))
ggplot(d2, aes(ChildAge, SignalCost)) + geom_jitter(aes(colour = Sex), width = .25) + facet_wrap(vars(AlloparentingFreq))

ggplot(d2, aes(ChildAge, TantrumFreqN)) + geom_jitter(aes(colour = Sex), width = .25) + facet_wrap(~factor(AlloparentingFreqN))
ggplot(d2, aes(ChildAge, CryFreqN)) + geom_jitter(aes(colour = Sex), width = .25) + facet_wrap(~factor(AlloparentingFreqN))
ggplot(d2, aes(ChildAge, SadFreqN)) + geom_jitter(aes(colour = Sex), width = .25) + facet_wrap(~factor(AlloparentingFreqN))
ggplot(d2, aes(ChildAge, RunawayFreqN)) + geom_jitter(aes(colour = Sex), width = .25) + facet_wrap(~factor(AlloparentingFreqN))

# Appendix 2: ordinal logistic regression of parental response

#+ message=F, warning=F, fig.height=9, fig.width=10

# note 6 instances of both positive and negative coded as NA
table(d2$PositiveResponse, d2$NegativeResponse)

mp7 <- polr(ParentalResponse ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + UnwantedTask + Punishment2 + Family2 + DiscomfortPainInjuryIllness, d2)
summary(mp7)

coeftest(mp7, vcov = vcovCL, type = "HC0", cluster = ~householdID)
plot_predictions(mp7, condition = c("Punishment2", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mp7, condition = c("DiscomfortPainInjuryIllness", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mp7, condition = c("Family2", "group"), type = "probs", vcov = ~householdID)

#' # Appendix 3: Models without health measures and with AlloparentingFreqN
# note some have interactions between ChildAge and sex that do not exist in the health models

#+ message=F, warning=F, fig.height=9, fig.width=10

# Signal Frequency

# partner status not significant unlike the signalfreqmax model
msf <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(msf)
plot(allEffects(msf))

# Maximum Signal Frequency

# OtherChildrenHH:LogIncome now marginally significant after changing from Neighborhood random effect to Household
mmsf <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH*LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(mmsf)
plot(allEffects(mmsf))

#Signal Cost
msc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(msc)
plot(allEffects(msc))

# Signal specific frequency
msad <- glmmTMB(SadFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(msad)
plot(allEffects(msad))

# ChildAge*Sex interaction (Marginal with AlloparentingFreqN)
mcry <- glmmTMB(CryFreqN ~ ChildAge*Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(mcry)
plot(allEffects(mcry))

mtantrum <- glmmTMB(TantrumFreqN ~ ChildAge+ Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(mtantrum)
plot(allEffects(mtantrum))

# conflict frequency
# AlloparentingFreqN seems to really hurt our power through reduced sample size

mconflict <- glmmTMB(ConflictFreqN ~ ChildAge+ Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + AlloparentingFreqN*Sex + EducationLevelYears + (1|householdID),data = d2, family = nbinom2)
summary(mconflict)
plot(allEffects(mconflict))

#' # Appendix 4: Current health sub-measures

#+ message=F, warning=F, fig.height=9, fig.width=10

# none are significant if all are included
# Currenthealth2 = not a significant predictor

# Currenthealth1: If an illness is 'going around', my child will get it.
mscH1 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth1 + (1|householdID),data = d2, family = nbinom2)
summary(mscH1)
plot(allEffects(mscH1))

msfH1 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth1 + (1|householdID),data = d2, family = nbinom2)
summary(msfH1)
plot(allEffects(msfH1))

mmsfH1 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH*LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth1 + (1|householdID),data = d2, family = nbinom2)
summary(mmsfH1)
plot(allEffects(mmsfH1))

# Currenthealth2: In general, my child is very susceptible to colds, flu, and other infectious diseases.

mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth2 + (1|householdID),data = d2, family = nbinom2)
summary(mscH2)
plot(allEffects(mscH2))

msfH2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth2 + (1|householdID),data = d2, family = nbinom2)
summary(msfH2)
plot(allEffects(msfH2))

mmsfH2 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH*LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth2 + (1|householdID),data = d2, family = nbinom2)
summary(mmsfH2)
plot(allEffects(mmsfH2))

# Currenthealth3: My child’s immune system protects him/her from most illnesses that other kids get.

mscH3 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth3 + (1|householdID),data = d2, family = nbinom2)
summary(mscH3)
plot(allEffects(mscH3))

msfH3 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth3 + (1|householdID),data = d2, family = nbinom2)
summary(msfH3)
plot(allEffects(msfH3))

mmsfH3 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH*LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth3 + (1|householdID),data = d2, family = nbinom2)
summary(mmsfH3)
plot(allEffects(mmsfH3))

# No Currenthealth variables predict tantrums

msadK1 <- glmmTMB(SadFreqN ~ ChildAge*Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth1 + (1|householdID),data = d2, family = nbinom2)
summary(msadK1)
plot(allEffects(msadK1))

# Marginally significant
mcryK1 <- glmmTMB(CryFreqN ~ ChildAge*Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth1 + (1|householdID),data = d2, family = nbinom2)
summary(mcryK1)
plot(allEffects(mcryK1))

msadK3 <- glmmTMB(SadFreqN ~ ChildAge*Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + Currenthealth3 + (1|householdID),data = d2, family = nbinom2)
summary(msadK3)
plot(allEffects(msadK3))

#' # Appendix 5: Runaway frequency model

#+ message=F, warning=F, fig.height=9, fig.width=10

mrunawayH <- glmmTMB(RunawayFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mrunawayH)
plot(allEffects(mrunawayH))

#' # Appendix 6: Ordinal logistic regression models for relative need and maternal investment

#+ message=F, warning=F, fig.height=9, fig.width=10

# Results with NeighborhoodF as a random effect and vc0 seem like outliers compared to all other methods

# Using preferred model for signal frequency (without partner status)
# signal Freq and Child age still significant after controlling for education years
# CurrentHealthMean has no predictive ability when added
mN1 <- polr(RelativeNeed3 ~ SignalFreq + ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults+ ConflictFreqN + AlloparentingFreqN + CurrentHealthMean, d2)
summary(mN1)
coeftest(mN1, vcov = vcovCL, type = "HC0", cluster = ~householdID)
plot_predictions(mN1, condition = c("ChildAge", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mN1, condition = c("SignalFreq", "group"), type = "probs", vcov = ~householdID)
plot_predictions(mN1, condition = c("CurrentHealthMean", "group"), type = "probs", vcov = ~householdID)

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

#' # New Analyses: Alloparenting within the family

#+ message=F, warning=F, fig.height=9, fig.width=12

# Adult alloparenting heavily weighted

#Signal Cost

mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + OtherAlloparentingN2b + (1|householdID),data = d2, family = nbinom2)
summary(mscH2)
plot(allEffects(mscH2))

mscH2i <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + OtherAlloparentingN2b*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mscH2i)
plot(allEffects(mscH2i))

#Signal Frequency

msf2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparentingN2b + (1|householdID),data = d2, family = nbinom2)
summary(msf2)
plot(allEffects(msf2))

msf2i <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparentingN2b*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(msf2i)
plot(allEffects(msf2i))

# Maximum Signal Frequency

mmsf2 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparentingN2b + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2)
plot(allEffects(mmsf2))

mmsf2i <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparentingN2b*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2i)
plot(allEffects(mmsf2i))

# Adult alloparenting lightly weighted

#Signal Cost

mscH3i <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + OtherAlloparentingN2a*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mscH3i)
plot(allEffects(mscH3i))

#Signal Frequency

msf3i <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparentingN2a*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(msf3i)
plot(allEffects(msf3i))

# Maximum Signal Frequency

mmsf3 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparentingN2a + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2)
plot(allEffects(mmsf2))

mmsf3i <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparentingN2a*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mmsf3i)
plot(allEffects(mmsf3i))

# Number of alloparents

hist(d2$OtherAlloparents) # Few beyond 4, most 3 and below

#Signal Cost

mscH2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + OtherAlloparents + (1|householdID),data = d2, family = nbinom2)
summary(mscH2)
plot(allEffects(mscH2))

mscH2i <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + OtherAlloparents*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mscH2i)
plot(allEffects(mscH2i))

#Signal Frequency

msf2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparents + (1|householdID),data = d2, family = nbinom2)
summary(msf2)
plot(allEffects(msf2))

msf2i <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparents*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(msf2i)
plot(allEffects(msf2i))

# Maximum Signal Frequency

mmsf2 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparents + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2)
plot(allEffects(mmsf2))

mmsf2i <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + OtherAlloparents*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mmsf2i)
plot(allEffects(mmsf2i))

#' # New Analyses: Child Age Ranges

#+ message=F, warning=F, fig.height=9, fig.width=9

# Age models with all ranges

# 2.5 - < 5 = different direction

mAllAgeSF <- glmmTMB(SignalFreq ~ ChildAge + Sex + HHOtherKidsAge0_2.5 + HHOtherKidsAge2.5_5 + HHOtherKidsAge5_10 + HHOtherKidsAge10Plus + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mAllAgeSF)
plot(allEffects(mAllAgeSF))

mAllAgeSC <- glmmTMB(SignalCost ~ ChildAge + Sex + HHOtherKidsAge0_2.5 + HHOtherKidsAge2.5_5 + HHOtherKidsAge5_10 + HHOtherKidsAge10Plus + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mAllAgeSC)
plot(allEffects(mAllAgeSC))

mAllAgeSFM <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + HHOtherKidsAge0_2.5 + HHOtherKidsAge2.5_5 + HHOtherKidsAge5_10 + HHOtherKidsAge10Plus + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mAllAgeSFM)
plot(allEffects(mAllAgeSFM))

# Models with just >= 2.5 - < 5 and total children
d2$OtherChildrenHH2 <- d2$OtherChildrenHH - d2$HHOtherKidsAge2.5_5

# Signal Frequency
# with AlloparentingFreq*HHOtherKidsAge2.5_5 interaction and AlloparentingFreqN*Sex interaction
ms2_5SFi <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH2 + HHOtherKidsAge2.5_5*AlloparentingFreqN + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(ms2_5SFi)
plot(allEffects(ms2_5SFi))

# with AlloparentingFreq*HHOtherKidsAge2.5_5 interaction
ms2_5SF <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH2 + HHOtherKidsAge2.5_5*AlloparentingFreqN + LogIncome + number_adults + ConflictFreqN + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(ms2_5SF)
plot(allEffects(ms2_5SF))

# SignalCost
# with AlloparentingFreq*HHOtherKidsAge2.5_5 interaction and AlloparentingFreqN*Sex interaction
ms2_5SCi <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH2 + HHOtherKidsAge2.5_5*AlloparentingFreqN + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(ms2_5SCi)
plot(allEffects(ms2_5SCi))

# with AlloparentingFreq*HHOtherKidsAge2.5_5 interaction
ms2_5SC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH2 + HHOtherKidsAge2.5_5*AlloparentingFreqN + LogIncome + number_adults + ConflictFreqN + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(ms2_5SC)
plot(allEffects(ms2_5SC))

# Max Signal Frequency
# with AlloparentingFreq*HHOtherKidsAge2.5_5 interaction and AlloparentingFreqN*Sex interaction
ms2_5SFMi <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH2 + HHOtherKidsAge2.5_5*AlloparentingFreqN + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(ms2_5SFMi)
plot(allEffects(ms2_5SFMi))

# with AlloparentingFreq*HHOtherKidsAge2.5_5 interaction
ms2_5SFM <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH2 + HHOtherKidsAge2.5_5*AlloparentingFreqN + LogIncome + number_adults + ConflictFreqN + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(ms2_5SFM)
plot(allEffects(ms2_5SFM))






# Previous model but for under 5
msfH5minus <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH2 + HHOtherKidsAge5Minus*AlloparentingFreqN + LogIncome + number_adults + ConflictFreqN + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(msfH5minus)
plot(allEffects(msfH5minus))

# AlloparentingFreqN*Sex back in
# Show


# Frequency of most common signal (running away excluded)

mmsfH <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + HHOtherKidsAge2.5_5 + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mmsfH)
plot(allEffects(mmsfH))

mmsfH <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + HHOtherKidsAge2.5_5*AlloparentingFreqN + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mmsfH)
plot(allEffects(mmsfH))

# This age range was not a significant predictor of specific signals. Marignal effect in some tantrum models.

# Notes:

# Aaron recommended his package for Z-scoring
# Also recommended not worrying about biomarkers

# Survey language as a way to get at differences in culture
# - Spanish: 322
# - English: 53
# - Nothing close to significant in our existing models.
# - Not a significant when the sole predictor.

# Sex of kids (not measured for the youngest kids)
# - no signaling data either so nothing to do here

# Revisit summed alloparenting contributions  (merged child alloparenting with adults)
# - OtherAlloparentingN2 marginally significant in some models.
# - Likely need to change its form: How to handle NAs and 0s?
# - Slightly different wording from previous questions is troublesome (is the caretaker who did the survey included?)

# Explore alloparenting*conflict
# - nothing close to an interaction here when replacing sex with conflict in the alloparenting interaction.

# interactions between conflict and kids and adults (couldn't read my writing here)

# more in the model predicting conflict

# "No significant relationships between Sidama infant signaling
# behaviors and the frequency of care exist. However, results indicate that signaling behaviors were
# positively associated with the size of infants' alloparental network."
# - Jessica Collins dissertation

# Want age categories similar to Helfrecht and Meehan 2015
# counts of sibs/children in specific age ranges? (seems like it makes sense to do this with relatedness)

# OtherAlloparentingN2 Models
# sex by OtherAlloparentingN2 marginalish


malloparenting <- glmmTMB(AlloparentingFreqN ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(malloparenting)
plot(allEffects(malloparenting))

# number of younger kids
# other child care

mscH <- glmmTMB(SignalCost ~ ChildAge*YoungerKids + OlderKids + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + (1|householdID),data = d2, family = nbinom2)
summary(mscH)
plot(allEffects(mscH))


#' # New Analyses: Anthropometrics

#+ message=F, warning=F, fig.height=9, fig.width=9

# NOBS = 72
mscHW <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + HeightZ + WeightZ + (1|householdID),data = d2, family = nbinom2)
summary(mscHW)
plot(allEffects(mscHW))

msfHW <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + HeightZ + WeightZ + (1|householdID),data = d2, family = nbinom2)
summary(msfHW)
plot(allEffects(msfHW))

mscBMI <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BMIZ + (1|householdID),data = d2, family = nbinom2)
summary(mscBMI)
plot(allEffects(mscBMI))

msfBMI <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BMIZ + (1|householdID),data = d2, family = nbinom2)
summary(msfBMI)
plot(allEffects(msfBMI))







msc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + GripR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msc)
plot(allEffects(msc))

# Tricep Residuals (N and N2)
msc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + TricepR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msc)
plot(allEffects(msc))



d2$BodyFat <- c(scale(d2$TricepR) + scale(d2$SubscapR))
msc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msc)
plot(allEffects(msc))

msc <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFat*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msc)
plot(allEffects(msc))





msc <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + SubscapR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msc)
plot(allEffects(msc))





msc <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + BodyFatPercentageR*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msc)
plot(allEffects(msc))

#' # New Analyses: Relatedness within the family

#+ message=F, warning=F, fig.height=9, fig.width=9

# Child relatedness
mCRC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanChildRelatedness + number_adults + (1|householdID),data = d2, family = nbinom2)
summary(mCRC)
plot(allEffects(mCRC))

mCRF <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanChildRelatedness + number_adults + (1|householdID),data = d2, family = nbinom2)
summary(mCRF)
plot(allEffects(mCRF))

# Sibling relatedness
mSRC <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + (1|householdID),data = d2, family = nbinom2)
summary(mSRC)
plot(allEffects(mSRC))

mSRF <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + (1|householdID),data = d2, family = nbinom2)
summary(mSRF)
plot(allEffects(mSRF))

mSRFM <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + (1|householdID),data = d2, family = nbinom2)
summary(mSRFM)
plot(allEffects(mSRFM))

#Sibling Relatedness with kids 2.5 to 5 years old
d2$OtherChildrenHH2 <- d2$OtherChildrenHH - d2$HHOtherKidsAge2.5_5

mSRC2 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH2 + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + HHOtherKidsAge2.5_5 + (1|householdID),data = d2, family = nbinom2)
summary(mSRC2)
plot(allEffects(mSRC2))

mSRF2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH2 + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + HHOtherKidsAge2.5_5 + (1|householdID),data = d2, family = nbinom2)
summary(mSRF2)
plot(allEffects(mSRF2))

mSRFM2 <- glmmTMB(SignalFreqMax ~ ChildAge + Sex + OtherChildrenHH2 + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + HHOtherKidsAge2.5_5 + (1|householdID),data = d2, family = nbinom2)
summary(mSRFM2)
plot(allEffects(mSRFM2))

# Models with OlderKids and Younger Kids instead of children 2.5 - 5

mSRC3 <- glmmTMB(SignalCost ~ ChildAge + Sex + OlderKids + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + YoungerKids + (1|householdID),data = d2, family = nbinom2)
summary(mSRC3)
plot(allEffects(mSRC3))

mSRF3 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OlderKids*ChildAge + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + YoungerKids + (1|householdID),data = d2, family = nbinom2)
summary(mSRF3)
plot(allEffects(mSRF3))

mSRF3 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OlderKids + LogIncome + number_adults + PartnerStatus + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + CurrentHealthMean + MeanSiblingRelatedness*number_adults + YoungerKids*ChildAge + (1|householdID),data = d2, family = nbinom2)
summary(mSRF3)
plot(allEffects(mSRF3))

#' # New Analyses: Adults who provide housework

# Huge difference compared to AlloparentingFreqN

msfhw1 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + AdultsHousework*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msfhw1)
plot(allEffects(msfhw1))

mschw1 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + AdultsHousework*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mschw1)
plot(allEffects(mschw1))

msfhw2 <- glmmTMB(SignalFreq ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + AdultsHousework*Sex + (1|householdID),data = d2, family = nbinom2)
summary(msfhw2)
plot(allEffects(msfhw2))

mschw1 <- glmmTMB(SignalCost ~ ChildAge + Sex + OtherChildrenHH + LogIncome + number_adults + ConflictFreqN + AlloparentingFreqN*Sex + EducationLevelYears + AdultsHousework*Sex + (1|householdID),data = d2, family = nbinom2)
summary(mschw1)
plot(allEffects(mschw1))









d5 <- d2 |>
  dplyr::filter(!is.na(CryFreq))
table(is.na(d5$AlloparentingFreq), d5$ChildAge, useNA = "a")
