library(utiladata2023) # data package
library(tidyverse)
library(labelled)
library(mgcv)
library(localgrowth)

source("recode.R")
source("dictionaries.R")

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

children2 <-
  children |>
  dplyr::filter(ChildID != "") |>
  dplyr::select(ChildID, householdID)

anthropometricMeans0 <-
  anthropometrics |>
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

anthropometricMeans <-
  anthropometricMeans0 |>
  dplyr::filter(
    measurements == 1 | year(Date.of.measurement) == 2023
  ) |>
  left_join(children2) # note that we filter out householdID based on it's column number

anthropometricMeans$WeightZ <- growthRef(AgeAtMeasurement, WeightMean_KG, Sex, anthropometricMeans, type = "Weight", pop = "CDC")
anthropometricMeans$HeightZ <- growthRef(AgeAtMeasurement, HeightMean, Sex, anthropometricMeans, type = "Height", pop = "CDC")
anthropometricMeans$BMIZ <- growthRef(AgeAtMeasurement, BMI, Sex, anthropometricMeans, type = "BMI", pop = "CDC")
anthropometricMeans$Sex2 <- factor(anthropometricMeans$Sex)

m <- mgcv::gam(GripStrengthMean ~ Sex+s(AgeAtMeasurement, by = Sex2, k = 4), data = anthropometricMeans)
anthropometricMeans$GripR <- residuals(m)

# Result for bodyfat predicitng crying freq changes if k is unspecified
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

d <-
  children |>
  remove_labels() |>
  #dplyr::filter(CompleteSurvey) %>%
  mutate(
    OnlyChild = n() == 1,
    YoungerKids = map_int(ChildAge, \(x) sum(x > ChildAge)),
    OlderKids = map_int(ChildAge, \(x) sum(x < ChildAge)),
    OlderGirls = map_int(ChildAge, \(x) sum(x < ChildAge & Sex == 'Female')),
    OlderBoys = map_int(ChildAge, \(x) sum(x < ChildAge & Sex == 'Male')),
    OneOlderGirl = OlderGirls > 0,
    OldestChild = ChildAge == max(ChildAge) & (!OnlyChild),
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
    NegativeResponseF = as.factor(NegativeResponseBinary)
  ) |>
  rowwise() |>
  mutate(
    SignalFreq = sum2(SadFreqN, CryFreqN, TantrumFreqN, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    SignalFreqMax = max2(SadFreqN, CryFreqN, TantrumFreqN, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    SignalCost = sum2(SadFreqN, CryFreqN * 2, TantrumFreqN * 3, na.rm = TRUE), # RunawayFreqN removed because interpretation not clear
    IllnessSusceptibilitySum = sum2(IllnessSusceptibility1, IllnessSusceptibility2, IllnessSusceptibility3, na.rm = TRUE),
    IllnessSusceptibilityMean = mean2(IllnessSusceptibility1, IllnessSusceptibility2, IllnessSusceptibility3, na.rm = TRUE),
    MedicalProblemsSum = sum2(SchoolAges1, SchoolAges2, na.rm = TRUE),
    MedicalProblemsMean = mean2(SchoolAges1, SchoolAges2, na.rm = TRUE) # 1 child has data for only 1
  ) |>
  ungroup() |>
  left_join(anthropometricMeans[-21], by = "ChildID") |>  # removing householdID column
  left_join(dplyr::select(caregivers, householdID, CPRatio, Neighborhood, ImmigrateUtila, IncomeCategory, IncomeCategoryN,
                          EducationLevel, EducationLevelYears, AdultsMoney, AdultsHousework, AdultsChildcare, AdultsNoChildcare, CaregiverAge,
                          number_children2, number_adults, UserLanguage, CurrentJob, contains("CaregiverMarital"), NeighborhoodQuality, HouseQuality,
                          contains("SocialSupport"), contains("FoodSecurity"), contains("HomeInstability_"), contains("Reality_")), by = "householdID") |>
  left_join(dplyr::select(causes, -householdID, -childHHid, -ChildID, -Sad1), by = "uniqueID")

# prepare household variables for analyses
d2 <-
  d |>
  mutate(
    LogIncome = log10(IncomeCategoryN),
    CaregiverAge = ifelse(is.na(CaregiverAge), mean(CaregiverAge, na.rm = TRUE), CaregiverAge), # imputing age for 1 caretaker where it is missing (affects 2 children). Also imputes caregiver ages for 7 other children. But these children do not have signaling data.
    StayAtHomeMom = ifelse(str_detect(CurrentJob, "home"), 1, 0),
    StayAtHomeMomF = as.factor(case_when(
      StayAtHomeMom == 0 ~ "No",
      StayAtHomeMom == 1 ~ "Yes")),
    BodyFat = c(scale(TricepR) + scale(SubscapR)),
    RecentSignalSeverity = com1[as.character(Sad1)],
    RecentSignalCause = com2[as.character(Sad1)],
    AdultsHousework = as.numeric(AdultsHousework),
    #AdultsChildcare = as.numeric(AdultsChildcare),
    Neighborhood2 = as.numeric(Neighborhood == "Camponado/Campolancho"),
    Neighborhood2F = as.factor(case_when(
      Neighborhood2 == 0 ~ "No",
      Neighborhood2 == 1 ~ "Yes")),
    NeighborhoodF = as.factor(Neighborhood),
    OtherChildrenHH = number_children2 - 1,
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
    HomeInstabilityMean = mean2(HomeInstability_1, HomeInstability_2, HomeInstability_3, HomeInstability_4)
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
modeldf <-
  d2 |>
  dplyr::filter(!(is.na(CryFreqN) & is.na(SadFreqN) & is.na(TantrumFreqN))) |>
  dplyr::select(
    householdID,
    childHHid,
    CryFreqN,
    SadFreqN,
    TantrumFreqN,
    SignalFreq,
    SignalCost,
    SignalFreqMax,
    ChildAge,
    OtherChildrenHH,
    YoungerKids,
    OtherChildAlloparentingFreqN,
    OlderGirls,
    OldestChild,
    LogIncome,
    number_adults,
    ConflictFreqN,
    AlloparentingFreqN,
    EducationLevelYears,
    IllnessSusceptibility1,
    IllnessSusceptibility2,
    IllnessSusceptibility3,
    MedicalProblemsMean,
    IllnessSusceptibilityMean,
    AdultsChildcare,
    AdultsHousework,
    SadFreqOF,
    CryFreqOF,
    TantrumFreqOF,
    ConflictFreqOF,
    AlloparentingFreq0,
    NeighborhoodQuality,
    HouseQuality,
    NumberOfChildren,
    MeanChildRelatedness,
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
    FlexedR,
    CaregiverAge,
    FoodSecurity,
    Neighborhood2,
    Sex,
    UserLanguage,
    ImmigrateUtila,
    PartnerStatus,
    RelativeNeed,
    RelativeMaternalInvestment,
    Neighborhood2F,
    StayAtHomeMomF,
    HomeInstability_1,
    HomeInstability_2,
    HomeInstability_3,
    HomeInstability_4,
    LifestyleReality_1,
    LifestyleReality_2,
    LifestyleReality_3,
    LifestyleReality_4,
    LifestyleReality_5,
    LifestyleReality_6,
    LifestyleReality_7,
    LifestyleReality_8
  )


# dataframes for analysis -------------------------------------------------

# new vars:
# HomeInstabilityMean
# StayAtHomeMom

# vars for consideration:
# Lifestyle perception/reality

# x <- tribble(
#   ~Var,         ~Rationale,
#   "ChildAge",   "Younger children have more needs they cannot address on their own and face lower reputational costs to crying. Older children engage in costlier signals."
# )

signal_vars_temp <-
  "ChildAge:Younger children have more needs they cannot address on their own and face lower reputational costs to crying. Older children engage in costlier signals.X
  Sex:The sexes differ in the adversity they face, the forms of competetion they will engage in with peers, their ability to bargain for increased support outside of signals of need, and the reputational costs of signaling need. These differences may be small or nonexistent in young children.X
  IllnessSusceptibilityMean:Children who are sick more require more energy to make up for the costs of the illness*immune response interaction. They may also benefit from increased investment if this has the potential to increase immune function.X
  MedicalProblemsMean:Children who are sick more require more energy to make up for the costs of the illness*immune response interaction. They may also benefit from increased investment if this has the potential to improve their health or increase immune function. No pathogenic health issues may limit the ability of the child to satisfy their own needs, something which can lead to greater payoffs for signaling need.X
  AlloparentingFreqN:Children might prefer to spend time outside of the family (e.g., to invest in their embodied capital or to find cooperative partners and maintain those relationships) and signal to reduce caregiver pressure to alloparent.X
  CaregiverAge:Overtime, caretakers learn how to better satsify child needs, resist child bargaining, communicate the reasons for their actions. Older moms also have lower residual reproductive value, somethign which may result in greater valuation of investment in existing children compared to younger moms.X
  StayAtHomeMom:Stay at home mom's are capable of providing more investment to children and see them more often.X
  EducationLevelYears:Greater education leads to differences in caretaking behavior and mental models of caretaking while also increasing caregiver earning potential.X
  AdultsChildcare:More adults helping with child care can lead to greater potential for investment. This might lead to decreased payoffs to singaling need in some circumstances. However, it may also increase signaling behaviors in others due to the presence of more signal targerts.X
  number_adults:More adults equals more targets for signaling.X
  AdultsHousework:More adults helping around the house likely leads to less pressure from caretakers for children engage in household labor.X
  NumberOfChildren:More children equals more behavioral conflict over investment.X
  OtherChildAlloparentingFreqN:More alloparenting effort from other children leads to more care for younger children and less pressure for the focal child to alloparent for other chlildren of alloparenting age.X
  LogIncome:The pool of availible resources influences the benefits to signaling.X
  HouseQuality:Children might determine relative status based in part on house quality. This alters the payoffs to signaling.X
  LifestyleReality_1:X
  LifestyleReality_2:X
  LifestyleReality_3:X
  LifestyleReality_4:X
  LifestyleReality_5:X
  LifestyleReality_6:X
  LifestyleReality_7:X
  LifestyleReality_8:X
  HomeInstability_1:X
  HomeInstability_2:X
  HomeInstability_3:X
  HomeInstability_4:X
  NeighborhoodQuality:Parental investment can help children succeed in dangerous environments through teaching and increased child status. This affects the payoff to signaling. Children also pick up on their relative status, which may motivate increased signaling to better compete with wealthier kids.X
  Neighborhood2:Campanado is percieved to be a lower quality environment for children and likely is based on demographics, population density, and differences in status and discrimination. Parental investment can help children succeed in dangerous environments through teaching and increased child status. This affects the payoff to signaling. Children also pick up on their relative status, which may motivate increased signaling to better compete with wealthier kids.X
  FoodSecurity:Children may be more likely to signal for more food when it is hard to come by.X
  ImmigrateUtila:There may be cultural differences in child signaling or parental responsiveness. Differences in resources by group influences the payoffs to signaling.X
  UserLanguage:There may be cultural differences in child signaling or parental responsiveness. Differences in resources by group influences the payoffs to signaling."

#HomeInstabilityMean:Children in unstable homes likely face more adversity and find themselves having to signal to new people or face the challenge of previous caretakers leaving or moving away from friends.X

signalvars <- data.frame(signal_vars_temp) |>
  separate_rows(signal_vars_temp, sep = "X") |>
  separate(col = signal_vars_temp, sep = ":", into = c("Var", "Rationale")) |>
  mutate(
    Var = str_trim(Var)
  )

SignalVars <- d2 |>
  dplyr::filter(!(is.na(CryFreqN) & is.na(SadFreqN) & is.na(TantrumFreqN) & is.na(ConflictFreqN))) |>
  dplyr::select(
    householdID,
    childHHid,
    SadFreqN, CryFreqN, TantrumFreqN, SignalFreq, SignalCost, ConflictFreqN, MeanChildRelatedness, # Fullsibs, Halfsibs, Stepsibs,
    all_of(signalvars$Var)
    # OnlyChild, Only children do not have to compete for attention or other forms of investment with existing children. They still may be motivated to signal for more investment which parents might prefer to devote to future children.X
    # OldestChild, Does not seem to add much beyond age + number of children.X
    # PartnerStatus, Does this add much beyond the number of adults in the household without better data on adult-child relatedness?
    # HouseQuality, Does this add much beyond neighborhood variables?
    ) |>
  mutate(
    across(starts_with("Lifestyle"), \(x) ifelse(x == "No", 0, 1)),
    Sex = ifelse(Sex == "Female", 0, 1),
    UserLanguage = ifelse(UserLanguage == "EN", 0, 1),
    ImmigrateUtila = ifelse(ImmigrateUtila == "No", 0, 1),
    # PartnerStatus = ifelse(PartnerStatus == "Unpartnered", 0, 1),
    AlloparentingXsex = AlloparentingFreqN * Sex,
  ) |>
  na.omit()

SignalVarsAnthro0 <-
  anthropometricMeans0 |>
  dplyr::filter(
    measurements == 2
  ) |>
  pivot_wider(names_from = "Year", values_from = BicepMean:BMI, id_cols = ChildID) |>
  left_join(d2[c("householdID", "childHHid", "ChildID")], by = 'ChildID') |>
  dplyr::filter(!is.na(householdID) & !is.na(childHHid))

SignalVarsAnthro <-
  SignalVars |>
  left_join(SignalVarsAnthro0) |>
  dplyr::filter(!is.na(HeightMean_2024))

# out <- skim(anthropometricMeans)
# nms <- out$skim_variable[out$skim_type == 'numeric' & out$complete_rate > 0.9]

tmp <- left_join(anthropometricMeans, d2[c("childHHid", "ChildID")], by = 'ChildID')

SignalVarsAnthro2 <-
  left_join(SignalVars, tmp) |>
  dplyr::select(
    -c(IllnessSusceptibilityMean:AlloparentingXsex),
    -MeanChildRelatedness,
    -AgeAtMeasurement,
    -Year,
    -measurements,
    -Date.of.measurement,
    -Sex2,
    -ChildID,
    -FlexedRb,
    -FlexedMean,
    -TricepMean,
    -SubscapMean,
    -HeightMean,
    -WeightMean_KG,
    -BMI,
    -BodyFatPercentageMean,
    -GripStrengthMean
    ) |>
  na.omit()

conflict_vars_temp <-
  "ChildAge:Younger children have more needs they cannot address on their own, something which can lead to conflict. Older children are more likely to have conflicts with parents over investment within vs. outside the family.X
  Sex:The sexes differ in the adversity they face and their ability to bargain for increased support outside of signals of need. These differences can lead to differences in conflict frequency. These differences may be small or nonexistent in young children.X
  AlloparentingFreqN:Children might prefer to spend time outside of the family, while parents might prefer more investment within the family.X
  OnlyChild:Only children do not have to compete for attention or other forms of investment with existing children, something which results in conflict in families with more children.X
  CaregiverAge:Caretakers learn how to better satsify child needs and communicate the reasons for their actions. Older moms also have lower residual reproductive value, resulting in greater valuation in investment compared to younger moms, all else being equal.X
  StayAtHomeMom:Stay at home mom's are capable of providing more investment to children and see them more often.X
  EducationLevelYears:Greater education leads to differences in caretaking behavior and increases caregiver earning potential, something which can affect conflict frequency.X
  AdultsChildcare:More adults helping with child care equals more investment for children and may diffuse conflict throughout the family.X
  number_adults:More adults equals more targets for signaling, which parents can see as conflict.X
  AdultsHousework:More adults helping around the house leads to less pressure from caretakers for children to alloparent.X
  NumberOfChildren:More children equals more conflict over investment.X
  OtherChildAlloparentingFreqN:More alloparenting effort from other children leads to more care for younger children and less pressure for the focal child to alloparent for older chlildren.X
  LogIncome:The pool of availible resources influences the benefits to signaling, which parents can see as conflict.X
  HomeInstabilityMean:Children in unstable homes likely face more adversity and find themselves having to signal to new people or face the challenge of previous caretakers leaving or moving away from friends. All of these can lead to conflict and signaling, which may be perceived as conflict.X
  NeighborhoodQuality:Parental investment can help children succeed in dangerous environments through teaching and increased child status.X
  Neighborhood2:Neighborhood can affect the strategies children use to compete among peers. Some of these strategies might lead to conflict with parents.X
  ImmigrateUtila:There may be cultural differences in child signaling or parental responsiveness. Differences in resources by group influences the payoffs to signaling. All of these might affect the frequency of conflict.X
  UserLanguage:There may be cultural differences in child signaling or parental responsiveness. Differences in resources by group influences the payoffs to signaling. All of these might affect the frequency of conflict."

# conflictvars <- data.frame(conflict_vars_temp) |>
#   separate_rows(conflict_vars_temp, sep = "X") |>
#   separate(col = conflict_vars_temp, sep = ":", into = c("Var", "Rationale")) |>
#   mutate(Var = str_trim(Var))
#
# # one option is to just use signaling vars as parents might see signaling as conflict
# ConflictVars <- d2 |>
#   dplyr::filter(!(is.na(CryFreqN) & is.na(SadFreqN) & is.na(TantrumFreqN))) |>
#   dplyr::select(
#     householdID,
#     childHHid,
#     ChildAge,
#     Sex,
#     IllnessSusceptibilityMean, # Children who are sick more require more energy to make up for the costs of the illness*immune response interaction. They may also benefit from increased investment if this has the potential to increase immune function. This can favor signaling, which parents see as conflict.
#     MedicalProblemsMean, # Children who are sick more require more energy to make up for the costs of the illness*immune response interaction. They may also benefit from increased investment if this has the potential to improve their health or increase immune function. This can favor signaling, which parents see as conflict.
#     AlloparentingFreqN,
#     # OnlyChild, Only children do not have to compete for attention or other forms of investment with existing children. They still may be motivated to signal for more investment which parents might prefer to devote to future children.
#     # OldestChild, Does not seem to add much beyond age + number of children.
#     CaregiverAge,
#     StayAtHomeMom,
#     # PartnerStatus, Does this add much beyond the number of adults in the househld?
#     EducationLevelYears,
#     AdultsChildcare,
#     number_adults,
#     NumberOfChildren,
#     OtherChildAlloparentingFreqN,
#     LogIncome,
#     HomeInstabilityMean,
#     NeighborhoodQuality,
#     Neighborhood2,
#     # HouseQuality, Does this add much?X
#     # FoodSecurity Children may be more likely to signal for more food when it is hard to come by. This can be seen as conflict by caretakers.
#     ImmigrateUtila,
#     UserLanguage
#   )

