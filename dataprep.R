
library(utiladata2023) # data package

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

d <-
  children |>
  remove_labels() |>
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
  left_join(dplyr::select(remove_labels(caregivers), householdID, CPRatio, Neighborhood, ImmigrateUtila, IncomeCategory, IncomeCategoryN,
                          EducationLevel, EducationLevelYears, AdultsMoney, AdultsHousework, AdultsChildcare, AdultsNoChildcare, CaregiverAge,
                          number_children2, number_adults, UserLanguage, CurrentJob, contains("CaregiverMarital"), NeighborhoodQuality, HouseQuality,
                          contains("SocialSupport"), contains("FoodSecurity"), contains("HomeInstability_"), contains("Reality_")), by = "householdID") |>
  left_join(dplyr::select(causes, -householdID, -childHHid, -text), by = "uniqueID")

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
    ConflictFamily2 = ifelse(ConflictFamily == "Maybe", 0, ConflictFamily),
    ConflictOutsideFamily2 = ifelse(ConflictOutsideFamily == "Maybe", 0, ConflictOutsideFamily),
    # WantNeed2 = ifelse(WantNeed == "Want" | WantNeed == "Need", WantNeed, NA),
    # WantNeedBinary = case_when(
    #   WantNeed2 == "Want" ~ 0,
    #   WantNeed2 == "Need" ~ 1),
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
utila_df <-
  d2 |>
  dplyr::filter(!(is.na(CryFreqN) & is.na(SadFreqN) & is.na(TantrumFreqN))) |>
  dplyr::select(
    uniqueID,
    householdID,
    childHHid,
    CryFreqN,
    SadFreqN,
    TantrumFreqN,
    SignalFreq,
    SignalCost,
    SignalFreqMax,
    CaregiverResponse,
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
    CaregiverAge,
    FoodSecurity,
    Neighborhood2,
    Sex,
    UserLanguage,
    ImmigrateUtila,
    PartnerStatus,
    RelativeNeed,
    RelativeMaternalInvestment,
    RelativeMaternalInvestment2,
    Neighborhood2F,
    StayAtHomeMom,
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

#HomeInstabilityMean:Children in unstable homes likely face more adversity and find themselves having to signal to new people or face the challenge of previous caretakers leaving or moving away from friends.X

householdIDs <- sample(utila_df$householdID, length(utila_df$householdID))
caregiverSex <- sort(caregivers$CaregiverSex[caregivers$householdID %in% utila_df$householdID], na.last = T)
