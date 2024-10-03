# long form dictionary (for use in summary tables)
longform_dict <-
  c(
    ChildAge = "Age (years)",
    SadFreqN = "Sad frequency (times per month)",
    CryFreqN = "Crying frequency (times per month)",
    TantrumFreqN = "Temper tantrum frequency (times per month)",
    AlloparentingFreqN = "Frequency of child's own alloparenting effort (times per month)",
    ConflictFreqN = "Frequency of conflict with focal caretaker (times per month)",
    SignalFreq = "Frequency of signals summed (times per month)",
    SignalCost = "Signal cost",
    SignalFreqMax = "Frequency of most common signal (times per month)",
    IllnessSusceptibilityMean = "Illness avoidance index",
    HeightZ = "Height Z-score (CDC)",
    WeightZ = "Weight Z-score (CDC)",
    BMIZ = "BMI Z-score (CDC)",
    SubscapMean = "Subscapular skinfold thickness (mm)",
    SubscapR = "Subscapular skinfold residuals",
    TricepMean = "Tricep skinfold thickness (mm)",
    TricepR = "Tricep skinfold residuals",
    BodyFat =  "Body fat composite", # child summary stats end here (1:18)
    CaregiverAge = "Caregiver age (years)",
    EducationLevelYears = "Caregiver education level (years)",
    IncomeCategoryN = "Monthly household income (Lempira)",
    number_adults = "Adults in household",
    AdultsChildcare = "Household adults who provide childcare" # adults summary table (19:23)
  )

# short form dictionary (for use in regression tables and correlation plots)
shortform_dict <-
  c(
    ChildAge = "Age (years)",
    Sex = "Sex",
    SexMale = "Sex (male)",
    OtherChildrenHH = "Other children",
    LogIncome = "Income (log)",
    number_adults = "Adults in home",
    PartnerStatusUnpartnered = "Partner status (single)",
    EducationLevelYears = "Parent ed. (years)",
    ConflictFreqN = "Conflict freq.",
    SignalFreq = "Signal freq.",
    SignalCost = "Signal cost",
    SignalFreqMax = "Max signal freq.",
    IllnessSusceptibilityMean = "Perceived illness susceptibility",
    SadFreqN = "Sadness freq.",
    CryFreqN = "Crying freq.",
    TantrumFreqN = "Tantrum freq.",
    AlloparentingFreqN = "Child allocare",
    TricepR = "Tricep skinfold (residuals)",
    SubscapR = "Subscapular skinfold (residuals)",
    BodyFat = "Body fat index",
    RelativeNeed3 = "Relative need (w/in household)",
    UnwantedTask = "Unwanted task",
    Punishment21 = "Punishment",
    Family21 = "Family conflict",
    DiscomfortPainInjuryIllness = "Discomfort/Pain/Injury/Illness",
    CaregiverResponse = "Caregiver Response",
    OlderKids = "Older children",
    YoungerKids = "Younger children",
    AdultsNoChildcare = "Adults that do not provide childcare",
    AdultsChildcare = "Adults that provide childcare",
    OtherChildAlloparentingFreqN = "Other household child alloparenting freq.",
    HeightMean_2024 = "2024 height (cm)",
    HeightMean_2023 = "2023 height (cm)",
    WeightMean_KG_2024 = "2024 weight (kg)",
    WeightMean_KG_2023 = "2023 weight (kg)"
  )

# shortform_dict[c(1,2,9,10,11,12,14,15,16)]

# kable(custom.summarize(modeldf, c(dict1[1,3], dict2[1,3,5]), facvar = "Sex")) %>% kable_styling()

# kable(custom.summarize(modeldf, vars_c[c(1,4,5)], facvar = "Sex")) %>% kable_styling()
