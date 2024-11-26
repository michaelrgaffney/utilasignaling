

# Data sharing after collection: With field data of the sort to be collected in
# this project, a balance must be struck between open science and privacy
# protection of participants, particularly those who are minors. Because Utila
# is a small island, many types of data might be easily triangulated to become
# individually identifiable. For this reason, we will publicly share data
# necessary for replication of primary analyses, such as growth data, hormone
# levels, and non-identifiable survey items, as well as composite measures such
# as factor loadings for energy availability and mortality risk.

# Identifiable items such as birth dates, household locations, and family
# composition will be shared publicly only in aggregated tables and not at the
# level of individual data. Additionally, once data is collected we will need to
# assess which other items might be considered identifying. For example, if only
# a few individuals have parents absent due to death, then parental death might
# be an identifiable trait, and we would therefore need to be cautious about
# sharing that data publicly in a way that could link it to individual
# participants.


source('dataprep.R')
library(sdcMicro)

AddedVars <-
  d2 |>
  dplyr::filter(!(is.na(CryFreqN) & is.na(SadFreqN) & is.na(TantrumFreqN) & is.na(ConflictFreqN))) |>
  dplyr::select(
    householdID,
    childHHid,
    RelativeNeed3,
    RelativeMaternalInvestment2,
    ConflictFamily:StatusConcerns,
    CaregiverResponse
  )

StudyVars <- SignalVars |>
  dplyr::select(-AlloparentingXsex) |>
  left_join(AddedVars, by = c('householdID', 'childHHid'))

# Definitely
# ChildAge
# Sex
# CaregiverAge
# NumberOfChildren
# number_adults
# Neighborhood2
# ImmigrateUtila
# UserLanguage

# Maybe
# StayAtHomeMom
# EducationLevelYears
# Lifestyle Reality vars
# FoodSecurity

# Probably not
# MeanChildRelatedness
# MedicalProblemsMean
# IllnessSusceptibilityMean
# HouseQuality
# LogIncome (but category based?)

sdcObj <- createSdcObj(
  StudyVars,
  keyVars = c("Sex", "Neighborhood2", "ImmigrateUtila", "UserLanguage"),
  numVars = c("ChildAge", "NumberOfChildren", "number_adults")
)
measure_risk(sdcObj)

sdcObj <- createSdcObj(
  StudyVars,
  keyVars = c("Sex", "Neighborhood2", "ImmigrateUtila"),
  numVars = c("ChildAge", "NumberOfChildren", "number_adults")
)
measure_risk(sdcObj)

sdcObj <- createSdcObj(
  StudyVars,
  keyVars = c("Sex", "Neighborhood2", "NumberOfChildren", "number_adults"),
  numVars = c("ChildAge")
)
measure_risk(sdcObj)


out <- measure_risk(
  as.data.frame(StudyVars[c('householdID', "Sex", "Neighborhood2", "ImmigrateUtila", "UserLanguage")]),
  keyVars = c("Sex", "Neighborhood2", "ImmigrateUtila", "UserLanguage"),
  hid = 'householdID'
)
out

out <- measure_risk(
  as.data.frame(StudyVars[c('householdID', "Sex", "Neighborhood2")]),
  keyVars = c("Sex", "Neighborhood2"),
  hid = 'householdID'
)
out

out <- measure_risk(
  as.data.frame(StudyVars),
  keyVars = c("Sex", "Neighborhood2", "ImmigrateUtila", "NumberOfChildren", "number_adults"),
  hid = 'householdID'
)
out

# Using microaggregation
sdcObj <- createSdcObj(
  StudyVars,
  keyVars = c("Sex", "Neighborhood2", "NumberOfChildren", "number_adults"),
  numVars = c("ChildAge")
)
out <- microaggregation(sdcObj)
