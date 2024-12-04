

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

# Run in a fresh R session

#+ message=F, warning=F
source("recode.R")
source("dictionaries.R")
source("dataprep.R")

library(sdcMicro)

# Cluster ages
discretize <- function(v, n){
  q <- quantile(v, seq(0, 1, 1/n))[-1]
  map_int(v, \(x) min(q[x<=q]))
}

# The householdID variable now has a unique random # for
# each child, not each household
new_hid <- sample(nrow(utila_df))
names(new_hid) <- utila_df$householdID

utila_df <-
  utila_df |>
  mutate(
    householdID = new_hid,
    childHHid = householdID,
    uniqueID = householdID,
    Neighborhood2 = ifelse(is.na(Neighborhood2), c(0,1), Neighborhood2),
    ChildAge = discretize(ChildAge, 4) # Child ages clustered into quartiles
  ) |>
  arrange(householdID) |>
  # These variables are randomized within the two neighborhoods
  mutate(
    NumberOfChildren = sample(NumberOfChildren, n()),
    number_adults = sample(number_adults, n()),
    UserLanguage = sample(UserLanguage, n()),
    .by = Neighborhood2
  )

# householdIDs are set to the new ones generated above
causes <-
  causes |>
  mutate(
    householdID = new_hid[householdID],
    childHHid = householdID,
    uniqueID = householdID
  )

# Assess disclosure risk
# Aim for k-anonymity >= 4
#+ message=T
sdcObj <- createSdcObj(
  utila_df,
  keyVars = c("ChildAge", "Sex", "Neighborhood2", "ImmigrateUtila")
)
measure_risk(sdcObj)

# Delete values to achieve k-anonymity >= 4
out <- localSuppression(sdcObj, k=4)
out
plot(out)

utila_df$ChildAge <- out@manipKeyVars$ChildAge

save(utila_df, causes, householdIDs, caregiverSex, file = "data/utilasignalingPublicData.rda")

