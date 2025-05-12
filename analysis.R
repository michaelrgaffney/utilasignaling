library(ggplot2)
library(forcats)
library(patchwork)
library(modelsummary)
library(ggcorrplot)
library(pvclust)
library(ggalluvial)
library(skimr)
library(glmnet)
library(tidymodels)
library(poissonreg)
library(marginaleffects)
library(easybgm)
library(tidygraph)
library(ggraph)
library(ordinalNet)
library(glmmTMB)
library(effects)

source("pca_plot_functions.R")

utila_FULLSIG <-
  utila_df |>
  dplyr::filter(!is.na(CryFreqN) & !is.na(SadFreqN) & !is.na(TantrumFreqN))

signal_vars <- c(
  ChildAge = "Younger children have more needs they cannot address on their own and face lower reputational costs to crying. Older children engage in costlier signals.",
  Sex = "The sexes differ in the adversity they face, the forms of competetion they will engage in with peers, their ability to bargain for increased support outside of signals of need, and the reputational costs of signaling need. These differences may be small or nonexistent in young children.",
  IllnessSusceptibilityMean = "Children who are sick more require more energy to make up for the costs of the illness*immune response interaction. They may also benefit from increased investment if this has the potential to increase immune function.",
  MedicalProblemsMean = "Children who are sick more require more energy to make up for the costs of the illness*immune response interaction. They may also benefit from increased investment if this has the potential to improve their health or increase immune function. No pathogenic health issues may limit the ability of the child to satisfy their own needs, something which can lead to greater payoffs for signaling need.",
  AlloparentingFreqN = "Children might prefer to spend time outside of the family (e.g., to invest in their embodied capital or to find cooperative partners and maintain those relationships) and signal to reduce caregiver pressure to alloparent.",
  CaregiverAge = "Overtime, caretakers learn how to better satsify child needs, resist child bargaining, communicate the reasons for their actions. Older moms also have lower residual reproductive value, somethign which may result in greater valuation of investment in existing children compared to younger moms.",
  StayAtHomeMom = "Stay at home mom's are capable of providing more investment to children and see them more often.",
  EducationLevelYears = "Greater education leads to differences in caretaking behavior and mental models of caretaking while also increasing caregiver earning potential.",
  AdultsChildcare = "More adults helping with child care can lead to greater potential for investment. This might lead to decreased payoffs to singaling need in some circumstances. However, it may also increase signaling behaviors in others due to the presence of more signal targerts.",
  number_adults = "More adults equals more targets for signaling.",
  AdultsHousework = "More adults helping around the house likely leads to less pressure from caretakers for children engage in household labor.",
  NumberOfChildren = "More children equals more behavioral conflict over investment.",
  OtherChildAlloparentingFreqN = "More alloparenting effort from other children leads to more care for younger children and less pressure for the focal child to alloparent for other chlildren of alloparenting age.",
  LogIncome = "The pool of availible resources influences the benefits to signaling.",
  HouseQuality = "Children might determine relative status based in part on house quality. This alters the payoffs to signaling.",
  LifestyleReality_1 = "",
  LifestyleReality_2 = "",
  LifestyleReality_3 = "",
  LifestyleReality_4 = "",
  LifestyleReality_5 = "",
  LifestyleReality_6 = "",
  LifestyleReality_7 = "",
  LifestyleReality_8 = "",
  HomeInstability_1 = "",
  HomeInstability_2 = "",
  HomeInstability_3 = "",
  HomeInstability_4 = "",
  NeighborhoodQuality = "Parental investment can help children succeed in dangerous environments through teaching and increased child status. This affects the payoff to signaling. Children also pick up on their relative status, which may motivate increased signaling to better compete with wealthier kids.",
  Neighborhood2 = "Campanado is percieved to be a lower quality environment for children and likely is based on demographics, population density, and differences in status and discrimination. Parental investment can help children succeed in dangerous environments through teaching and increased child status. This affects the payoff to signaling. Children also pick up on their relative status, which may motivate increased signaling to better compete with wealthier kids.",
  FoodSecurity = "Children may be more likely to signal for more food when it is hard to come by.",
  ImmigrateUtila = "There may be cultural differences in child signaling or parental responsiveness. Differences in resources by group influences the payoffs to signaling.",
  UserLanguage = "There may be cultural differences in child signaling or parental responsiveness. Differences in resources by group influences the payoffs to signaling."
)

SignalVars <-
  utila_df |>
  dplyr::filter(
    !(is.na(CryFreqN) &
    is.na(SadFreqN) &
    is.na(TantrumFreqN) &
    is.na(ConflictFreqN))
  ) |>
  dplyr::select(
    householdID,
    childHHid,
    SadFreqN,
    CryFreqN,
    TantrumFreqN,
    SignalFreq,
    SignalCost,
    ConflictFreqN,
    MeanChildRelatedness, # Fullsibs, Halfsibs, Stepsibs,
    all_of(names(signal_vars))
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
  na.omit() |>
  mutate(
    OtherKidsSignalCost = sum(SignalCost) - SignalCost,
    OtherKidsConflict = sum(ConflictFreqN) - ConflictFreqN,
    .by = householdID
  )

## main paper plots --------------------------------------------------------

signalfreqdf <-
  utila_df |>
  dplyr::select(
    householdID,
    childHHid,
    # uniqueID,
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
    across(
      Sadness:Tantrums,
      \(x) factor(c(freq_short[x]), levels = c(freq_short))
    )
  ) |>
  na.omit()

signaldf_long <-
  sfdf |>
  pivot_longer(
    cols = Sadness:Tantrums,
    names_to = "Signal",
    values_to = "Freq"
  ) |>
  mutate(
    Signal = factor(Signal, levels = c("Sadness", "Crying", "Tantrums"))
  ) |>
  na.omit()

barplot_SignalFreq <-
  signaldf_long |>
  group_by(Signal, Freq) |>
  summarise(N = n()) |>
  ggplot(aes(N, Signal, fill = fct_rev(Freq))) +
  geom_col(position = "stack", color = "grey") +
  scale_fill_viridis_d(option = "B", direction = -1) +
  labs(x = "Number of children", y = "") +
  guides(
    fill = guide_legend(
      "Times/month:",
      nrow = 1,
      reverse = TRUE,
      position = "top"
    )
  ) +
  theme_minimal(15)
barplot_SignalFreq
ggsave("Figures/barplot_SignalFreq.png", barplot_SignalFreq)


signal_subset <- utila_df[c(
  "ConflictFreqN",
  "Sex",
  "ChildAge",
  "SadFreqN",
  "CryFreqN",
  "TantrumFreqN",
  "SignalFreq",
  "SignalCost",
  "AlloparentingFreqN",
  "NeighborhoodQuality"
)]
names(signal_subset) <- shortform_dict[names(signal_subset)]

signal_corrplot <-
  signal_subset |>
  mutate(
    `Child sex` = as.numeric(`Child sex`)
  ) |>
  cor(use = "pairwise.complete.obs") |>
  ggcorrplot(
    type = "upper",
    hc.order = TRUE,
    hc.method = "ward.D",
    lab = TRUE,
    lab_col = "black",
    lab_size = 4.5,
    tl.cex = 16
  ) +
  scico::scale_fill_scico(
    palette = "vik",
    midpoint = 0,
    begin = .1,
    end = .9,
    limits = c(-1, 1)
  ) +
  guides(fill = guide_colorbar(title = "Correlation coefficients"))
signal_corrplot
ggsave("Figures/signal_corrplot.pdf", width = 12, height = 9)
ggsave("Figures/signal_corrplot.png", width = 12, height = 9)

barplot_conflict <-
  utila_df |>
  mutate(
    ConflictFreq2 = factor(freq_short[ConflictFreqOF], levels = c(freq_short))
  ) |>
  filter_at(vars(ConflictFreqOF), all_vars(!is.na(.))) |>
  ggplot(aes(x = ConflictFreq2, fill = Sex)) +
  geom_bar(position = "stack") +
  ylim(0, 150) +
  labs(title = "Conflict") +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 0.7) +
  labs(x = "Frequency of conflict (per month)", y = "Number of children") +
  theme_minimal(15)
barplot_conflict

barplot_alloparenting <-
  utila_df |>
  mutate(
    AlloparentingFreq02 = factor(
      freq_short[AlloparentingFreq0],
      levels = c(freq_short)
    )
  ) |>
  filter_at(vars(AlloparentingFreq0), all_vars(!is.na(.))) |>
  ggplot(aes(x = AlloparentingFreq02, fill = Sex)) +
  geom_bar(position = "stack") +
  ylim(0, 150) +
  labs(title = "Alloparenting") +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 0.7) +
  labs(
    x = "Frequency of child alloparenting (per month)",
    y = "Number of children"
  ) +
  theme_minimal(15)
barplot_alloparenting

plot_conflict_allo <-
  (barplot_conflict + barplot_alloparenting) +
  plot_layout(axes = 'collect_y', axis_titles = "collect_y", guides = 'collect')
plot_conflict_allo
ggsave("Figures/plot_conflict_allo.pdf", width = 12, height = 9)
ggsave("Figures/plot_conflict_allo.png", width = 12, height = 9)

cause_cluster_analysis <-
  causes |>
  mutate(
    across(ConflictFamily:StatusConcerns, as.numeric)
  ) |>
  dplyr::select(ConflictFamily:StatusConcerns) |>
  na.omit() |>
  pvclust(method.hclust = "ward.D2")

x <- cutree(cause_cluster_analysis$hclust, k = 3)

causes2 <-
  causes |>
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

cause_response0 <-
  utila_df |>
  left_join(causes2[c("uniqueID", "CauseType")]) |>
  dplyr::filter(
    !is.na(CauseType),
    !is.na(CaregiverResponse),
    !is.na(ChildAge)
  ) |>
  mutate(
    AgeCategory = ifelse(ChildAge <= 10, "Younger", "Older"), # Median age in this subsample
    AgeCategory = factor(AgeCategory, c("Younger", "Older")),
    CaregiverResponse = case_when(
      CaregiverResponse == -1 ~ "Negative",
      CaregiverResponse == 0 ~ "Neutral",
      CaregiverResponse == 1 ~ "Positive"
    ),
    CaregiverResponse = factor(
      CaregiverResponse,
      levels = c("Positive", "Neutral", "Negative")
    )
  )

barplot_CaregiverResponse <-
  cause_response0 |>
  summarize(
    N = n(),
    .by = c(CauseType, CaregiverResponse)
  ) |>
  mutate(
    CauseType = factor(
      CauseType,
      levels = c("Family Conflict", "Adversity", "Transgression", "Unknown")
    )
  ) |>
  ggplot(aes(N, CauseType, fill = CaregiverResponse)) +
  geom_col(position = "stack", color = "grey") +
  scale_fill_viridis_d(option = "B", direction = -1) +
  scale_y_discrete(limits = rev) +
  labs(x = "Number of children", y = "") +
  guides(
    fill = guide_legend(
      "Caregiver response:",
      nrow = 1,
      reverse = T,
      position = "top"
    )
  ) +
  theme_minimal(15)
barplot_CaregiverResponse
ggsave("Figures/barplot_CaregiverResponse.pdf", width = 8, height = 8)
ggsave("Figures/barplot_CaregiverResponse.png", width = 8, height = 8)

cause_response0 |>
  summarize(
    N = n(),
    .by = c(CauseType, CaregiverResponse)
  ) |>
  ggplot(aes(N, CaregiverResponse, fill = CauseType)) +
  geom_col(position = "stack", color = "grey") +
  scale_fill_viridis_d(option = "B", direction = -1) +
  scale_y_discrete(limits = rev) +
  labs(x = "Number of children", y = "") +
  guides(
    fill = guide_legend(
      "Caregiver response:",
      nrow = 1,
      reverse = T,
      position = "top"
    )
  ) +
  theme_minimal(15)

sfdfsum <-
  sfdf |>
  mutate(
    across(
      Sadness:Tantrums,
      \(x) factor(c(freq_short2[x]), levels = rev(unique(freq_short2)))
    )
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

signal_alluvial_plot <-
  ggplot(
    sfdfsum,
    aes(axis1 = Sadness, axis2 = Crying, axis3 = Tantrums, y = Freq)
  ) +
  geom_alluvium(aes(fill = Sadness), color = "black", show.legend = F) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(
    limits = c("Sadness", "Crying", "Tantrums"),
    expand = c(.2, .05)
  ) +
  ylab("Number of\nchildren") +
  theme_minimal(20) +
  theme(axis.title.y = element_text(angle = 0, hjust = 1)) +
  scale_fill_viridis_d()
signal_alluvial_plot

# Frequency by signal type ------------------------------------------------

e <-
  utila_df |>
  dplyr::filter(!is.na(CryFreqN) & !is.na(SadFreqN) & !is.na(TantrumFreqN)) |>
  pivot_longer(CryFreqN:TantrumFreqN, names_to = 'Signal', values_to = 'Frequency') |>
  mutate(
    Signal = factor(Signal, levels = c("SadFreqN", "CryFreqN", "TantrumFreqN"))
    )

m_freq_signal <- glmmTMB(Frequency ~ Signal * ChildAge + Signal * NeighborhoodQuality + Sex * ChildAge + (1|householdID/uniqueID), family = nbinom2, e)
# summary(m_freq_signal)
# plot(allEffects(m_freq_signal))
# plot(Effect(c("Signal", "ChildAge"), mod=m_freq_signal), cols = 3)
# plot(Effect(c("Signal", "NeighborhoodQuality"), mod=m_freq_signal), cols = 3)
# plot(Effect(c("ChildAge", "Sex"), mod=m_freq_signal), cols = 3)
