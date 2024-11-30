library(utiladata2023) # data package
library(tidyverse)
library(modelsummary)
library(labelled)
library(ggcorrplot)
library(pvclust)
library(ggalluvial)
library(hagenutils)
library(knitr)
library(kableExtra)

source("dataprep.R")

## main paper plots --------------------------------------------------------

signalfreqdf <-
  d2 |>
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

signaldf_long <-
  sfdf |>
  pivot_longer(cols = Sadness:Tantrums, names_to = "Signal", values_to = "Freq")  |>
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
  guides(fill = guide_legend("Times per month:", nrow = 1, reverse = TRUE, position = "top")) +
  theme_minimal(15)
barplot_SignalFreq


signal_subset <- d2[c("ConflictFreqN", "Sex", "ChildAge", "SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalCost", "AlloparentingFreqN", "NeighborhoodQuality")]
names(signal_subset) <- shortform_dict[names(signal_subset)]

signal_corrplot <-
  signal_subset |>
  mutate(
    `Child sex` = as.numeric(`Child sex`)
  ) |>
  cor( use = "pairwise.complete.obs") |>
  ggcorrplot(hc.order = TRUE, hc.method = "ward.D",lab = TRUE, lab_col = "black", lab_size = 4.5) +
  scico::scale_fill_scico(palette = "vik", midpoint = 0, begin = .1, end = .9, limits = c(-1, 1)) +
  guides(fill = guide_colorbar(title = "Correlation coefficients"))
signal_corrplot


barplot_conflict <-
  modeldf |>
  mutate(
    ConflictFreq2 = factor(freq_short[ConflictFreqOF], levels = c(freq_short))
  )|>
  filter_at(vars(ConflictFreqOF),all_vars(!is.na(.))) |>
  ggplot(aes(x=ConflictFreq2, fill= Sex)) +
  geom_bar(position = "stack") +
  labs(title = "Conflict") +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 0.7) +
  labs(x = "Frequency of conflict (per month)", y = "Number of children") +
  theme_minimal(15)
barplot_conflict

barplot_alloparenting <-
  modeldf |>
  mutate(
    AlloparentingFreq02 = factor(freq_short[AlloparentingFreq0], levels = c(freq_short))
  )|>
  filter_at(vars(AlloparentingFreq0),all_vars(!is.na(.))) |>
  ggplot(aes(x= AlloparentingFreq02, fill = Sex)) +
  geom_bar(position = "stack") +
  labs(title = "Alloparenting") +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 0.7) +
  labs(x = "Frequency of child alloparenting (per month)", y = "Number of children") +
  theme_minimal(15)
barplot_alloparenting

cause_cluster_analysis <-
  causes |>
  mutate(
    across(ConflictFamily:StatusConcerns, as.numeric)
  ) |>
  dplyr::select(ConflictFamily:StatusConcerns) |>
  na.omit() |>
  pvclust(method.hclust = "ward.D2")

x <- cutree(cause_cluster_analysis $hclust, k = 3)

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
  d2 |>
  left_join(causes2[c("uniqueID", "CauseType")]) |>
  dplyr::filter(!is.na(CauseType), !is.na(CaregiverResponse), !is.na(ChildAge)) |>
  mutate(
    AgeCategory = ifelse(ChildAge <= 10, "Younger", "Older"), # Median age in this subsample
    AgeCategory = factor(AgeCategory, c("Younger", "Older")),
    CaregiverResponse = case_when(
      CaregiverResponse == -1 ~ "Negative",
      CaregiverResponse ==  0 ~ "Neutral",
      CaregiverResponse ==  1 ~ "Positive"
    ),
    CaregiverResponse = factor(CaregiverResponse, levels = c("Positive", "Neutral", "Negative"))
  )

barplot_CaregiverResponse <-
  cause_response0 |>
  summarize(
    N = n(),
    .by = c(CauseType, CaregiverResponse)
  ) |>
  ggplot(aes(N, CauseType, fill = CaregiverResponse)) +
  geom_col(position = "stack", color = "grey") +
  scale_fill_viridis_d(option = "B", direction = -1) +
  scale_y_discrete(limits = rev) +
  labs(x = "Number of children", y = "") +
  guides(fill = guide_legend("Caregiver response:", nrow = 1, reverse = T, position = "top")) +
  theme_minimal(15)
barplot_CaregiverResponse

sfdfsum <-
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

signal_alluvial_plot <-
  ggplot(sfdfsum, aes(axis1 = Sadness, axis2 = Crying, axis3 = Tantrums, y = Freq)) +
  geom_alluvium(aes(fill = Sadness), color = "black", show.legend = F) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Sadness", "Crying", "Tantrums"), expand = c(.2, .05)) +
  ylab("Number of\nchildren") +
  theme_minimal(20) +
  theme(axis.title.y = element_text(angle = 0, hjust = 1)) +
  scale_fill_viridis_d()
signal_alluvial_plot

# PCA ---------------------------------------------------------------------
e <- mdf2 |>
  dplyr::select(
    - SignalFreq,
    - SignalCost,
    - SignalFreqMax,
    - OldestChild,
    - OlderGirls,
    - PartnerStatus,
    - CaregiverAge,
    - YoungerKids,
    - contains("Xsex")
  )

e <- set_names(e, shortform_dict[names(e)])
# names(e) <- shortform_dict[names(e)]

m <- prcomp(e, scale. = T) # Remove composite signaling vars
plot_loadings <- pca_loadings_plot(m, 1:2) # + theme(legend.position = 'top', legend.title = element_blank())
plot_biplot <- pca_biplot(m) + theme_minimal(15)

plot_pca <- plot_loadings + plot_biplot + plot_layout(widths = c(1,2)) +
  plot_annotation(tag_levels = "A")
ggsave("Figures/plot_pca.svg", plot_pca)
ggsave("Figures/plot_pca.pdf", plot_pca)

