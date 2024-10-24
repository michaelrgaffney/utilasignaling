
library(glmnet)
library(glmnetUtils)
library(hagenutils)
library(patchwork)
library(skimr)
library(pvclust)
library(ggcorrplot)
library(tidymodels)
library(poissonreg)
library(marginaleffects)
library(easybgm)
library(tidygraph)
library(ggraph)
library(BDgraph)

# Regularized regression --------------------------------------------------

glmnet2 <- function(d, outcome, indices, alpha = 1, fam = 'quasipoisson'){

  reg <- ifelse(fam == 'quasipoisson', poisson_reg, linear_reg)
  fam <- ifelse(fam == 'quasipoisson', quasipoisson, gaussian)

    # Scale numeric vars by 2 SD, leave binary vars alone
  # d[indices] <- apply(d[indices], MARGIN = 2, FUN = \(x) if (length(unique(x)) == 2) return (x) else return(c(scale(x)/2)))
  d[indices] <- scale(d[indices]) # Scale all vars
  cvm <- list()
  for (i in 1:20){
    m <- cv.glmnet(
      as.matrix(d[indices]),
      d[[outcome]],
      alpha = alpha, # 0: ridge; 1: lasso
      relax = F,
      standardize = F,
      family = fam # MASS::negative.binomial(0.4) #quasipoisson()
    )
    cvm <- c(cvm, list(m$cvm))
  }

  cvm.mean <- rowMeans(matrix(list_c(cvm), nrow = length(cvm[[1]])))
  # return(tibble(lambda = m$lambda, cvm = rowMeans(cvm.mat)))

  lambda.min <- m$lambda[which.min(cvm.mean)]

  m2 <-
    reg(mixture = alpha, penalty = lambda.min) |>
    set_engine("glmnet", family = fam) |> # MASS::negative.binomial(0.4)
    fit_xy(d[indices], d[[outcome]])

  coefs <- tidy(m2)
  lmbda <- coefs$lambda[which.min(abs(coefs$lambda - lambda.min))]
  # print(lmbda)
  coefs <- coefs |> dplyr::filter(lambda == lmbda)
  names(coefs$estimate) <- shortform_dict[coefs$term]

  p <-
    hagenutils::ggdotchart(coefs$estimate[-1]) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    xlab("Coefficients") +
    ggtitle(shortform_dict[outcome])

  return(list(data = d, model = m2, coefs = coefs, coefplot = p))
}

signalparams <- expand_grid(Outcome = c("SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalCost"), alpha = c(0, 1))
names(signalparams$Outcome) <- str_c(signalparams$Outcome, signalparams$alpha)
signalparams$out <- map2(signalparams$Outcome, signalparams$alpha, \(x, y) glmnet2(SignalVars, x, 9:ncol(SignalVars), alpha = y), .progress = T)

conflictparams <- expand_grid(Outcome = c("ConflictFreqN"), alpha = c(0, 1))
names(conflictparams$Outcome) <- str_c(conflictparams$Outcome, conflictparams$alpha)
conflictparams$out <- map2(conflictparams$Outcome, conflictparams$alpha, \(x, y) glmnet2(SignalVars, x, 9:ncol(SignalVars), alpha = y), .progress = T)

plot_lasso <-
  signalparams$out$SadFreqN1$coefplot + signalparams$out$CryFreqN1$coefplot +
  signalparams$out$TantrumFreqN1$coefplot + signalparams$out$SignalFreq1$coefplot +
  signalparams$out$SignalCost1$coefplot + conflictparams$out$ConflictFreqN1$coefplot +
  plot_layout(axis_titles = 'collect', ncol = 2)
plot_lasso
ggsave("Figures/plot_lasso.pdf", plot_lasso, width = 12, height = 12)
ggsave("Figures/plot_lasso.svg", plot_lasso, width = 12, height = 12)

plotpredictions2 <- function(params, cond, title = NULL, plot = T){

  d <-
    map(params$out, \(x) plot_predictions(x$model, condition = cond, newdata = x$data, draw = F)) |>
    list_rbind(names_to = 'Outcome') |>
    dplyr::filter(str_ends(Outcome, '1') & Outcome != 'SignalCost1') |>
    mutate(Outcome = shortform_dict[str_remove(Outcome, '1')])

  if (!plot) return(d)

  ggplot(d, aes_string(cond, 'estimate', color = 'Outcome')) +
    geom_line(linewidth = 2) +
    ylim(0, 40) +
    labs(x = shortform_dict[cond], y = "Frequency per month") +
    theme_minimal(15)
}

# Protective factors (negative coefficients)
signal_effects_plot <-
  plotpredictions2(signalparams, 'ChildAge', 'Adjusted signaling frequencies') +
  plotpredictions2(signalparams, 'NeighborhoodQuality', 'Adjusted signaling frequencies') +
  plotpredictions2(signalparams, 'AdultsChildcare', 'Adjusted signaling frequencies') +
  plotpredictions2(signalparams, 'LogIncome', 'Adjusted signaling frequencies') +
  plotpredictions2(signalparams, 'EducationLevelYears', 'Adjusted signaling frequencies') +
  plotpredictions2(signalparams, 'MedicalProblemsMean', 'Adjusted signaling frequencies') +
  plot_layout(ncol = 2, guides = 'collect', axes = 'collect_y') + plot_annotation(title = "Adjusted signaling frequencies") &
  scico::scale_color_scico_d(palette = 'roma', begin = 0.1, end = 0.8)

signal_effects_plot &
  # scale_colour_viridis_d(option = 'B', begin = 0.2, end = 0.8)
  scico::scale_color_scico_d(palette = 'roma', begin = 0.1, end = 0.8)

ggsave("Figures/signal_effects_plot.pdf", signal_effects_plot, width = 12, height = 9)
ggsave("Figures/signal_effects_plot.svg", signal_effects_plot, width = 12, height = 9)


# Conflict effects plot

conflictpredictors0 <- c("MedicalProblemsMean", "EducationLevelYears", "ChildAge", "NeighborhoodQuality")
# names(conflictpredictors0) <- shortform_dict[conflictpredictors0]
conflictpredictors <- map(conflictpredictors0, \(x) c(x, "Sex"))

conflict_df <-
  map(conflictpredictors, \(x) plot_predictions(conflictparams$out$ConflictFreqN1$model, condition = x, newdata = conflictparams$out$ConflictFreqN1$data, draw = F) |> dplyr::select(estimate, all_of(x))) |>
  list_rbind(names_to = 'Predictor') |>
  dplyr::select(Predictor, estimate, Sex, all_of(conflictpredictors0)) |>
  pivot_longer(all_of(conflictpredictors0), values_drop_na = T) |>
  mutate(
    Sex = ifelse(as.numeric(Sex) < 2, "Female", "Male"),
    name = shortform_dict[name]
    )

conflict_effects_plot <-
  ggplot(conflict_df, aes(value, estimate, colour = Sex)) +
  geom_line(linewidth = 2) +
  scale_color_binary() +
  guides(color = guide_legend(reverse = T)) +
  ylim(0, NA) +
  labs(title = "Conflict frequency predictors", x = "Standardized value of predictor", y = "Adjusted frequency per month") +
  facet_wrap(~name, ncol = 2) +
  theme_bw(15)
conflict_effects_plot

ggsave("Figures/conflict_effects_plot.pdf", conflict_effects_plot, width = 12, height = 9)
ggsave("Figures/conflict_effects_plot.svg", conflict_effects_plot, width = 12, height = 9)

# Need to plot this by sex
# out <- plotpredictions2(signalparams, c('AlloparentingXsex', 'Sex'), plot = F)
#
# ggplot(out, aes(AlloparentingXsex, estimate, color = Sex)) +
#   geom_line(linewidth = 2) +
#   # ylim(0, 40) +
#   labs(y = "Frequency per month") +
#   facet_wrap(~Outcome) +
#   theme_minimal(15)

# PCA ---------------------------------------------------------------------

sv2 <-
  SignalVars[-c(1,2)] |>
  mutate(`Possession score` = rowSums(pick(LifestyleReality_1:LifestyleReality_8))) |>
  dplyr::select(-starts_with("Lifestyle"), -AlloparentingXsex)
sv2 <- set_names(sv2, shortform_dict[names(sv2)])

pca1 <- prcomp(sv2, scale. = T)
plot_loadings1 <- pca_loadings_plot(pca1, 1:2, reverse = 1) + theme_bw(20)
plot_biplot1 <- pca_biplot(pca1)

plot_pca <- plot_loadings1 + plot_biplot1 + plot_layout(widths = c(1,2))
plot_pca

ggsave("Figures/plot_pca.pdf", plot_pca, width = 20, height = 12)
ggsave("Figures/plot_pca.svg", plot_pca, width = 20, height = 12)

# pca2 <- prcomp(mdf2[-c(4:6)], scale. = T) # Remove composite signaling vars
# plot_loadings2 <- pca_loadings_plot(pca2, 1:2) # + theme(legend.position = 'top', legend.title = element_blank())
# plot_biplot2 <- pca_biplot(pca2) + theme_minimal(15)

# Graphical models ---------------------------------------------------------

# possessions <- SignalVars |> dplyr::select(starts_with("Lifestyle"))
#
# m <- bdgraph(data = possessions, method = "gcgm", g.prior = 0.1, iter = 10000, burnin = 7000)
# msum <- summary(m)

discrete_vars0 <-
  SignalVars |>
  summarize(discrete = across(everything(), \(x) length(unique(x)) <= 3)) |>
  dplyr::select(discrete)
discrete_vars <- as.integer(discrete_vars0$discrete)
names(discrete_vars) <- names(discrete_vars0$discrete)

# m2 <- bdgraph(data = SignalVars[-c(1,2)], method = "gcgm", not.cont = x2[-c(1,2)], g.prior = 0.1, cores = 'all', iter = 50000, burnin = 7000)
# m2sum <- summary(m2)

m3 <- easybgm(SignalVars[-c(1,2)], type = 'mixed', not_cont = discrete_vars[-c(1,2)], iter = 50000, package = 'BDgraph', g.prior = 0.05)
plot_network(m3, layout = "spring", vsize = 4, label.cex = 1)
plot_structure(m3, vsize = 4, label.cex = 1)
plot_structure_probabilities(m3, as_BF = F)
# plot_parameterHDI(m3) # need to have save = T

g3nodes <-
  as_tbl_graph(m3$structure) |>
  activate(nodes) |>
  data.frame() |>
  mutate(name = shortform_dict[name])

g3struct <-
  as_tbl_graph(m3$structure) |>
  activate(edges) |>
  data.frame() |>
  mutate(id = paste(to, from)) |>
  dplyr::select(-weight)

g3params <-
  as_tbl_graph(as.matrix(m3$parameters)) |>
  activate(edges) |>
  data.frame() |>
  mutate(
    id = paste(to, from),
    Correlation = ifelse(weight < 0, 'Negative', 'Positive'),
    weight2 = abs(weight),
    weight = 1 - weight2 # distance
    ) |>
  dplyr::select(id, weight, Correlation, weight2)

g3edges <- left_join(g3struct, g3params, by = 'id')
g3edges$id <- NULL

g3 <- tbl_graph(nodes = g3nodes, edges = g3edges, directed = F)

graph_plot <- function(g, weights = NULL, layout = 'stress', title = NULL){

  p <- ggraph(g, layout = layout)

  if (is.null(weights)) p <- p + geom_edge_link(aes(color = Correlation), linewidth = 2)
  else p <- p + geom_edge_link(aes(color = Correlation, linewidth = .data[[weights]]))

  p +
    geom_node_point(size = 5) +
    geom_node_text(aes(label = name), repel = T, max.overlaps = Inf) +
    scale_edge_color_manual(values = viridisLite::magma(2, begin = 0.2, end = 0.8)) +
    ggtitle(title) +
    theme_graph(base_family="sans")
}

# Plot full graph
plot_full_graph <- graph_plot(g3, layout = 'stress', title = "Bayesian graphical model", weights = 'weight2')
plot_full_graph

# Three options for MST:
# algorithm: unweighted, weighted (prim), weights = 1
# weights = 1: rep(1, times = igraph::gsize(g3))

mst_weighted <- igraph::mst(g3, algorithm = 'prim')
plot_mst <- graph_plot(mst_weighted, weights = 'weight2', title = 'MST (weighted)')
plot_mst

ggsave('Figures/plot_mst.svg', plot_mst, width = 12, height = 9)
ggsave("Figures/plot_mst.pdf", plot_mst, width = 12, height = 9)

mst_unweighted <- igraph::mst(g3, algorithm = 'unweighted')
plot_mst_uw <- graph_plot(mst_unweighted, title = 'MST (unweighted)')
plot_mst_uw

ggsave('Figures/plot_mst_uw.svg', plot_mst_uw, width = 12, height = 9)
ggsave("Figures/plot_mst_uw.pdf", plot_mst_uw, width = 12, height = 9)

mst_unit <- igraph::mst(g3, algorithm = 'prim', weights = rep(1, times = igraph::gsize(g3)))
plot_mst_unit <- graph_plot(mst_unit, title = 'MST (unit weights)')
plot_mst_unit

ggsave('Figures/plot_mst_unit.svg', plot_mst_unit, width = 12, height = 9)
ggsave("Figures/plot_mst_unit.pdf", plot_mst_unit, width = 12, height = 9)


# annotate("rect", xmin = 12.4, xmax = 14.5, ymin = -0.75, ymax = 3.75, colour = 'red', fill = NA, linewidth = 2) +


# Anthropometrics ---------------------------------------------------------

m_height <- lm(HeightMean_2024 ~ HeightMean_2023 + Sex*ChildAge, SignalVarsAnthro)
SignalVarsAnthro$heightR <- m_height$residuals

m_weight <- lm(WeightMean_KG_2024 ~ WeightMean_KG_2023, SignalVarsAnthro)
SignalVarsAnthro$weightR <- m_weight$residuals

SignalVarsAnthro2 <- dplyr::select(SignalVarsAnthro, -ConflictFreqN)

anthroparams <- expand_grid(Outcome = c("heightR", "weightR"), alpha = c(0, 0.5, 1))
names(anthroparams$Outcome) <- str_c(anthroparams$Outcome, anthroparams$alpha)
anthroparams$out <- map2(anthroparams$Outcome, anthroparams$alpha, \(x, y) glmnet2(SignalVarsAnthro2, x, 3:41, alpha = y, fam = 'gaussian'), .progress = T)

# Predicting signals

# anthro2params <- expand_grid(Outcome = c("SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalCost"), alpha = c(0, 1))
# names(anthro2params$Outcome) <- str_c(anthro2params$Outcome, anthro2params$alpha)
# anthro2params$out <- map2(anthro2params$Outcome, anthro2params$alpha, \(x, y) glmnet2(SignalVarsAnthro2, x, 9:ncol(SignalVarsAnthro2), alpha = y), .progress = T)

# conflictAnthroparams <- expand_grid(Outcome = c("ConflictFreqN"), alpha = c(0, 0.5, 1))
# names(conflictAnthroparams$Outcome) <- str_c(conflictAnthroparams$Outcome, conflictAnthroparams$alpha)
# conflictAnthroparams$out <- map2(conflictAnthroparams$Outcome, conflictAnthroparams$alpha, \(x, y) glmnet2(SignalVarsAnthro2, x, 9:ncol(SignalVarsAnthro2), alpha = y), .progress = T)


# Ordinal regressions -----------------------------------------------------

library(ordinalNet)

# For binary predictors
ordinal_plot <- function(fit, predictor, data, title){
  x <- predict(fit, type = 'response')
  d <- data |>
    mutate(
      Negative = x[,1],
      Neutral = x[,2],
      Positive = x[,3],
      {{predictor}} := ifelse({{predictor}} < 0, 0, 1)
    ) |>
    dplyr::select(
      {{predictor}}, Negative:Positive
    ) |>
    pivot_longer(Negative:Positive)
  ggplot(d, aes({{predictor}}, value, colour = name)) +
    geom_count(position = position_dodge(width = 0.1)) +
    geom_smooth(method='lm', se = F, position = position_dodge(width = 0.1)) +
    scale_x_continuous(breaks = c(0, 1)) +
    scale_color_viridis_d(option = 'B', end = 0.8) +
    guides(color = guide_legend("Response", reverse = T)) +
    ylim(0, NA) +
    labs(title = title, x = shortform_dict[deparse(substitute(predictor))], y = "Probability")+
    theme_minimal(15)
}

# For continuous predictors
ordinal_plot2 <- function(fit, predictor, data, title){
  x <- predict(fit, type = 'response')
  d <- data |>
    mutate(
      Less = x[,1],
      Equal = x[,2],
      More = x[,3]
      # {{predictor}} := ifelse({{predictor}} < 0, 0, 1)
    ) |>
    summarise(
      Less = mean(Less),
      Equal = mean(Equal),
      More = mean(More),
      .by = {{predictor}}
      # {{predictor}} := ifelse({{predictor}} < 0, 0, 1)
    ) |>
    pivot_longer(Less:More) |>
    mutate(name = factor(name, levels = c("Less", "Equal", "More")))

  ggplot(d, aes({{predictor}}, value, colour = name)) +
    # geom_smooth(se = F, linewidth = 2) +
    geom_line(linewidth = 2) +
    scale_color_viridis_d(option = 'B', end = 0.8) +
    guides(color = guide_legend("Response", reverse = T)) +
    ylim(0, NA) +
    labs(title = title, x = shortform_dict[deparse(substitute(predictor))], y = "Probability")+
    theme_minimal(15)
}

# Perceived need

SignalVars4 <-
  left_join(SignalVars, d2[c("householdID", "childHHid", "RelativeNeed3")]) |>
  relocate(RelativeNeed3, .after = 'childHHid') |>
  mutate(
    across(-c(1:3), as.numeric),
    across(SadFreqN:AlloparentingXsex, \(x) c(scale(x)))
  ) |>
  na.omit()

out <- ordinalNetCV(
  as.matrix(SignalVars4[-c(1:3)]),
  SignalVars4[[3]],
  standardize = F,
  alpha = 1.0,
  family = "cumulative",
  link = "logit",
  lambdaMinRatio = 1e-04,
  printProgress = T
)
summary(out)
colMeans(summary(out))
# coef(out$fit, matrix = TRUE, whichLambda = 1)
plot_need_coefs <- ggdotchart(coef(out$fit)[-c(1:2)])
plot_need_coefs
plot_need_age <- ordinal_plot2(out$fit, ChildAge, SignalVars4, 'Relative need')
plot_need_sad <- ordinal_plot2(out$fit, SadFreqN, SignalVars4, 'Relative need')

plot_need_combined <- plot_need_age / plot_need_sad + ggtitle("") + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A")
ggsave("Figures/plot_need_combined.pdf", plot_need_combined, width = 12, height = 12)
ggsave("Figures/plot_need_combined.svg", plot_need_combined, width = 12, height = 12)

# Relative investment

SignalVars5 <-
  left_join(SignalVars, d2[c("householdID", "childHHid", "RelativeMaternalInvestment2", "RelativeNeed3")]) |>
  relocate(RelativeMaternalInvestment2, .after = 'childHHid') |>
  mutate(
    across(-c(1:3), as.numeric),
    across(SadFreqN:RelativeNeed3, \(x) c(scale(x)))
  ) |>
  na.omit()

out <- ordinalNetCV(
  as.matrix(SignalVars5[-c(1:3)]),
  SignalVars5[[3]],
  standardize = F,
  alpha = 1.0,
  family = "cumulative",
  link = "logit",
  lambdaMinRatio = 1e-04,
  printProgress = T
)
summary(out)
colMeans(summary(out))
# coef(out$fit, matrix = TRUE, whichLambda = 1)
plot_invest_coefs <- ggdotchart(coef(out$fit)[-c(1:2)])
plot_invest_coefs
plot_invest_need <- ordinal_plot2(out$fit, RelativeNeed3, SignalVars5, 'Relative investment')
ggsave("Figures/plot_invest_need.pdf", plot_invest_need, width = 12, height = 12)
ggsave("Figures/plot_invest_need.svg", plot_invest_need, width = 12, height = 12)

# Caregiver response

causematrix <-
  causes |>
  mutate(
    across(Family:StatusConcerns, \(x) ifelse(x == "Maybe" | x == 'maybe', "0", x)),
    across(Family:StatusConcerns, as.numeric)
  ) |>
  dplyr::select(householdID, childHHid, Family:StatusConcerns) |>
  na.omit()

cause_count <- map_int(causematrix[-c(1:5)], \(x) sum(x, na.rm = T))
cause_vars <- names(cause_count[cause_count > 4])

SignalVars3 <-
  left_join(SignalVars, d2[c("householdID", "childHHid", "CaregiverResponse")]) |>
  left_join(causematrix[c("householdID", "childHHid", cause_vars)]) |>
  relocate(CaregiverResponse, .after = childHHid) |>
  mutate(
    across(-c(1:3), as.numeric),
    across(SadFreqN:StatusConcerns, \(x) c(scale(x)))
    ) |>
  na.omit()

out <- ordinalNetCV(
  as.matrix(SignalVars3[-c(1:3)]),
  SignalVars3[[3]],
  standardize = F,
  alpha = 1.0,
  family = "cumulative",
  link = "logit",
  lambdaMinRatio = 1e-04,
  printProgress = T
)
summary(out)
colMeans(summary(out))
# coef(out$fit, matrix = TRUE, whichLambda = 1)
plot_caregiverresponse_coefs <- ggdotchart(coef(out$fit)[-c(1:2)])
plot_caregiver_pain <- ordinal_plot(out$fit, DiscomfortPainInjuryIllness, data = SignalVars3, title = 'Caregiver response')
plot_caregiver_punish <- ordinal_plot(out$fit, Punishment, data = SignalVars3, title = 'Caregiver response')
plot_caregiver_response_combined <- plot_caregiver_punish / plot_caregiver_pain + ggtitle("") + plot_layout(guides = 'collect')

ggsave("Figures/plot_caregiver_response_combined.pdf", plot_caregiver_response_combined, width = 12, height = 12)
ggsave("Figures/plot_caregiver_response_combined.svg", plot_caregiver_response_combined, width = 12, height = 12)

# newdat <- datagrid(NeighborhoodQuality = seq(-2.5, 1.5, 0.1), newdata = SignalVars3[-c(1:3)])

# sv3 <-
#   SignalVars3 |>
#   mutate(
#     P1 = x[,1],
#     P2 = x[,2],
#     P3 = x[,3]
#   ) |>
#   dplyr::select(
#     Punishment, P1:P3
#   ) |>
#   pivot_longer(P1:P3)
#
# ggplot(sv3, aes(Punishment, value, colour = name)) +
#   geom_count(position = position_dodge(width = 0.1)) +
#   geom_smooth(method='lm', se = F, position = position_dodge(width = 0.1)) +
#   theme_minimal(15)

# Hacking polr
# AIC.polr <- function(x, k = 2){
#   dev <- x$deviance
#   nparams <- x$edf
#   dev + k*nparams
# }
#
# library(MASS)
# m <- polr(CaregiverResponse ~ ., data = SignalVars3[-c(1,2)])
# m2 <- stepAIC(m)

# Old code ----------------------------------------------------------------

# space2newline <- function(v) str_replace_all(v, " ", "\n")
#
# plotpredictions <- function(params, outcome, condition, title = NULL){
#   model <- params$out[[outcome]]$model
#   data <- params$out[[outcome]]$data
#   if(is.null(title)) title <- outcome
#   if(length(condition) == 2){
#     if(condition[2] == 'Sex') lbls = c("Female", "Male") else lbls = c("No", "Yes")
#     scl <- scale_color_binary(labels = lbls, name = shortform_dict[condition[2]] |> space2newline())
#   } else {
#     scl <- scale_color_manual(values = 'black')
#   }
#   p <-
#     plot_predictions(model, condition = condition, newdata = data) +
#     scl +
#     ylim(0, 30) +
#     labs(title = title, x = shortform_dict[condition[1]]) +
#     theme_minimal(15)
#   p$layers[[1]]$aes_params$linewidth <- 2
#   return(p)
# }
#
# plotpredictions(signalparams, "CryFreqN1", c("ChildAge", "LifestyleReality_5"), "Crying")
# plotpredictions(conflictparams, "ConflictFreqN1", c("ChildAge", "Sex"))

# Other stuff

# out <- skim(modeldf)
# nms <- out$skim_variable[out$skim_type == 'numeric' & out$complete_rate > 0.9]
# nms <- c(nms, "Sex", "UserLanguage", "ImmigrateUtila", "PartnerStatus", "OldestChild")
# mdf2 <-
#   modeldf |>
#   dplyr::select(
#     all_of(nms),
#     # -householdID,
#     -childHHid,
#     -NegativeResponse,
#     -PositiveResponse,
#     -IncomeCategoryN,
#     -ConflictFreqN,
#     -OtherChildrenHH,
#     -contains('IllnessSusceptibility'),
#
#     # New omits
#     -OlderGirls,
#     -YoungerKids,
#     -FoodSecurity,
#     -UserLanguage
#
#   ) |>
#   na.omit() |>
#   mutate(
#     Sex = ifelse(Sex == "Female", 0, 1),
#     # UserLanguage = ifelse(UserLanguage == "EN", 0, 1),
#     ImmigrateUtila = ifelse(ImmigrateUtila == "No", 0, 1),
#     PartnerStatus = ifelse(PartnerStatus == "Unpartnered", 0, 1),
#     across(-c(1:6), \(x) c(scale(x))),
#     AlloparentingXsex = AlloparentingFreqN * Sex,
#     # ChildAgeXsex = ChildAge * Sex,
#     # OldestXsex = OldestChild * Sex
#   )

# Final models ------------------------------------------------------------

# frmla <- paste("CryFreqN ~", paste(names(mdf2)[c(8:10, 12:21)], collapse = ' + '), " + AlloparentingFreqN*Sex", " + (1|householdID)")
# m <- glmmTMB(as.formula(frmla), data = mdf2, family = nbinom2)
# summary(m)
#
# glmmTMBformulas <- c(std_formulas, bodyfat_formulas, conflict_formula)
# names(glmmTMBformulas) <- c(paste0(signals, "illness"), paste0(signals, "bodyfat"), "ConflictFreqN")
#
# glmmTMBmodels <- tibble(
#   Outcome = c(signals, signals, "ConflictFreqN"),
#   Model = map(glmmTMBformulas, \(f) glmmTMB(as.formula(f), family = nbinom2, data = d2)),
#   AIC = map_dbl(Model, AIC)
# )

# Cross validate both lambda and alpha
# nvar <- ncol(mdf2)
# m <- cva.glmnet(
#   as.matrix(mdf2[7:nvar]),
#   mdf2[["TantrumFreqN"]],
#   standardize = T,
#   family =  quasipoisson() # MASS::negative.binomial(3) #
# )
# minlossplot(m, cv.type = 'min')
# coefs <- coef(m$modlist[[1]], s = 'lambda.min')[-1,]
# ggdotchart(coefs) + geom_vline(xintercept = 0, linetype = 'dotted')


# # Sad: alpha = 0
# # Cry: alpha = 1
# # Tantrum: alpha = 0

# # Now fit with those alpha values

# glmnetSignal <- function(signal, alpha = 0, relax = T, s = "lambda.min"){
#   m <- cv.glmnet(
#     as.matrix(mdf2[8:nvar]),
#     mdf2[[signal]],
#     alpha = alpha, # alpha = 0: ridge; alpha = 1: lasso
#     relax = relax,
#     standardize = T,
#     family =  quasipoisson() # MASS::negative.binomial(3) #
#   )

#   coefs <- coef(m, s=s)[-1,]
#   ggdotchart(coefs) +
#     geom_vline(xintercept = 0, linetype = 'dotted') +
#     labs(title = signal, subtitle = str_glue("alpha = {alpha}; relax = {relax}; s = {s}"))
# }

# p1 <- glmnetSignal("SadFreqN", alpha = 0.0, relax = T, s = "lambda.min")
# p2 <- glmnetSignal("CryFreqN", alpha = 1.0, relax = F, s = "lambda.min")
# p3 <- glmnetSignal("TantrumFreqN", alpha = 0.0, relax = T, s = "lambda.min")

# p1/p2/p3 + plot_layout(axes = 'collect')

# Regularized regression --------------------------------------------------

# glmnetSignal <- function(signal, alpha = 0, relax = T, s = "lambda.min"){
#   m <- cv.glmnet(
#     as.matrix(mdf2[8:nvar]),
#     mdf2[[signal]],
#     alpha = alpha, # alpha = 0: ridge; alpha = 1: lasso
#     relax = relax,
#     standardize = T,
#     family =  quasipoisson() # MASS::negative.binomial(3) #
#   )
#
#   coefs <- coef(m, s=s)[-1,]
#   ggdotchart(coefs) +
#     geom_vline(xintercept = 0, linetype = 'dotted') +
#     labs(title = signal, subtitle = str_glue("alpha = {alpha}; relax = {relax}; s = {s}"))
# }
#
# p1 <- glmnetSignal("SadFreqN", alpha = 0.0, relax = T, s = "lambda.min")
# p2 <- glmnetSignal("CryFreqN", alpha = 1.0, relax = F, s = "lambda.min")
# p3 <- glmnetSignal("TantrumFreqN", alpha = 0.0, relax = T, s = "lambda.min")
#
# p1/p2/p3 + plot_layout(axes = 'collect')


# Plot predictions

# library(tidymodels)

# glmnet2 <- function(d, outcome, indices, alpha = 1){
#   # Scale numeric vars by 2 SD, leave binary vars alone
#   # d[indices] <- apply(d[indices], MARGIN = 2, FUN = \(x) if (length(unique(x)) == 2) return (x) else return(c(scale(x)/2)))
#   d[indices] <- scale(d[indices])
#   lambdas <- c()
#   for (i in 1:20){
#     m <- cv.glmnet(
#       as.matrix(d[indices]),
#       d[[outcome]],
#       alpha = alpha, # 0: ridge; 1: lasso
#       relax = F,
#       standardize = F,
#       family = quasipoisson() # MASS::negative.binomial(0.4) #quasipoisson()
#     )
#     lambdas <- c(lambdas, m$lambda.min)
#   }
#   lambda.mean <- mean(lambdas)
#   # print(lambda.mean)
#
#   m2 <-
#     poisson_reg(mixture = alpha, penalty = lambda.mean) %>%
#     set_engine("glmnet", family = quasipoisson) %>% # MASS::negative.binomial(0.4)
#     fit_xy(d[indices], d[[outcome]])
#
#   coefs <- tidy(m2)
#   lmbda <- coefs$lambda[which.min(abs(coefs$lambda - lambda.mean))]
#   # print(lmbda)
#   coefs <- coefs |> dplyr::filter(lambda == lmbda)
#   names(coefs$estimate) <- shortform_dict[coefs$term]
#
#   p <-
#     hagenutils::ggdotchart(coefs$estimate[-1]) +
#     geom_vline(xintercept = 0, linetype = 'dotted') +
#     ggtitle(outcome)
#
#   return(list(data = d, åmodel = m2, coefs = coefs, coefplot = p))
# }
#
# params <- expand_grid(Outcome = c("SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalCost"), alpha = c(0, 1))
# names(params$Outcome) <- str_c(params$Outcome, params$alpha)
# params$out <- map2(params$Outcome, params$alpha, \(x, y) glmnet2(SignalVars, x, 9:ncol(SignalVars), alpha = y), .progress = T)
#
# (params$out$SadFreqN1$coefplot + params$out$CryFreqN1$coefplot) / (params$out$TantrumFreqN1$coefplot + params$out$SignalFreq1$coefplot) / (params$out$SignalCost1$coefplot + plot_spacer())
# plot_predictions(params$out$CryFreqN1$model, condition = c('ChildAge', 'Sex'), newdata = SignalVars)
#
# SignalVars2 <-
#   SignalVars |>
#   mutate(
#     AlloparentingXage = c(scale(AlloparentingFreqN) * scale(ChildAge))
#   )
# params <- expand_grid(Outcome = c("ConflictFreqN"), alpha = c(0, 1))
# names(params$Outcome) <- str_c(params$Outcome, params$alpha)
# params$out <- map2(params$Outcome, params$alpha, \(x, y) glmnet2(SignalVars2, x, 9:ncol(SignalVars2), alpha = y), .progress = T)
# params$out$SadFreqN1$coefplot

# Clustering --------------------------------------------------------------

# mdf3 <-
#   mdf2 |>
#   dplyr::select(-matches("Older|Oldest|Xsex|Max|Kids"))
#
# out <- pvclust(mdf3, method.hclust = 'ward.D2', method.dist = 'abscor')
# plot(out)
#
# var_clusters <- cutree(out$hclust, k = 3)
# signal_vars <- names(var_clusters[var_clusters == 1])
# household_vars <- names(var_clusters[var_clusters == 2])
# neighborhood_vars <- names(var_clusters[var_clusters == 3])

# Theory based var groups

# signal_vars <- c("CryFreqN", "SadFreqN", "TantrumFreqN", "SignalFreq", "SignalCost")
# household_vars <- c(
#   "ChildAge",
#   "number_adults",
#   "ConflictFreqN",
#   "AlloparentingFreqN",
#   "NumberOfChildren",
#   "CaregiverAge",
#   "Sex",
#   "PartnerStatus",
#   "LogIncome",
#   "EducationLevelYears",
#   "FoodSecurity",
#   "HouseQuality",
#   "UserLanguage",
#   "ImmigrateUtila"
# )
#
# neighborhood_vars <- c("NeighborhoodQuality", "Neighborhood2")

# Correlation matrices ----------------------------------------------------

# ggcorrplot(
#   cor(
#     mdf2 |> dplyr::select(-matches("Xsex|Old|Kids|Max")
#     ),
#     use = 'pairwise.complete.obs'),
#   hc.order = T,
#   hc.method = 'centroid'
# ) +
#   scico::scale_fill_scico(palette = 'vik', midpoint = 0)

# experimental

# sv3 <-
#   SignalVars[-c(1,2)] |>
#   mutate(`Possession score` = rowSums(pick(LifestyleReality_1:LifestyleReality_8))) |>
#   dplyr::select(-starts_with("Lifestyle"))
#
# signalparams <- expand_grid(Outcome = c("SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalCost"), alpha = c(0, 1))
# names(signalparams$Outcome) <- str_c(signalparams$Outcome, signalparams$alpha)
# signalparams$out <- map2(signalparams$Outcome, signalparams$alpha, \(x, y) glmnet2(sv3, x, 7:ncol(sv3), alpha = y), .progress = T)
#
# conflictparams <- expand_grid(Outcome = c("ConflictFreqN"), alpha = c(0, 1))
# names(conflictparams$Outcome) <- str_c(conflictparams$Outcome, conflictparams$alpha)
# conflictparams$out <- map2(conflictparams$Outcome, conflictparams$alpha, \(x, y) glmnet2(sv3, x, 7:ncol(sv3), alpha = y), .progress = T)

