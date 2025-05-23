
# Regularized regression --------------------------------------------------

glmnet2 <- function(d, outcome, indices, alpha = 1, fam = 'quasipoisson'){

  reg <- ifelse(fam == 'quasipoisson', poisson_reg, linear_reg)
  fam <- ifelse(fam == 'quasipoisson', quasipoisson, gaussian)

  d[indices] <- scale(d[indices]) # Scale all predictor vars
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
    ggdotchart(coefs$estimate[-1]) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    xlab("Coefficients") +
    ggtitle(shortform_dict[outcome])

  return(list(data = d, lambda.min = lambda.min, model = m2, coefs = coefs, coefplot = p))
}

sv <-
  SignalVars |>
  dplyr::select(where(\(x) sd(x, na.rm = T) > 0))

signalparams <- expand_grid(Outcome = c("SadFreqN", "CryFreqN", "TantrumFreqN", "SignalFreq", "SignalCost"), alpha = c(0, 1))
names(signalparams$Outcome) <- str_c(signalparams$Outcome, signalparams$alpha)
signalparams$out <- map2(signalparams$Outcome, signalparams$alpha, \(x, y) glmnet2(sv, x, 9:ncol(sv), alpha = y), .progress = T)

conflictparams <- expand_grid(Outcome = c("ConflictFreqN"), alpha = c(0, 1))
names(conflictparams$Outcome) <- str_c(conflictparams$Outcome, conflictparams$alpha)
conflictparams$out <- map2(conflictparams$Outcome, conflictparams$alpha, \(x, y) glmnet2(sv, x, 9:ncol(sv), alpha = y), .progress = T)

plot_lasso <-
  signalparams$out$SadFreqN1$coefplot + signalparams$out$CryFreqN1$coefplot +
  signalparams$out$TantrumFreqN1$coefplot + signalparams$out$SignalFreq1$coefplot +
  signalparams$out$SignalCost1$coefplot + conflictparams$out$ConflictFreqN1$coefplot +
  plot_layout(axis_titles = 'collect', ncol = 2)
plot_lasso
ggsave("Figures/plot_lasso.pdf", plot_lasso, width = 12, height = 12)
ggsave("Figures/plot_lasso.png", plot_lasso, width = 12, height = 12)

plotpredictions2 <- function(params, cond, title = NULL, plot = T){

  d <-
    map(params$out, \(x) plot_predictions(x$model, condition = cond, newdata = x$data, draw = F)) |>
    list_rbind(names_to = 'Outcome') |>
    dplyr::filter(str_ends(Outcome, '1') & Outcome != 'SignalCost1') |>
    mutate(
      Outcome = shortform_dict[str_remove(Outcome, '1')],
      Outcome = factor(Outcome, levels = c("Signal freq.", "Tantrum freq.", "Crying freq.", "Sadness freq."))
      )

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
ggsave("Figures/signal_effects_plot.png", signal_effects_plot, width = 12, height = 9)


# Alloparenting special case (SignalCost)

e <- SignalVars |>
  dplyr::select(-householdID, -childHHid, -c(SadFreqN:SignalFreq), -ConflictFreqN, -AlloparentingXsex) |>
  dplyr::select(where(\(x) sd(x, na.rm = T) > 0)) |>
  relocate(AlloparentingFreqN, .after = SignalCost)
e[-1] <- scale(e[-1])

f <- paste(names(e[-c(1,2)]), collapse = " + ")
f <- as.formula(paste("SignalCost ~ ", "AlloparentingFreqN*Sex +", f))

m_alloparent <-
  poisson_reg(mixture = 1, penalty = signalparams$out$SignalCost1$lambda.min) |>
  set_engine("glmnet", family = quasipoisson) |>
  fit(f, data = e)

plot_alloparenting_cost <-
  plot_predictions(m_alloparent, condition = c("AlloparentingFreqN", "Sex"), newdata = e) +
  scale_color_binary(labels = c("Female", "Male")) +
  ylim(0, NA) +
  labs(x = "Alloparenting Frequency (standardized)", y = "Signal cost") +
  theme_minimal(15)

plot_alloparenting_cost$layers[[1]]$aes_params$linewidth <- 2
plot_alloparenting_cost

# Alloparenting special case (TantrumFreq)

e2 <- SignalVars |>
  dplyr::select(-householdID, -childHHid, -c(SadFreqN:CryFreqN), -c(SignalFreq:SignalCost), -ConflictFreqN, -AlloparentingXsex) |>
  dplyr::select(where(\(x) sd(x, na.rm = T) > 0)) |>
  relocate(AlloparentingFreqN, .after = TantrumFreqN)
e2[-1] <- scale(e2[-1])

f2 <- paste(names(e2[-c(1,2)]), collapse = " + ")
f2 <- as.formula(paste("TantrumFreqN ~ ", "AlloparentingFreqN*Sex +", f2))

m_alloparent2 <-
  poisson_reg(mixture = 1, penalty = signalparams$out$TantrumFreqN1$lambda.min) |>
  set_engine("glmnet", family = quasipoisson) |>
  fit(f2, data = e2)

plot_alloparenting_cost_tantrum <-
  plot_predictions(m_alloparent2, condition = c("AlloparentingFreqN", "Sex"), newdata = e2) +
  scale_color_binary(labels = c("Female", "Male")) +
  ylim(0, NA) +
  labs(x = "Alloparenting Frequency (standardized)", y = "Tantrums (per month)") +
  theme_minimal(15)

plot_alloparenting_cost_tantrum$layers[[1]]$aes_params$linewidth <- 2
plot_alloparenting_cost_tantrum

plot_alloparentingXsex <- plot_alloparenting_cost / plot_alloparenting_cost_tantrum +
  plot_layout(guides = "collect", axes = "collect") &
  theme_bw(15)
plot_alloparentingXsex

ggsave("Figures/plot_alloparentingXsex.pdf", plot_alloparentingXsex, width = 12, height = 9)
ggsave("Figures/plot_alloparentingXsex.png", plot_alloparentingXsex, width = 12, height = 9)

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
ggsave("Figures/conflict_effects_plot.png", conflict_effects_plot, width = 12, height = 9)

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
  SignalVars |>
  dplyr::select(-householdID, -childHHid) |>
  mutate(`Possession score` = rowSums(pick(LifestyleReality_1:LifestyleReality_8))) |>
  dplyr::select(-starts_with("Lifestyle"), -AlloparentingXsex) |>
  dplyr::select(where(\(x) sd(x, na.rm = T) > 0)) # remove vars set to mean for anonymity
sv2 <- set_names(sv2, shortform_dict[names(sv2)])

pca1 <- prcomp(sv2, scale. = T)
plot_loadings1 <- pca_loadings_plot(pca1, 1:2, reverse = 1) + theme_bw(16) + theme(legend.position = 'none')
plot_biplot1 <- pca_biplot(pca1, threshold = 0.19)

plot_pca <- plot_loadings1 + plot_biplot1 + plot_layout(widths = c(1,2))
plot_pca

ggsave("Figures/plot_pca.pdf", plot_pca, width = 20, height = 9)
ggsave("Figures/plot_pca.svg", plot_pca, width = 20, height = 10)
ggsave("Figures/plot_pca.png", plot_pca, width = 20, height = 12)

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

# Very long run time
# set.seed(456)
# m3 <- easybgm(SignalVars[-c(1,2,3,4,5)], type = 'mixed', not_cont = discrete_vars[-c(1,2,3,4,5)], iter = 1000000, package = 'BDgraph', g.prior = 0.05)
# save(m3, file = 'data/m3_1e6_cost_freq.rda')
load("data/m3_1e6_cost_freq.rda")

plot_network(m3, layout = "spring", vsize = 4, label.cex = 1, exc_prob = 0.9)
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

ggsave('Figures/plot_mst.png', plot_mst, width = 12, height = 9)
ggsave("Figures/plot_mst.pdf", plot_mst, width = 12, height = 9)

mst_unweighted <- igraph::mst(g3, algorithm = 'unweighted')
plot_mst_uw <- graph_plot(mst_unweighted, title = 'MST (unweighted)')
plot_mst_uw

ggsave('Figures/plot_mst_uw.png', plot_mst_uw, width = 12, height = 9)
ggsave("Figures/plot_mst_uw.pdf", plot_mst_uw, width = 12, height = 9)

mst_unit <- igraph::mst(g3, algorithm = 'prim', weights = rep(1, times = igraph::gsize(g3)))
plot_mst_unit <- graph_plot(mst_unit, title = 'MST (unit weights)')
plot_mst_unit

ggsave('Figures/plot_mst_unit.png', plot_mst_unit, width = 12, height = 9)
ggsave("Figures/plot_mst_unit.pdf", plot_mst_unit, width = 12, height = 9)


# annotate("rect", xmin = 12.4, xmax = 14.5, ymin = -0.75, ymax = 3.75, colour = 'red', fill = NA, linewidth = 2) +


# Ordinal regressions -----------------------------------------------------

# For binary predictors
ordinal_plot <- function(fit, predictor, data, title, ylabel = "Probability"){
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

  dsum <- d |>
    summarise(mean = mean(value), .by = c(name, {{predictor}}))
  print(dsum)

  ggplot(d, aes({{predictor}}, value, colour = name)) +
    geom_count(position = position_dodge(width = 0.1)) +
    geom_smooth(method='lm', se = F, position = position_dodge(width = 0.1)) +
    scale_x_continuous(breaks = c(0, 1)) +
    scale_color_viridis_d(option = 'B', end = 0.8) +
    guides(color = guide_legend("Response", reverse = T)) +
    ylim(0, NA) +
    labs(title = title, x = shortform_dict[deparse(substitute(predictor))], y = ylabel)+
    theme_minimal(15)
}

# For continuous predictors
ordinal_plot2 <- function(fit, predictor, data, title, ylabel="Probability"){
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
    labs(title = title, x = shortform_dict[deparse(substitute(predictor))], y = ylabel)+
    theme_minimal(15)
}

# Perceived need

SignalVars4 <-
  left_join(SignalVars, utila_df[c("householdID", "childHHid", "RelativeNeed3")]) |>
  dplyr::select(where(\(x) !all(duplicated(x)[-1L]))) |>
  relocate(RelativeNeed3, .after = 'childHHid') |>
  mutate(
    across(-c(1:3), as.numeric),
    across(SadFreqN:OtherKidsConflict, \(x) c(scale(x)))
  ) |>
  na.omit()

m_ordinalcv1 <- ordinalNetCV(
  as.matrix(SignalVars4[-c(1:3)]),
  SignalVars4[[3]],
  standardize = F,
  alpha = 1.0,
  family = "cumulative",
  link = "logit",
  lambdaMinRatio = 1e-04,
  printProgress = T
)
summary(m_ordinalcv1)
colMeans(summary(m_ordinalcv1))
# coef(m_ordinalcv1$fit, matrix = TRUE, whichLambda = 1)
m_ordinalcv1coefs <- coef(m_ordinalcv1$fit)[-c(1:2)]
names(m_ordinalcv1coefs) <- shortform_dict[names(m_ordinalcv1coefs)]
plot_need_coefs <- ggdotchart(m_ordinalcv1coefs)
plot_need_coefs
ggsave("Figures/plot_need_coefs.pdf", plot_need_coefs, width = 12, height = 12)
ggsave("Figures/plot_need_coefs.png", plot_need_coefs, width = 12, height = 12)

plot_need_age <- ordinal_plot2(m_ordinalcv1$fit, ChildAge, SignalVars4, 'Relative need', ylabel = "Probability of relative need")
plot_need_sad <- ordinal_plot2(m_ordinalcv1$fit, SadFreqN, SignalVars4, 'Relative need', ylabel = "Probability of relative need")

plot_need_combined <- plot_need_age / plot_need_sad + ggtitle("") + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a")
ggsave("Figures/plot_need_combined.pdf", plot_need_combined, width = 12, height = 12)
ggsave("Figures/plot_need_combined.png", plot_need_combined, width = 12, height = 12)

# Relative investment

SignalVars5 <-
  left_join(SignalVars, utila_df[c("householdID", "childHHid", "RelativeMaternalInvestment2", "RelativeNeed3")]) |>
    dplyr::select(where(\(x) !all(duplicated(x)[-1L]))) |>
  relocate(RelativeMaternalInvestment2, .after = 'childHHid') |>
  mutate(
    across(-c(1:3), as.numeric),
    across(SadFreqN:RelativeNeed3, \(x) c(scale(x)))
  ) |>
  na.omit()

m_ordinalcv2 <- ordinalNetCV(
  as.matrix(SignalVars5[-c(1:3)]),
  SignalVars5[[3]],
  standardize = F,
  alpha = 1.0,
  family = "cumulative",
  link = "logit",
  lambdaMinRatio = 1e-04,
  printProgress = T
)
summary(m_ordinalcv2)
colMeans(summary(m_ordinalcv2))
# coef(m_ordinalcv2$fit, matrix = TRUE, whichLambda = 1)
m_ordinalcv2coefs <- coef(m_ordinalcv2$fit)[-c(1:2)]
names(m_ordinalcv2coefs) <- shortform_dict[names(m_ordinalcv2coefs)]
plot_invest_coefs <- ggdotchart(m_ordinalcv2coefs)
plot_invest_coefs
ggsave("Figures/plot_invest_coefs.pdf", plot_invest_coefs, width = 12, height = 12)
ggsave("Figures/plot_invest_coefs.png", plot_invest_coefs, width = 12, height = 12)

plot_invest_need <- ordinal_plot2(m_ordinalcv2$fit, RelativeNeed3, SignalVars5, 'Relative investment', ylabel = "Probability of relative investment")
ggsave("Figures/plot_invest_need.pdf", plot_invest_need, width = 12, height = 12)
ggsave("Figures/plot_invest_need.png", plot_invest_need, width = 12, height = 12)

# Caregiver response

causematrix <-
  causes |>
  mutate(
    across(ConflictFamily:StatusConcerns, \(x) ifelse(x == "Maybe" | x == 'maybe', "0", x)),
    across(ConflictFamily:StatusConcerns, as.numeric)
  ) |>
  dplyr::select(householdID, childHHid, ConflictFamily:StatusConcerns) |>
  na.omit()
cause_count <- map_int(causematrix[-c(1:2)], \(x) sum(x, na.rm = T))
cause_vars <- names(cause_count[cause_count > 5])

SignalVars3 <-
  left_join(SignalVars, utila_df[c("householdID", "childHHid", "CaregiverResponse")]) |>
  left_join(causematrix[c("householdID", "childHHid", cause_vars)]) |>
  dplyr::select(where(\(x) !all(duplicated(x)[-1L]))) |>
  relocate(CaregiverResponse, .after = childHHid) |>
  mutate(
    across(-c(1:3), as.numeric),
    across(SadFreqN:StatusConcerns, \(x) c(scale(x))) # changed StatusConcerns to SeparationAttentionSeeking
    ) |>
  na.omit()

m_ordinalcv3 <- ordinalNetCV(
  as.matrix(SignalVars3[-c(1:3)]),
  SignalVars3[[3]],
  standardize = F,
  alpha = 1.0,
  family = "cumulative",
  link = "logit",
  lambdaMinRatio = 1e-04,
  printProgress = T
)
summary(m_ordinalcv3)
colMeans(summary(m_ordinalcv3))
# coef(m_ordinalcv3$fit, matrix = TRUE, whichLambda = 1)
m_ordinalcv3coefs <- coef(m_ordinalcv3$fit)[-c(1:2)]
names(m_ordinalcv3coefs) <- shortform_dict[names(m_ordinalcv3coefs)]
plot_caregiverresponse_coefs <- ggdotchart(m_ordinalcv3coefs)
plot_caregiverresponse_coefs
ggsave("Figures/plot_caregiverresponse_coefs.pdf", plot_caregiverresponse_coefs, width = 12, height = 12)
ggsave("Figures/plot_caregiverresponse_coefs.png", plot_caregiverresponse_coefs, width = 12, height = 12)

plot_caregiver_familyconflict <- ordinal_plot(m_ordinalcv3$fit, ConflictFamily, data = SignalVars3, title = 'Caregiver response', ylabel = "Probability of caregiver response")
plot_caregiver_trangsression <- ordinal_plot(m_ordinalcv3$fit, TransgressionMade, data = SignalVars3, title = 'Caregiver response', ylabel = "Probability of caregiver response")
plot_caregiver_loss <- ordinal_plot(m_ordinalcv3$fit, LossOfPrivlegesOrItem, data = SignalVars3, title = 'Caregiver response', ylabel = "Probability of caregiver response")
plot_caregiver_pain <- ordinal_plot(m_ordinalcv3$fit, DiscomfortPainInjuryIllness, data = SignalVars3, title = 'Caregiver response', ylabel = "Probability of caregiver response")

plot_caregiver_response_combined <- plot_caregiver_familyconflict + plot_caregiver_loss + ggtitle("") + plot_caregiver_trangsression + ggtitle("") + plot_caregiver_pain + ggtitle("") + plot_layout(guides = 'collect')
ggsave("Figures/plot_caregiver_response_combined.pdf", plot_caregiver_response_combined, width = 12, height = 12)
ggsave("Figures/plot_caregiver_response_combined.png", plot_caregiver_response_combined, width = 12, height = 12)

plot_caregiver_punish <- ordinal_plot(m_ordinalcv3$fit, Punishment, data = SignalVars3, title = 'Caregiver response')
ggsave("Figures/plot_caregiver_punish.pdf", plot_caregiver_punish, width = 12, height = 12)
ggsave("Figures/plot_caregiver_punish.png", plot_caregiver_punish, width = 12, height = 12)

# want to add relativeneed and relativeinvestment but code does not work

SignalVars6 <- SignalVars5 |>
  mutate(RelativeMaternalInvestment2 = as.numeric(RelativeMaternalInvestment2)) |>
  dplyr::select(-matches("Xsex|Old|Kids|Max|childHH|householdID"))

# code breaks when inserting 550:552 in here
plot_fullcorrmat <- ggcorrplot(
  cor(
    SignalVars6 |>
      set_names(nm = shortform_dict[names(SignalVars6)])
    ,
    use = 'pairwise.complete.obs'),
  hc.order = T,
  hc.method = 'ward.D2'
) +
  scico::scale_fill_scico(palette = 'vik', midpoint = 0)
plot_fullcorrmat

ggsave("Figures/plot_fullcorrmat.pdf", plot_fullcorrmat, width = 12, height = 12)
ggsave("Figures/plot_fullcorrmat.png", plot_fullcorrmat, width = 12, height = 12)
ggsave("Figures/plot_fullcorrmat.svg", plot_fullcorrmat, width = 12, height = 12)
