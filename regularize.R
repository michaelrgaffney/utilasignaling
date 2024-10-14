
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


# Regularized regression --------------------------------------------------

glmnet2 <- function(d, outcome, indices, alpha = 1){
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
      family = quasipoisson() # MASS::negative.binomial(0.4) #quasipoisson()
    )
    cvm <- c(cvm, list(m$cvm))
  }
  
  cvm.mean <- rowMeans(matrix(list_c(cvm), nrow = length(cvm[[1]])))
  # return(tibble(lambda = m$lambda, cvm = rowMeans(cvm.mat)))
  
  lambda.min <- m$lambda[which.min(cvm.mean)]

  m2 <-
    poisson_reg(mixture = alpha, penalty = lambda.min) |> 
    set_engine("glmnet", family = quasipoisson) |> # MASS::negative.binomial(0.4)
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


space2newline <- function(v) str_replace_all(v, " ", "\n")

plotpredictions <- function(params, outcome, condition, title = NULL){
  model <- params$out[[outcome]]$model
  data <- params$out[[outcome]]$data
  if(is.null(title)) title <- outcome
  if(length(condition) == 2){
    if(condition[2] == 'Sex') lbls = c("Female", "Male") else lbls = c("No", "Yes")
    scl <- scale_color_binary(labels = lbls, name = shortform_dict[condition[2]] |> space2newline())
  } else {
    scl <- scale_color_manual(values = 'black')
  }
  p <- 
    plot_predictions(model, condition = condition, newdata = data) +
    scl +
    ylim(0, 30) +
    labs(title = title, x = shortform_dict[condition[1]]) +
    theme_minimal(15)
  p$layers[[1]]$aes_params$linewidth <- 2
  return(p)
}

plotpredictions(signalparams, "CryFreqN1", c("ChildAge", "LifestyleReality_5"), "Crying")
plotpredictions(conflictparams, "ConflictFreqN1", c("ChildAge", "Sex"))

# plot_predictions(signalparams$out$CryFreqN1$model, condition = c('ChildAge', 'Sex'), newdata = signalparams$out$CryFreqN1$data)
# plot_predictions(conflictparams$out$ConflictFreqN1$model, condition = c('ChildAge', 'Sex'), newdata = conflictparams$out$ConflictFreqN1$data)

# SignalVars2 <-
#   SignalVars |>
#   mutate(
#     AlloparentingXage = c(scale(AlloparentingFreqN) * scale(ChildAge))
#   )

# Other stuff

out <- skim(modeldf)
nms <- out$skim_variable[out$skim_type == 'numeric' & out$complete_rate > 0.9]
nms <- c(nms, "Sex", "UserLanguage", "ImmigrateUtila", "PartnerStatus", "OldestChild")
mdf2 <-
  modeldf |>
  dplyr::select(
    all_of(nms),
    # -householdID,
    -childHHid,
    -NegativeResponse,
    -PositiveResponse,
    -IncomeCategoryN,
    -ConflictFreqN,
    -OtherChildrenHH,
    -contains('IllnessSusceptibility'),

    # New omits
    -OlderGirls,
    -YoungerKids,
    -FoodSecurity,
    -UserLanguage

    ) |>
  na.omit() |>
  mutate(
    Sex = ifelse(Sex == "Female", 0, 1),
    # UserLanguage = ifelse(UserLanguage == "EN", 0, 1),
    ImmigrateUtila = ifelse(ImmigrateUtila == "No", 0, 1),
    PartnerStatus = ifelse(PartnerStatus == "Unpartnered", 0, 1),
    across(-c(1:6), \(x) c(scale(x))),
    AlloparentingXsex = AlloparentingFreqN * Sex,
    # ChildAgeXsex = ChildAge * Sex,
    # OldestXsex = OldestChild * Sex
  )

# PCA ---------------------------------------------------------------------

m <- prcomp(mdf2[-c(4:6)], scale. = T) # Remove composite signaling vars
plot_loadings <- pca_loadings_plot(m, 1:2) # + theme(legend.position = 'top', legend.title = element_blank())
plot_biplot <- pca_biplot(m) + theme_minimal(15)

plot_pca <- plot_loadings + plot_biplot + plot_layout(widths = c(1,2))
plot_pca

# Regularized regression --------------------------------------------------


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


# Clustering --------------------------------------------------------------

mdf3 <-
  mdf2 |>
  dplyr::select(-matches("Older|Oldest|Xsex|Max|Kids"))

out <- pvclust(mdf3, method.hclust = 'ward.D2', method.dist = 'abscor')
plot(out)

var_clusters <- cutree(out$hclust, k = 3)
signal_vars <- names(var_clusters[var_clusters == 1])
household_vars <- names(var_clusters[var_clusters == 2])
neighborhood_vars <- names(var_clusters[var_clusters == 3])

# Theory based var groups

signal_vars <- c("CryFreqN", "SadFreqN", "TantrumFreqN", "SignalFreq", "SignalCost")
household_vars <- c(
  "ChildAge",
  "number_adults",
  "ConflictFreqN",
  "AlloparentingFreqN",
  "NumberOfChildren",
  "CaregiverAge",
  "Sex",
  "PartnerStatus",
  "LogIncome",
  "EducationLevelYears",
  "FoodSecurity",
  "HouseQuality",
  "UserLanguage",
  "ImmigrateUtila"
  )

neighborhood_vars <- c("NeighborhoodQuality", "Neighborhood2")

# Correlation matrices ----------------------------------------------------

ggcorrplot(
  cor(
    mdf2 |> dplyr::select(-matches("Xsex|Old|Kids|Max")
    ),
    use = 'pairwise.complete.obs'),
  hc.order = T,
  hc.method = 'centroid'
) +
  scico::scale_fill_scico(palette = 'vik', midpoint = 0)


# Final models ------------------------------------------------------------

frmla <- paste("CryFreqN ~", paste(names(mdf2)[c(8:10, 12:21)], collapse = ' + '), " + AlloparentingFreqN*Sex", " + (1|householdID)")
m <- glmmTMB(as.formula(frmla), data = mdf2, family = nbinom2)
summary(m)

glmmTMBformulas <- c(std_formulas, bodyfat_formulas, conflict_formula)
names(glmmTMBformulas) <- c(paste0(signals, "illness"), paste0(signals, "bodyfat"), "ConflictFreqN")

glmmTMBmodels <- tibble(
  Outcome = c(signals, signals, "ConflictFreqN"),
  Model = map(glmmTMBformulas, \(f) glmmTMB(as.formula(f), family = nbinom2, data = d2)),
  AIC = map_dbl(Model, AIC)
)
