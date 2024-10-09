
library(glmnet)
library(glmnetUtils)
library(hagenutils)
library(patchwork)
library(skimr)
library(pvclust)
library(ggcorrplot)

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
nvar <- ncol(mdf2)
m <- cva.glmnet(
  as.matrix(mdf2[7:nvar]),
  mdf2[["TantrumFreqN"]],
  standardize = T,
  family =  quasipoisson() # MASS::negative.binomial(3) #
)
minlossplot(m, cv.type = 'min')
coefs <- coef(m$modlist[[1]], s = 'lambda.min')[-1,]
ggdotchart(coefs) + geom_vline(xintercept = 0, linetype = 'dotted')


# Sad: alpha = 0
# Cry: alpha = 1
# Tantrum: alpha = 0

# Now fit with those alpha values

glmnetSignal <- function(signal, alpha = 0, relax = T, s = "lambda.min"){
  m <- cv.glmnet(
    as.matrix(mdf2[8:nvar]),
    mdf2[[signal]],
    alpha = alpha, # alpha = 0: ridge; alpha = 1: lasso
    relax = relax,
    standardize = T,
    family =  quasipoisson() # MASS::negative.binomial(3) #
  )

  coefs <- coef(m, s=s)[-1,]
  ggdotchart(coefs) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    labs(title = signal, subtitle = str_glue("alpha = {alpha}; relax = {relax}; s = {s}"))
}

p1 <- glmnetSignal("SadFreqN", alpha = 0.0, relax = T, s = "lambda.min")
p2 <- glmnetSignal("CryFreqN", alpha = 1.0, relax = F, s = "lambda.min")
p2p3 <- glmnetSignal("TantrumFreqN", alpha = 0.0, relax = T, s = "lambda.min")

p1/p2/p3 + plot_layout(axes = 'collect')


# Plot predictions

m <- cv.glmnet(
  as.matrix(mdf2[8:nvar]),
  mdf2[["CryFreqN"]],
  alpha = 1, # alpha = 0: ridge; alpha = 1: lasso
  relax = F,
  standardize = T,
  family =  quasipoisson() # MASS::negative.binomial(3) #
)


m2 <- glmnet(
  as.matrix(mdf2[8:nvar]),
  mdf2[["CryFreqN"]],
  alpha = 1, # alpha = 0: ridge; alpha = 1: lasso
  lambda = m$lambda.min,
  relax = F,
  standardize = T,
  family =  quasipoisson() # MASS::negative.binomial(3) #
)

# library(tidymodels)

fit <-
  poisson_reg(mixture = 1, penalty = m$lambda.min) %>%
  set_engine("glmnet", family = quasipoisson) %>%
  fit(CryFreqN ~ ., data = mdf2[c(2, 8:22)])

tdy <- tidy(fit) |> dplyr::filter(lambda == m$lambda.min)
est <- tdy$estimate[-1]
names(est) <- tdy$term[-1]
ggdotchart(est)

nd <- datagrid(ChildAge = seq(-1.6, 2.5, 0.01), Sex = rep(c(-1, 1), nrow(mdf2)), newdata=mdf2[8:22])
out <- predictions(fit, type = "numeric", newdata = nd)
ggplot(out, aes(ChildAge, estimate, color = factor(Sex))) + geom_line()


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
