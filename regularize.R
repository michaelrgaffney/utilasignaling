
library(glmnet)
library(hagenutils)
library(patchwork)
library(skimr)

out <- skim(modeldf)
nms <- out$skim_variable[out$skim_type == 'numeric' & out$complete_rate > 0.9]
nms <- c(nms, "Sex", "UserLanguage", "ImmigrateUtila", "PartnerStatus")
mdf2 <-
  modeldf |>
  dplyr::select(
    all_of(nms),
    -householdID,
    -childHHid,
    -NegativeResponse,
    -PositiveResponse,
    -IncomeCategoryN,
    -ConflictFreqN,
    -OtherChildrenHH, # omit or keep?
    -contains('IllnessSusceptibility')
    ) |>
  na.omit() |>
  mutate(
    Sex = ifelse(Sex == "Female", 0, 1),
    UserLanguage = ifelse(UserLanguage == "EN", 0, 1),
    ImmigrateUtila = ifelse(ImmigrateUtila == "No", 0, 1),
    PartnerStatus = ifelse(PartnerStatus == "Unpartnered", 0, 1),
    across(-c(1:6), \(x) c(scale(x))),
    AlloparentingXsex = AlloparentingFreqN * Sex,
    ChildAgeXsex = ChildAge * Sex
  )

# PCA ---------------------------------------------------------------------

m <- prcomp(mdf2[-c(4:6)], scale. = T) # Remove composite signaling vars
plot_loadings <- pca_loadings_plot(m, 1:2) # + theme(legend.position = 'top', legend.title = element_blank())
plot_biplot <- pca_biplot(m) + theme_minimal(15)

plot_pca <- plot_loadings + plot_biplot + plot_layout(widths = c(1,2))
plot_pca

# Regularized regression --------------------------------------------------

nvar <- ncol(mdf2)

glmnetSignal <- function(signal, alpha = 0, relax = T, s = "lambda.min"){
  m <- cv.glmnet(
    as.matrix(mdf2[7:nvar]),
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

p1 <- glmnetSignal("SadFreqN", alpha = 0.0, relax = F, s = "lambda.min")
p2 <- glmnetSignal("CryFreqN", alpha = 0.0, relax = F, s = "lambda.min")
p3 <- glmnetSignal("TantrumFreqN", alpha = 0.0, relax = F, s = "lambda.min")

p1/p2/p3 + plot_layout(axes = 'collect')
