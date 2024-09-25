
library(mgcv)
library(lmerTest)
library(marginaleffects)
library(gratia)

signals <-
  d2 |>
  dplyr::select(uniqueID, ChildAge, SadFreqN, CryFreqN, TantrumFreqN) |>
  pivot_longer(SadFreqN:TantrumFreqN, names_to = 'Type', values_to = 'Freq') |>
  mutate(uniqueID = factor(uniqueID), Type = factor(Type)) |>
  na.omit()

m <- gam(Freq ~ s(ChildAge, by = Type) + s(uniqueID, bs = 're'), family = nb, data = signals)
summary(m)
draw(m, select = 1:3, fun = exp) # Correct?

plot_predictions(m, condition = c('ChildAge', 'Type'), type = 'link', transform = exp)
plot_predictions(m, condition = c('ChildAge', 'Type'), type = 'response')

m2 <- lmer()
