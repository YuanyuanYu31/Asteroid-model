library(dplyr)
library(tidyverse)
library(highcharter)
library(knitr)
library(kableExtra)
library(ggcorrplot)
library(brms)
library(car)
library(bayesplot)
library(tidybayes)
library(ggplot2)
library(rjags)
library(runjags)
library(MCMCpack)
library(BANOVA)

impact <- read.csv("./impacts.csv", header = TRUE) # impact data
dat <- read.csv("./merged.csv", header = TRUE) # merged impact and orbit data
head(dat)
summary(dat)
nrow(dat)
nrow(impact)
unique(dat$Epoch..TDB.)
print(unique(dat$Object.Classification))
# exploratory data analysis
impacts2 <- read.csv("./impacts.csv", header = TRUE) %>%
  mutate(Period.Length = Period.End - Period.Start, .after = Period.End)
%>%
  select(Period.Length, Possible.Impacts, Asteroid.Velocity,
         Asteroid.Magnitude, Asteroid.Diameter..km., Cumulative.Palermo.Scale)
impacts2 %>%
  keep(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
impacts3 <- read.csv("./impacts.csv", header = TRUE) %>%
  mutate(Log.Cumulative.Impact.Probability =
           log(Cumulative.Impact.Probability), .after =
           Cumulative.Impact.Probability) %>%
  select(Cumulative.Impact.Probability, Log.Cumulative.Impact.Probability)
impacts3 %>%
  keep(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
dat2 <- read.csv("./merged.csv", header = TRUE) %>%
  select(Perihelion.Distance..AU., Orbit.Eccentricity)
dat2 %>%
  keep(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
# question 1 analysis
impact %>%
  mutate(Period_Length = Period_End - Period_Start, .after = Period_End)
#check the multicollinearity
fit_lp <- lm(log(Cumulative_Impact_Probability) ~ Period_Length +
               Asteroid_Velocity + Asteroid_Magnitude + Asteroid_Diameter +
               Cumulative_Palermo_Scale, data = impact)
vif(fit_lp)
#fit the Bayesian regression model
fit_bp <- brm(log(Cumulative_Impact_Probability) ~ Period_Length +
                Asteroid_Velocity + Asteroid_Magnitude + Asteroid_Diameter +
                Cumulative_Palermo_Scale, data = impact)
get_variables(fit_bp)
plot(fit_bp)
pp_check(fit_bp)
#drop Asteroid_Diameter
fit_bp1 <- brm(log(Cumulative_Impact_Probability) ~ Period_Length +
                 Asteroid_Velocity + Asteroid_Magnitude + Cumulative_Palermo_Scale, 
               data = impact)
plot(fit_bp1)
pp_check(fit_bp1)
#prediction
y1 <- predict(fit_bp1)
prob1 <- data.frame(exp(y1))
prob1[which.max(prob1$Estimate), ]
impacts[569,]
prior_summary(fit_bp1)
# question 2 analysis
data = impact
# create factor variable
data$Impact.Amount = as.factor(ifelse(data$Possible.Impacts <
                                        7,"Few",ifelse(data$Possible.Impacts > 44,"Many","Moderate")))
data$Impact.Amount = factor(data$Impact.Amount,levels = c("Few",
                                                          "Moderate", "Many"))
data$Period.Length = data$Period.End-data$Period.Start
data$id = 1:nrow(data) # id variable for BANOVA

data$logAsteroid.Velocity = log(data$Asteroid.Velocity+1) # transform
velocity and diameter
data$logAsteroid.Diameter = log(data$Asteroid.Diameter..km.)
data$Class <- as.factor(gsub("\\s*\\([^\\)]+\
\)","",as.character(data$Object.Classification)))
# Plots of relevant variables
ggplot(data = data,aes(x = Asteroid.Magnitude,fill = Impact.Amount)) +
  geom_density(alpha=.7) +
  labs(x = "Magnitude")
ggplot(data = data,aes(x = log(Asteroid.Diameter..km.),fill =
                         Impact.Amount)) +
  geom_density(alpha=.7) +
  labs(x = "log(Diameter)")
ggplot(data = data,aes(x = log(Asteroid.Velocity+1),fill = Impact.Amount))
+
  geom_density(alpha=0.7) +
  labs(x = "log(Velocity + 1)")
ggplot(data = data,aes(x = Cumulative.Palermo.Scale,fill = Impact.Amount))
+
  geom_density(alpha=0.7) +
  labs(x = "Cumulative Palermo Scale")

# Velocity
fit_vel = BANOVA.Normal(logAsteroid.Velocity ~ Impact.Amount, l1_hyper =
                          c(1,1,0.0001), data = data, id = data$id, burnin
                        = 5000,
                        sample = 1000, thin = 10)
## Convergence Tests
conv.diag(fit_vel)
## Posterior means
table.predictions(fit_vel)
## Coefficients and Bayesian p-values
fit_vel$coef.tables
## Effect Sizes and credible intervals
BAnova(fit_vel)
## Trace plots
trace.plot(fit_vel)
# Magnitude
fit_mag = BANOVA.Normal(Asteroid.Magnitude ~ Impact.Amount, l1_hyper =
                          c(1,1,0.0001), data = data, id = impact$id,
                        burnin = 5000,
                        sample = 1000, thin = 10)
## Convergence Tests
conv.diag(fit_mag)
## Posterior means
table.predictions(fit_mag)
## Coefficients and Bayesian p-values
fit_mag$coef.tables
## Effect Sizes and credible intervals
BAnova(fit_mag)
## Trace plots
trace.plot(fit_mag)

# Diameter
fit_diam = BANOVA.Normal(logAsteroid.Diameter ~ Impact.Amount, l1_hyper =
                           c(1,1,0.0001), data = data, id = data$id,
                         burnin = 5000,
                         sample = 1000, thin = 10)
## Convergence Tests
conv.diag(fit_diam)
## Posterior means
table.predictions(fit_diam)
## Coefficients and Bayesian p-values
fit_diam$coef.tables
## Effect Sizes and credible intervals
BAnova(fit_diam)
## Trace plots
trace.plot(fit_diam)

# Cumulative Palermo Scale
fit_cps = BANOVA.Normal(Cumulative.Palermo.Scale ~ Impact.Amount, l1_hyper
                        =
                          c(1,1,0.0001), data = data, id = impact$id,
                        burnin = 5000,
                        sample = 1000, thin = 10)
## Convergence Tests
conv.diag(fit_cps)
## Posterior means
table.predictions(fit_cps)
## Coefficients and Bayesian p-values
fit_cps$coef.tables
## Effect Sizes and credible intervals
BAnova(fit_cps)
## Trace plots
trace.plot(fit_cps)

# question 3 analysis
asteroid_merge <- dat
head(asteroid_merge)
table(asteroid_merge$Object.Classification)
# The palermo scale for different object classification
ggplot(data = asteroid_merge, aes(x =Cumulative.Palermo.Scale, fill=
                                    Object.Classification))
+ geom_density(alpha=0.7)
+ labs(x = "Cumulative Palermo Scale")
#fit the model
fit <- brm(log(Cumulative.Impact.Probability) ~ 1 + (1
                                                     +Cumulative.Palermo.Scale|Object.Classification)
           + (1 + Perihelion.Distance..AU. +Asteroid.Magnitude+
                Orbit.Eccentricity|Epoch..TDB.)
           + Orbit.Eccentricity + Asteroid.Magnitude +
             Cumulative.Palermo.Scale
           +Perihelion.Distance..AU., data = asteroid_merge, control =
             list(adapt_delta = 0.98),
           chains = 2, iter = 2000, warmup = 1000)
summary(fit)

plot(fit)
pp_check(fit)
#prediction
y_fit <- predict(fit)
y_fit
prob_fit <- data.frame(exp(y_fit))
prob_fit[which.max(prob_fit$Estimate),]
asteroid_merge[567,]

