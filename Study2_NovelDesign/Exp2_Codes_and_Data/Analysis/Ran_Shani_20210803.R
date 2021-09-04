

# Model predicting outcome with "Theta", only for 100% contingency (w. subset)
DM_dat_long = DM_dat_long %>% 
  mutate(theta_10 = theta + rnorm(n = nrow(DM_dat_long), sd = .10),
         theta_11 = theta + rnorm(n = nrow(DM_dat_long), sd = .09),
         theta_12 = theta + rnorm(n = nrow(DM_dat_long), sd = .08),
         theta_13 = theta + rnorm(n = nrow(DM_dat_long), sd = .07),
         theta_14 = theta + rnorm(n = nrow(DM_dat_long), sd = .06),
         theta_14 = theta + rnorm(n = nrow(DM_dat_long), sd = .05),
         )

# Our plan - each time predict choices with a slightly different Theta, and plot how much better the AIC gets
AIC(glmer(Outcome ~ 1 + theta_10 + (1 + theta_10|subjectID), data = subset(DM_dat_long, Contingency==1) , family = "binomial"))
AIC(glmer(Outcome ~ 1 + theta_11 + (1 + theta_11|subjectID), data = subset(DM_dat_long, Contingency==1) , family = "binomial"))
AIC(glmer(Outcome ~ 1 + theta_12 + (1 + theta_12|subjectID), data = subset(DM_dat_long, Contingency==1) , family = "binomial"))
AIC(glmer(Outcome ~ 1 + theta_13 + (1 + theta_13|subjectID), data = subset(DM_dat_long, Contingency==1) , family = "binomial"))
AIC(glmer(Outcome ~ 1 + theta + (1 + theta|subjectID), data = subset(DM_dat_long, Contingency==1) , family = "binomial"))

m2 = glmer(Outcome ~ 1 + theta + (1 + theta|subjectID), data = subset(DM_dat_long, Contingency==1) , family = "binomial")

# Incorrect ()
m3 = glmer(Outcome ~ 1 + theta + (1 + theta|subjectID), data = subset(DM_dat_long, Contingency==1))
m3 = glmer(Outcome ~ 1 + rand + (1 + rand|subjectID), data = subset(DM_dat_long, Contingency==1) , family = "binomial")
summary(m2)
summary(m3)
AIC(m2)
AIC(m3)

# Expit function (log-odds to p) , Odds = p/(1-p)
expit = function(x){exp(x)/(1+exp(x))}

