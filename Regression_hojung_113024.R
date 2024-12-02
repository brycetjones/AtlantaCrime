library(tidyverse)
library(ggplot2)
library(jtools)
library(GGally)
library(stargazer)
library(car)
library(lmtest)
library(sandwich)
library(fastDummies)
library(broom.mixed)
library(sjPlot)
library(leaps)
library(lmtest)
library(MASS)
library(psych)
library(pscl)

getwd()
setwd("C:/Users/Hojung Yu/Documents/GitHub/AtlantaCrime/")

df <- read.csv("data_final.csv")

head(df)
nrow(df)
ncol(df)
summary(df)


colnames(df)
df["nonwhite_ratio"] = 1 - df["white_ratio"]
df <- df %>% relocate("nonwhite_ratio", .before = "median_incomeE")
reg_variables = c("pop_den", "white_ratio", "nonwhite_ratio", "median_incomeE",
                  "less_than_hs_ratio", "Commercial", "HighdensityResidential","Industrial",
                  "Institutional", "LowdensityResidential","ResidentialCommercial",
                  "min_station_dist", "vio_crimerate","nonvio_crimerate")


ggpairs(df[,reg_variables])

#####################################################Violent crime model


cor_matrix <- cor(df[reg_variables])
round(cor_matrix,2)

#1. Normal OLS
ols_1 <- lm(I(vio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100) 
            + median_incomeE + I(less_than_hs_ratio*100) + Commercial + HighdensityResidential + Industrial + Institutional +
              + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df)
summary(ols_1)
plot(ols_1)

cd_ols1 <- cooks.distance(ols_1)
df$cd_ols1 <- cd_ols1

plot(cd_ols1, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean (cd_ols1, na.rm = T), col = "red")
text(x = 1: length(cd_ols1) + 5,
     y = cd_ols1,
     col = "red",
     labels = ifelse(cd_ols1 > 4 * mean(cd_ols1, na.rm = T),
                     names(cd_ols1),
                     ""))

#Remove influential observations
df_noout_vio <- df[df$cd_ols1 < 4 * mean(cd_ols1, na.rm = T), ] # exclude hih-influence points not using 4 * mean
ols_1_noout <- lm(I(vio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100)
                  + median_incomeE + I(less_than_hs_ratio*100) + Commercial + HighdensityResidential + Industrial + Institutional +
                    + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout_vio)

summary(ols_1_noout)
# Compare no outliers to outliers using scaled coefficient plots
plot_summs(ols_1, ols_1_noout, scale = TRUE)
AIC(ols_1)
AIC(ols_1_noout)

# Not conducting log-regression since it is not significant with these datasets.
# Using df_noout_vio without outliers since it has smaller AIC score.
# Variables Selection

##Stepwise
null.model <- lm(I(vio_crimerate*100) ~ 1, data = df_noout_vio)
full.model <- lm(I(vio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100)
                 + median_incomeE + I(less_than_hs_ratio*100) + Commercial + HighdensityResidential + Industrial + Institutional +
                   + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout_vio)

step.model.for <- step(null.model,
                       scope = formula(full.model),
                       direction = "forward",
                       trace = 0)

step.model.back <- step(full.model,
                        direction = "backward",
                        trace = 0)

step.model.both <- step(null.model,
                        scope = formula(full.model),
                        direction = "both",
                        trace = 0)

stargazer(step.model.for, step.model.back, step.model.both,
          type = "text",
          add.lines = list(c("AIC", round(AIC(step.model.for),1), round(AIC(step.model.back),1), round(AIC(step.model.both),1))),
          column.labels = c("Forward", "Backward", "Both"))

subset.model <- regsubsets(formula(full.model),
                           data = df,
                           nvmax = 11,
                           method = "exhaustive")

##Subset models
reg.summary <- summary(subset.model)
reg.summary

par(mfrow = c(2,2)) # This code specifies that, instead of showing only one graph in the plots window, you want 4 (2 x 2) graphs.

plot(reg.summary$rsq, # This code plots the number of variables on x-axis and r-squared values on y-axis
     xlab = "Number of variables", 
     ylab = "R-Squared", 
     type = "l") # type = "l" means that you want to display the graph with a line, not points.

plot2 <- plot(reg.summary$adjr2, # Adjusted r-squared values on y-axis
              xlab = "Number of variables", 
              ylab = "Adjusted R-Sqaured", 
              type = "l")

points(which.max(reg.summary$adjr2), # This code inserts a red dot on the line graph. The red dot denotes the point where the adjusted R-squared is at its highest.
       reg.summary$adjr2[which.max(reg.summary$adjr2)], 
       col="red",
       cex=2,
       pch=20)

plot(reg.summary$rss, # sum of squares residual values on y-axis
     xlab = "Number of variables", 
     ylab = "RSS", 
     type = "l")

plot(reg.summary$bic, # BIC values on y-axis. BIC is very similar to AIC except that it puts more penalty to the number of variables included.
     xlab = "Number of variables", 
     ylab = "BIC", 
     type = "l")

points(which.min(reg.summary$bic), reg.summary$bic[which.min(reg.summary$bic)], col="red",cex=2,pch=20)

######################optimal regression in violent crime
#BIC lowest model
vio_ols_v1 <- lm(I(vio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100) +  I(less_than_hs_ratio * 100) 
                 + Industrial + Institutional + LowdensityResidential + min_station_dist, data = df_noout_vio)
#Adjusted R-Squared highest model
vio_ols_v2 <- lm(I(vio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100) +  I(less_than_hs_ratio * 100) 
                 + Industrial + Institutional + min_station_dist
                 + median_incomeE + Commercial + HighdensityResidential + ResidentialCommercial, data = df_noout_vio)

summary(vio_ols_v1)
summary(vio_ols_v2)

AIC(vio_ols_v1)
AIC(vio_ols_v2)

##Check Multicollinearity
dev.off()
vif_v1 <- vif(vio_ols_v1)
barplot(vif_v1, main = "VIF Values", horiz = TRUE, col = "steelblue") # seems pretty good

##Check Heteroskedasticity
plot(predict(vio_ols_v1), 
     vio_ols_v1$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(vio_ols_v1) 
#Chisquare = 108.5503, Shows strong heteroskedasticity!!! So, I generalized vio_crimerate

vio_ols_v1_log <- lm(I(log(vio_crimerate*100 + 0.01)) ~ pop_den + I(nonwhite_ratio * 100) +  I(less_than_hs_ratio * 100) 
                     + Industrial + Institutional + LowdensityResidential + min_station_dist, data = df_noout_vio)
summary(vio_ols_v1_log)

plot(predict(vio_ols_v1_log), 
     vio_ols_v1_log$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(vio_ols_v1_log) 
#Chisquare 89.95286 even bigger

#Still not fixing!!

#Considering interactive term(Yanfu)
#for non_white ratio and population density, there might be an interaction in how they 
#influence the crime rate.

#check pop_den and non-white for interaction
range(df_noout_vio$nonwhite_ratio)
range(df_noout_vio$pop_den)

#low and high value subsets for non_white
nonwhite_low <- df_noout_vio[df_noout_vio$nonwhite_ratio <= 0.2, ]
nonwhite_high <- df_noout_vio[df_noout_vio$nonwhite_ratio >= 0.8, ]
nw_low_mod <- lm(I(log(vio_crimerate * 100 + 0.001)) ~ pop_den, data = nonwhite_low)
summary(nw_low_mod)
nw_high_mod <- lm(I(log(vio_crimerate * 100 + 0.001)) ~ pop_den, data = nonwhite_high)
summary(nw_high_mod)

par(mfrow = c(1,1))
plot(x = nonwhite_low$pop_den, y = log(nonwhite_low$vio_crimerate) * 100 + 0.001, 
     pch = 19, xlab = "population density", ylab = "log_crimerate", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25)) +
  abline(nw_low_mod, col = "blue", lwd = 2) +
  points(x = nonwhite_high$pop_den, log(nonwhite_high$vio_crimerate) * 100 + 0.001, 
         col = rgb(red = 1, green = 0, blue = 0, alpha = 0.25), pch = 19) +
  abline(nw_high_mod, col = "red", lwd = 2)
#from this plot we can see the slopes of low non white and high non white are different. There is
#an interaction between non white and pop den.

#We can now add the interaction term to model
vio_ols_v1_log_interaction <- lm(I(log(vio_crimerate*100 + 0.001)) ~ pop_den + I(nonwhite_ratio * 100) + pop_den:I(nonwhite_ratio * 100)
                                 + I(less_than_hs_ratio * 100) + Industrial + Institutional
                                 + LowdensityResidential + min_station_dist, data = df_noout_vio)

summary(vio_ols_v1_log_interaction)

plot(predict(vio_ols_v1_log_interaction), 
     vio_ols_v1_log_interaction$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)

#Check heteroskedasticity
ncvTest(vio_ols_v1_log_interaction)
#Chisquare 170 even bigger

#The interaction term is statistically significant.

# Considering glm function (WIP_HOJUNG)
# create binary column for violent crime rate
median_crimerate <- median(df_noout_vio$vio_crimerate, na.rm = TRUE)
df_noout_vio$binary_vio <- ifelse(df_noout_vio$vio_crimerate > median_crimerate, 1, 0)

median_crimerate <- median(df$vio_crimerate, na.rm = TRUE)
df$binary_vio <- ifelse(df$vio_crimerate > median_crimerate, 1, 0)

# binary logical regression
vio_glm_v1 <- glm(binary_vio ~ I(pop_den * 1000) + I(nonwhite_ratio * 100) + 
                    I(less_than_hs_ratio * 100) + Industrial + Institutional + 
                    LowdensityResidential + min_station_dist, data = df_noout_vio, family = "binomial")

# binary logical regression
vio_glm_v2 <- glm(binary_vio ~ I(pop_den * 1000) + I(nonwhite_ratio * 100) + 
                    I(less_than_hs_ratio * 100) + + Institutional +
                    LowdensityResidential, data = df_noout_vio, family = "binomial")
summary(vio_glm_v1)
summary(vio_glm_v2)
vif(vio_glm_v1)
vif(vio_glm_v2)
# creating odds ratio
# Regression result (with odds ratio conversion)
round( # Rounds the numbers up to 3 digits
  cbind( # Column-bind Odds Ratio to the regerssion output
    "Odds Ratio" = exp(vio_glm$coefficients),
    summary(vio_glm)$coefficients
  ),3)

# pseudo R-squared 0.4361733
pR2(vio_glm_v1)[4]
pR2(vio_glm_v2)[4]

# plotting the logistic curve

# population density with violent crime
ggplot(df_noout_vio, aes(I(pop_den * 1000), binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Population density") +
  ylab("Probability of violent crime")

# nonwhite_ratio with violent crime
ggplot(df_noout_vio, aes(I(nonwhite_ratio * 100), binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Non white ratio") +
  ylab("Probability of violent crime")


summary(df$nonvio_crimerate)

##Below 5 or 10, Since all VIF values are well below 5 or 10, 
#there is no evidence of significant multicollinearity in your model based on these results.



###############################################################
#############################################################NONViolent crime model
# 1. Normal OLS
outliers_regression_non <- lm(I(nonvio_crimerate*100) ~ pop_den + I(nonwhite_ratio*100)
                          + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                            + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df)
summary(outliers_regression_non)


cd_ols2 <- cooks.distance(outliers_regression_non)
df$cd_ols2 <- cd_ols2

plot(cd_ols2, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean (cd_ols2, na.rm = T), col = "red")
text(x = 1: length(cd_ols2) + 5,
     y = cd_ols2,
     col = "red",
     labels = ifelse(cd_ols2 > 4 * mean(cd_ols2, na.rm = T),
                     names(cd_ols2),
                     ""))

#Remove influential observations
df_noout_nonvio <- subset(df[df$cd_ols2 < 4 * mean(cd_ols2, na.rm = T), ])
outliers_regression_noout_non <- lm(I(nonvio_crimerate*100) ~ pop_den + I(nonwhite_ratio*100)
                                + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                                  + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)
summary(outliers_regression_noout_non)

#Plot the difference

plot_summs(outliers_regression_non, outliers_regression_noout_non, scale = TRUE)
AIC(outliers_regression_non)
AIC(outliers_regression_noout_non)

# Not conducting log-regression since it is not significant with these datasets.
# Using df_noout_NONvio without outliers since it has smaller AIC score.
# Variables Selection


##Stepwise
null.model <- lm(I(nonvio_crimerate*100) ~ 1, data = df_noout_nonvio)
full.model <- lm(I(nonvio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100)
                 + median_incomeE + I(less_than_hs_ratio*100) + Commercial + HighdensityResidential + Industrial + Institutional +
                   + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)

step.model.for <- step(null.model,
                       scope = formula(full.model),
                       direction = "forward",
                       trace = 0)

step.model.back <- step(full.model,
                        direction = "backward",
                        trace = 0)

step.model.both <- step(null.model,
                        scope = formula(full.model),
                        direction = "both",
                        trace = 0)

stargazer(step.model.for, step.model.back, step.model.both,
          type = "text",
          add.lines = list(c("AIC", round(AIC(step.model.for),1), round(AIC(step.model.back),1), round(AIC(step.model.both),1))),
          column.labels = c("Forward", "Backward", "Both"))

subset.model <- regsubsets(formula(full.model),
                           data = df,
                           nvmax = 11,
                           method = "exhaustive")


reg.summary <- summary(subset.model)
reg.summary

par(mfrow = c(2,2)) # This code specifies that, instead of showing only one graph in the plots window, you want 4 (2 x 2) graphs.

plot(reg.summary$rsq, # This code plots the number of variables on x-axis and r-squared values on y-axis
     xlab = "Number of variables", 
     ylab = "R-Squared", 
     type = "l") # type = "l" means that you want to display the graph with a line, not points.

plot2 <- plot(reg.summary$adjr2, # Adjusted r-squared values on y-axis
              xlab = "Number of variables", 
              ylab = "Adjusted R-Sqaured", 
              type = "l")

points(which.max(reg.summary$adjr2), # This code inserts a red dot on the line graph. The red dot denotes the point where the adjusted R-squared is at its highest.
       reg.summary$adjr2[which.max(reg.summary$adjr2)], 
       col="red",
       cex=2,
       pch=20)

plot(reg.summary$rss, # sum of squares residual values on y-axis
     xlab = "Number of variables", 
     ylab = "RSS", 
     type = "l")

plot(reg.summary$bic, # BIC values on y-axis. BIC is very similar to AIC except that it puts more penalty to the number of variables included.
     xlab = "Number of variables", 
     ylab = "BIC", 
     type = "l")

points(which.min(reg.summary$bic), reg.summary$bic[which.min(reg.summary$bic)], col="red",cex=2,pch=20)



######################optimal regression in nonviolent model
#BIC lowest model (7 variables)
nonvio_ols_v1 <- lm(I(nonvio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100) 
                 + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)
#Adjusted R-Squared highest model (8 variables)
nonvio_ols_v2 <- lm(I(nonvio_crimerate*100) ~ pop_den + I(nonwhite_ratio * 100) +  I(less_than_hs_ratio * 100) 
                    + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)
AIC(nonvio_ols_v1)
AIC(nonvio_ols_v2)
BIC(nonvio_ols_v1)
BIC(nonvio_ols_v2)

summary(nonvio_ols_v1)
summary(nonvio_ols_v2)
#AIC and R2 are both better on model v2


##Check Heteroskedasticity
dev.off()
plot(predict(nonvio_ols_v2), 
     nonvio_ols_v1$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(nonvio_ols_v1) 
#Chisquare = 113.1714, Shows strong heteroskedasticity!!! So, I generalized nonvio_crimerate


nonvio_ols_v2_log <- lm(I(log(nonvio_crimerate*100 + 0.001)) ~ pop_den + I(nonwhite_ratio * 100) +  I(less_than_hs_ratio * 100) 
                    + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)

summary(nonvio_ols_v2_log)

plot(predict(nonvio_ols_v2_log), 
     nonvio_ols_v2_log$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(nonvio_ols_v2_log) 
# After generalizing the dependent variable it looks good!

##Check Multicollinearity
vif_v2_log <- vif(nonvio_ols_v2_log)
barplot(vif_v2_log, main = "VIF Values", horiz = TRUE, col = "steelblue") 
 
#All values below 2, looks good!