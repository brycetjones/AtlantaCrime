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

getwd()
setwd("C:/Users/Hojung Yu/Documents/GitHub/AtlantaCrime/")

df <- read.csv("data_final.csv")

head(df)
nrow(df)
ncol(df)
summary(df)

colnames(df)

ggpairs(df[,reg_variables])

#####################################################Violent crime model
reg_variables = c("pop_den", "black_ratio", "white_ratio","other_ratio", "median_incomeE", "Commercial", 
              "HighdensityResidential", "Industrial","Institutional",
              "LowdensityResidential", "ResidentialCommercial",
              "min_station_dist", "violent_ratio", "nonviolent_ratio")

cor_matrix <- cor(df[reg_variables])
round(cor_matrix,2)

#1. without logs
outliers_regression <- lm(I(violent_ratio*1000 + 0.001) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                          + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                            + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df)
summary(outliers_regression)

cooks_distance_nolog <- cooks.distance(outliers_regression)
df$cooks_distance_nolog <- cooks.distance(outliers_regression)

plot(cooks_distance_nolog, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean (cooks_distance_nolog, na.rm = T), col = "red")
text(x = 1: length(cooks_distance_nolog) + 5,
     y = cooks_distance_nolog,
     col = "red",
     labels = ifelse(cooks_distance_nolog > 4 * mean(cooks_distance_nolog, na.rm = T),
                     names(cooks_distance_nolog),
                     ""))

#Remove influential observations
noout <- subset(df[df$cooks_distance_nolog < .20, ])
outliers_regression_noout <- lm(I(violent_ratio*1000 + 0.001) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                                + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                                  + LowdensityResidential + ResidentialCommercial + min_station_dist, data = noout)
summary(outliers_regression_noout)

#Plot the difference
plot_summs(outliers_regression, outliers_regression_noout, scale = TRUE)

#2. with logs
outlier.reg <- lm(I(log(violent_ratio*1000 + 0.001)) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                  + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                    + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df)
summary(outlier.reg)

cook.dist <- cooks.distance(outlier.reg)
df$cook.dist <- cooks.distance(outlier.reg)
plot(cook.dist, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean (cook.dist, na.rm = T), col = "red")
text(x = 1: length(cook.dist) + 5,
     y = cook.dist,
     col = "red",
     labels = ifelse(cook.dist > 4 * mean(cook.dist, na.rm = T),
                     names(cook.dist),
                     ""))

#If we want to remove outliers...
df_noout <- subset(df[df$cook.dist <.20, ])
outlier.reg_noout <- lm(I(log(violent_ratio*1000 + 0.001)) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                        + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                          + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout)
summary(outlier.reg_noout)

#Compare no outliers to outliers using scaled coefficient plots
#Will only let me run if scale = FALSE...not if it = TRUE
plot_summs(outlier.reg, outlier.reg_noout, scale = FALSE)


# #3. Real Regression_with logs / Generalized violent ratio and black_ratio, lowdensity
# model_1 <- lm(I(log(violent_ratio*1000 + 0.001)) ~ I(black_ratio*100) +
#               LowdensityResidential, data = df)
# summary(model_1)
# ncvTest(model_1)

##Stepwise
null.model <- lm(I(log(violent_ratio*1000 + 0.001)) ~ 1, data = df)
full.model <- lm(I(log(violent_ratio*1000 + 0.001)) ~ pop_den + median_incomeE + I(black_ratio*100) + median_incomeE + Commercial + 
                   LowdensityResidential + min_station_dist + I(less_than_hs_ratio*100), data = df)

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
                           nvmax = 5,
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
subsetted.model_v1 <- lm(I(log(violent_ratio*1000 + 0.001)) ~ I(black_ratio * 100) + Commercial + LowdensityResidential + I(less_than_hs_ratio * 100), data = df)
summary(subsetted.model_v1)
dev.off()

##Check Heteroskedasticity
plot(resid(subsetted.model_v1), main = "Residual Plot", xlab = "Values",ylab = "Residuals")
abline(h = 0, col = "red")
ncvTest(subsetted.model_v1)

##Check Multicollinearity
vif_values_v1 <- vif(subsetted.model_v1)
barplot(vif_values_v1, main = "VIF Values", horiz = TRUE, col = "steelblue") 
##Below 5 or 10, Since all VIF values are well below 5 or 10, 
#there is no evidence of significant multicollinearity in your model based on these results.

#############################################################NONViolent crime model
#1. without logs
outliers_regression_non <- lm(I(nonviolent_ratio*1000 + 0.001) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                          + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                            + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df)
summary(outliers_regression_non)

cooks_distance_nolog_non <- cooks.distance(outliers_regression_non)
df$cooks_distance_nolog_non <- cooks.distance(outliers_regression_non)

plot(cooks_distance_nolog_non, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean (cooks_distance_nolog_non, na.rm = T), col = "red")
text(x = 1: length(cooks_distance_nolog_non) + 5,
     y = cooks_distance_nolog_non,
     col = "red",
     labels = ifelse(cooks_distance_nolog_non > 4 * mean(cooks_distance_nolog_non, na.rm = T),
                     names(cooks_distance_nolog_non),
                     ""))

#Remove influential observations
noout_non <- subset(df[df$cooks_distance_nolog_non < .20, ])
outliers_regression_noout_non <- lm(I(nonviolent_ratio*1000 + 0.001) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                                + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                                  + LowdensityResidential + ResidentialCommercial + min_station_dist, data = noout_non)
summary(outliers_regression_noout_non)

#Plot the difference
plot_summs(outliers_regression_non, outliers_regression_noout_non, scale = TRUE)

#2. with logs
outlier.reg.non <- lm(I(log(nonviolent_ratio*1000 + 0.001)) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                  + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                    + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df)
summary(outlier.reg.non)

cook.dist.non <- cooks.distance(outlier.reg.non)
df$cook.dist.non <- cooks.distance(outlier.reg.non)
plot(cook.dist.non, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean (cook.dist.non, na.rm = T), col = "red")
text(x = 1: length(cook.dist.non) + 5,
     y = cook.dist.non,
     col = "red",
     labels = ifelse(cook.dist.non > 4 * mean(cook.dist, na.rm = T),
                     names(cook.dist.non),
                     ""))

#If we want to remove outliers...
df_noout_non <- subset(df[df$cook.dist.non <.20, ])
outlier.reg_noout_non <- lm(I(log(nonviolent_ratio*1000 + 0.001)) ~ pop_den + I(black_ratio*100) + I(white_ratio*100) + I(other_ratio*100)
                        + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                          + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout_non)
summary(outlier.reg_noout_non)

#Compare no outliers to outliers using scaled coefficient plots
#Will only let me run if scale = FALSE...not if it = TRUE
plot_summs(outlier.reg.non, outlier.reg_noout_non, scale = FALSE)

##Stepwise
null.model <- lm(I(log(nonviolent_ratio*1000 + 0.001)) ~ 1, data = df)
full.model <- lm(I(log(nonviolent_ratio*1000 + 0.001)) ~ pop_den + median_incomeE + I(black_ratio*100) +
                   Commercial + LowdensityResidential + ResidentialCommercial + less_than_hs_ratio + min_station_dist, data = df)

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
                           nvmax = 5,
                           method = "exhaustive")



##Best subset
full.model_n <- lm(I(log(nonviolent_ratio*1000 + 0.001)) ~ pop_den + median_incomeE + I(black_ratio*100) +
                     Commercial + LowdensityResidential + ResidentialCommercial + less_than_hs_ratio + min_station_dist, data = df)

subset.model_n <- regsubsets(formula(full.model_n),
                             data = df,
                             nvmax = 6,
                             method = "exhaustive")

reg.summary_n <- summary(subset.model_n)
reg.summary_n

par(mfrow = c(2,2)) # This code specifies that, instead of showing only one graph in the plots window, you want 4 (2 x 2) graphs.

plot(reg.summary_n$rsq, # This code plots the number of variables on x-axis and r-squared values on y-axis
     xlab = "Number of variables", 
     ylab = "R-Squared", 
     type = "l") # type = "l" means that you want to display the graph with a line, not points.

plot2 <- plot(reg.summary_n$adjr2, # Adjusted r-squared values on y-axis
              xlab = "Number of variables", 
              ylab = "Adjusted R-Sqaured", 
              type = "l")

points(which.max(reg.summary_n$adjr2), # This code inserts a red dot on the line graph. The red dot denotes the point where the adjusted R-squared is at its highest.
       reg.summary_n$adjr2[which.max(reg.summary_n$adjr2)], 
       col="red",
       cex=2,
       pch=20)

plot(reg.summary_n$rss, # sum of squares residual values on y-axis
     xlab = "Number of variables", 
     ylab = "RSS", 
     type = "l")

plot(reg.summary_n$bic, # BIC values on y-axis. BIC is very similar to AIC except that it puts more penalty to the number of variables included.
     xlab = "Number of variables", 
     ylab = "BIC", 
     type = "l")

points(which.min(reg.summary_n$bic), reg.summary_n$bic[which.min(reg.summary_n$bic)], col="red",cex=2,pch=20)

######################optimal regression in nonviolent model
subsetted.model_v1_n <- lm(I(log(nonviolent_ratio*1000 + 0.001)) ~ pop_den + Commercial + LowdensityResidential
                           + ResidentialCommercial + less_than_hs_ratio, data = df)
summary(subsetted.model_v1_n)
######################Check
##Check Heteroskedasticity
dev.off()
plot(resid(subsetted.model_v1_n), main = "Residual Plot", xlab = "Values",ylab = "Residuals")
abline(h = 0, col = "red")
ncvTest(subsetted.model_v1_n)

##Check Multicollinearity
vif_values_v1_n <- vif(subsetted.model_v1_n)
barplot(vif_values_v1_n, main = "VIF Values", horiz = TRUE, col = "steelblue") 
##Below 5 or 10, Since all VIF values are well below 5 or 10, 
#there is no evidence of significant multicollinearity in your model based on these results.
