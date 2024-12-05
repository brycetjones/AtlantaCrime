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
library(gridExtra)



getwd()
setwd("C:/Users/Hojung Yu/Documents/GitHub/AtlantaCrime/")

df_org <- read.csv("data_final.csv")

head(df_org)
nrow(df_org)
ncol(df_org)
summary(df_org)
colnames(df_org)

df <- df_org

#Normalize Variables between 0 and 100.
summary(df_org)
df["pop_den"] = df_org["pop_den"] * 1000
df["nonwhite_ratio"] = (df_org['black_ratio'] + df_org['other_ratio'])*100
df["median_incomeE"] = df_org["median_incomeE"]*0.0001
df["less_than_hs_ratio"] = df_org["less_than_hs_ratio"]*100
df["vio_crimerate"] = df_org["vio_crimerate"] * 100
df["nonvio_crimerate"] = df_org["nonvio_crimerate"]*100
df["min_station_dist"] = df_org["min_station_dist"] * 0.01

summary(df)

df <- df %>% relocate("nonwhite_ratio", .before = "median_incomeE")
reg_variables = c("pop_den", "nonwhite_ratio", "median_incomeE",
                  "less_than_hs_ratio", "Commercial", "HighdensityResidential","Industrial",
                  "Institutional", "LowdensityResidential","ResidentialCommercial",
                  "min_station_dist", "vio_crimerate","nonvio_crimerate")

ggpairs(df[,reg_variables])

plots <- list()
for (var in reg_variables) {
  # Create a histogram for each variable
  p <- ggplot(df, aes_string(x = var)) +
    geom_histogram(fill = "#D79A1E", color = "black", alpha = 0.7) +
    ggtitle(paste("Histogram of", var)) +
    xlab(var) +
    ylab("Frequency") +
    theme_minimal()
  
  # Store the plot in the list
  plots[[var]] <- p
}

# View the plots (e.g., display the plot for the first variable)
print(plots[[reg_variables[1]]])
grid.arrange(grobs = plots, ncol = 4)


#####################################################Violent crime model

cor_matrix <- cor(df[reg_variables])
round(cor_matrix,2)

#1. Normal OLS
ols_1 <- lm(vio_crimerate ~ pop_den + nonwhite_ratio +
              median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
              LowdensityResidential + ResidentialCommercial + min_station_dist, data = df)
summary(ols_1)

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
ols_1_noout <- lm(vio_crimerate ~ pop_den + nonwhite_ratio +
                    median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                    LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout_vio)

summary(ols_1_noout)

# Compare no outliers to outliers using scaled coefficient plots
dev.off()
plot_summs(ols_1, ols_1_noout, scale = TRUE) +
  ggtitle("Comparison of Regression Models",
          subtitle = "With Outliers vs Without Outliers") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  theme(plot.subtitle = element_text(size = 12, hjust = 0.5, color = "blue"))

AIC(ols_1, ols_1_noout)

# Not conducting log-regression since it is not significant with these datasets.
# Using df_noout_vio without outliers since it has smaller AIC score.
# Variables Selection

##Stepwise
null.model <- lm(vio_crimerate ~ 1, data = df_noout_vio)
full.model <- lm(vio_crimerate ~ pop_den + nonwhite_ratio
                 + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
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
vio_ols_v1 <- lm(vio_crimerate ~ pop_den + nonwhite_ratio +  less_than_hs_ratio 
                 + Industrial + Institutional + LowdensityResidential + min_station_dist, data = df_noout_vio)
#Adjusted R-Squared highest model
vio_ols_v2 <- lm(vio_crimerate ~ pop_den + nonwhite_ratio +  less_than_hs_ratio 
                 + Industrial + Institutional + min_station_dist
                 + median_incomeE + Commercial + HighdensityResidential + ResidentialCommercial, data = df_noout_vio)

summary(vio_ols_v1)
summary(vio_ols_v2)

AIC(vio_ols_v1, vio_ols_v2)

#Since Vio_ols_v1 AIC is similar to v2 and R-squared is similar to each other
#vio_ols_v1 is more simple. Accept vio_ols_v1
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

###1. Log
vio_ols_v1_log <- lm(log(vio_crimerate+1) ~ pop_den + nonwhite_ratio +  less_than_hs_ratio 
                     + Industrial + Institutional + LowdensityResidential + min_station_dist, data = df_noout_vio)

summary(vio_ols_v1_log) 
#Dependent Variable log formation's regression model is different from the original one
#Saying that less_than_hs_ratio is insignificant. 

plot(predict(vio_ols_v1_log), 
     vio_ols_v1_log$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(vio_ols_v1_log) 

###2. Higher term

#Normal Term
vio_ols_v1 
#Quadratic terms
vio_ols_v2 <- lm(vio_crimerate ~ 
                   pop_den + pop_den^2 + 
                   nonwhite_ratio + nonwhite_ratio^2 + 
                   less_than_hs_ratio+ less_than_hs_ratio^2 +
                   Industrial + I(Industrial^2) +
                   Institutional + I(Institutional^2) +
                   LowdensityResidential + I(LowdensityResidential^2) +
                    min_station_dist + I((min_station_dist)^2), 
                 data = df_noout_vio)

#Cubic Terms
vio_ols_v3 <- lm(vio_crimerate ~ 
                   pop_den + pop_den^2 + + pop_den^3 +
                   nonwhite_ratio + nonwhite_ratio^2 + nonwhite_ratio^3 + 
                   less_than_hs_ratio + less_than_hs_ratio^2 + less_than_hs_ratio^3 + 
                   Industrial + I(Industrial^2) + I(Industrial^3) + 
                   Institutional + I(Institutional^2) + I(Institutional^3) + 
                   LowdensityResidential + I(LowdensityResidential^2) + I(LowdensityResidential^3) + 
                   I(min_station_dist) + I((min_station_dist)^2) + I((min_station_dist)^3), 
                 data = df_noout_vio)

AIC(vio_ols_v1,vio_ols_v2,vio_ols_v3)
plot_summs(vio_ols_v1, vio_ols_v2,vio_ols_v3, scale = TRUE)

plot_summs(vio_ols_v1, vio_ols_v2,vio_ols_v3, scale = TRUE) +
  ggtitle("Comparison of Regression Models",
          subtitle = "Normal vs. Quadratic vs. Cubic") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  theme(plot.subtitle = element_text(size = 12, hjust = 0.5, color = "blue"))

 

# # Stepwise selection based on AIC
# step_model <- stepAIC(vio_ols_v2, direction = "both")
# AIC(vio_ols_v2,step_model)
# 
# # Check the summary of the reduced model
# summary(step_model)
# summary(vio_ols_v1)
# vif(step_model) # Multicollinearity on pop_den
# 
# # Centering variables
# df_noout_vio$pop_den_centered <- scale(df_noout_vio$pop_den, center = TRUE, scale = FALSE)
# 
# # Use centered variable in the model
# step_model <- lm(formula = vio_crimerate ~  pop_den_centered + I(pop_den_centered^2) + 
#      nonwhite_ratio^2 + less_than_hs_ratio + 
#      Industrial + Institutional + LowdensityResidential + I(min_station_dist^2), 
#    data = df_noout_vio)
# vif(step_model) # Multicollinearity on pop_den
# 

dev.off()
plot(predict(vio_ols_v2), 
     vio_ols_v2$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(vio_ols_v2) 
ncvTest(vio_ols_v1) #decrease very small


#Considering interactive term(Yanfu)
#for non_white ratio and population density, there might be an interaction in how they 
#influence the crime rate.

#check pop_den and non-white for interaction
range(df_noout_vio$nonwhite_ratio)
range(df_noout_vio$pop_den)

#low and high value subsets for non_white
nonwhite_low <- df_noout_vio[df_noout_vio$nonwhite_ratio <= 20, ]
nonwhite_high <- df_noout_vio[df_noout_vio$nonwhite_ratio >= 80, ]
nw_low_mod <- lm(I(log(vio_crimerate + 1)) ~ pop_den, data = nonwhite_low)
summary(nw_low_mod)
nw_high_mod <- lm(I(log(vio_crimerate + 1)) ~ pop_den, data = nonwhite_high)
summary(nw_high_mod)

par(mfrow = c(1,1))
plot(x = nonwhite_low$pop_den, y = log(nonwhite_low$vio_crimerate+1), 
     pch = 19, xlab = "population density", ylab = "log_crimerate", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25)) +
  abline(nw_low_mod, col = "blue", lwd = 2) +
  points(x = nonwhite_high$pop_den, log(nonwhite_high$vio_crimerate+1), 
         col = rgb(red = 1, green = 0, blue = 0, alpha = 0.25), pch = 19) +
  abline(nw_high_mod, col = "red", lwd = 2)
#from this plot we can see the slopes of low non white and high non white are different. There is
#an interaction between non white and pop den.

#We can now add the interaction term to model
vio_ols_v1_log_interaction <- lm(I(log(vio_crimerate + 1)) ~ pop_den + nonwhite_ratio + pop_den:nonwhite_ratio
                                 + less_than_hs_ratio + Industrial + Institutional
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
vio_glm_v1 <- glm(binary_vio ~ pop_den + nonwhite_ratio + 
                    less_than_hs_ratio + Industrial + Institutional + 
                    LowdensityResidential + min_station_dist, data = df_noout_vio, family = "binomial"(link = "logit"))

# binary logical regression
vio_glm_v2 <- glm(binary_vio ~ pop_den + nonwhite_ratio + 
                    less_than_hs_ratio + Institutional +
                    LowdensityResidential, data = df_noout_vio, family = "binomial"(link = "logit"))
summary(vio_glm_v1)
summary(vio_glm_v2)

vif(vio_glm_v1)
vif(vio_glm_v2)
# creating odds ratio
# Regression result (with odds ratio conversion)


round( # Rounds the numbers up to 3 digits
  cbind( # Column-bind Odds Ratio to the regerssion output
    "Odds Ratio" = exp(vio_glm_v2$coefficients),
    summary(vio_glm_v2)$coefficients
  ),3)

# pseudo R-squared 0.4361733
pR2(vio_glm_v1)[4]
pR2(vio_glm_v2)[4]

# Wes hould select vio_glm_v2

# plotting the logistic curve

# population density with violent crime
ggplot(df_noout_vio, aes(pop_den, binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Population density") +
  ylab("Probability of violent crime")

ggplot(df_noout_vio, aes(nonwhite_ratio, binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Non white ratio") +
  ylab("Probability of violent crime")

ggplot(df_noout_vio, aes(less_than_hs_ratio, binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Less than HIghschool's Degree") +
  ylab("Probability of violent crime")

ggplot(df_noout_vio, aes(Industrial, binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Industrial Land Use") +
  ylab("Probability of violent crime")

ggplot(df_noout_vio, aes(Institutional, binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Institutional Land Use") +
  ylab("Probability of violent crime")

ggplot(df_noout_vio, aes(LowdensityResidential, binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Lowdensity-Residential Land Use") +
  ylab("Probability of violent crime")


ggplot(df_noout_vio, aes(min_station_dist, binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Minimum Distance to Police Station") +
  ylab("Probability of violent crime")

AIC(vio_glm_v1,vio_ols_v1_log,vio_ols_v1,vio_ols_v1_log_interaction)

##Below 5 or 10, Since all VIF values are well below 5 or 10, 
#there is no evidence of significant multicollinearity in your model based on these results.

###############################################################
#############################################################NONViolent crime model
# 1. Normal OLS
outliers_regression_non <- lm(nonvio_crimerate ~ pop_den + nonwhite_ratio
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
outliers_regression_noout_non <- lm(nonvio_crimerate ~ pop_den + nonwhite_ratio
                                + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
                                  + LowdensityResidential + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)
summary(outliers_regression_noout_non)

#Plot the difference

plot_summs(outliers_regression_non, outliers_regression_noout_non, scale = TRUE) +
  ggtitle("Comparison of Regression Models",
          subtitle = "With Outliers vs Without Outliers") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  theme(plot.subtitle = element_text(size = 12, hjust = 0.5, color = "blue"))

AIC(outliers_regression_non)
AIC(outliers_regression_noout_non)

# Not conducting log-regression since it is not significant with these datasets.
# Using df_noout_NONvio without outliers since it has smaller AIC score.
# Variables Selection


##Stepwise
null.model <- lm(nonvio_crimerate ~ 1, data = df_noout_nonvio)
full.model <- lm(nonvio_crimerate ~ pop_den + nonwhite_ratio
                 + median_incomeE + less_than_hs_ratio + Commercial + HighdensityResidential + Industrial + Institutional +
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
nonvio_ols_v1 <- lm(nonvio_crimerate ~ pop_den + nonwhite_ratio 
                 + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)
#Adjusted R-Squared highest model (8 variables)
nonvio_ols_v2 <- lm(nonvio_crimerate ~ pop_den + nonwhite_ratio +  less_than_hs_ratio 
                    + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)
AIC(nonvio_ols_v1, nonvio_ols_v2)

summary(nonvio_ols_v1)
summary(nonvio_ols_v2)
#AIC and R2 are both better, but simple model v1


##Check Heteroskedasticity
dev.off()
plot(predict(nonvio_ols_v1), 
     nonvio_ols_v1$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(nonvio_ols_v1) 
#Chisquare = 122.9289, Shows strong heteroskedasticity!!! So, I generalized nonvio_crimerate

nonvio_ols_v1_log <- lm(I(log(nonvio_crimerate + 1)) ~ pop_den + nonwhite_ratio +
                    + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)

summary(nonvio_ols_v1_log)

plot(predict(nonvio_ols_v1_log), 
     nonvio_ols_v1_log$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)
ncvTest(nonvio_ols_v1_log) 
# After generalizing the dependent variable it looks good!

##Check Multicollinearity
vif_v1_log <- vif(nonvio_ols_v1_log)
barplot(vif_v1_log, main = "VIF Values", horiz = TRUE, col = "steelblue") 
 
#All values below 2, looks good!

#Considering interactive term(Yanfu)
#for non_white ratio and population density, there might be an interaction in how they 
#influence the crime rate.

#check pop_den and non-white for interaction
range(df_noout_nonvio$nonwhite_ratio)
range(df_noout_nonvio$pop_den)

#low and high value subsets for non_white
nonwhite_low_1 <- df_noout_nonvio[df_noout_nonvio$nonwhite_ratio <= 20, ]
nonwhite_high_1 <- df_noout_nonvio[df_noout_nonvio$nonwhite_ratio >= 80, ]
nw_low_mod_1 <- lm(I(log(vio_crimerate + 1)) ~ pop_den, data = nonwhite_low_1)
summary(nw_low_mod_1)
nw_high_mod_1 <- lm(I(log(vio_crimerate + 1)) ~ pop_den, data = nonwhite_high_1)
summary(nonwhite_high_1)

par(mfrow = c(1,1))
plot(x = nonwhite_low_1$pop_den, y = log(nonwhite_low_1$vio_crimerate+1), 
     pch = 19, xlab = "population density", ylab = "log_crimerate", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25)) +
  abline(nw_low_mod_1, col = "blue", lwd = 2) +
  points(x = nonwhite_high_1$pop_den, log(nonwhite_high_1$vio_crimerate+1), 
         col = rgb(red = 1, green = 0, blue = 0, alpha = 0.25), pch = 19) +
  abline(nw_high_mod_1, col = "red", lwd = 2)
#from this plot we can see the slopes of low non white and high non white are different. There is
#an interaction between non white and pop den.

#We can now add the interaction term to model
lm(nonvio_crimerate ~ pop_den + nonwhite_ratio 
   + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)
nonvio_ols_v1_log_interaction <- lm(I(log(nonvio_crimerate + 1)) ~ pop_den + nonwhite_ratio + pop_den:nonwhite_ratio
                                    + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio)

summary(nonvio_ols_v1_log_interaction)

plot(predict(nonvio_ols_v1_log_interaction), 
     nonvio_ols_v1_log_interaction$residuals, 
     xlab = "Predicted values",
     ylab = "Residuals",
     main = "Residuals vs. Predicted")
abline(h = 0, col = 'red', lty = 2)

AIC(nonvio_ols_v1_log_interaction,nonvio_ols_v1_log)


# Considering glm function for non_vio
# create binary column for violent crime rate
median_noncrimerate <- median(df_noout_nonvio$nonvio_crimerate, na.rm = TRUE)
df_noout_nonvio$binary_nonvio <- ifelse(df_noout_nonvio$nonvio_crimerate > median_noncrimerate, 1, 0)

# binary logical regression
non_vio_glm_v1 <- glm(binary_nonvio ~ pop_den + nonwhite_ratio +
                      + Commercial + HighdensityResidential + Industrial + ResidentialCommercial + min_station_dist, data = df_noout_nonvio, family = "binomial"(link = "logit"))



summary(non_vio_glm_v1)

vif(non_vio_glm_v1)
# creating odds ratio
# Regression result (with odds ratio conversion)
round( # Rounds the numbers up to 3 digits
  cbind( # Column-bind Odds Ratio to the regerssion output
    "Odds Ratio" = exp(non_vio_glm_v1$coefficients),
    summary(non_vio_glm_v1)$coefficients
  ),3)

# pseudo R-squared 0.4361733
pR2(vio_glm_v1)[4]

# plotting the logistic curve

# population density with violent crime
ggplot(df_noout_nonvio, aes(pop_den, binary_nonvio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Population density") +
  ylab("Probability of nonviolent crime")

ggplot(df_noout_nonvio, aes(I(nonwhite_ratio), binary_vio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Non white ratio") +
  ylab("Probability of nonviolent crime")


ggplot(df_noout_nonvio, aes(nonwhite_ratio, binary_nonvio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Non white ratio") +
  ylab("Probability of nonviolent crime")

ggplot(df_noout_nonvio, aes(Commercial, binary_nonvio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Commercial Land Use") +
  ylab("Probability of nonviolent crime")

ggplot(df_noout_nonvio, aes(HighdensityResidential, binary_nonvio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("HighdensityResidential Land Use") +
  ylab("Probability of nonviolent crime")

ggplot(df_noout_nonvio, aes(Industrial, binary_nonvio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Industrial Land Use") +
  ylab("Probability of nonviolent crime")

ggplot(df_noout_nonvio, aes(ResidentialCommercial, binary_nonvio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("ResidentialCommercial Land Use") +
  ylab("Probability of nonviolent crime")

ggplot(df_noout_nonvio, aes(min_station_dist, binary_nonvio)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("Minimum Distance to Police Station") +
  ylab("Probability of nonviolent crime")

AIC(non_vio_glm_v1, nonvio_ols_v1_log, nonvio_ols_v1)