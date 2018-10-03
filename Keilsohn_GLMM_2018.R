### William Keilsohn
### (Generalized) Mixed Linear Models

### Load Packages
library(nlme)
library(dplyr)

BACD$concentration<-as.factor(BACD$concentration)

# Degregation of medium as a response to bacterial strain

dotchart(BACD$degradation, groups = factor(BACD$string),
         ylab = "Strain", xlab = "Degradation (ml/L/day)",
         main = "Cleveland dotplot")

boxplot(BACD$degradation~BACD$string, 
        ylab = "Degradation (mg/L/day)",
        xlab = "Strain",
        main = "Bar chart")



pairs(BACD)

# Bacterial consentration in medium

M1<-lm(degradation~concentration * string,
       data = BACD)
hist(resid(M1), ylab = "Residuals")

#op<-par(mfrow = c(1,2))
plot(resid(M1)~BACD$concentration, type = "b",
     xlab = "Concentration (bacteria/ul)", ylab = "Residuals")
plot(resid(M1)~BACD$string, type = "b",
     xlab = "Strain", ylab = "Residuals")
#par(op)

plot(x=fitted(M1), y = resid(M1), xlab = "Fitted Values",
     ylab = "Residuals", main = "Model 1")

M2<-gls(degradation~concentration * string,
        weights = varIdent(form = ~1|concentration*string),
        data = BACD)
op2<-par(mfrow = c(2,2), mar = c(5,4,2,2))
hist(resid(M2), ylab = "Residuals")
plot(resid(M2, type = c("pearson"))~BACD$concentration,
     xlab = "Concentration (bacteria/ul)", ylab = "Residuals")
plot(resid(M2, type = c("pearson"))~BACD$string,
     xlab = "Strain", ylab = "Residuals")
plot(x=fitted(M2), y = resid(M2, type = c("pearson")), xlab = "Fitted Values",
     ylab = "Residuals", main = "Model 2")
par(op2)

anova(M2)


# Bacterial responce to temperature and location

SND$Site<-as.factor(SND$Site)
SND$Density<-as.integer(SND$Density)

M3<-lme(Density~Temperature, data = SND, random = ~1|Site)
summary(M3)

hist(resid(M3), ylab = "Rediduals", xlab = "Temperature (degrees C)")


op3<-par(mfrow = c(1,2))
boxplot(resid(M3)~factor(SND$Temperature), type = "b",
     xlab = "Temperature (degrees C)", ylab = "Residuals")
plot(resid(M3)~SND$Site, type = "b",
     xlab = "Site", ylab = "Residuals")
par(op3)

plot(x=fitted(M3), y = resid(M3), xlab = "Fitted Values",
     ylab = "Residuals", main = "Model 3")


M4<-gls(Density~Temperature, data = SND,
        weights = varIdent(form = ~1|Site))
summary(M4)

op2<-par(mfrow = c(2,2), mar = c(5,4,2,2))
hist(resid(M4), ylab = "Residuals")
plot(resid(M4, type = c("pearson"))~SND$Temperature,
     xlab = "Temperature (Degrees C)", ylab = "Residuals")
plot(resid(M4, type = c("pearson"))~SND$Site,
     xlab = "Site", ylab = "Residuals")
plot(x=fitted(M4), y = resid(M4, type = c("pearson")), xlab = "Fitted Values",
     ylab = "Residuals", main = "Model 4")
par(op2)

anova(M3)
anova(M3)

### Second attempt at previous model

M5<-lme(Density~Temperature, data = SND, random = ~1|Site,
        weights = varIdent(form = ~1|Site))
summary(M5)

op2<-par(mfrow = c(2,2), mar = c(5,4,2,2))
hist(resid(M5), ylab = "Residuals")
plot(resid(M5, type = c("pearson"))~SND$Temperature,
     xlab = "Temperature (Degrees C)", ylab = "Residuals")
plot(resid(M5, type = c("pearson"))~SND$Site,
     xlab = "Site", ylab = "Residuals")
plot(x=fitted(M5), y = resid(M5, type = c("pearson")), xlab = "Fitted Values",
     ylab = "Residuals", main = "Model 5")
par(op2)


anova(M5)
summary(M5)