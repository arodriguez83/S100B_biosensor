# S100B_biosensor
Thisnotebook analyses the data from the biosensor experiment lead by Alexander Rodríguez
This is an R notebook

#Read data AuEs-PBS
PBS <- read.csv("PBS.csv") 
Ag_PBS <- as.factor(PBS$Niveles)
dRct_PBS<- PBS$dRct

# Normalization
PBS$y = log(dRct_PBS)
PBS$x = Ag_PBS

boxplot(PBS$y ~ PBS$x, col=456, xlab = "[S100B](pg/mL)", ylab= " logdRct", main = "S100B test in PBS (AuEs)")

require(MASS)
AjustePBS<-fitdistr(PBS$y, "normal")
AjustePBS
###Kolmogorov-Smirnov AuEs-PBS
KsnPBS<- ks.test(PBS$y, "pnorm", mean =AjustePBS$estimate[1], sd= AjustePBS$estimate[2])
KsnPBS
##Shapiro AuEs-PBS
SwnPBS<-shapiro.test(PBS$y)
SwnPBS

#Exploratory AuEs-PBS
par(mfrow=c(1,3))
hist(PBS$y, xlab = "log_dRct", ylab = "Frequency", las=1, 
     main = "log_dRct (PBS & AuEs)", ylim = c(0,12), col = "gray")
boxplot(PBS$y, col = "gray", ylab="log_dRct")
plot(density(PBS$y), xlab = "log_dRct", ylab = "Density", las=1, main = "")

#Welch ANOVA AuEs-PBS
library(rstatix)
ANOVA_PBS <- PBS %>% welch_anova_test(y ~ x)
ANOVA_PBS
summary(ANOVA_PBS)

#Post hoc analysis
#Games Howell AuEs-PBS
games_howell_test(PBS, y ~ Niveles, conf.level = 0.95, detailed = FALSE)

Modelo1_lm <- lm(y ~ x, data = PBS)
plot(Modelo1_lm, which = c(1:3))

#Read data AuEs-Plasma
Plasma <- read.csv("Plasma.csv") 
Ag_Plasma <- as.factor(Plasma$Niveles)
dRct_Plasma<- Plasma$dRct

# Normalization
Plasma$y = log(dRct_Plasma)
Plasma$x = Ag_Plasma

par(mfrow=c(1,1))
boxplot(Plasma$y ~ Plasma$x, col=456, xlab = "[S100B](pg/mL)", ylab= " logdRct", main = "S100B test in Plasma (AuEs)")

##EstimaciÃ³n los parÃ¡metros de la distribuciÃ³n hipotÃ©tica (normal)
AjustePlasma<-fitdistr(Plasma$y, "normal")
AjustePlasma
###Kolmogorov-Smirnov Plasma
KsnPlasma<- ks.test(Plasma$y, "pnorm", mean =AjustePlasma$estimate[1], sd= AjustePlasma$estimate[2])
KsnPlasma
##Shapiro Plasma
SwnPlasma<-shapiro.test(Plasma$y)
SwnPlasma

#Exploratory AuEs-Plasma
par(mfrow=c(1,3))
hist(Plasma$y, xlab = "log_dRct", ylab = "Frequency", las=1, 
     main = "log_dRct (Plasma & AuEs)", ylim = c(0,10), col = "gray")
boxplot(Plasma$y, col = "gray", ylab="log_dRct")
plot(density(Plasma$y), xlab = "log_dRct", ylab = "Density", las=1, main = "")

#Welch ANOVA Plasma
library(rstatix)
ANOVA_Plasma <- Plasma %>% welch_anova_test(y ~ x)
ANOVA_Plasma
summary(ANOVA_Plasma)

#Post hoc analysis AuEs-plasma
##Games Howell Plasma
games_howell_test(Plasma, y ~ Niveles, conf.level = 0.95, detailed = FALSE)

Modelo2_lm <- lm(y ~ x, data = Plasma)
plot(Modelo2_lm, which = c(1:3))


#Read data IDEs 
IDE2 <- read.csv("IDEdRct.csv") 
Ag_IDE2 <- as.factor(IDE2$Niveles)
dRct_IDE2<- IDE2$drct

par(mfrow=c(1,1))
boxplot(dRct_IDE2 ~ Ag_IDE2, col=456, xlab = "[S100B](pg/mL)", ylab= "dRct", main = "S100B test in Plasma (AuIDEs)")

##EstimaciÃ³n los parÃ¡metros de la distribuciÃ³n hipotÃ©tica (normal)
require(MASS)
AjusteIDE2<-fitdistr(IDE2$drct, "normal")
AjusteIDE2
###Kolmogorov-Smirnov IDEs 
KsnIDE2<- ks.test(IDE2$drct, "pnorm", mean =AjusteIDE2$estimate[1], sd= AjusteIDE2$estimate[2])
KsnIDE2

##Shapiro IDEs 
SwnIDE2<-shapiro.test(dRct_IDE2)
SwnIDE2

#Exploratory IDEs dRCt
par(mfrow=c(1,3))
hist(IDE2$drct, xlab = "log_dRct", ylab = "Frequency", las=1, 
     main = "log_dRct (Plasma & IDEs)", ylim = c(0,12), col = "gray")
boxplot(IDE2$drct, col = "gray", ylab="log_dRct")
plot(density(IDE2$drct), xlab = "log_dRct", ylab = "Density", las=1, main = "")

#Homoscedasticity IDEs 
library(car)
leveneTest(dRct_IDE2 ~ Ag_IDE2, data = IDE2, center = "median")

bartlett.test(dRct_IDE2 ~ Ag_IDE2, data = IDE2 )

##ANOVA IDEs 
anova_IDE<- aov(dRct_IDE2 ~ Ag_IDE2)
summary(anova_IDE)

##Tukey IDEs
TukeyHSD(anova_IDE)
par(mfrow=c(1,1))
plot(TukeyHSD(anova_IDE))

#Welch ANOVA IDEs
library(rstatix)
ANOVAw <- PBS %>% welch_anova_test(dRct_IDE2 ~ Ag_IDE2)
ANOVAw
summary(ANOVAw)

#Prueba Post hoc
##Games Howell IDEs
games_howell_test(IDE2, dRct_IDE2 ~ Niveles, conf.level = 0.95, detailed = FALSE)

#Q-Q plot IDEs
library(ggplot2)
ggplot(IDE2, aes(sample = dRct_IDE2, col = Niveles)) + stat_qq()
+ stat_qq_line() + facet_grid(.~ Niveles)

##cuantiles teÃ³ricos vs cuantiles muestrales IDEs
qqPlot(dRct_IDE2, xlab="Cuantiles teÃ³ricos", ylab="Cuantiles muestrales", 
       las=1,main="Q-Q plot dRct_IDEs")

Modelo_IDE_lm <- lm(dRct_IDE2 ~ Ag_IDE2, data = IDE2)
plot(Modelo_IDE_lm, which = c(1:3))
##DispersiÃ³n
g1 <- ggplot(IDE2, aes(x=Niveles, y = dRct_IDE2, color =  Ag_IDE2)) 
g1 + geom_point() + geom_smooth(method="lm")

#Read data controles
ctr <- read.csv("controles.csv") 
condition <- as.factor(ctr$Niveles)
dRct <- ctr$dRct
Rct <- ctr$Rct


##AnÃ¡lisis para Rct
boxplot(Rct ~ condition, col=456, xlab = "Controls", ylab= "Rct", main = "Control tests for AuEs")

require(MASS)
Ajuste_ctr<-fitdistr(ctr$Rct, "normal")
Ajuste_ctr
###Kolmogorov-Smirnov controles
Ksn_ctr<- ks.test(Rct, "pnorm", mean = Ajuste_ctr$estimate[1], 
                  sd= Ajuste_ctr$estimate[2])
Ksn_ctr

#Homoscedasticity controls
library(car)
leveneTest(Rct ~ condition, data = ctr, center = "median")

bartlett.test(Rct ~ condition, data = ctr )


#ANOVA ctr
anova_ctr<- aov(Rct ~ condition)
summary(anova_ctr)

##Tukey ctr
TukeyHSD(anova_ctr)
par(mfrow=c(1,1))
plot(TukeyHSD(anova_ctr))

#Read data controls-IDE
ctrIDE <- read.csv("controlesIDE.csv") 
condition_IDE <- as.factor(ctrIDE$Niveles)
dRct_IDE <- ctrIDE$dRct
Rct_IDE <- ctrIDE$Rct


##Analysis controls-IDE
boxplot(Rct_IDE ~ condition_IDE, col=456, xlab = "Controls", ylab= "Rct", main = "Control tests for AuIDEs")

require(MASS)
Ajuste_ctrIDE<-fitdistr(ctrIDE$Rct, "normal")
Ajuste_ctrIDE
###Kolmogorov-Smirnov controles IDE
Ksn_ctrIDE<- ks.test(Rct_IDE, "pnorm", mean = Ajuste_ctrIDE$estimate[1], 
                  sd= Ajuste_IDE$estimate[2])
Ksn_ctrIDE

#Homoscedasticity IDEs
library(car)
leveneTest(Rct_IDE ~ condition_IDE, data = ctrIDE, center = "median")

bartlett.test(Rct_IDE ~ condition_IDE, data = ctrIDE )


#ANOVA ctr
anova_ctrIDE<- aov(Rct_IDE ~ condition_IDE)
summary(anova_ctrIDE)

##Tukey ctr
TukeyHSD(anova_ctrIDE)
par(mfrow=c(1,1))
plot(TukeyHSD(anova_ctrIDE))
