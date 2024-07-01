#### Paper Sub-Antarctic Marine Isoscapes - Riccialdelli et al. 2024 - Progresss in Oceanography -under review ####

######### Generalized Least Square (GLS) models to compare between SIA values. Large, Regional and Local approaches ######
## Last modification June 2024 - Sami Dodino ####

rm(list=ls())
ls()

library(lme4)
library(lattice)
library(nlme)
library(geepack)
library(MuMIn)
library(sjPlot)
library(dplyr)
library(rstatix)


#########################################   LARGE SCALE  #########################################

table<-read.csv("baseline_data_GIS_8Junio2023_isoscapes_plusClima2.csv",header=T)

#######  Only zoo ######
large_zoo <- table %>% 
  group_by(d13Cs,d15N) %>% 
  filter(Subgroup=="ZP"& Season=="autumn")

str(large_zoo)
names(large_zoo)

#Evaluation of normality

# For  nitrogen
shapiro.test(large_zoo$d15N[large_zoo$MarineArea == "1_BC"])#normal
hist(large_zoo$d15N[large_zoo$MarineArea == "1_BC"])
shapiro.test(large_zoo$d15N[large_zoo$MarineArea == "2_CA"])
hist(large_zoo$d15N[large_zoo$MarineArea == "2_CA"])
shapiro.test(large_zoo$d15N[large_zoo$MarineArea == "3_BB"])#normal
shapiro.test(large_zoo$d15N[large_zoo$MarineArea == "4_PFZ"])#normal
shapiro.test(large_zoo$d15N[large_zoo$MarineArea == "5_PW"])#normal

large_zoo_d15N<-lm(d15N~MarineArea,data=large_zoo, na.action=na.fail)
large_zoo_d15N
plot(large_zoo_d15N)#qqplot normal

#For Carbon
shapiro.test(large_zoo$d13Cs[large_zoo$MarineArea == "1_BC"])#normal
hist(large_zoo$d13Cs[large_zoo$MarineArea == "1_BC"])
shapiro.test(large_zoo$d13Cs[large_zoo$MarineArea == "2_CA"])
hist(large_zoo$d13Cs[large_zoo$MarineArea == "2_CA"])
shapiro.test(large_zoo$d13Cs[large_zoo$MarineArea == "3_BB"])#normal
shapiro.test(large_zoo$d13Cs[large_zoo$MarineArea == "4_PFZ"])#normal
shapiro.test(large_zoo$d13Cs[large_zoo$MarineArea == "5_PW"])#normal

large_zoo_d13C<-lm(d13Cs~MarineArea,data=large_zoo, na.action=na.fail)
large_zoo_d13C
plot(large_zoo_d13C)#qqplot normal


# We need to model heteroscedasticity, we know that the groups are unbalanced
# So, we used generalized least square models (GLS) with varIdent variance function 

vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls1<-gls(d13Cs~MarineArea,weights = vf2, data=large_zoo, method="ML", na.action=na.fail)
summary(M.gls1)
dredge(M.gls1)
confint(M.gls1)

anova(M.gls1, test = "F")

range(large_zoo$d13Cs[large_zoo$MarineArea == "1_BC"])
range(large_zoo$d13Cs[large_zoo$MarineArea == "2_CA"])
range(large_zoo$d13Cs[large_zoo$MarineArea == "3_BB"])
range(large_zoo$d13Cs[large_zoo$MarineArea == "4_PFZ"])
range(large_zoo$d13Cs[large_zoo$MarineArea == "5_PW"])

sum(large_zoo$d13Cs & large_zoo$MarineArea == "1_BC")
sum(large_zoo$d13Cs & large_zoo$MarineArea == "2_CA")
sum(large_zoo$d13Cs & large_zoo$MarineArea == "3_BB")
sum(large_zoo$d13Cs & large_zoo$MarineArea == "4_PFZ")
sum(large_zoo$d13Cs & large_zoo$MarineArea == "5_PW")


## posthoc test for multiple comparisons 

large_zoo <- as.data.frame(large_zoo)

games_howell_test(data=large_zoo, d13Cs~MarineArea, conf.level = 0.95, detailed = FALSE)

# We now evaluate the variables that may be explaining these isotopic differences:
# Excluding latitude (correlation with SST) with mosaic depth

#evaluate correlations

cor.test(large_zoo$Lat,large_zoo$Chla_sat_2)# no correlac
cor.test(large_zoo$Long,large_zoo$Chla_sat_2)# no correlac

cor.test(large_zoo$Lat,large_zoo$SST_sat_2)# correlac, r=0.91
cor.test(large_zoo$Long,large_zoo$SST_sat_2)# no correlac

cor.test(large_zoo$SST_sat_2,large_zoo$Chla_sat_2)# no correlac

cor.test(large_zoo$Lat,large_zoo$depth_mosaico)# no correlac
cor.test(large_zoo$Long,large_zoo$depth_mosaico)# no correlac
cor.test(large_zoo$Chla_sat_2,large_zoo$depth_mosaico)# no correlac
cor.test(large_zoo$SST_sat_2,large_zoo$depth_mosaico)# no correlac

# There are rows with some NAs in chla and sst satellite,
# We create a new object that does not contain NAs and includes all the columns that interest us

large_zoo_no_na <- na.omit(large_zoo[, c("d15N","d13Cs", "SST_sat_2", "Long", "depth_mosaico", "Chla_sat_2", "MarineArea")])

M.gls1.2 <- gls(d13Cs ~ SST_sat_2 + Long + Chla_sat_2 + depth_mosaico, weights = vf2, data = large_zoo_no_na, method = "ML",na.action = na.fail)
summary(M.gls1.2)

dredge(M.gls1.2)
coef(M.gls1.2)
confint(M.gls1.2)

tab_model(M.gls1.2)

ModelSel.NB<-model.sel(M.gls112,M.gls111,rank="AICc")
ModelSel.NB
MMI.NB<-model.avg(M.gls112,M.gls111,beta=F,revised.var=F)
summary(MMI.NB)
coef(MMI.NB)
confint(MMI.NB)

#######  Only Eupha ######

large_EUFA <- table %>%
  filter(SP_red == "EUFA" & Season=="autumn")

#Nitrogen
shapiro.test(large_EUFA$d15N[large_EUFA$MarineArea == "1_BC"])#normal
shapiro.test(large_EUFA$d15N[large_EUFA$MarineArea == "2_CA"])#normal
shapiro.test(large_EUFA$d15N[large_EUFA$MarineArea == "3_BB"])#normal
shapiro.test(large_EUFA$d15N[large_EUFA$MarineArea == "4_PFZ"])#normal
shapiro.test(large_EUFA$d15N[large_EUFA$MarineArea == "5_PW"])#normal

large_eufa_d15N<-lm(d15N~MarineArea,data=large_EUFA, na.action=na.fail)
plot(large_eufa_d15N)#qqplot accept normality


#Carbon

shapiro.test(large_EUFA$d13Cs[large_EUFA$MarineArea == "1_BC"])#normal
shapiro.test(large_EUFA$d13Cs[large_EUFA$MarineArea == "2_CA"])#normal
shapiro.test(large_EUFA$d13Cs[large_EUFA$MarineArea == "3_BB"])#normal
shapiro.test(large_EUFA$d13Cs[large_EUFA$MarineArea == "4_PFZ"])#normal
shapiro.test(large_EUFA$d13Cs[large_EUFA$MarineArea == "5_PW"])#no normal, see qqplot

large_eufa_d13Cs<-lm(d13Cs~MarineArea,data=large_EUFA, na.action=na.fail)
plot(large_eufa_d13Cs)#qqplot accept normality


#### unbalanced groups, we model heterocedacea ###
vf2 <- varIdent(form= ~ 1 | MarineArea)

#Nitrogen  

M.gls2.1<-gls(d15N~MarineArea,weights = vf2, data=large_EUFA, method="ML", na.action=na.fail)
summary(M.gls2.1)
anova(M.gls2.1, test = "F")

range(large_EUFA$d15N[large_EUFA$MarineArea == "1_BC"])
range(large_EUFA$d15N[large_EUFA$MarineArea == "2_CA"])
range(large_EUFA$d15N[large_EUFA$MarineArea == "3_BB"])
range(large_EUFA$d15N[large_EUFA$MarineArea == "4_PFZ"])
range(large_EUFA$d15N[large_EUFA$MarineArea == "5_PW"])

games_howell_test(data = large_EUFA, d15N ~ MarineArea, detailed = FALSE)


#Carbon 

M.gls1.1<-gls(d13Cs~MarineArea,weights = vf2, data=large_EUFA, method="ML", na.action=na.fail)
summary(M.gls1.1)
anova(M.gls1.1, test = "F")


range(large_EUFA$d13Cs[large_EUFA$MarineArea == "1_BC"])
sum(large_EUFA$d13Cs & large_EUFA$MarineArea == "1_BC")
range(large_EUFA$d13Cs[large_EUFA$MarineArea == "2_CA"])
sum(large_EUFA$d13Cs & large_EUFA$MarineArea == "2_CA")
range(large_EUFA$d13Cs[large_EUFA$MarineArea == "3_BB"])
sum(large_EUFA$d13Cs & large_EUFA$MarineArea == "3_BB")
range(large_EUFA$d13Cs[large_EUFA$MarineArea == "4_PFZ"])
sum(large_EUFA$d13Cs & large_EUFA$MarineArea == "4_PFZ")
range(large_EUFA$d13Cs[large_EUFA$MarineArea == "5_PW"])
sum(large_EUFA$d13Cs & large_EUFA$MarineArea == "5_PW")


large_EUFA <- as.data.frame(large_EUFA)

games_howell_test(data = large_EUFA, d13Cs ~ MarineArea, detailed = FALSE)

#### Models evaluating explanatory variables #####
## There are rows with some NAs in chla and sst satellite,
#We create a new object that does not contain NAs and includes all the columns that interest us

large_EUFA_no_na <- na.omit(large_EUFA[, c("d15N","d13Cs", "SST_sat_2","Lat", "Long", "depth_mosaico", "Chla_sat_2", "MarineArea")])
View(large_EUFA_no_na)

#evaluate correlations

cor.test(large_EUFA_no_na$Lat,large_EUFA_no_na$Chla_sat_2)# no correlac
cor.test(large_EUFA_no_na$Long,large_EUFA_no_na$Chla_sat_2)# no correlac
cor.test(large_EUFA_no_na$Lat,large_EUFA_no_na$SST_sat_2)# correlac, r=0.89
cor.test(large_EUFA_no_na$Long,large_EUFA_no_na$SST_sat_2)# no correlac
cor.test(large_EUFA_no_na$SST_sat_2,large_EUFA_no_na$Chla_sat_2)# no correlac
cor.test(large_EUFA_no_na$Lat,large_EUFA_no_na$depth_mosaico)# no correlac
cor.test(large_EUFA_no_na$Long,large_EUFA_no_na$depth_mosaico)# no correlac
cor.test(large_EUFA_no_na$Chla_sat_2,large_EUFA_no_na$depth_mosaico)# no correlac
cor.test(large_EUFA_no_na$SST_sat_2,large_EUFA_no_na$depth_mosaico)# no correlac

#Nitrogen
M.gls1.2.2 <- gls(d15N ~ SST_sat_2 + Long + Chla_sat_2 + depth_mosaico, weights = vf2, data = large_EUFA_no_na, method = "ML",na.action = na.fail)
summary(M.gls1.2.2)

dredge(M.gls1.2.2)
coef(M.gls1.2.2)
confint(M.gls1.2.2)

tab_model(M.gls1.2.2)

#Carbon
M.gls1.1.2 <- gls(d13Cs ~ SST_sat_2 + Long + Chla_sat_2 + depth_mosaico, weights = vf2, data = large_EUFA_no_na, method = "ML",na.action = na.fail)
summary(M.gls1.1.2)

dredge(M.gls1.1.2)
coef(M.gls1.1.2)
confint(M.gls1.1.2)

tab_model(M.gls1.1.2)

#######  Only Amphipod ######

large_Amphi<- table %>%
  filter(SP_red %in% c("Amphipod")& Season=="autumn")

#Nitrogen
shapiro.test(large_Amphi$d15N)#normality
large_amphi_d15N<-lm(d15N~MarineArea,data=large_Amphi, na.action=na.fail)
plot(large_amphi_d15N)

#Carbon
shapiro.test(large_Amphi$d13Cs)#normality
large_amphi_d13Cs<-lm(d13Cs~MarineArea,data=large_Amphi, na.action=na.fail)
plot(large_amphi_d13Cs)

# unbalanced groups, we model heterocedacea
vf2 <- varIdent(form= ~ 1 | MarineArea)

#Nitrogen 
M.gls2.1.1<-gls(d15N~MarineArea,weights = vf2, data=large_Amphi, method="ML", na.action=na.fail)
summary(M.gls2.1.1)
dredge(M.gls2.1.1)
confint(M.gls2.1.1)
anova(M.gls2.1.1, test = "F")

range(large_Amphi$d15N[large_Amphi$MarineArea == "1_BC"])
range(large_Amphi$d15N[large_Amphi$MarineArea == "2_CA"])
range(large_Amphi$d15N[large_Amphi$MarineArea == "3_BB"])
range(large_Amphi$d15N[large_Amphi$MarineArea == "4_PFZ"])
range(large_Amphi$d15N[large_Amphi$MarineArea == "5_PW"])

games_howell_test(data = large_Amphi, d15N ~ MarineArea, detailed = FALSE)

#Carbon
M.gls1.1.1<-gls(d13Cs~MarineArea,weights = vf2, data=large_Amphi, method="ML", na.action=na.fail)
summary(M.gls1.1.1)
dredge(M.gls1.1.1)
confint(M.gls1.1.1)

anova(M.gls1.1.1, test = "F")

large_Amphi<- as.data.frame(large_Amphi)

games_howell_test(data = large_Amphi, d13Cs ~ MarineArea, detailed = FALSE)

range(large_Amphi$d13Cs[large_Amphi$MarineArea == "1_BC"])
sum(large_Amphi$d13Cs & large_Amphi$MarineArea == "1_BC")
range(large_Amphi$d13Cs[large_Amphi$MarineArea == "2_CA"])
sum(large_Amphi$d13Cs & large_Amphi$MarineArea == "2_CA")
range(large_Amphi$d13Cs[large_Amphi$MarineArea == "3_BB"])
sum(large_Amphi$d13Cs & large_Amphi$MarineArea == "3_BB")
range(large_Amphi$d13Cs[large_Amphi$MarineArea == "4_PFZ"])
sum(large_Amphi$d13Cs & large_Amphi$MarineArea == "4_PFZ")
range(large_Amphi$d13Cs[large_Amphi$MarineArea == "5_PW"])
sum(large_Amphi$d13Cs & large_Amphi$MarineArea == "5_PW")

#### Models evaluating explanatory variables #####
## There are rows with some NAs in chla and sst satellite,
#We create a new object that does not contain NAs and includes all the columns that interest us
large_Am_no_na <- na.omit(large_Amphi[, c("d15N","d13Cs", "SST_sat_2","Lat", "Long", "depth_mosaico", "Chla_sat_2", "MarineArea")])

#Correlations 
cor.test(large_Am_no_na$Lat,large_Am_no_na$Chla_sat_2)# no correlac
cor.test(large_Am_no_na$Long,large_Am_no_na$Chla_sat_2)# no correlac
cor.test(large_Am_no_na$Lat,large_Am_no_na$SST_sat_2)# correlac, r=0.93
cor.test(large_Am_no_na$Long,large_Am_no_na$SST_sat_2)# no correlac
cor.test(large_Am_no_na$SST_sat_2,large_Am_no_na$Chla_sat_2)# no correlac
cor.test(large_Am_no_na$Lat,large_Am_no_na$depth_mosaico)# no correlac
cor.test(large_Am_no_na$Long,large_Am_no_na$depth_mosaico)# no correlac
cor.test(large_Am_no_na$Chla_sat_2,large_Am_no_na$depth_mosaico)# no correlac
cor.test(large_Am_no_na$SST_sat_2,large_Am_no_na$depth_mosaico)# no correlac

#Nitrogen
M.gls1.3 <- gls(d15N ~ SST_sat_2 + Long + Chla_sat_2 + depth_mosaico, weights = vf2, data = large_Am_no_na, method = "ML",na.action = na.fail)
summary(M.gls1.3)
dredge(M.gls1.3)
coef(M.gls1.3)
confint(M.gls1.3)
tab_model(M.gls1.3)
#Carbon
M.gls1.1.1.2 <- gls(d13Cs ~ SST_sat_2 + Long + Chla_sat_2 + depth_mosaico, weights = vf2, data = large_Am_no_na, method = "ML",na.action = na.fail)
summary(M.gls1.1.1.2)
dredge(M.gls1.1.1.2)
coef(M.gls1.1.1.2)
confint(M.gls1.1.1.2)
tab_model(M.gls1.1.1.2)



#########################################   REGIONAL SCALE #########################################

## Only Phyto

regional_FP <- table %>%
  filter(Subgroup == "FP" & MarineArea %in% c("1_BC", "2_CA", "3_BB"))

#Nitrogen - Autumn

regional_FP_d15N_autumn <- regional_FP %>%
  select("d15N","Season","MarineArea")%>%
  filter(Season == "autumn")%>%
  na.omit()

#Normality
shapiro.test(regional_FP_d15N_autumn$d15N)#normality
model_regional_FP_d15N_autumn<-lm(d15N~MarineArea,data=regional_FP_d15N_autumn, na.action=na.fail)
plot(model_regional_FP_d15N_autumn)

#Heterocedasticity

vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls3.2.0<-gls(d15N~MarineArea,weights = vf2, data=regional_FP_d15N_autumn, na.action=na.fail)
summary(M.gls3.2.0)
dredge(M.gls3.2.0)
confint(M.gls3.2.0)
anova(M.gls3.2.0, test = "F")

games_howell_test(data = regional_FP_d15N_autumn, d15N ~ MarineArea, detailed = FALSE)

range(regional_FP_d15N_autumna$d15N[regional_FP_d15N_autumn$MarineArea == "1_BC"])
sum(regional_FP_d15N_autumn$d15N & regional_FP_d15N_autumn$MarineArea == "1_BC")
range(regional_FP_d15N_autumn$d15N[regional_FP_d15N_autumn$MarineArea == "2_CA"])
sum(regional_FP_d15N_autumn$d15N & regional_FP_d15N_autumn$MarineArea == "2_CA")
range(regional_FP_d15N_autumn$d15N[regional_FP_d15N_autumn$MarineArea == "3_BB"])
sum(regional_FP_d15N_autumn$d15N & regional_FP_d15N_autumn$MarineArea == "3_BB")

#Carbon - Autumn

regional_FP_d13C_autumn <- regional_FP %>%
  select("d13Cs","Season","MarineArea")%>%
  filter(Season == "autumn")
#Normality
shapiro.test(regional_FP_d13C_autumn$d13Cs)#normality
model_regional_FP_d13C_autumn<-lm(d13Cs~MarineArea,data=regional_FP_d13C_autumn)
plot(model_regional_FP_d13C_autumn)

#Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls3.1.0<-gls(d13Cs~MarineArea,weights = vf2, data=regional_FP_d13C_autumn, na.action=na.fail)
summary(M.gls3.1.0)
dredge(M.gls3.1.0)
anova(M.gls3.1.0, test = "F")

games_howell_test(data = regional_FP_d13C_autumn, d13Cs ~ MarineArea, detailed = FALSE)

range(regional_FP_d13C_autumn$d13Cs[regional_FP_d13C_autumn$MarineArea == "1_BC"])
sum(regional_FP_d13C_autumn$d13Cs & regional_FP_d13C_autumn$MarineArea == "1_BC")
range(regional_FP_d13C_autumn$d13Cs[regional_FP_d13C_autumn$MarineArea == "2_CA"])
sum(regional_FP_d13C_autumn$d13Cs & regional_FP_d13C_autumn$MarineArea == "2_CA")
range(regional_FP_d13C_autumn$d13Cs[regional_FP_d13C_autumn$MarineArea == "3_BB"])
sum(regional_FP_d13C_autumn$d13Cs & regional_FP_d13C_autumn$MarineArea == "3_BB")


#Nitrogen - Spring

regional_FP_d15N_spring <- regional_FP %>%
  select("d15N","Season","MarineArea")%>%
  filter(Season == "spring")%>%
  na.omit()
#Normality
shapiro.test(regional_FP_d15N_spring$d15N)#check qqplot
model_regional_FP_d15N_spring<-lm(d15N~MarineArea,data=regional_FP_d15N_spring)
plot(model_regional_FP_d15N_spring)

#Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls3.2.00<-gls(d15N~MarineArea,weights = vf2, data=regional_FP_d15N_spring, na.action=na.fail)
summary(M.gls3.2.00)
anova(M.gls3.2.00, test = "F")

range(regional_FP_d15N_spring$d15N[regional_FP_d15N_spring$MarineArea == "1_BC"])
range(regional_FP_d15N_spring$d15N[regional_FP_d15N_spring$MarineArea == "2_CA"])
range(regional_FP_d15N_spring$d15N[regional_FP_d15N_spring$MarineArea == "3_BB"])
sum(regional_FP_d15N_spring$d15N & regional_FP_d15N_spring$MarineArea == "1_BC")
sum(regional_FP_d15N_spring$d15N & regional_FP_d15N_spring$MarineArea == "2_CA")
sum(regional_FP_d15N_spring$d15N & regional_FP_d15N_spring$MarineArea == "3_BB")

games_howell_test(data = regional_FP_d15N_spring, d15N~ MarineArea, detailed = FALSE)

# Carbon - Spring

regional_FP_d13C_spring <- regional_FP %>%
  select("d13Cs","Season","MarineArea")%>%
  filter(Season == "spring")%>%
  na.omit()

#Normality
shapiro.test(regional_FP_d13C_spring$d13Cs)#check qqplot
model_regional_FP_d13C_spring<-lm(d13Cs~MarineArea,data=regional_FP_d13C_spring)
plot(model_regional_FP_d13C_spring)# qqplot accept normality


#Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls3.1.00<-gls(d13Cs~MarineArea,weights = vf2, data=regional_FP_d13C_spring, na.action=na.fail)
summary(M.gls3.1.00)
anova(M.gls3.1.00, test = "F")


games_howell_test(data = regional_FP_d13C_spring, d13Cs ~ MarineArea, detailed = FALSE)

range(regional_FP_d13C_spring$d13Cs[regional_FP_d13C_spring$MarineArea == "1_BC"])
sum(regional_FP_d13C_spring$d13Cs & regional_FP_d13C_spring$MarineArea == "1_BC")
range(regional_FP_d13C_spring$d13Cs[regional_FP_d13C_spring$MarineArea == "2_CA"])
sum(regional_FP_d13C_spring$d13Cs & regional_FP_d13C_spring$MarineArea == "2_CA")
range(regional_FP_d13C_spring$d13Cs[regional_FP_d13C_spring$MarineArea == "3_BB"])
sum(regional_FP_d13C_spring$d13Cs & regional_FP_d13C_spring$MarineArea == "3_BB")

#remove NAs from satellital data

regional_FP_no_na <- na.omit(regional_FP[, c("d15N","d13Cs", "Temp", "Long", "depth_mosaico", "Sal", "MarineArea", "Season")])

#Correlations

cor.test(regional_FP$Temp, regional_FP$Long) #no correlac r=-0.5
cor.test(regional_FP$Temp,regional_FP$depth_mosaico)# no correlac r=0.2
cor.test(regional_FP$Temp,regional_FP$Sal) #no correlac r=-0.41
cor.test(regional_FP$Long,regional_FP$depth_mosaico)#no correlac r=-0.4
cor.test(regional_FP$Long,regional_FP$Sal)# correlac r=0.8
cor.test(regional_FP$Sal,regional_FP$depth_mosaico)# no correlac r=-0.3

#Models for Nitrogen

M.gls3.2.0<-gls(d15N~MarineArea,weights = vf2, data=regional_FP_no_na, na.action=na.fail)
summary(M.gls3.2.0)
dredge(M.gls3.2.0)
confint(M.gls3.2.0)
tab_model(M.gls3.2.0)
anova(M.gls3.2.0, test = "F")

range(regional_FP_no_na$d15N[regional_FP_no_na$MarineArea == "1_BC"])
sum(regional_FP_no_na$d15N & regional_FP_no_na$MarineArea == "1_BC")
range(regional_FP_no_na$d15N[regional_FP_no_na$MarineArea == "2_CA"])
sum(regional_FP_no_na$d15N & regional_FP_no_na$MarineArea == "2_CA")
range(regional_FP_no_na$d15N[regional_FP_no_na$MarineArea == "3_BB"])
sum(regional_FP_no_na$d15N & regional_FP_no_na$MarineArea == "3_BB")


regional_FP_no_na <- as.data.frame(regional_FP_no_na)

games_howell_test(data = regional_FP_no_na, d15N ~ MarineArea, detailed = FALSE)

M.gls3.2<-gls(d15N~Temp+depth_mosaico+Sal+Season,weights = vf2, data=regional_FP_no_na, na.action=na.fail)
summary(M.gls3.2)
dredge(M.gls3.2)
tab_model(M.gls3.2)


#Models for Carbon

M.gls3.1<-gls(d13Cs~Temp+depth_mosaico+Sal+Season,weights = vf2, data=regional_FP_no_na, na.action=na.fail)
summary(M.gls3.1)
dredge(M.gls3.1)
tab_model(M.gls3.1)


## Only SPOM CB vs BB autumn

regional_SPOM <- table %>%
  filter(Subgroup == "SPOM" & MarineArea %in% c("1_BC", "2_CA", "3_BB"))

# Autumn, cb vs bb

regional_SPOM_autumn<- regional_SPOM %>%
  filter(Season=="autumn" & MarineArea %in% c("1_BC", "3_BB"))

# Nitrogen
# Normality
shapiro.test(regional_SPOM_autumn$d15N)#normality
model_regional_SPOM_d15N_autumn<-lm(d15N~MarineArea,data=regional_SPOM_autumn)
plot(model_regional_SPOM_d15N_autumn)# qqplot accept normality

#Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls4.2.00<-gls(d15N~MarineArea,weights = vf2, data=regional_SPOM_autumn, na.action=na.fail)
summary(M.gls4.2.00)
anova(M.gls4.2.00, test = "F")

range(regional_SPOM_autumn$d15N[regional_SPOM_autumn$MarineArea == "1_BC"])
sum(regional_SPOM_autumn$d15N & regional_SPOM_autumn$MarineArea == "1_BC")
range(regional_SPOM_autumn$d15N[regional_SPOM_autumn$MarineArea == "3_BB"])
sum(regional_SPOM_autumn$d15N & regional_SPOM_autumn$MarineArea == "3_BB")

games_howell_test(data = regional_SPOM_autumn, d15N ~ MarineArea, detailed = FALSE)

# Carbon
# Normality
shapiro.test(regional_SPOM_autumn$d13Cs)#check qqplot
model_regional_SPOM_d13C_autumn<-lm(d13Cs~MarineArea,data=regional_SPOM_autumn)
plot(model_regional_SPOM_d13C_autumn)# qqplot accept normality

#Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)
M.gls4.1.0<-gls(d13Cs~MarineArea,weights = vf2, data=regional_SPOM_autumn, na.action=na.fail)
summary(M.gls4.1.0)
anova(M.gls4.1.0, test = "F")

range(regional_SPOM_autumn$d13Cs[regional_SPOM_autumn$MarineArea == "1_BC"])
sum(regional_SPOM_autumn$d13Cs & regional_SPOM_autumn$MarineArea == "1_BC")
range(regional_SPOM_autumn$d13Cs[regional_SPOM_autumn$MarineArea == "3_BB"])
sum(regional_SPOM_autumn$d13Cs & regional_SPOM_autumn$MarineArea == "3_BB")

games_howell_test(data = regional_SPOM_autumn, d13Cs ~ MarineArea, detailed = FALSE)


# Correlations
cor.test(regional_SPOM_autumn$Temp, regional_SPOM_autumn$Long) # correlac r=-0.8
cor.test(regional_SPOM_autumn$Temp,regional_SPOM_autumn$depth_mosaico)# no correlac r=0.09
cor.test(regional_SPOM_autumn$Temp,regional_SPOM_autumn$Sal) # correlac r=-0.66
cor.test(regional_SPOM_autumn$Long,regional_SPOM_autumn$depth_mosaico)#no correlac r=-0.5
cor.test(regional_SPOM_autumn$Long,regional_SPOM_autumn$Sal)# correlac r=0.6

# Remove NAs in satellital data

regional_SPOM_no_na <- na.omit(regional_SPOM[, c("d15N","d13Cs", "Temp", "depth_mosaico", "Sal", "Season", "Long","MarineArea")])

# Models for Nitrogen
M.gls4.2<-gls(d15N~Temp+depth_mosaico,weights = vf2, data=regional_SPOM_no_na, na.action=na.fail)
summary(M.gls4.2)
dredge(M.gls4.2)
tab_model(M.gls4.2)

#Models for Carbon
M.gls4.1<-gls(d13Cs~Temp+depth_mosaico,weights = vf2, data=regional_SPOM.0.1, na.action=na.fail)
summary(M.gls4.1)
dredge(M.gls4.1)
tab_model(M.gls4.1)


## Only SPOM CB, AC BB - Spring

regional_SPOM_spring <- table %>%
  filter(Subgroup == "SPOM" & MarineArea %in% c("1_BC", "2_CA", "3_BB")& Season=="spring")

# Nitrogen
# Normality
shapiro.test(regional_SPOM_spring$d15N)#check qqplot
model_regional_SPOM_d15N_spring<-lm(d15N~MarineArea,data=regional_SPOM_spring)
plot(model_regional_SPOM_d15_spring)# qqplot accept normality

#Heterocedasticity

vf2 <- varIdent(form= ~ 1 | MarineArea)

#one value NA in d15N column.
regional_SPOM_spring_d15N <- regional_SPOM_spring %>%
  select("d15N","MarineArea")%>%
  na.omit()

M.gls4.3<-gls(d15N~MarineArea,weights = vf2, data=regional_SPOM_spring_d15N, na.action=na.fail)
summary(M.gls4.3)
anova(M.gls4.3, test = "F")

range(regional_SPOM_spring_d15N$d15N[regional_SPOM_spring_d15N$MarineArea == "1_BC"])
sum(regional_SPOM_spring_d15N$d15N & regional_SPOM_spring_d15N$MarineArea == "1_BC")
range(regional_SPOM_spring_d15N$d15N[regional_SPOM_spring_d15N$MarineArea == "2_CA"])
sum(regional_SPOM_spring_d15N$d15N & regional_SPOM_spring_d15N$MarineArea == "2_CA")
range(regional_SPOM_spring_d15N$d15N[regional_SPOM_spring_d15N$MarineArea == "3_BB"])
sum(regional_SPOM_spring_d15N$d15N & regional_SPOM_spring_d15N$MarineArea == "3_BB")

games_howell_test(data = regional_SPOM_spring_d15N, d15N ~ MarineArea, detailed = FALSE)

# Carbon
# Normality
shapiro.test(regional_SPOM_spring$d13Cs)#check qqplot
model_regional_SPOM_d13C_spring<-lm(d13Cs~MarineArea,data=regional_SPOM_spring)
plot(model_regional_SPOM_d13C_spring)# qqplot accept normality

#Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls4.3.1<-gls(d13Cs~MarineArea,weights = vf2, data=regional_SPOM_spring, na.action=na.fail)
summary(M.gls4.3.1)
anova(M.gls4.3.1, test = "F")

range(regional_SPOM_spring$d13Cs[regional_SPOM_spring$MarineArea == "1_BC"])
sum(regional_SPOM_spring$d13Cs & regional_SPOM_spring$MarineArea == "1_BC")
range(regional_SPOM_spring$d13Cs[regional_SPOM_spring$MarineArea == "2_CA"])
sum(regional_SPOM_spring$d13Cs & regional_SPOM_spring$MarineArea == "2_CA")
range(regional_SPOM_spring$d13Cs[regional_SPOM_spring$MarineArea == "3_BB"])
sum(regional_SPOM_spring$d13Cs & regional_SPOM_spring$MarineArea == "3_BB")

games_howell_test(data = regional_SPOM_spring, d13Cs ~ MarineArea, detailed = FALSE)

# Models

#correlations

cor.test(regional_SPOM_spring$Temp,regional_SPOM_spring$Long) # correlac r=-0.9
cor.test(regional_SPOM_spring$Temp,regional_SPOM_spring$depth_mosaico)# no correlac r=0.4
cor.test(regional_SPOM_spring$Temp,regional_SPOM_spring$Sal) # correlac r=-0.95
cor.test(regional_SPOM_spring$Long,regional_SPOM_spring$depth_mosaico)#no correlac r=-0.4
cor.test(regional_SPOM_spring$Long,regional_SPOM_spring$Sal)# correlac r=0.9


## Remove NAs satellital data
regional_SPOM_spring_no_na <- na.omit(regional_SPOM_spring[, c("d15N","d13Cs", "Temp", "depth_mosaico", "Sal", "Season", "Long","MarineArea")])

# For nitrogen
M.gls4.2.2<-gls(d15N~Temp+depth_mosaico,weights = vf2, data=regional_SPOM_spring_no_na, na.action=na.fail)
summary(M.gls4.2.2)
dredge(M.gls4.2.2)
tab_model(M.gls4.2.2)

# For Carbon
M.gls4.1.2<-gls(d13Cs~Temp+depth_mosaico,weights = vf2, data=regional_SPOM_spring_no_na, na.action=na.fail)
summary(M.gls4.1.2)
dredge(M.gls4.1.2)
tab_model(M.gls4.1.2)


## Only COPE - Autumn

regional_COPE <- table %>%
  filter(SP_red == "COPE" & MarineArea %in% c("1_BC", "2_CA", "3_BB"))

# Nitrogen
# Normality
shapiro.test(regional_COPE$d15N)# check qqplot
model_regional_COPE_d15N<-lm(d15N~MarineArea,data=regional_COPE)
plot(model_regional_COPE_d15N)#qqplot normal

# Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls5.2.0<-gls(d15N~MarineArea,weights = vf2, data=regional_COPE, na.action=na.fail)
summary(M.gls5.2.0)
anova(M.gls5.2.0, test = "F")


range(regional_COPE$d15N[regional_COPE$MarineArea == "1_BC"])
sum(regional_COPE$d15N & regional_COPE$MarineArea == "1_BC")
range(regional_COPE$d15N[regional_COPE$MarineArea == "2_CA"])
sum(regional_COPE$d15N & regional_COPE$MarineArea == "2_CA")
range(regional_COPE$d15N[regional_COPE$MarineArea == "3_BB"])
sum(regional_COPE$d15N & regional_COPE$MarineArea == "3_BB")

games_howell_test(data = regional_COPE, d15N ~ MarineArea, detailed = FALSE)

#Carbon
shapiro.test(regional_COPE$d13Cs)# normality p=0.8
model_regional_COPE_d13C<-lm(d13Cs~MarineArea,data=regional_COPE)
plot(model_regional_COPE_d13C)#qqplot normal.

# Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls5.1.0<-gls(d13Cs~MarineArea,weights = vf2, data=regional_COPE, na.action=na.fail)
summary(M.gls5.1.0)
anova(M.gls5.1.0, test = "F")

games_howell_test(data = regional_COPE, d13Cs ~ MarineArea, detailed = FALSE)

range(regional_COPE$d13Cs[regional_COPE$MarineArea == "1_BC"])
sum(regional_COPE$d13Cs & regional_COPE$MarineArea == "1_BC")
range(regional_COPE$d13Cs[regional_COPE$MarineArea == "2_CA"])
sum(regional_COPE$d13Cs & regional_COPE$MarineArea == "2_CA")
range(regional_COPE$d13Cs[regional_COPE$MarineArea == "3_BB"])
sum(regional_COPE$d13Cs & regional_COPE$MarineArea == "3_BB")


# Correlations

cor.test(regional_COPE$Temp, regional_COPE$Long) # correlac r=-0.9
cor.test(regional_COPE$Temp,regional_COPE$depth_mosaico)# no correlac r=0.4
cor.test(regional_COPE$Temp,regional_COPE$Sal) # correlac r=-0.6
cor.test(regional_COPE$Long,regional_COPE$depth_mosaico)#no correlac r=-0.4
cor.test(regional_COPE$Long,regional_COPE$Sal)# correlac r=0.8

# Remove NAs satellital data

regional_COPE_no_na <- na.omit(regional_COPE[, c("d15N","d13Cs", "Temp", "depth_mosaico", "Sal","Long","MarineArea")])

# Models
# For Nitrogen
M.gls5.2<-gls(d15N~Temp+depth_mosaico,weights = vf2, data=regional_COPE_no_na, na.action=na.fail)
summary(M.gls5.2)
dredge(M.gls5.2)
tab_model(M.gls5.2)
# For Carbon
M.gls5.1<-gls(d13Cs~Temp+depth_mosaico,weights = vf2, data=regional_COPE_no_na, na.action=na.fail)
summary(M.gls5.1)
dredge(M.gls5.1)
tab_model(M.gls5.1)

#########################################   LOCAL SCALE #########################################

#### ONLY FP, autumn vs sring

local_FP <- table %>%
  filter(Subgroup == "FP" & MarineArea %in% c("1_BC", "2_CA"))

# Nitrogen
# Normality
shapiro.test(local_FP$d15N)# check qqlot
model_local_FP_d15N<-lm(d15N~MarineArea,data=local_FP)
plot(model_local_FP_d15N) #qqplot ok

# Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls6.1<-gls(d15N~MarineArea+Season,weights = vf2, data=local_FP, method="ML", na.action=na.fail)
summary(M.gls6.1)
dredge(M.gls6.1)
tab_model(M.gls6.1)

# Carbon
shapiro.test(local_FP$d13Cs)# check qqplot
model_local_FP_d13C<-lm(d13Cs~MarineArea,data=local_FP)
plot(model_local_FP_d13C) #qqplot ok
# Heterocedasticity
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls6.2<-gls(d13Cs~MarineArea+Season,weights = vf2, data=local_FP_no_na, method="ML", na.action=na.fail)
summary(M.gls6.2)
dredge(M.gls6.2)
coef(M.gls6.2)
confint(M.gls6.2)
tab_model(M.gls6.2)


## Only FP, autumn

local_FP_autumn <- table %>%
  filter(Subgroup == "FP" & MarineArea %in% c("1_BC", "2_CA") & Season=="autumn")

# Nitrogen
# Normality
shapiro.test(local_FP_autumn$d15N)# check qqplot
model_local_FP_autumn_d15N<-lm(d15N~MarineArea,data=local_FP_autumn)
plot(model_local_FP_autumn_d15N) #qqplot ok

# Carbon
# Normality
shapiro.test(local_FP_autumn$d13Cs)#normal
model_local_FP_autumn_d13C<-lm(d13Cs~MarineArea,data=local_FP_autumn)
plot(model_local_FP_autumn_d13C) #qqplot ok


## Remove NAs
local_FP_autumn_no_na <- na.omit(local_FP_autumn[, c("d15N","d13Cs","MarineArea","NH4",
                                                     "NO3", "ChlaTot", "Temp", "Si","PO4")])

#Correlations

cor.test(local_FP_autumn_no_na$ChlaTot,local_FP_autumn_no_na$NO3)# no correlac
cor.test(local_FP_autumn_no_na$ChlaTot,local_FP_autumn_no_na$PO4)# no correlac
cor.test(local_FP_autumn_no_na$ChlaTot,local_FP_autumn_no_na$Si)# no correlac
cor.test(local_FP_autumn_no_na$ChlaTot,local_FP_autumn_no_na$NH4)# no correlac
cor.test(local_FP_autumn_no_na$ChlaTot,local_FP_autumn_no_na$Temp)# no correlac
cor.test(local_FP_autumn_no_na$NO3,local_FP_autumn_no_na$PO4)# no correlac
cor.test(local_FP_autumn_no_na$NO3,local_FP_autumn_no_na$Si)# correlac r=-0.6
cor.test(local_FP_autumn_no_na$NO3,local_FP_autumn_no_na$NH4)# correlac r=-0.9
cor.test(local_FP_autumn_no_na$NO3,local_FP_autumn_no_na$Temp)# correlac r=-0.8
cor.test(local_FP_autumn_no_na$PO4,local_FP_autumn_no_na$Si)# no correlac
cor.test(local_FP_autumn_no_na$PO4,local_FP_autumn_no_na$NH4)# no correlac
cor.test(local_FP_autumn_no_na$PO4,local_FP_autumn_no_na$Temp)# no correlac
cor.test(local_FP_autumn_no_na$Si,local_FP_autumn_no_na$NH4)# correlac r=-0.7
cor.test(local_FP_autumn_no_na$Si,local_FP_autumn_no_na$Temp)# no correlac
cor.test(local_FP_autumn_no_na$NH4,local_FP_autumn_no_na$Temp)# correlac r=0.9

#Models 
vf2 <- varIdent(form= ~ 1 | MarineArea)

# For carbon 

M.gls7.1<-gls(d13Cs~ChlaTot+PO4+NO3,weights = vf2, data=local_FP_autumn_no_na, method="ML", na.action=na.fail)
summary(M.gls7.1)
dredge(M.gls7.1)
coef(M.gls7.1)
confint(M.gls7.1)
tab_model(M.gls7.1)

# For nitrogen 

M.gls7.2<-gls(d15N~ChlaTot+PO4+NO3,weights = vf2, data=local_FP_autumn_no_na, method="ML", na.action=na.fail)
summary(M.gls7.2)
dredge(M.gls7.2)
coef(M.gls7.2)
confint(M.gls7.2)
tab_model(M.gls7.2)


## Only ZP, autumn

local_ZP_autumn <- table %>%
  filter(Subgroup == "ZP" & MarineArea %in% c("1_BC", "2_CA") & Season=="autumn")

# Nitrogen
shapiro.test(local_ZP_autumn$d15N)# normal
model_local_ZP_autumn_d15N<-lm(d15N~MarineArea,data=local_ZP_autumn)
plot(model_local_ZP_autumn_d15N) #qqplot ok

# Carbon
# Normality
shapiro.test(local_ZP_autumn$d13Cs)#normal
model_local_ZP_autumn_d13C<-lm(d13Cs~MarineArea,data=local_ZP_autumn)
plot(model_local_ZP_autumn_d13C) #qqplot ok

## Remove NAs
local_ZP_autumn_no_na <- na.omit(local_ZP_autumn[, c("d15N","d13Cs","MarineArea","NH4",
                                                     "NO3", "ChlaTot", "Temp", "Si","PO4")])

#Correlations

cor.test(local_ZP_autumn_no_na$ChlaTot,local_ZP_autumn_no_na$NO3)# correlac r=-0.6
cor.test(local_ZP_autumn_no_na$ChlaTot,local_ZP_autumn_no_na$PO4)# no correlac
cor.test(local_ZP_autumn_no_na$ChlaTot,local_ZP_autumn_no_na$Si)# no correlac
cor.test(local_ZP_autumn_no_na$ChlaTot,local_ZP_autumn_no_na$NH4)# no correlac
cor.test(local_ZP_autumn_no_na$ChlaTot,local_ZP_autumn_no_na$Temp)# correlac? r=0.6

cor.test(local_ZP_autumn_no_na$NO3,local_ZP_autumn_no_na$PO4)# no correlac
cor.test(local_ZP_autumn_no_na$NO3,local_ZP_autumn_no_na$Si)# no correlac 
cor.test(local_ZP_autumn_no_na$NO3,local_ZP_autumn_no_na$NH4)# correlac r=-0.8
cor.test(local_ZP_autumn_no_na$NO3,local_ZP_autumn_no_na$Temp)# correlac r=-0.7

cor.test(local_ZP_autumn_no_na$PO4,local_ZP_autumn_no_na$Si)# no correlac
cor.test(local_ZP_autumn_no_na$PO4,local_ZP_autumn_no_na$NH4)# no correlac
cor.test(local_ZP_autumn_no_na$PO4,local_ZP_autumn_no_na$Temp)# no correlac

cor.test(local_ZP_autumn_no_na$Si,local_ZP_autumn_no_na$NH4)# no correlac
cor.test(local_ZP_autumn_no_na$Si,local_ZP_autumn_no_na$Temp)# no correlac

cor.test(local_ZP_autumn_no_na$NH4,local_ZP_autumn_no_na$Temp)# correlac r=0.9

# Models
# For Nitrogen
vf2 <- varIdent(form= ~ 1 | MarineArea)

M.gls8.2<-gls(d15N~ChlaTot+NH4+Si+PO4,weights = vf2, data=local_ZP_autumn_no_na, method="ML", na.action=na.fail)
summary(M.gls8.2)
dredge(M.gls8.2)
coef(M.gls8.2)
confint(M.gls8.2)
tab_model(M.gls8.2)

# For Carbon

M.gls8.1<-gls(d13Cs~ChlaTot+NH4+Si+PO4,weights = vf2, data=local_ZP_autumn_no_na, method="ML", na.action=na.fail)
summary(M.gls8.1)
dredge(M.gls8.1)
coef(M.gls8.1)
confint(M.gls8.1)
tab_model(M.gls8.1)




