library(unmarked)
library(AICcmodavg)
library(dplyr)
library(ggplot2)
library(ggpubr)

## ----import and check data----
data <- read.csv("Data/occupancy_data.csv") 
str(data)
df <- na.omit(data) ##omit rows with NA values 
df$elionurus <- as.factor(df$elionurus)

### Calculate interval between reference time (9:00 am) and time of survey
df$ref_time
tref <- strptime(df$ref_time, format = "%H:%M:%S")
t1 <- strptime(df$time_j1, format = "%H:%M:%S")
t2 <- strptime(df$time_j2, format = "%H:%M:%S")
t3 <- strptime(df$time_j3, format = "%H:%M:%S")
int1 <- as.numeric(difftime(t1, tref, units = "mins")) 
int2 <- as.numeric(difftime(t2, tref, units = "mins"))
int3 <- as.numeric(difftime(t3, tref, units = "mins"))

#bind the numeric values of intervals with df
df <- cbind(df, int1, int2, int3)

#check correlation between survey covariates
atemp_cor <- stack(df, select = c("atemp_j1", "atemp_j2","atemp_j3"))
colnames(atemp_cor) <- c("atemp", "survey")
hum_cor <- stack(df, select = c("hum_j1", "hum_j2","hum_j3"))
colnames(hum_cor) <- c("hum", "survey")
wind_cor <- stack(df, select = c("wind_j1", "wind_j2","wind_j3"))
colnames(wind_cor) <- c("wind", "survey")
stemp_cor <- stack(df, select = c("stemp_j1", "stemp_j2","stemp_j3"))
colnames(stemp_cor) <- c("stemp", "survey")
int_cor <- stack(df, select = c("int1", "int2","int3"))
colnames(int_cor) <- c("int", "survey")
df_cor <- cbind(atemp_cor$atemp, hum_cor$hum, wind_cor$wind, stemp_cor$stemp, int_cor$int)
colnames(df_cor) <- c("atemp","hum","wind","stemp","int")
round(cor(df_cor), 2)

det_cor <- stack(df, select = c("j1","j2","j3"))
names(det_cor) <- c("det", "survey")

##check correlation between site covariates
df_sitecovs <- df %>% select(grain_size, t_veg, grass_suit, grass_non_suit, t_grass,
                             shrubs, trees, compaction, size)
cor_df_sitecovs <- round(cor(df_sitecovs), 2)
write.csv(cor_df_sitecovs, "cor_df_sitecovs.csv")

##preliminary analysis of survey covs with detection data-----
plot(atemp_cor$atemp, det_cor$det, col = "red", xlab = "Air temperature", ylab = "Detection (0/1)") 
plot(hum_cor$hum, det_cor$det, col = "red")
plot(stemp_cor$stemp, det_cor$det, col = "red")
plot(int_cor$int, det_cor$det, col = "red", xlab = "Time of survey (since 9 am)", ylab = "Detection (0/1)")
plot(wind_cor$wind, det_cor$det, col = "red")

## ----preparing data for occupancy analysis----
ynew <- as.matrix(df[c("j1", "j2", "j3")])
atemp <- as.matrix(df[c("atemp_j1", "atemp_j2", "atemp_j3")])
hum <- as.matrix(df[c("hum_j1", "hum_j2", "hum_j3")])
wind <- as.matrix(df[c("wind_j1", "wind_j2", "wind_j3")])
stemp <- as.matrix(df[c("stemp_j1", "stemp_j2", "stemp_j3")])
interval <- as.matrix(df[c("int1", "int2", "int3")])
elionurus <- as.factor(df[, "elionurus"])
t_veg.orig <- df[,"t_veg"]
dung.orig <- df[,"dung"]
grass_suit.orig <- df[,"grass_suit"]
grass_non_suit.orig <- df[,"grass_non_suit"]
t_grass.orig <- df[,"t_grass"]
t_herbs.orig <- df[,"t_herbs"]
shrubs.orig <- df[,"shrubs"]
trees.orig <- df[,"trees"]
comp.orig <- df[,"compaction"]
size.orig <- df[,"size"]
grain_size.orig <- df[,"grain_size"]
area <- df[,"area"]

##scaling the continuous covariates such that mean = 0 and sd = 1, i.e x[i]-mean(x)/sd(x)
t_veg.s <- scale(t_veg.orig, scale = T, center = T) 
grass_suit.s <- scale(grass_suit.orig, scale = T, center = T)
grass_non_suit.s <- scale(grass_non_suit.orig, scale = T, center = T)
t_grass.s <- scale(t_grass.orig, scale = T, center = T)
shrubs.s <- scale(shrubs.orig, scale = T, center = T)
trees.s <- scale(trees.orig, scale = T, center = T)
comp.s <- scale(comp.orig, scale = T, center = T)
size.s <- scale(size.orig, scale = T, center = T)
grain_size.s <- scale(grain_size.orig, scale = T, center = T)
atemp.s <- scale(atemp, scale = T, center = T)
hum.s <- scale(hum, scale = T, center = T)
wind.s <- scale(wind, scale = T, center = T)
stemp.s <- scale(stemp, scale = T, center = T)
interval.s <- scale(interval, scale = T, center = T)

modocc <- unmarkedFrameOccu(y = ynew, siteCovs = data.frame(t_veg.s = t_veg.s, grass_suit.s = grass_suit.s, grass_non_suit.s = grass_non_suit.s, t_grass.s = t_grass.s,
                                                            shrubs.s = shrubs.s, trees.s = trees.s, comp.s = comp.s, size.s = size.s, grain.s = grain_size.s, elionurus = elionurus, area = area), 
                            obsCovs = list(atemp.s = atemp.s, hum.s = hum.s, wind.s = wind.s, stemp.s = stemp.s, interval.s = interval.s))
summary(modocc)

## ----detection models and AIC-based model selection----
m1 <- occu(~1 ~1, modocc)
summary(m1)
m2 <- occu(~t_veg.s ~1, modocc) #no clear cut effect of total vegetation
summary(m2)
confint(m2, type = "det")
m3 <- occu(~trees.s ~1, modocc) #no clear cut effect of trees
summary(m3)
confint(m3, type = "det")
m4 <- occu(~atemp.s ~1, modocc) #clearly air temperature increases detection
summary(m4)
confint(m4, type = "det")
m5 <- occu(~wind.s ~1, modocc)
summary(m5)
confint(m5, type = "det")
m6 <- occu(~stemp.s ~1, modocc)
summary(m6)
confint(m6, type = "det")
m6.a <- occu(~stemp.s + wind.s ~1, modocc)
summary(m6.a)
confint(m6.a, type = "det")
m7 <- occu(~atemp.s + I(atemp.s^2) ~1, modocc)
summary(m7)
confint(m7, type = "det")
m8 <- occu(~stemp.s + I(stemp.s^2) ~1, modocc) #stemp doesnt have clear effects
summary(m8)
confint(m8, type = "det")
m8.a <- occu(~stemp.s + I(stemp.s^2) + wind.s ~1, modocc) #stemp doesnt have clear effects
summary(m8.a)
confint(m8.a, type = "det")
m9 <- occu(~atemp.s + wind.s ~1, modocc)
summary(m9)
confint(m9, type = "det")
m10 <- occu(~atemp.s + I(atemp.s^2) + wind.s ~1, modocc)
summary(m10)
confint(m10, type = "det")
m11 <- occu(~interval.s ~1, modocc)
summary(m11)
confint(m11, type = "det")
m12 <- occu(~interval.s + I(interval.s^2) ~1, modocc)
summary(m12)
confint(m12, type = "det")
m13 <- occu(~interval.s + wind.s ~1, modocc)
summary(m13)
confint(m13, type = "det")
m14 <- occu(~interval.s + + I(interval.s^2) + wind.s ~1, modocc)
summary(m14)
confint(m14, type = "det")

##AIC selection for detection models
detmodels <- fitList('psi(.)p(.)' = m1,
                   'psi(.)p(tveg)' = m2,
                   'psi(.)p(trees)' = m3,
                   'psi(.)p(atemp)'= m4,
                   'psi(.)p(wind)'= m5,
                   'psi(.)p(stemp)' = m6,
                   'psi(.)p(stemp+wind)' = m6.a,
                   'psi(.)p(atemp+atemp^2)' = m7,
                   'psi(.)p(stemp+stemp^2)' = m8,
                   'psi(.)p(stemp+stemp^2+wind)' = m8.a,
                   'psi(.)p(atemp+wind)' = m9,
                   'psi(.)p(atemp+atemp^2+wind)' = m10,
                   'psi(.)p(time)' = m11,
                   'psi(.)p(time+time^2)' = m12,
                   'psi(.)p(time+wind)' = m13,
                   'psi(.)p(time+time^2+wind)' = m14
)
detmodels.result <- modSel(detmodels)
detmodels.result
detection_aic <- detmodels.result@Full

write.csv(detection_aic, "results_tables/detection_aic_results.csv")

## ----occupancy models and AIC-based model selection----
m15 <- occu(~atemp.s + wind.s ~elionurus, modocc)
summary(m15)
confint(m15, type = "state")
predict(m15, type = "state") ## checking effect of presence of elionurus sp.
m16 <- occu(~atemp.s + wind.s ~grass_suit.s, modocc)
summary(m16)
confint(m16, type = "state")
m17 <- occu(~atemp.s + wind.s ~trees.s, modocc)
summary(m17)
confint(m17, type = "state")
m18 <- occu(~atemp.s + wind.s ~shrubs.s, modocc) 
summary(m18)
confint(m18, type = "state")
m19 <- occu(~atemp.s + wind.s ~grain.s, modocc)
summary(m19)
confint(m19, type = "state")
m20 <- occu(~atemp.s + wind.s ~grain.s + I(grain.s^2), modocc)
summary(m20)
confint(m20, type = "state")
m21 <- occu(~atemp.s + wind.s ~size.s, modocc)
summary(m21)
confint(m21, type = "state")
m22 <- occu(~atemp.s + wind.s ~comp.s, modocc)
summary(m22)
confint(m22, type = "state")
m23 <- occu(~atemp.s + wind.s ~elionurus + grass_suit.s, modocc)
summary(m23)
confint(m23, type = "state")
m24 <- occu(~atemp.s + wind.s ~elionurus + trees.s, modocc)
summary(m24)
confint(m24, type = "state")
m25 <- occu(~atemp.s + wind.s ~elionurus + grain.s, modocc)
summary(m25)
confint(m25, type = "state")
m26 <- occu(~atemp.s + wind.s ~elionurus + trees.s + grain.s, modocc)
summary(m26)
confint(m26, type = "state")
m27 <- occu(~atemp.s + wind.s ~elionurus + trees.s + grain.s + I(grain.s^2), modocc)
summary(m27)
confint(m27, type = "state")
m28 <- occu(~atemp.s + wind.s ~elionurus + trees.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m28)
confint(m28, type = "state")
m29 <- occu(~atemp.s + wind.s ~elionurus + grain.s + I(grain.s^2), modocc)
summary(m29)
confint(m29, type = "state")
m30 <- occu(~atemp.s + wind.s ~elionurus + grass_suit.s + grain.s + I(grain.s^2), modocc)
summary(m30)
confint(m30, type = "state")
m31 <- occu(~atemp.s + wind.s ~elionurus + comp.s, modocc)
summary(m31)
confint(m31, type = "state")
m32 <- occu(~atemp.s + wind.s ~elionurus + comp.s + grain.s + I(grain.s^2), modocc)
summary(m32)
confint(m32, type = "state")
m33 <- occu(~atemp.s + wind.s ~elionurus + shrubs.s, modocc)
summary(m33)
confint(m33, type = "state")
m34 <- occu(~atemp.s + wind.s ~comp.s + grain.s, modocc)
summary(m34)
confint(m34, type = "state")
m35 <- occu(~atemp.s + wind.s ~comp.s + grain.s + I(grain.s^2), modocc)
summary(m35)
confint(m35, type = "state")
m36 <- occu(~atemp.s + wind.s ~trees.s + grain.s, modocc)
summary(m36)
confint(m36, type = "state")
m37 <- occu(~atemp.s + wind.s ~trees.s + grain.s + I(grain.s^2), modocc)
summary(m37)
confint(m37, type = "state")
m38 <- occu(~atemp.s + wind.s ~trees.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m38)
confint(m38,type = "state")
m39 <- occu(~atemp.s + wind.s ~grass_suit.s + grain.s + I(grain.s^2), modocc)
summary(m39)
confint(m39,type = "state")
m40 <- occu(~atemp.s + wind.s ~trees.s + grass_suit.s, modocc)
summary(m40)
confint(m40,type = "state")
m41 <- occu(~atemp.s + wind.s ~trees.s + grass_suit.s + grain.s + I(grain.s^2), modocc)
summary(m41)
confint(m41,type = "state")

##global model 
mglobal <- occu(~atemp.s + I(atemp.s^2) + wind.s + t_veg.s + trees.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(mglobal)
confint(mglobal, type = "state")
confint(mglobal, type = "det")

##occupancy model list
occumodels <- fitList('psi(.)p(atemp+wind)' = m9, ##null model for occupancy
                      'psi(elionurus)p(atemp+wind)' = m15,
                      'psi(grass_suit)p(atemp+wind)' = m16,
                      'psi(trees)p(atemp+wind)' = m17,
                      'psi(shrub)p(atemp+wind)' = m18,
                      'psi(grain)p(atemp+wind)' = m19,
                      'psi(grain+grain^2)p(atemp+wind)' = m20,
                      'psi(area)p(atemp+wind)' = m21,
                      'psi(comp)p(atemp+wind)' = m22,
                      'psi(elionurus+grass_suit)p(atemp+wind)' = m23,
                      'psi(elionurus+trees)p(atemp+wind)' = m24,
                      'psi(elionurus+grain)p(atemp+wind)' = m25,
                      'psi(elionurus+trees+grain)p(atemp+wind)' = m26,
                      'psi(elionurus+trees+grain+grain^2)p(atemp+wind)' = m27,
                      'psi(elionurus+trees+area+grain+grain^2)p(atemp+wind)' = m28,
                      'psi(elionurus+grain+grain^2)p(atemp+wind)' = m29,
                      'psi(elionurus+grass_suit+grain+grain^2)p(atemp+wind)' = m30,
                      'psi(elionurus+comp)p(atemp+wind)' = m31,
                      'psi(elionurus+comp+grain+grain^2)p(atemp+wind)' = m32,
                      'psi(elionurus+shrub)p(atemp+wind)' = m33,
                      'psi(comp+grain)p(atemp+wind)' = m34,
                      'psi(comp+grain+grain^2)p(atemp+wind)' = m35,
                      'psi(trees+grain)p(atemp+wind)' = m36,
                      'psi(trees+grain+grain^2)p(atemp+wind)' = m37,
                      'psi(trees+area+grain+grain^2)p(atemp+wind)' = m38,
                      'psi(grass_suit+grain+grain^2)p(atemp+wind)' = m39,
                      'psi(trees+grass_suit)p(atemp+wind)' = m40,
                      'psi(trees+grass_suit+grain+grain^2)p(atemp+wind)' = m41
)

##AIC-based model selection for occupancy
modsel.occ <- modSel(occumodels)
occupancy_aic <- modsel.occ@Full

write.csv(occupancy_aic, "results_tables/occupancy_aic_results.csv")

##GOF Bootstrap for the global model (Mackenzie & Bailey, 2004)
gof.boot1 <- mb.gof.test(mglobal, nsim = 10000)
gof.boot1

## ----QAIC-based model selection with c-hat = 1.31 and occupancy fixed to full list of covariates for det----
m1.qaic <- occu(~1 ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m1.qaic)
m2.qaic <- occu(~t_veg.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc) #no clear cut effect of total vegetation
summary(m2.qaic)
confint(m2.qaic, type = "det")
m3.qaic <- occu(~trees.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc) #no clear cut effect of trees
summary(m3.qaic)
confint(m3.qaic, type = "det")
m4.qaic <- occu(~atemp.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc) #clearly air temperature increases detection
summary(m4.qaic)
confint(m4.qaic, type = "det")
m5.qaic <- occu(~wind.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m5.qaic)
confint(m5.qaic, type = "det")
m6.qaic <- occu(~stemp.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m6.qaic)
confint(m6.qaic, type = "det")
m6.a.qaic <- occu(~stemp.s + wind.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m6.a.qaic)
confint(m6.a.qaic, type = "det")
m7.qaic <- occu(~atemp.s + I(atemp.s^2) ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m7.qaic)
confint(m7.qaic, type = "det")
m8.qaic <- occu(~stemp.s + I(stemp.s^2) ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc) #stemp doesnt have clear effects
summary(m8.qaic)
confint(m8.qaic, type = "det")
m8.a.qaic <- occu(~stemp.s + I(stemp.s^2) + wind.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc) #stemp doesnt have clear effects
summary(m8.a.qaic)
confint(m8.a.qaic, type = "det")
m9.qaic <- occu(~atemp.s + wind.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m9.qaic)
confint(m9.qaic, type = "det")
m10.qaic <- occu(~atemp.s + I(atemp.s^2) + wind.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m10.qaic)
confint(m10.qaic, type = "det")
m11.qaic <- occu(~interval.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m11.qaic)
confint(m11.qaic, type = "det")
m12.qaic <- occu(~interval.s + I(interval.s^2) ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m12.qaic)
confint(m12.qaic, type = "det")
m13.qaic <- occu(~interval.s + wind.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m13.qaic)
confint(m13.qaic, type = "det")
m14.qaic <- occu(~interval.s + + I(interval.s^2) + wind.s ~elionurus + grass_suit.s + trees.s + shrubs.s + comp.s + size.s + grain.s + I(grain.s^2), modocc)
summary(m14.qaic)
confint(m14.qaic, type = "det")

detmodels.qaic <- list('psi(*)p(.)' = m1.qaic,
                     'psi(*)p(tveg)' = m2.qaic,
                     'psi(*)p(trees)' = m3.qaic,
                     'psi(*)p(atemp)'= m4.qaic,
                     'psi(*)p(wind)'= m5.qaic,
                     'psi(*)p(stemp)' = m6.qaic,
                     'psi(*)p(stemp+wind)' = m6.a.qaic,
                     'psi(*)p(atemp+atemp^2)' = m7.qaic,
                     'psi(*)p(stemp+stemp^2)' = m8.qaic,
                     'psi(*)p(stemp+stemp^2+wind)' = m8.a.qaic,
                     'psi(*)p(atemp+wind)' = m9.qaic,
                     'psi(*)p(atemp+atemp^2+wind)' = m10.qaic,
                     'psi(*)p(time)' = m11.qaic,
                     'psi(*)p(time+time^2)' = m12.qaic,
                     'psi(*)p(time+wind)' = m13.qaic,
                     'psi(*)p(time+time^2+wind)' = m14.qaic
)

qaictable.det <- data.frame(aictab(detmodels.qaic, second.ord = FALSE, c.hat = 1.31))
write.csv(qaictable.det, "results_tables/detection_qaic_results.csv")

occumodels.qaic <- list('psi(.)p(atemp+wind)' = m9, ##null model for occupancy
                      'psi(elionurus)p(atemp+wind)' = m15,
                      'psi(grass_suit)p(atemp+wind)' = m16,
                      'psi(trees)p(atemp+wind)' = m17,
                      'psi(shrub)p(atemp+wind)' = m18,
                      'psi(grain)p(atemp+wind)' = m19,
                      'psi(grain+grain^2)p(atemp+wind)' = m20,
                      'psi(area)p(atemp+wind)' = m21,
                      'psi(comp)p(atemp+wind)' = m22,
                      'psi(elionurus+grass_suit)p(atemp+wind)' = m23,
                      'psi(elionurus+trees)p(atemp+wind)' = m24,
                      'psi(elionurus+grain)p(atemp+wind)' = m25,
                      'psi(elionurus+trees+grain)p(atemp+wind)' = m26,
                      'psi(elionurus+trees+grain+grain^2)p(atemp+wind)' = m27,
                      'psi(elionurus+trees+area+grain+grain^2)p(atemp+wind)' = m28,
                      'psi(elionurus+grain+grain^2)p(atemp+wind)' = m29,
                      'psi(elionurus+grass_suit+grain+grain^2)p(atemp+wind)' = m30,
                      'psi(elionurus+comp)p(atemp+wind)' = m31,
                      'psi(elionurus+comp+grain+grain^2)p(atemp+wind)' = m32,
                      'psi(elionurus+shrub)p(atemp+wind)' = m33,
                      'psi(comp+grain)p(atemp+wind)' = m34,
                      'psi(comp+grain+grain^2)p(atemp+wind)' = m35,
                      'psi(trees+grain)p(atemp+wind)' = m36,
                      'psi(trees+grain+grain^2)p(atemp+wind)' = m37,
                      'psi(trees+area+grain+grain^2)p(atemp+wind)' = m38,
                      'psi(grass_suit+grain+grain^2)p(atemp+wind)' = m39,
                      'psi(trees+grass_suit)p(atemp+wind)' = m40,
                      'psi(trees+grass_suit+grain+grain^2)p(atemp+wind)' = m41,
                      'psi(*)p(*)' = mglobal
)

qaictable.occ <- data.frame(aictab(occumodels.qaic, second.ord = FALSE, c.hat = 1.31))
write.csv(qaictable.occ, "results_tables/occupancy_qaic_results.csv")

##coefficient estimates for the top-ranked occupancy model - m27, corrected for overdispersion c-hat = 1.31
summary_m27 <- data.frame(summaryOD(m27, c.hat = 1.31)$outMat)
write.csv(summary_m27, "results_tables/coef_top-ranked_model.csv", row.names = T)

##creating new data for predicting occupancy and detection based on top-ranked model, corrected for c-hat = 1.31
newgrain_size <- seq(min(df$grain_size), max(df$grain_size), length.out = 100) 
newtrees <- seq(min(df$trees), max(df$trees), length.out = 100)
newatemp <- seq(min(atemp_cor$atemp), max(atemp_cor$atemp), length.out = 100)
newwind <- seq(min(wind_cor$wind), max(wind_cor$wind), length.out = 100)

##scaling the newly created covariates
newgrain_size.s <- scale(newgrain_size, scale = T, center = T)
newtrees.s <- scale(newtrees, scale = T, center = T)
newatemp.s <- scale(newatemp, scale = T, center = T)
newwind.s <- scale(newwind, scale = T, center = T)

##predicted occupancy as a function of sand grain size on predicted occupancy, in presence and absence of Elionurus sp.  
#all other covariates set to mean value
df_grain_pred_woelio <- data.frame(grain = newgrain_size, grain.s = newgrain_size.s, 
                                   trees.s = 0, elionurus = 0,
                                   atemp.s = 0, wind.s = 0)
predicted_grain_woelio <- modavgPred(cand.set = list(m27 = m27), newdata = df_grain_pred_woelio, parm.type = "psi", c.hat = 1.31)
predicted_grain_woelio <- cbind(df_grain_pred_woelio, predicted_grain_woelio$mod.avg.pred, predicted_grain_woelio$uncond.se, predicted_grain_woelio$lower.CL, predicted_grain_woelio$upper.CL) 

colnames(predicted_grain_woelio) <- c("grain", "grain.s", "trees.s", 
                                      "elionurus", "atemp.s", "wind.s",
                                      "mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
df_grain_pred_elio <- data.frame(grain = newgrain_size, grain.s = newgrain_size.s, 
                                 trees.s = 0, elionurus = 1,
                                 atemp.s = 0, wind.s = 0)
predicted_grain_elio <- modavgPred(cand.set = list(m27 = m27), newdata = df_grain_pred_elio, parm.type = "psi", c.hat = 1.31)
predicted_grain_elio <- cbind(df_grain_pred_elio, predicted_grain_elio$mod.avg.pred, predicted_grain_elio$uncond.se, predicted_grain_elio$lower.CL, predicted_grain_elio$upper.CL) 
colnames(predicted_grain_elio) <- c("grain", "grain.s", "trees.s", 
                                    "elionurus", "atemp.s", "wind.s",
                                    "mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")

predicted_grain_df <- rbind(predicted_grain_woelio, predicted_grain_elio)

##predicted occupancy as a function of trees, in presence and absence of Elionurus sp.  
#all other covariates set to mean value
df_trees_pred_woelio <- data.frame(grain.s = 0, 
                                   trees = newtrees,
                                   trees.s = newtrees.s, elionurus = 0,
                                   atemp.s = 0, wind.s = 0)
predicted_trees_woelio <- modavgPred(cand.set = list(m27 = m27), newdata = df_trees_pred_woelio, parm.type = "psi", c.hat = 1.31)
predicted_trees_woelio <- cbind(df_trees_pred_woelio, predicted_trees_woelio$mod.avg.pred, predicted_trees_woelio$uncond.se, predicted_trees_woelio$lower.CL, predicted_trees_woelio$upper.CL) 

colnames(predicted_trees_woelio) <- c("grain", "trees", "trees.s", 
                                      "elionurus", "atemp.s", "wind.s",
                                      "mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
df_trees_pred_elio <- data.frame(grain.s = 0, 
                                 trees = newtrees,
                                 trees.s = newtrees.s, elionurus = 1,
                                 atemp.s = 0, wind.s = 0)
predicted_trees_elio <- modavgPred(cand.set = list(m27 = m27), newdata = df_trees_pred_elio, parm.type = "psi", c.hat = 1.31)
predicted_trees_elio <- cbind(df_trees_pred_elio, predicted_trees_elio$mod.avg.pred, predicted_trees_elio$uncond.se, predicted_trees_elio$lower.CL, predicted_trees_elio$upper.CL) 

colnames(predicted_trees_elio) <- c("grain", "trees", "trees.s", 
                                    "elionurus", "atemp.s", "wind.s",
                                    "mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")

predicted_trees_df <- rbind(predicted_trees_woelio, predicted_trees_elio)

##predicted detection as a function of air temperature 
df_atemp_pred <- data.frame(grain.s = 0, 
                            trees.s = 0, elionurus = 0,
                            atemp = newatemp, atemp.s = newatemp.s, wind.s = 0)
predicted_atemp <- modavgPred(cand.set = list(m27 = m27), newdata = df_atemp_pred, parm.type = "detect", c.hat = 1.31)
predicted_atemp <- cbind(df_atemp_pred, predicted_atemp$mod.avg.pred, predicted_atemp$uncond.se, predicted_atemp$lower.CL, predicted_atemp$upper.CL) 

colnames(predicted_atemp) <- c("grain", "trees", "elionurus", "atemp", "atemp.s", "wind.s",
                                "mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")

##predicted detection as a function of wind 
df_wind_pred <- data.frame(grain.s = 0, 
                           trees.s = 0, elionurus = 0,
                           atemp.s = 0, wind = newwind, wind.s = newwind.s)
predicted_wind <- modavgPred(cand.set = list(m27 = m27), newdata = df_wind_pred, parm.type = "detect", c.hat = 1.31)
predicted_wind <- cbind(df_wind_pred, predicted_wind$mod.avg.pred, predicted_wind$uncond.se, predicted_wind$lower.CL, predicted_wind$upper.CL) 

colnames(predicted_wind) <- c("grain", "trees", "elionurus", "atemp", "wind", "wind.s",
                               "mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")

## predicted occupancy per study area
predicted_occu_area <- modavgPred(cand.set = list(m27 = m27), c.hat = 1.31,
                                  parm.type = "psi", 
                                  newdata = modocc@siteCovs)[c("mod.avg.pred","uncond.se","lower.CL","upper.CL")]
predicted_occu_area_df <- data.frame(predicted_occu = predicted_occu_area$mod.avg.pred,
                                     se = predicted_occu_area$uncond.se,
                                     lower.CL = predicted_occu_area$lower.CL,
                                     upper.CL = predicted_occu_area$upper.CL,
                                     modocc@siteCovs)
psych::describe(predicted_occu_area_df$predicted_occu)
psych::describeBy(predicted_occu_area_df$predicted_occu, group = predicted_occu_area_df$area)

## ----graphs for results----
## predicted occupancy and detection against covariates
my_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 18),
  legend.title = element_blank(), 
  legend.text = element_text(size = 18),
  legend.position = "bottom")

fig_pred_grain <- ggplot(predicted_grain_df, aes(x = grain, y = mod.avg.pred,
                                                 group = factor(elionurus),
                                                 color = factor(elionurus),
                                                 fill = factor(elionurus))) + 
                   geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.2, colour = NA) +
                   geom_line(size = 0.6) +
                   scale_x_continuous(breaks = seq(300, 800, 100)) +
                   scale_y_continuous(breaks = seq(0, 1, 0.2)) +
                   scale_fill_manual(values = c("1" = "blue","0" = "red"),
                                     guide = "none") +
                   scale_colour_manual(values = c("1" = "blue","0" = "red"),
                                       labels = c("Elionurus sp. absent", "Elionurus sp. present")) +
                   xlab("Sand grain size (μm)") + 
                   ylab("Occupancy probability")+
                   my_theme

fig_pred_trees <- ggplot(predicted_trees_df, aes(x = trees, y = mod.avg.pred,
                                                 group = factor(elionurus),
                                                 color = factor(elionurus),
                                                 fill = factor(elionurus))) + 
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.2, colour = NA) +
  geom_line(size = 0.6) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c("1" = "blue","0" = "red"),
                    guide = "none") +
  scale_colour_manual(values = c("1" = "blue","0" = "red"),
                      labels = c("Elionurus sp. absent", "Elionurus sp. present")) +
  xlab("Tree cover (%)") + 
  ylab("Occupancy probability")+
  my_theme

pred_occ_graph <- ggpubr::ggarrange(fig_pred_grain, fig_pred_trees, common.legend = T, legend = "bottom")

ggsave(path = "figs/main",
       filename = "predicted_occu_graph.png",
       plot = pred_occ_graph, 
       width = 40, height = 15, units = "cm",dpi = 600)

fig_pred_atemp <- ggplot(predicted_atemp, aes(x = atemp, y = mod.avg.pred)) + 
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.2, colour = NA) +
  geom_line(size = 0.6) +
  scale_x_continuous(breaks = seq(20, 45, 5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Air temperature (⁰C) ") + 
  ylab("Detection probability")+
  my_theme

fig_pred_wind <- ggplot(predicted_wind, aes(x = wind, y = mod.avg.pred)) + 
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.2, colour = NA) +
  geom_line(size = 0.6) +
  scale_x_continuous(breaks = seq(0, 8, 2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("Wind (m/s) ") + 
  ylab("Detection probability")+
  my_theme

pred_det_graph <- ggarrange(fig_pred_atemp, fig_pred_wind)

ggsave(path = "figs/main",
       filename = "predicted_det_graph.png",
       plot = pred_det_graph, 
       width = 40, height = 15, units = "cm",dpi = 600)

fig_pred_occu_area <- ggplot(data = predicted_occu_area_df, aes(x = area, y = predicted_occu))+
  geom_boxplot() +
  stat_summary(fun = mean, size = 1.5, shape = 8, color = "red", geom = "point") +
  labs(x = "Area",y = "Occupancy probability")+
  my_theme

ggsave(path = "figs/main",
       filename = "predicted_occ_area_graph.png",
       plot = fig_pred_occu_area, 
       width = 20, height = 15, units = "cm",dpi = 600)

## ----supplementary: summary graphs of site covariates----
area_elionurus <- ggplot(df, aes(x = area, fill = elionurus)) + 
  geom_bar(stat = "count",  alpha = 0.6) +
  scale_fill_manual(values = c("1" = "blue",
                               "0" = "red"),  
                    labels = c("Elionurus sp. absent", "Elionurus sp. present")) +
  scale_y_continuous(breaks = seq(0, 14, by = 2)) +
  labs(x = "Area", y = "Number of sites") +
  my_theme
ggsave(path = "figs/supp",
       filename = "area_elionurus.png",
       plot = area_elionurus, 
       width = 20, height = 13, units = "cm",dpi = 600)

area_trees <- ggplot(df, aes(x = area, y = trees))+ 
  geom_boxplot()+
  xlab("Area") + 
  ylab("Trees (%)") +
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  my_theme
ggsave(path = "figs/supp", 
       filename = "area_trees.png",
       plot = area_trees, width = 20, height = 13, units = "cm",dpi = 600)

area_sandgrain <- ggplot(df, aes(x = area, y = grain_size))+ 
  geom_boxplot()+
  xlab("Area") + 
  ylab("Sand grain size (μm)") +
  scale_y_continuous(breaks = seq(200, 800, by = 100)) +
  my_theme
ggsave(path = "figs/supp", 
       filename = "area_sandgrain.png",
       plot = area_sandgrain, width = 20, height = 13, units = "cm",dpi = 600)

area_bunchgrass <- ggplot(df, aes(x = area, y = grass_suit))+ 
  geom_boxplot()+
  xlab("Area") + 
  ylab("Bunch grass (%)") +
  my_theme
ggsave(path = "figs/supp", 
       filename = "area_bunchgrass.png",
       plot = area_bunchgrass, width = 20, height = 13, units = "cm",dpi = 600)

area_shrubs <- ggplot(df, aes(x = area, y = shrubs))+ 
  geom_boxplot()+
  xlab("Area") + 
  ylab("Shrubs (%)") +
  my_theme
ggsave(path = "figs/supp", 
       filename = "area_shrubs.png",
       plot = area_shrubs, width = 20, height = 13, units = "cm",dpi = 600)

area_compaction <- ggplot(df, aes(x = area, y = compaction))+ 
  geom_boxplot()+
  xlab("Area") + 
  ylab("Soil compaction (tons/ft2)") +
  my_theme
ggsave(path = "figs/supp", 
       filename = "area_compaction.png",
       plot = area_compaction, width = 20, height = 13, units = "cm",dpi = 600)

area_size <- ggplot(df, aes(x = area, y = size))+ 
  geom_boxplot()+
  xlab("Area") + 
  ylab("Area of site (m2)") +
  my_theme
ggsave(path = "figs/supp", 
       filename = "area_size.png",
       plot = area_size, width = 20, height = 13, units = "cm",dpi = 600)

summary_graphs_area <- ggarrange(area_elionurus, area_trees,
                                 area_sandgrain, area_bunchgrass,
                                 area_shrubs, area_compaction,
                                 area_size, nrow = 4, ncol = 2)
