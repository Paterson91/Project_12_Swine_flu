#Myeloid analysis 
#check Chanice's data - number of fields per pig and pig numbers 

library(readr)
library(dplyr)
library(plotrix) #good for calculating standard error for plotting
library(ggplot2)
library(lme4) #make sure data is a dataframe 
library(emmeans) #replacement to lsmeans(); note to use plim argument within pwpp() the latest version needs to be downloaded 
#from github, until it's updated on the cran 
library(cowplot)
MYELOID_MASTERSHEET_AT <- read_csv("MYELOID_MASTERSHEET_AT.csv")
df <- MYELOID_MASTERSHEET_AT

#check pig id's
table(df$Pig) #64 pigs, not all pigs have ten fields 
#match pig id's to master metadata - check they join and there aren't typo's!
metadata_swiv_timecourse <- read_csv("metadata-swiv-timecourse.csv")
metadata <- metadata_swiv_timecourse
names(metadata)[names(metadata) == "pig"] <- "Pig"

#use inner join to return all matching rows in metadata
id_check <- semi_join(metadata, df, "Pig") #correct - pig id's match up 
table(df$Replicate)

#plot figures for single molecules in lamina propria: use t-cell analysis 'nasal_septum_figures.R' as template script
#Log.Total.LP.MHCII
#Log.Total.LP.CD172
#Log.Total.LP.CD16
#Log.Total.LP.CD11R1

#calculate mean of log MHCII for each pig across replicate fields on each day
#calculate standard error of the mean and use this for plotting confidence intervals, use package plotrix to calculate standard error 
#MHCII
geo_mean_mhcii <- df %>%
  group_by(Pig, Day, Replicate, Infected) %>%
  summarise(geomean_mhcii = mean(Log.Total.LP.MHCII), se_geomean_cmhcii = std.error(Log.Total.LP.MHCII))

table(geo_mean_mhcii$Replicate)
#write.csv(geo_mean_mhcii, "geometric_mean_LP_mhcii.csv")

#CD172
geo_mean_cd172 <- df %>%
  group_by(Pig, Day, Replicate, Infected) %>%
  summarise(geomean_cd172 = mean(Log.Total.LP.CD172), se_geomean_cd172 = std.error(Log.Total.LP.CD172))
#write.csv(geo_mean_cd172, "geometric_mean_LP_cd172.csv")

#CD16
geo_mean_cd16 <- df %>%
  group_by(Pig, Day, Replicate, Infected) %>%
  summarise(geomean_cd16 = mean(Log.Total.LP.CD16), se_geomean_cd16 = std.error(Log.Total.LP.CD16))
#write.csv(geo_mean_cd16, "geometric_mean_LP_cd16.csv")

#CD11R1
geo_mean_cd11r1 <- df %>%
  group_by(Pig, Day, Replicate, Infected) %>%
  summarise(geomean_cd11r1 = mean(Log.Total.LP.CD11R1), se_geomean_cd11r1 = std.error(Log.Total.LP.CD11R1))
#write.csv(geo_mean_cd11r1, "geometric_mean_LP_cd11r1.csv")


#-------------------------------------------------------------------------------------------

#FIGURES
#multipanel figure 
nasal_septum_multipanel_allreps <- read_csv("nasal_septum_multipanel_allreps_myeloid.csv")
multi_df2 <- nasal_septum_multipanel_allreps
multi_df2$molecule <- as.factor(multi_df2$molecule)
multi_df2$Infected <- as.factor(multi_df2$Infected)
multi_df2$Replicate <- as.factor(multi_df2$Replicate)
class(multi_df2$Day)
class(df)
df <- as.data.frame(df)

#re-order 'molecule' so that order of panels in multipanel figure reflects: viral load, cd8, mhcI, cd4, mhcII
multi_df2$molecule = factor(multi_df2$molecule, levels = c("Viral load", "MHCII", "CD172", "CD16", "CD11R1"))
levels(multi_df2$molecule)

#when using facet wrap do not use $ to call to object!!

myColors <- c("blue", "red", "darkgray")

y_title <- expression(paste("log" ["10"], " total CD8 alpha"))
x_title <- expression(paste("log" ["10"], " REU/ml"))
y_title_cd4 <- expression(paste("log" ["10"], " total CD4"))
y_title_mhci <- expression(paste("log" ["10"], " total MHCI"))
y_title_mhcii <- expression(paste("log" ["10"], " total MHCII"))
y_title_multi <- expression(paste("log" ["10"]))
y_title_mhciicd16 <- expression(paste("log" ["10"], " total MHCIICD16"))
y_title_mhciicd172 <- expression(paste("log" ["10"], " total MHCIICD172"))

#multipanel figure showing standard error of mean across fields for each pig, different shape for each replicate
tiff("multi_panel_nasalseptum_allreps_myeloid.tiff", width = 7, height = 9, units = 'in', res = 300, compression = 'lzw')
multi_panelv2 <- ggplot(multi_df2, aes(x = Day, y = geo_mean, group = Infected, colour = Infected)) + 
  labs(x = "Day", y = y_title_multi) + geom_point(aes(shape = Replicate)) + scale_color_manual(labels = c("Control", "Infected", "In-contact"), values = myColors) +
  scale_x_continuous(breaks=c(-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13), labels=c("", "-1", "0", "1", "2","3","4","5", "6", "7","8","9", "10", "11", "12", "13")) +
  theme(axis.title =element_text(family = "Arial", color = "black", size =14)) + theme(strip.text.x = element_text(size = 12)) +   labs(title = "Lamina Propria\n", colour = "Infection Status") +
  geom_errorbar(aes(ymin=geo_mean-se_geomean, ymax=geo_mean+se_geomean, width = 0.3))
multi_panelfacetv2 <- multi_panelv2 + facet_grid(molecule~., scales = "free")
multi_panelfacetv2
dev.off()

#multipanel figure showing mean of pigs at each day 'grand mean'

#calculate mean of log mhcii/cd172/cd16/cd11r1 at each day across the four biological replicate groups
grand_geo_mean_mhcii <- df %>%
  group_by(Day, Infected) %>%
  summarise(geomean_mhcii = mean(Log.Total.LP.MHCII), se_geomean_mhcii = std.error(Log.Total.LP.MHCII))
#write.csv(grand_geo_mean_mhcii, "grand_geometric_mean_LP_mhcii.csv")

grand_geo_mean_cd172 <- df %>%
  group_by(Day, Infected) %>%
  summarise(geomean_cd172 = mean(Log.Total.LP.CD172), se_geomean_cd172 = std.error(Log.Total.LP.CD172))
#write.csv(grand_geo_mean_cd172, "grand_geometric_mean_LP_cd172.csv")

grand_geo_mean_cd16 <- df %>%
  group_by(Day, Infected) %>%
  summarise(geomean_cd16 = mean(Log.Total.LP.CD16), se_geomean_cd16 = std.error(Log.Total.LP.CD16))
#write.csv(grand_geo_mean_cd16, "grand_geometric_mean_LP_cd16.csv")

grand_geo_mean_cd11r1 <- df %>%
  group_by(Day, Infected) %>%
  summarise(geomean_cd11r1 = mean(Log.Total.LP.CD11R1), se_geomean_cd11r1 = std.error(Log.Total.LP.CD11R1))
#write.csv(grand_geo_mean_cd11r1, "grand_geometric_mean_LP_cd11r1.csv")

#create csv in excel of grand_geometric_mean for plotting and and save as multipanel_fig__nasalseptum_myeloid.csv

multipanel_fig_nasalseptum <- read_csv("multipanel_fig__nasalseptum_myeloid.csv")
multi_df <- multipanel_fig_nasalseptum
multi_df$molecule <- as.factor(multi_df$molecule)
multi_df$Infected <- as.factor(multi_df$Infected)
class(multi_df$Day)

#re-order 'molecule' so that order of panels in multipanel figure reflects: viral load, cd8, mhcI, cd4, mhcII
multi_df$molecule = factor(multi_df$molecule, levels = c("Viral load", "MHCII", "CD172", "CD16", "CD11R1"))
levels(multi_df$molecule)

tiff("multi_panel_nasalseptum_myeloid_grandmean.tiff", width = 7, height = 9, units = 'in', res = 300, compression = 'lzw')
multi_panel <- ggplot(multi_df, aes(x = Day, y = geo_mean, group = Infected, colour = Infected)) + 
  labs(x = "Day", y = y_title_multi) + geom_point() + scale_color_manual(labels = c("Control", "Infected", "In-contact"), values = myColors) +
  scale_x_continuous(breaks=c(-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13), labels=c("", "-1", "0", "1", "2","3","4","5", "6", "7","8","9", "10", "11", "12", "13")) +
  theme(axis.title =element_text(family = "Arial", color = "black", size =14)) + theme(strip.text.x = element_text(size = 12)) +   labs(title = "Nasal Septum, Lamina Propria\n", colour = "Infection Status") +
  geom_errorbar(aes(ymin=geo_mean-se_geomean, ymax=geo_mean+se_geomean, width = 0.3))
multi_panel
multi_panelfacet <- multi_panel + facet_grid(molecule~., scales = "free")
multi_panelfacet
dev.off()

#looks like there is an association with viral load and CD11R1 - but this trend is less apparent when looking on a per-replicate-basis due to variation
#let's see what the stat's say 

#-------------------------------------------------------------------------------------------

#STATISTICS 
#Mixed effect models - lsmeans - plot w/ se 
#phase of infection = fixed effect, pig as random effect, mean estimates adjusted for by replicate. 
#I think that the data supplied to lme4 lme() has to be a dataframe, otherwise the lsmeans() won't work on S4 objects with null pointers (whatever that means?)
df$Pig <- as.factor(df$Pig)
df$field <- as.factor(df$field)
class(df$Log.Total.LP.MHCII)
df$Infected <- as.factor(df$Infected)
df$Phase2 <- as.factor(df$Phase2)
df$Replicate <- as.factor(df$Replicate)
class(df$Pig)
class(df$field)
class(df$Infected)
class(df$Phase2)

#First, check random effects

#include response with no fixed effect and random effects only
#look for change in residual deviance and degrees of freedom between models compared to chi-sqaured distribution 

#nested random effect of field within pig 
mhcii_both <- lmer(Log.Total.LP.MHCII ~ 1 + (1|Pig/field), data = df, REML = FALSE)
# Error: number of levels of each grouping factor must be < number of observations

#random effect for field
mhcii_field <- lmer(Log.Total.LP.MHCII ~ 1 + (1|field), data = df, REML = FALSE)
# singular fit
summary(mhcii_field)

#complex random effects so each subject (pig) has a different slope and a different intercept 
mhcii_pig_si <- lmer(Log.Total.LP.MHCII ~ 1 + (Log.Total.LP.MHCII|Pig), data = df, REML = FALSE)
summary(mhcii_pig_si)
#singular fit

#simple random effects so each subject (pig) has the same slope but a different intercept 
mhcii_pig <- lmer(Log.Total.LP.MHCII ~ 1 + (1|Pig), data = df, REML = FALSE)
summary(mhcii_pig)

#evidence to keep random effect of pig with varying intercept and common slope... satisfied with random effects in the model
#evidence to keep because no other random effects fitted, if they did, then a LRTS could be calculated for competing models 

#is MHCII associated with phase of infection? 

library(nlme)
mhcii_null2 <- lme(Log.Total.LP.MHCII ~ 1, random = ~1|Pig, data = df, method = "ML")
summary(mhcii_null2) #AIC = 43.014

mhcii_rep <- lme(Log.Total.LP.MHCII ~ Replicate, random = ~1|Pig, data = df, method = "ML")
summary(mhcii_rep)#AIC = 19.5876

mhcii_phase2 <- lme(Log.Total.LP.MHCII ~ Phase2 + Replicate, random = ~ 1|Pig, data = df, method = "ML")
summary(mhcii_phase2) #AIC = 16.883

anova(mhcii_null2,mhcii_phase2) #p = <.0001
anova(mhcii_rep,mhcii_phase2) # p = 0.0263
#? when comparing two nested models for the overall p-value for the effect of Phase, this can be done using the likelihood ratio test,
# compare model with replicate and phase to model with replicate. 

#obtain least square means and pairwise contrasts for Phase, averaged over replicate, adjust for multiple comparisons (default Tukey)

#plot ls-mean mhcii
mhcii_ytitle <- expression(paste("log" ["10"], " MHCII"))

#estimated marginal means - for vignette see: https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html
mhcii_em2 <- emmeans(mhcii_phase2, "Phase2")
mhcii_lsm2_lsmeans <- as.data.frame(mhcii_em2)
mhcii_lsm2_contrasts <- as.data.frame(pairs(mhcii_em2))
#pairwise comparisons averaged over the levels of replicate, tukey adjustment
em_plot_mhcii <- plot(mhcii_em2, comparisons = TRUE, ylab= "Phase", xlab = mhcii_ytitle)
plot_emm_mhcii <- pwpp(mhcii_em2, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1)) #plot adj p-values and marginal means
plot_emm_mhcii

#save results
mhcii_em2_contrasts <- as.data.frame(pairs(mhcii_em2))
mhcii_em2_emeans <- as.data.frame(emmeans(mhcii_phase2, "Phase2"))

plot_mhcii_mean_bar <- ggplot(mhcii_lsm2_lsmeans, aes(x=Phase2, y = emmean, fill = Phase2)) + 
  geom_bar(stat = "identity", position = "identity") + geom_errorbar(aes(ymin=mhcii_lsm2_lsmeans$lower.CL, ymax = mhcii_lsm2_lsmeans$upper.CL)) +
  scale_fill_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(mhcii_ytitle) + xlab("Phase") + labs(fill = "Phase")
plot_mhcii_mean_bar

plot_mhcii_mean_point <- ggplot(mhcii_lsm2_lsmeans, aes(x=Phase2, y = emmean, color = Phase2)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=mhcii_lsm2_lsmeans$lower.CL, ymax = mhcii_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(mhcii_ytitle) + xlab("Phase") + labs(fill = "Phase") +labs(colour= "Phase")
plot_mhcii_mean_point


#is #CD172 associated with phase of infection? 
cd172_null2 <- lme(Log.Total.LP.CD172 ~ 1, random = ~1|Pig, data = df, method = "ML")
summary(cd172_null2) #AIC = 649.8904

cd172_rep <- lme(Log.Total.LP.CD172 ~ Replicate, random = ~1|Pig, data = df, method = "ML")
summary(cd172_rep) #AIC = 626.9545

cd172_phase2 <- lme(Log.Total.LP.CD172 ~ Phase2 + Replicate, random = ~ 1|Pig, data = df, method = "ML")
summary(cd172_phase2) #AIC = 628.8642

anova(cd172_null2,cd172_phase2) #p = <.0001
anova(cd172_rep,cd172_phase2) # p = 0.1513

#obtain least square means and pairwise contrasts for Phase, averaged over replicate
#plot ls-mean cd172
cd172_ytitle <- expression(paste("log" ["10"], " CD172"))

cd172_em2 <- emmeans(cd172_phase2, "Phase2")
cd172_lsm2_lsmeans <- as.data.frame(cd172_em2)
cd172_lsm2_contrasts <- as.data.frame(pairs(cd172_em2))

#use emmeans 
cd172_em2 <- emmeans(cd172_phase2, "Phase2")
pairs(cd172_em2) #pairwise comparisons averaged over the levels of replicate, tukey adjustment
em_plot_cd172 <- plot(cd172_em2, comparisons = TRUE, ylab= "Phase", xlab = cd172_ytitle)
plot_emm_cd172 <- pwpp(cd172_em2, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001, 1)) #plot adj p-values and marginal means
plot_emm_cd172

#save results
cd172_em2_contrasts <- as.data.frame(pairs(cd172_em2))
cd172_em2_emeans <- as.data.frame(emmeans(cd172_phase2, "Phase2"))

plot_cd172_mean_bar <- ggplot(cd172_lsm2_lsmeans, aes(x=Phase2, y = emmean, fill = Phase2)) + 
  geom_bar(stat = "identity", position = "identity") + geom_errorbar(aes(ymin=cd172_lsm2_lsmeans$lower.CL, ymax = cd172_lsm2_lsmeans$upper.CL)) +
  scale_fill_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(cd172_ytitle) + xlab("Phase") + labs(fill = "Phase")
plot_cd172_mean_bar

plot_cd172_mean_point <- ggplot(cd172_lsm2_lsmeans, aes(x=Phase2, y = emmean, color = Phase2)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=cd172_lsm2_lsmeans$lower.CL, ymax = cd172_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(cd172_ytitle) + xlab("Phase") + labs(fill = "Phase")
plot_cd172_mean_point


#is #CD16 associated with phase of infection? 
cd16_null2 <- lme(Log.Total.LP.CD16 ~ 1, random = ~1|Pig, data = df, method = "ML")
summary(cd16_null2) #AIC = 1014.032

cd16_rep <- lme(Log.Total.LP.CD16 ~ Replicate, random = ~1|Pig, data = df, method = "ML")
summary(cd16_rep) #AIC = 1015.417 

cd16_phase2 <- lme(Log.Total.LP.CD16 ~ Phase2 + Replicate, random = ~ 1|Pig, data = df, method = "ML")
summary(cd16_phase2) #AIC = 1008.543 

anova(cd16_null2,cd16_phase2) #p = 0.006
anova(cd16_rep,cd16_phase2) # p =  0.0047

#obtain least square means and pairwise contrasts for Phase, averaged over replicate
cd16_em2 <- emmeans(cd16_phase2, "Phase2")
cd16_lsm2_lsmeans <- as.data.frame(cd16_em2)
cd16_lsm2_contrasts <- as.data.frame(pairs(cd16_em2))

#plot ls-mean cd16
cd16_ytitle <- expression(paste("log" ["10"], " CD16"))

#use emmeans 
cd16_em2 <- emmeans(cd16_phase2, "Phase2")
pairs(cd16_em2) #pairwise comparisons averaged over the levels of replicate, tukey adjustment
em_plot_cd16 <- plot(cd16_em2, comparisons = TRUE, ylab= "Phase", xlab = cd16_ytitle)
plot_emm_cd16 <- pwpp(cd16_em2, method = "pairwise", ylab = "Phase", sort = "F", values = T, add.space = 1, plim = c(0.001,1))
#plot adj p-values and marginal means
plot_emm_cd16

#save results
cd16_em2_contrasts <- as.data.frame(pairs(cd16_em2))
cd16_em2_emeans <- as.data.frame(emmeans(cd16_phase2, "Phase2"))

plot_cd16_mean_bar <- ggplot(cd16_lsm2_lsmeans, aes(x=Phase2, y = emmean, fill = Phase2)) + 
  geom_bar(stat = "identity", position = "identity") + geom_errorbar(aes(ymin=cd16_lsm2_lsmeans$lower.CL, ymax = cd16_lsm2_lsmeans$upper.CL)) +
  scale_fill_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(cd16_ytitle) + xlab("Phase") + labs(fill = "Phase")
plot_cd16_mean_bar

plot_cd16_mean_point <- ggplot(cd16_lsm2_lsmeans, aes(x=Phase2, y = emmean, color = Phase2)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=cd16_lsm2_lsmeans$lower.CL, ymax = cd16_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(cd16_ytitle) + xlab("Phase") + labs(fill = "Phase")
plot_cd16_mean_point

#----------------
#is CD11R1 associated with phase of infection? 
cd11r1_null2 <- lme(Log.Total.LP.CD11R1 ~ 1, random = ~1|Pig, data = df, method = "ML")
summary(cd11r1_null2) #AIC = 1415.282

cd11r1_rep <- lme(Log.Total.LP.CD11R1 ~ Replicate, random = ~1|Pig, data = df, method = "ML")
summary(cd11r1_rep) #AIC = 1383.951 

cd11r1_phase2 <- lme(Log.Total.LP.CD11R1 ~ Phase2 + Replicate, random = ~ 1|Pig, data = df, method = "ML")
summary(cd11r1_phase2) #AIC = 1384.491 

anova(cd11r1_null2,cd11r1_phase2) #p = <.0001
anova(cd11r1_rep,cd11r1_phase2) #p =  0.092

#obtain least square means and pairwise contrasts for Phase, averaged over replicate
cd11r1_em2 <- emmeans(cd11r1_phase2, "Phase2")
cd11r1_lsm2_lsmeans <- as.data.frame(cd11r1_em2)
cd11r1_lsm2_contrasts <- as.data.frame(pairs(cd11r1_em2))

#plot ls-mean cd11r1
cd11r1_ytitle <- expression(paste("log" ["10"], " CD11R1"))

#use emmeans 
cd11r1_em2 <- emmeans(cd11r1_phase2, "Phase2")
pairs(cd11r1_em2) #pairwise comparisons averaged over the levels of replicate, tukey adjustment
em_plot_cd11r1 <- plot(cd11r1_em2, comparisons = TRUE, ylab= "Phase", xlab = cd11r1_ytitle)
#plot adj p-values and marginal means
plot_emm_cd11r1 <- pwpp(cd11r1_em2, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_emm_cd11r1

#save results
cd11r1_em2_contrasts <- as.data.frame(pairs(cd11r1_em2))
cd11r1_em2_emeans <- as.data.frame(emmeans(cd16_phase2, "Phase2"))

plot_cd11r1_mean_bar <- ggplot(cd11r1_lsm2_lsmeans, aes(x=Phase2, y = emmean, fill = Phase2)) + 
  geom_bar(stat = "identity", position = "identity") + geom_errorbar(aes(ymin=cd11r1_lsm2_lsmeans$lower.CL, ymax = cd11r1_lsm2_lsmeans$upper.CL)) +
  scale_fill_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(cd11r1_ytitle) + xlab("Phase") + labs(fill = "Phase")
plot_cd11r1_mean_bar

plot_cd11r1_mean_point <- ggplot(cd11r1_lsm2_lsmeans, aes(x=Phase2, y = emmean, color = Phase2)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=cd11r1_lsm2_lsmeans$lower.CL, ymax = cd11r1_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(cd11r1_ytitle) + xlab("Phase") + labs(fill = "Phase")
plot_cd11r1_mean_point

#-------------
#multi-panel figures
#plot means + sem
plot_grid(plot_mhcii_mean_bar, plot_cd172_mean_bar, plot_cd16_mean_bar, plot_cd11r1_mean_bar)

multi_plot_mean_sem <- plot_grid(plot_mhcii_mean_point + theme(legend.position="none"), plot_cd172_mean_point + theme(legend.position="none"), 
                                 plot_cd16_mean_point + theme(legend.position="none"), plot_cd11r1_mean_point + theme(legend.position="none"))
multi_plot_mean_sem
#get shared legend and add to figure 
legend_a <- get_legend(plot_mhcii_mean_point + guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom"))
multi_plot_mean_sem2 <- plot_grid(multi_plot_mean_sem, legend_a, ncol = 1, rel_heights = c(1,0.1))
multi_plot_mean_sem2
#vertical arrangement 
multi_plot_mean_sem_vertical <- plot_grid(plot_mhcii_mean_point + theme(legend.position="none"), plot_cd172_mean_point + theme(legend.position="none"), 
                                 plot_cd16_mean_point + theme(legend.position="none"), plot_cd11r1_mean_point + theme(legend.position="none"), ncol = 1, nrow = 4)
multi_plot_mean_sem_vertical <- plot_grid(multi_plot_mean_sem_vertical, legend_a, ncol = 1, rel_heights = c(1,0.1))

#could improve figure by including day in x-axis labels, indicating which days fall into which phases. 
#bars start at y=0, so appear upside down.  

#plot adjusted p-values
multi_plot_adj_pvalues <- plot_grid(plot_emm_mhcii,plot_emm_cd172,plot_emm_cd16,plot_emm_cd11r1, labels = c("MHCII", "CD172", "CD16", "CD11R1"))
#plot vertical vertical arrangement
multi_plot_adj_pvalues_vertical <- plot_grid(plot_emm_mhcii,plot_emm_cd172,plot_emm_cd16,plot_emm_cd11r1, ncol = 1, nrow = 5, rel_heights = c(1,1,1,1,0.1))

#plot marginal means and adj-pvalues side-by-side
multi_emmeans_pvalues <- plot_grid(multi_plot_mean_sem_vertical, multi_plot_adj_pvalues_vertical, align = "v", ncol = 2, axis = "b")
multi_emmeans_pvalues

multi_stats <- plot_grid(plot_mhcii_mean_point + theme(legend.position="none"),plot_emm_mhcii, plot_cd172_mean_point + theme(legend.position="none"), plot_emm_cd172,
                         plot_cd16_mean_point + theme(legend.position="none"),plot_emm_cd16, plot_cd11r1_mean_point + theme(legend.position="none"), plot_emm_cd11r1,
                         legend_a, ncol = 2, nrow = 5, rel_heights = c(1,1,1,1,0.1), axis = "b")
multi_stats

#align the panels on the x-axis label 
r1 <- plot_grid(plot_cd172_mean_point + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_emm_cd172 + theme(axis.title.x = element_blank()), align = "h", axis = "b")

r2 <- plot_grid(plot_cd16_mean_point + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_emm_cd16 + theme(axis.title.x = element_blank()), align = "h", axis = "b")

r3 <- plot_grid(plot_cd11r1_mean_point + theme(legend.position="none"),plot_emm_cd11r1, align = "h", axis = "b")

r4 <- plot_grid(plot_mhcii_mean_point + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_emm_mhcii + theme(axis.title.x = element_blank()), align = "h", axis = "b")


multi_plot_final <- plot_grid(r1,r2,r3,r4,legend_a, nrow = 5, rel_heights = c(1,1,1,1,0.1))
multi_plot_final

save_plot("myeloid_multi_mean_adjp_FV.tiff", nrow = 5, ncol = 2, multi_plot_final, base_width = 5.2, dpi = 300)
save_plot("myeloid_multi_mean_adjp_FV_lowres.tiff", nrow = 5, ncol = 2, multi_plot_final, base_width = 5.2, dpi = 100)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#PCA - useful tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/ 

library(factoextra) #ggplot2 based PCA graphics
colnames(df[,45:59]) #for multivariate use log.prop.lp of 15 molecules 
log.prop.lp_forpca <- print(list(colnames(df[,45:59])))
pca_df <- df[,45:59]
pca_df_pigrownames <- df[,c(1:4,6,45:59)]
pca_df_pigrownames$Pig_Rep_Inf_Day_Field <- paste(pca_df_pigrownames$Pig, pca_df_pigrownames$Replicate, pca_df_pigrownames$Infected, pca_df_pigrownames$Day, pca_df_pigrownames$field, sep="_")
pca_df_pigrownames <- pca_df_pigrownames[,6:21]
pca_df_pigrownames <- data.frame(pca_df_pigrownames, row.names = 16)
pca_df_pigrownames2 <- cbind(df$Infected, pca_df_pigrownames)
colnames(pca_df_pigrownames2)[1] <- "Infected"
pca_df_pigrownames_phase <- cbind(df$Phase2, pca_df_pigrownames)
colnames(pca_df_pigrownames_phase)[1] <- "Phase"
pca_df_pigrownames_phase_rep <- cbind(df$Phase2, df$Replicate, pca_df_pigrownames)
colnames(pca_df_pigrownames_phase_rep)[c(1,2)] <- c("Phase", "Replicate")
pca_df_pigrownames_phase_rep_pig <- cbind(df$Phase2, df$Replicate, df$Pig, pca_df_pigrownames)
colnames(pca_df_pigrownames_phase_rep_pig)[c(1,2,3)] <- c("Phase", "Replicate", "Pig")

#compute PCA using prcomp()
res.pca <- prcomp(df[,45:59], scale = T)
res.pca2 <- prcomp(pca_df_pigrownames, scale = T)

#visualise eigenvalues (scree plot)
fviz_eig(res.pca)
fviz_eig(res.pca2) #uses data with rownames set

#graph of individuals 
fviz_pca_ind(res.pca2,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)    # Avoid text overlapping

#graph of variables - positive correlations point to the same side of the plot, 
# negative correlations point to opposite sides of the plot
fviz_pca_var(res.pca2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Avoid text overlapping

#eigenvalues
eig.val <- get_eigenvalue(res.pca2)

# Results for Variables
res.var <- get_pca_var(res.pca2)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

#Contribution of variable (n=15) to a given principle component (in percentage) 
#var.coord = loadings * the component standard deviations
#var.cos2 = var.coord^2
#var.contrib. The contribution of a variable to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component)

#Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- res.pca$rotation
variable.loadings.df <- as.data.frame(loadings)
write.csv(variable.loadings.df, "myeloid_variable_loadings.csv")
sdev <- res.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) #multiplies the loadings by the standard deviation exaplined by each principle component 
head(var.coord[, 1:4])

# Compute Cos2 (variances/square of loadings*standard deviation)
#::::::::::::::::::::::::::::::::::::::::
var.cos2 <- var.coord^2
head(var.cos2[, 1:4])

# Compute contributions
#::::::::::::::::::::::::::::::::::::::::
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])
print(var.contrib[, 1:4])

#export variable contributions as table 
var.contrib.df <- as.data.frame(var.contrib)
write.csv(var.contrib.df, "PCA_myeloid_variable_contributions.csv")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#LINEAR MODEL - ASSESSING ASSOCIATION OF PRINCIPLE COMPONENTS TO PHASE OF INFECTION 
#linear model on principle components 
#is phase of infection associated with PC1, PC2, PC3 or PCR4?
#separate models (not multivariate response) as each principle component is independent
pca_df_pigrownames #pca performed on this dataset 
pca_df_pigrownames_phase #this dataset contains phase 
res.pca3 <- prcomp(pca_df_pigrownames_phase[,2:16], scale = T) #PCA performed on dataset which includes phase, but pca performed on immune molecules only 
summary(res.pca3)
res.pca3$x
res.pca$rotation

#linear model 
dim(res.pca3$rotation) #want the 'x' component from prcomp, not rotation, x provides the matrix of scores
dim(res.pca3$x)
pc1_lm2 <- lm(res.pca3$x[,1] ~ Phase, data = pca_df_pigrownames_phase) #pc1 as dependent variable
summary(pc1_lm) #3.52e-08
pc2_lm <- lm(res.pca3$x[,2] ~ Phase, data = pca_df_pigrownames_phase) #pc2 as dependent variable
summary(pc2_lm) #6.142e-15
pc3_lm <- lm(res.pca3$x[,3] ~ Phase, data = pca_df_pigrownames_phase) #pc3 as dependent variable
summary(pc3_lm) #0.003404
pc4_lm <- lm(res.pca3$x[,4] ~ Phase, data = pca_df_pigrownames_phase) #pc4 as dependent variable
summary(pc4_lm) #0.0003683

library(emmeans)
#pairwise comparisons 
pc1_lm_em <- emmeans(pc1_lm, "Phase")
pairs(pc1_lm_em)

pc2_lm_em <- emmeans(pc2_lm, "Phase")
pairs(pc2_lm_em)

pc3_lm_em <- emmeans(pc3_lm, "Phase")
pairs(pc3_lm_em)

pc2_lm_em <- emmeans(pc2_lm, "Phase")
pairs(pc2_lm_em)

#for comparison, compare prcomp with princomp, following Emily's code for PCA published in www.jimmunol.org/cgi/doi/10.4049/jimmunol.1502632
local({
  .PC <- princomp(pca_df_pigrownames_phase[,2:16], cor = TRUE)
cat("\nComponent loadings:\n")
print(unclass(loadings(.PC))) #compare to prcomp res.pca$rotation
cat("\nComponent variances:\n")
print(.PC$sd^2)
cat("\n")
print(summary(.PC))
screeplot(.PC)
complete_dataset <<- within(pca_df_pigrownames_phase, {
PC15 <- .PC$scores[,15]
PC14 <- .PC$scores[,14]
PC13 <- .PC$scores[,13]
PC12 <- .PC$scores[,12]
PC11 <- .PC$scores[,11]
PC10 <- .PC$scores[,10]
PC9 <- .PC$scores[,9]
PC8 <- .PC$scores[,8]
PC7 <- .PC$scores[,7]
PC6 <- .PC$scores[,6]
PC5 <- .PC$scores[,5]
PC4 <- .PC$scores[,4]
PC3 <- .PC$scores[,3]
PC2 <- .PC$scores[,2]
PC1 <- .PC$scores[,1]
})
})

lm.1 <- lm(PC1 ~ Phase, data = complete_dataset)
summary(lm.1)

#great - lm results are the same using princomp w/ rmcdr or prcomp 
#in addition, when the PC scores are added to the dataframe for each observation, they match the scores accessed using $x from the prcomp() summary
#lesson learnt, is that to access the score for lm then need to take 'x' from prcomp model summary 

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#RE-EVALUATE MODELS WITH PHASE AND REPLICATE AS PREDICTORS, IN LINE WITH MODELS ELSEWHERE TO ADJUST FOR REPLICATE EFFECT
#linear model on principle components 
#is phase of infection associated with PC1, PC2, PC3 or PCR4 (adjusting for replicate)?
#separate models (not multivariate response) as each principle component is independent

pca_df_pigrownames_phase_rep #dataframe with phase and replicate - perform pca on this dataset 
res.pca.pr <- prcomp(pca_df_pigrownames_phase_rep[,3:17], scale = T) #PCA performed on dataset which includes phase, but pca performed on immune molecules only 
summary(res.pca.pr)
res.pca.pr$x #scores
res.pca.pr$rotation #loadings 

#linear model 
dim(res.pca.pr$rotation) #want the 'x' component from prcomp, not rotation, x provides the matrix of scores
dim(res.pca.pr$x)
pc1_lm_pr <- lm(res.pca.pr$x[,1] ~ Phase + Replicate, data = pca_df_pigrownames_phase_rep) #pc1 as dependent variable
summary(pc1_lm_pr) #p = < 2.2e-16   (p = 3.52e-08 with just phase)
pc2_lm_pr <- lm(res.pca.pr$x[,2] ~ Phase + Replicate, data = pca_df_pigrownames_phase_rep) #pc2 as dependent variable
summary(pc2_lm_pr) #p = < 2.2e-16 (p = 6.142e-15 with just phase)
pc3_lm_pr <- lm(res.pca.pr$x[,3] ~ Phase + Replicate, data = pca_df_pigrownames_phase_rep) #pc3 as dependent variable
summary(pc3_lm_pr) #p = 1.198e-15 (p = 0.003404 with just phase)
pc4_lm_pr <- lm(res.pca3$x[,4] ~ Phase + Replicate, data = pca_df_pigrownames_phase_rep) #pc4 as dependent variable
summary(pc4_lm_pr) #p = 1.175e-15 (p = 0.0003683 with just phase)

library(emmeans)
#pairwise comparisons 
pc1_lm_em <- emmeans(pc1_lm_pr, "Phase")
pairs(pc1_lm_em)

pc2_lm_em <- emmeans(pc2_lm_pr, "Phase")
pairs(pc2_lm_em)

pc3_lm_em <- emmeans(pc3_lm_pr, "Phase")
pairs(pc3_lm_em)

pc4_lm_em <- emmeans(pc4_lm_pr, "Phase")
pairs(pc4_lm_em)

#Plot estimated marginal means and adjusted p-values in multipanel figure for linear models for PC1 to PC4
#pairwise comparisons averaged over the levels of replicate, tukey adjustment
#PC1
pc1_em2 <- emmeans(pc1_lm_pr, "Phase")
pc1_lsm2_lsmeans <- as.data.frame(pc1_em2)
pc1_lsm2_contrasts <- as.data.frame(pairs(pc1_em2))

#plot adj p-values 
plot_emm_pc1 <- pwpp(pc1_lm_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_emm_pc1

#plot estimated means 
plot_pc1_mean_point <- ggplot(pc1_lsm2_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc1_lsm2_lsmeans$lower.CL, ymax = pc1_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC1") + xlab("Phase") + labs(fill = "Phase")
plot_pc1_mean_point

#PC2
pc2_em2 <- emmeans(pc2_lm_pr, "Phase")
pc2_lsm2_lsmeans <- as.data.frame(pc2_em2)
pc2_lsm2_contrasts <- as.data.frame(pairs(pc2_em2))

#plot adj p-values 
plot_emm_pc2 <- pwpp(pc2_lm_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_emm_pc2

#plot estimated means 
plot_pc2_mean_point <- ggplot(pc2_lsm2_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc2_lsm2_lsmeans$lower.CL, ymax = pc2_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC2") + xlab("Phase") + labs(fill = "Phase")
plot_pc2_mean_point

#PC3
pc3_em2 <- emmeans(pc3_lm_pr, "Phase")
pc3_lsm2_lsmeans <- as.data.frame(pc3_em2)
pc3_lsm2_contrasts <- as.data.frame(pairs(pc3_em2))

#plot adj p-values 
plot_emm_pc3 <- pwpp(pc3_lm_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_emm_pc3

#plot estimated means 
plot_pc3_mean_point <- ggplot(pc3_lsm2_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc3_lsm2_lsmeans$lower.CL, ymax = pc3_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC3") + xlab("Phase") + labs(fill = "Phase")
plot_pc3_mean_point

#PC4
pc4_em2 <- emmeans(pc4_lm_pr, "Phase")
pc4_lsm2_lsmeans <- as.data.frame(pc4_em2)
pc4_lsm2_contrasts <- as.data.frame(pairs(pc4_em2))

#plot adj p-values 
plot_emm_pc4 <- pwpp(pc4_lm_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_emm_pc4

#plot estimated means 
plot_pc4_mean_point <- ggplot(pc4_lsm2_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc4_lsm2_lsmeans$lower.CL, ymax = pc4_lsm2_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC4") + xlab("Phase") + labs(fill = "Phase")
plot_pc4_mean_point

#multipanel figure
legend_b <- get_legend(plot_pc1_mean_point + guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom"))

pc_r1 <- plot_grid(plot_pc1_mean_point + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_emm_pc1 + theme(axis.title.x = element_blank()), align = "h", axis = "b")

pc_r2 <- plot_grid(plot_pc2_mean_point + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_emm_pc2 + theme(axis.title.x = element_blank()), align = "h", axis = "b")

pc_r3 <- plot_grid(plot_pc3_mean_point + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_emm_pc3 + theme(axis.title.x = element_blank()), align = "h", axis = "b")

pc_r4 <- plot_grid(plot_pc4_mean_point + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_emm_pc4 + theme(axis.title.x = element_blank()), align = "h", axis = "b")


multi_plot_final_PCA <- plot_grid(pc_r1,pc_r2,pc_r3,pc_r4,legend_b, nrow = 5, rel_heights = c(1,1,1,1,0.1))
multi_plot_final_PCA

save_plot("myeloid_multi_mean_adjp_PCA_FV.tiff", nrow = 5, ncol = 2, multi_plot_final_PCA, base_width = 5.2, dpi = 300)
save_plot("myeloid_multi_mean_adjp_PCA_FV_lowres.tiff", nrow = 5, ncol = 2, multi_plot_final_PCA, base_width = 5.2, dpi = 100)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#RE-EVALUATE MODELS WITH PHASE AND REPLICATE AS PREDICTORS, RANDOM EFFECT OF PIG, IN LINE WITH MODELS ELSEWHERE TO ADJUST FOR REPLICATE EFFECT AND BETWEEN PIG VARIATION 
#linear model on principle components 
#is phase of infection associated with PC1, PC2, PC3 or PCR4 (adjusting for replicate)?
#separate models (not multivariate response) as each principle component is independent

pca_df_pigrownames_phase_rep_pig #dataframe with phase and replicate - perform pca on this dataset 
res.pca.prp <- prcomp(pca_df_pigrownames_phase_rep_pig[,4:18], scale = T) #PCA performed on dataset which includes phase, but pca performed on immune molecules only 
summary(res.pca.prp)
res.pca.prp$x #scores
res.pca.prp$rotation #loadings 
lme_pca_df <- cbind(pca_df_pigrownames_phase_rep_pig[,1:3],res.pca.prp$x)

#linear mixed effect model - random effect (intercept) for pig 
dim(res.pca.prp$rotation) #want the 'x' component from prcomp, not rotation, x provides the matrix of scores
dim(res.pca.prp$x)
pc1_lme_prp <- lme(PC1 ~ Phase + Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML") #pc1 as dependent variable
summary(pc1_lme_prp) 
pc1_lme_un <- lme(PC1 ~ Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML")
#overall significance of phase 
anova(pc1_lme_prp, pc1_lme_un) #p = 0.0386 (p = < 2.2e-16 with phase + rep) #much more sensible!!

pc2_lme_prp <- lme(PC2 ~ Phase + Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML") #pc2 as dependent variable
summary(pc2_lme_prp) 
pc2_lme_un <- lme(PC2 ~ Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML")
#overall significance of phase 
anova(pc2_lme_prp, pc2_lme_un) #p = 0.0335 (p = < 2.2e-16 with phase + rep), again much better

pc3_lme_prp <- lme(PC3 ~ Phase + Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML") #pc3 as dependent variable
summary(pc3_lme_prp) 
pc3_lme_un <- lme(PC3 ~ Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML")
#overall significance of phase 
anova(pc3_lme_prp, pc3_lme_un) #p = 0.5992 (p = 1.198e-15  with phase + rep)

pc4_lme_prp <- lme(PC4 ~ Phase + Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML") #pc4 as dependent variable
summary(pc4_lme_prp) 
pc4_lme_un <- lme(PC4 ~ Replicate, random = ~ 1|Pig, data = lme_pca_df, method = "ML")
#overall significance of phase 
anova(pc4_lme_prp, pc4_lme_un) #p = 0.4896  (p = 1.175e-15 with phase + rep)

library(emmeans)
#pairwise comparisons 
pc1_lme_prp_em <- emmeans(pc1_lme_prp, "Phase")
pairs(pc1_lme_prp_em , adjust = "sidak")
pairs(pc1_lme_prp_em)

pc2_lme_prp_em <- emmeans(pc2_lme_prp, "Phase")
pairs(pc2_lme_prp_em, adjust = "fdr")

pc3_lme_prp_em <- emmeans(pc3_lme_prp, "Phase")
pairs(pc3_lme_prp_em)

pc4_lme_prp_em <- emmeans(pc4_lme_prp, "Phase")
pairs(pc4_lme_prp_em)

#plot estimates marginal means and adjusted p-values for 

#plot adj p-values 
plot_pc1_lme_prp_em <- pwpp(pc1_lme_prp_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_pc1_lme_prp_em
plot_pc2_lme_prp_em <- pwpp(pc2_lme_prp_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_pc2_lme_prp_em
plot_pc3_lme_prp_em <- pwpp(pc3_lme_prp_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_pc3_lme_prp_em
plot_pc4_lme_prp_em <- pwpp(pc4_lme_prp_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))
plot_pc4_lme_prp_em

#plot estimated means
#PC1
pc1_lme_prp_lsmeans <- as.data.frame(pc1_lme_prp_em)
plot_pc1_mean_point <- ggplot(pc1_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc1_lme_prp_lsmeans$lower.CL, ymax = pc1_lme_prp_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC1") + xlab("Phase") + labs(fill = "Phase")
plot_pc1_mean_point

plot_pc1_mean_point_sem <- ggplot(pc1_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC1") + xlab("Phase") + labs(fill = "Phase")
plot_pc1_mean_point_sem

#PC2
pc2_lme_prp_lsmeans <- as.data.frame(pc2_lme_prp_em)
plot_pc2_mean_point <- ggplot(pc2_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc2_lme_prp_lsmeans$lower.CL, ymax = pc2_lme_prp_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC2") + xlab("Phase") + labs(fill = "Phase")
plot_pc2_mean_point

plot_pc2_mean_point_sem <- ggplot(pc2_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC2") + xlab("Phase") + labs(fill = "Phase")
plot_pc2_mean_point_sem

#PC3
pc3_lme_prp_lsmeans <- as.data.frame(pc3_lme_prp_em)
plot_pc3_mean_point <- ggplot(pc3_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc3_lme_prp_lsmeans$lower.CL, ymax = pc3_lme_prp_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC3") + xlab("Phase") + labs(fill = "Phase")
plot_pc3_mean_point

plot_pc3_mean_point_sem <- ggplot(pc3_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC3") + xlab("Phase") + labs(fill = "Phase")
plot_pc3_mean_point_sem

#PC4
pc4_lme_prp_lsmeans <- as.data.frame(pc4_lme_prp_em)
plot_pc4_mean_point <- ggplot(pc4_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=pc4_lme_prp_lsmeans$lower.CL, ymax = pc4_lme_prp_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC4") + xlab("Phase") + labs(fill = "Phase")
plot_pc4_mean_point

plot_pc4_mean_point_sem <- ggplot(pc4_lme_prp_lsmeans, aes(x=Phase, y = emmean, color = Phase)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab("PC4") + xlab("Phase") + labs(fill = "Phase")
plot_pc4_mean_point_sem

#multipanel figure
legend_b <- get_legend(plot_pc1_mean_point + guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom"))

pc_r1_lme <- plot_grid(plot_pc1_mean_point_sem + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_pc1_lme_prp_em + theme(axis.title.x = element_blank()), align = "h", axis = "b")
pc_r1_lme

pc_r2_lme <- plot_grid(plot_pc2_mean_point_sem + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_pc2_lme_prp_em + theme(axis.title.x = element_blank()), align = "h", axis = "b")

pc_r3_lme <- plot_grid(plot_pc3_mean_point_sem + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_pc3_lme_prp_em + theme(axis.title.x = element_blank()), align = "h", axis = "b")

pc_r4_lme <- plot_grid(plot_pc4_mean_point_sem + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_pc4_lme_prp_em + theme(axis.title.x = element_blank()), align = "h", axis = "b")


multi_plot_final_PCA_lme <- plot_grid(pc_r1_lme,pc_r2_lme,pc_r3_lme,pc_r4_lme,legend_b, nrow = 5, rel_heights = c(1,1,1,1,0.1))
multi_plot_final_PCA_lme

save_plot("myeloid_multi_mean_adjp_PCA_lme_FV.tiff", nrow = 5, ncol = 2, multi_plot_final_PCA_lme, base_width = 5.2, dpi = 300)
save_plot("myeloid_multi_mean_adjp_PCA_FV_lme_lowres.tiff", nrow = 5, ncol = 2, multi_plot_final_PCA_lme, base_width = 5.2, dpi = 100)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Interactions w/ MHCII and CD16; MHCII and CD172 
#Mixed effect models 

mhciicd16_lme_prp <- lme(Log.Total.Prop.LP.MHCII.CD16 ~ Phase2 + Replicate, random = ~ 1|Pig, data = df, method = "ML") #pc1 as dependent variable
summary(mhciicd16_lme_prp) 
mhciicd16_lme_un <- lme(Log.Total.Prop.LP.MHCII.CD16 ~ Replicate, random = ~ 1|Pig, data = df, method = "ML")
#overall significance of phase 
anova(mhciicd16_lme_prp, mhciicd16_lme_un) #p = 0.0018

mhciicd172_lme_prp <- lme(Log.Total.Prop.MHCII.CD172 ~ Phase2 + Replicate, random = ~ 1|Pig, data = df, method = "ML") #pc1 as dependent variable
summary(mhciicd172_lme_prp) 
mhciicd172_lme_un <- lme(Log.Total.Prop.MHCII.CD172 ~ Replicate, random = ~ 1|Pig, data = df, method = "ML")
#overall significance of phase 
anova(mhciicd172_lme_prp, mhciicd172_lme_un) #p = 0.1257

#estimated marginal means
mhciicd16_lme_prp_em <- emmeans(mhciicd16_lme_prp, "Phase2")
pairs(mhciicd16_lme_prp_em)

mhciicd172_lme_prp_em <- emmeans(mhciicd172_lme_prp, "Phase2")
pairs(mhciicd172_lme_prp_em)

mhciicd16_lme_p_lsmeans <- as.data.frame(mhciicd16_lme_prp_em)
plot_mhciicd16_mean_point_sem <- ggplot(mhciicd16_lme_p_lsmeans, aes(x=Phase2, y = emmean, color = Phase2)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(y_title_mhciicd16) + xlab("Phase") + labs(fill = "Phase2")
plot_mhciicd16_mean_point_sem

mhciicd172_lme_p_lsmeans <- as.data.frame(mhciicd172_lme_prp_em)
plot_mhciicd172_mean_point_sem <- ggplot(mhciicd172_lme_p_lsmeans, aes(x=Phase2, y = emmean, color = Phase2)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(y_title_mhciicd172) + xlab("Phase") + labs(fill = "Phase2")
plot_mhciicd172_mean_point_sem

plot_mhciicd172_mean_point <- ggplot(mhciicd172_lme_p_lsmeans, aes(x=Phase2, y = emmean, color = Phase2)) + geom_point(aes(y=emmean)) + geom_errorbar(aes(ymin=mhciicd172_lme_p_lsmeans$lower.CL, ymax = mhciicd172_lme_p_lsmeans$upper.CL), width = .6) +
  scale_color_manual(values = c("blue", "red", "red", "red", "red", "darkgray")) + ylab(y_title_mhciicd172) + xlab("Phase") + labs(fill = "Phase")
plot_mhciicd172_mean_point

#plot pairwise comparisons
plot_mhciicd16_lme_prp_em <- pwpp(mhciicd16_lme_prp_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))

plot_mhciicd172_lme_prp_em <- pwpp(mhciicd172_lme_prp_em, method = "pairwise", ylab = "Phase", sort = "F", add.space = 1, plim = c(0.001,1))


#multipanel figure
legend_combinations <- get_legend(plot_mhciicd16_mean_point_sem + guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom"))

mhciicd16_lme_plot <- plot_grid(plot_mhciicd16_mean_point_sem + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_mhciicd16_lme_prp_em + theme(axis.title.x = element_blank()), align = "h", axis = "b")

mhciicd172_lme_plot <- plot_grid(plot_mhciicd172_mean_point_sem + theme(legend.position="none") + theme(axis.title.x = element_blank()),plot_mhciicd172_lme_prp_em + theme(axis.title.x = element_blank()), align = "h", axis = "b")

multi_plot_final_myeloidcombinations_lme <- plot_grid(mhciicd16_lme_plot,mhciicd172_lme_plot,legend_b, nrow = 3, rel_heights = c(1,1,0.1))
multi_plot_final_myeloidcombinations_lme

save_plot("myeloidcombinations_multi_mean_adjp_lme_FV.tiff", nrow = 3, ncol = 2, multi_plot_final_myeloidcombinations_lme, base_width = 5.2, dpi = 300)
save_plot("myeloidcombinations_multi_mean_adjp_lme_FV_lowres.tiff", nrow = 3, ncol = 2, multi_plot_final_myeloidcombinations_lme, base_width = 5.2, dpi = 100)


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#PCA results for individuals - full data including value for each field for each molecule/day 

# Results for individuals
res.ind <- get_pca_ind(res.pca2)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

#colour individuals by groups 
groups <- as.factor(pca_df_pigrownames2$Infected)
myColors <- c("blue", "red", "darkgray")

pca_ind_group <- fviz_pca_ind(res.pca,
                              col.ind = groups, # color by groups
                              palette = myColors,
                              addEllipses = TRUE, # Concentration ellipses
                              ellipse.type = "confidence",
                              legend.title = "Infection Status",
                              repel = TRUE)
pca_ind_group

#::::::::::::::::::::::::::::
#For visualising PCA on individuals, take mean of fields as too many points for plotting when each field is considered 
library(plotrix)
pca_df_meanoffields <- df[,c(1:6,45:59)]
pca_df_meanoffields2 <- pca_df_meanoffields %>%
  group_by(Pig, Day, Replicate, Infected, Phase2) %>%
  summarise(MHCII = mean(Log.Prop.LP.MHCII), CD172 = mean(Log.Prop.LP.CD172),
            MHCIICD172 = mean(Log.Prop.LP.MHCIICD172), CD16 = mean(Log.Prop.LP.CD16),
            MHCIICD16 = mean(Log.Prop.LP.MHCIICD16), CD172CD16 = mean(Log.Prop.LP.CD172CD16),
            MHCIICD172CD16 = mean(Log.Prop.LP.MHCIICD172CD16), CD11R1 = mean(Log.Prop.LP.CD11R1),
            MHCIICD11R1 = mean(Log.Prop.LP.MHCIICD11R1), CD172CD11R1 = mean(Log.Prop.LP.CD172CD11R1),
            MHCIICD172CD11R1 = mean(Log.Prop.LP.MHCIICD172CD11R1), CD16CD11R1 = mean(Log.Prop.LP.CD16CD11R1),
            MHCIICD16CD11R1 = mean(Log.Prop.LP.MHCIICD16CD11R1), CD172CD16CD11R1 = mean(Log.Prop.LP.CD172CD16CD11R1),
            MHCIICD172CD16CD11R1= mean(Log.Prop.LP.MHCIICD172CD16CD11R1))
pca_df_meanoffields3 <- pca_df_meanoffields2
pca_df_meanoffields3$Pig_Day_Rep_Infected <- paste(pca_df_meanoffields3$Pig, pca_df_meanoffields3$Day, pca_df_meanoffields3$Replicate, pca_df_meanoffields3$Infected, sep="_")
pca_df_meanoffields3 <- pca_df_meanoffields3[,6:21]
pca_df_meanoffields3_rownames <- data.frame(pca_df_meanoffields3, row.names = 16)
pca_df_meanoffields3_rownames_group <- cbind(pca_df_meanoffields2$Infected, pca_df_meanoffields3_rownames)
colnames(pca_df_meanoffields3_rownames_group)[1] <- "Infected"


#perform PCA on mean of fields for plotting of individuals (n=64 vs n=627)
df.pca.fieldmean <- pca_df_meanoffields3_rownames_group[,2:16]
res.pca.fieldmean <- prcomp(df.pca.fieldmean, scale = T)

#colour individuals by groups
#colour by infection status 
groups_fieldmean <- as.factor(pca_df_meanoffields3_rownames_group$Infected)

pca_ind_fieldmean_group <- fviz_pca_ind(res.pca.fieldmean,
                              col.ind = groups_fieldmean, # color by groups
                              palette = myColors,
                              addEllipses = TRUE, # Concentration ellipses
                              ellipse.type = "confidence",
                              legend.title = "Infection Status",
                              repel = TRUE, label = "none")
pca_ind_fieldmean_group 
#saved as myeloid_variables_pca_plot_indiv_64obs.tiff

#colour by phase 

pca_df_meanoffields3_rownames_phase <- cbind(pca_df_meanoffields2$Phase2, pca_df_meanoffields3_rownames)
colnames(pca_df_meanoffields3_rownames_phase)[1] <- "Phase"
phase_fieldmean <- as.factor(pca_df_meanoffields3_rownames_phase$Phase)
mycolors_phase <- c("blue", "yellow", "tan1", "tomato3", "red4", "darkgray")
pca_ind_fieldmean_phase <- fviz_pca_ind(res.pca.fieldmean,
                                        col.ind = phase_fieldmean, # color by groups
                                        palette = mycolors_phase,
                                        addEllipses = TRUE, # Concentration ellipses
                                        ellipse.type = "confidence",
                                        legend.title = "Infection Status",
                                        repel = TRUE, label = "none")
pca_ind_fieldmean_phase
#saved as myeloid_variables_pca_plot_indiv_64obs_phase.tiff


