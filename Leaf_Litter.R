############
#Augusta Creek leaf litter weights
############
#Packages
library(dplyr)
library(data.table)
library(ggplot2)
library(plyr)
library(vegan)
library(MASS)
library(indicspecies)
library(funrar)
library(doBy)
library(pastecs)
library(e1071) 
library(rstatix)
library(ggpubr)
library(car)
library(emmeans)
library(tidyverse)
library(lme4)
library(reshape)
library(ecole)
library(devtools)

#Functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

normalit<-function(m){
  ((m)/max(m))*100
}

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

#function parameter
#@model: a lm or glm fitted model
#@name_f: a character string with the name of the categorical (factor) variable
#@name_x: a character string with the name of the interacting variable, by default the intercept
add_se <- function(model,name_f,name_x="Intercept"){
  #grab the standard error of the coefficients
  se_vec <- summary(model)$coefficients[,2]
  if(name_x=="Intercept"){
    #the stabdard error of the intercept
    se_x <- se_vec[1]
    #get the level-specific standard errors
    se_f <- se_vec[grep(name_f,names(se_vec))]
    se_f <- se_f[-grep(":",names(se_f))]
    #get the covariance between the intercept and the level-specific parameters
    vcov_f <- vcov(model)[grep(name_f,rownames(vcov(model))),grep(name_x,colnames(vcov(model)))]
    vcov_f <- vcov_f[-grep(":",names(vcov_f))]
    #the estimated average value at each level
    coef_f <- coef(model)[1]+coef(model)[names(vcov_f)]
  }
  else{
    #similar code for the case of another variable than the intercept
    se_x <- se_vec[name_x]
    se_f <- se_vec[grep(name_f,names(se_vec))]
    se_f <- se_f[grep(":",names(se_f))]
    vcov_f <- vcov(model)[grep(name_f,rownames(vcov(model))),grep(name_x,colnames(vcov(model)))][,1]
    vcov_f <- vcov_f[grep(":",names(vcov_f))]
    coef_f <- coef(model)[name_x]+coef(model)[names(vcov_f)]
  }
  #compute the summed SE
  se_f <- sqrt(se_x**2+se_f**2+2*vcov_f)
  #create the output dataframe
  out <- data.frame(Coef=coef_f,SE=se_f)
  return(out)
}

#Color vectors
leaftaxacolvec<-c("#7fc97f","#beaed4","#fdc086","#386cb0")
leaftaxacolvec_nc<-c("#7fc97f","#beaed4","#386cb0")
macrocolvec<-c("#1b9e77", "#d95f02","#7570b3")
reach_col_vec<-c("#1b9e77", "#d95f02", "#7570b3")
twoleaftaxacolvec<-c("#7fc97f","#386cb0")
ffgcolvec<-c("#1b9e77","#7570b3","#66a61e","#e7298a","#d95f02")

#Create water chemistry dataset
H2OCh<-read.csv("AC_H2O_Chem.csv", sep = ",", header = T )
anova(lm(Water_Temp_C~Reach+Date, data=H2OCh))
#date significant
anova(lm(DO_perc~Reach+Date, data=H2OCh))
#nothing significant
anova(lm(ORP~Reach+Date, data=H2OCh))
#date significant
anova(lm(pH~Reach+Date, data=H2OCh))
#Reach significant


Dens<-read.csv("Densiometer_AC.csv", sep = ",", header = T )
Dens$Sum<-rowSums(Dens[,3:ncol(Dens)])
Dens$PercShade<-Dens$Sum*0.26
anova(lm(PercShade~Site..Reach.+Date, data=Dens))
#site and reach significantly different
Dens %>%
  group_by(Site..Reach.) %>%
  get_summary_stats(PercShade, type = "mean_se")

H2OChmeans<-aggregate(. ~ Reach, H2OCh, function(x) c(mean = mean(x)))
H2OChse<-aggregate(. ~ Reach, H2OCh, function(x) c(stder = sd(x)/sqrt(sum(!is.na(x)))))

#Upload AFDM datasheet
AFDM<-read.csv("AFDM_Data_Sheet.csv", sep = ",", header = T )
#Create column for leaf pack dry mass
AFDM$TotalWetM<-AFDM$D.TotalWetM-AFDM$DishMass
AFDM$PartialWM<-AFDM$D.PartialWM-AFDM$DishMass
AFDM$PDryM<-AFDM$D.PDryM-AFDM$DishMass
AFDM$ProportionLeaf<-AFDM$PDryM/AFDM$PartialWM
AFDM$LPDM<-AFDM$TotalWetM*AFDM$ProportionLeaf
#calculate afdm
#A
AFDM$ADM<-AFDM$AP.DM-AFDM$APanM
AFDM$AAM<-AFDM$AP.AshM-AFDM$APanM
AFDM$APercAsh<-(AFDM$AAM/AFDM$ADM)
AFDM$APercOrg<-1-AFDM$APercAsh
AFDM$AAFDM<-AFDM$LPDM*AFDM$APercOrg
#B
AFDM$BDM<-AFDM$BP.DM-AFDM$BPM
AFDM$BAM<-AFDM$BP.AM-AFDM$BPM
AFDM$BPercAsh<-(AFDM$BAM/AFDM$BDM)
AFDM$BPercOrg<-1-AFDM$BPercAsh
AFDM$BAFDM<-AFDM$LPDM*AFDM$BPercOrg
#AverageAFDM
AFDM$AFDM<-rowMeans(AFDM[c('AAFDM', 'BAFDM')], na.rm=TRUE)
#Upload metadata
LPMetadata<-read.csv("Field_Exp_Sample_Metadata.csv", sep = ",", header = T)
LPMetadata$Reach<-factor(LPMetadata$Reach, levels=c("US","Gap","DS"))

#Combine metadata with dataset
LeafPackM<-merge(AFDM, LPMetadata, by="Pack_ID")
#find initial AFDM
AFDMI<-subset(LeafPackM, Time_Point=="0")
AFDMIs<-subset(AFDMI, select = c(Leaf_Type,Reach,AFDM) )
AFDMIsag<-aggregate(AFDMIs$AFDM, by=list("Reach"=AFDMIs$Reach,"Leaf_Type"=AFDMIs$Leaf_Type),
                    FUN=mean, na.rm=TRUE)
colnames(AFDMIsag)[colnames(AFDMIsag)=="x"] <- "initialAFDM"
LeafPackMi<-merge(LeafPackM, AFDMIsag, by=c("Leaf_Type","Reach"))
#convert AFDM for each leaf pack to % AFDM remaining
LeafPackMi[is.nan(LeafPackMi)] <- 0
LeafPackMi$percAFDMremain<-(LeafPackMi$AFDM/LeafPackMi$initialAFDM)*100

#determine if there are any outliers
LeafPackMi %>%
  group_by(Reach, Leaf_Type, Days_Exposure) %>%
  identify_outliers(percAFDMremain)
#no outliers

#regress the natural log of mean % AFDM remaining on days exposure for each leaf type and reach
LeafPackMi$lnpercAFDMremain<-log(LeafPackMi$percAFDMremain+1)
afdm.lm<-lm(lnpercAFDMremain~Days_Exposure*Leaf_Type*Reach, data=LeafPackMi)
summary(afdm.lm)
#Intercept, Days and Days:Cotton interaction significant
#try again with intercept set as constant
afdm.lm.ni<-lm(lnpercAFDMremain~Days_Exposure+Days_Exposure:Leaf_Type+Days_Exposure:Reach+
                 Days_Exposure:Reach:Leaf_Type, data=LeafPackMi)
summary(afdm.lm.ni)
#intercept, time, time:buckthorn and time:cotton significant
#output to excel file
#use variance covariance matrix to calculate SE for -k estimates
vcovafdm<-vcov(afdm.lm.ni)
write.csv(vcovafdm,'vcovafdm.csv')
#create k value and ses .csv based on calculations
#or use this code
anova(afdm.lm.ni)
#only leaf type significant
afdm.lst.ni<- lstrends(afdm.lm.ni, "Leaf_Type", var="Days_Exposure")
pairs(afdm.lst.ni)

#Calculate mean and sd for each leaf type and sampling day
AFDMsummary<-summarySE(LeafPackMi, measurevar=c("percAFDMremain"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
AFDMsummary<-subset(AFDMsummary, percAFDMremain!="NaN")
AFDMsummary$Time_Point[AFDMsummary$Time_Point == 1] <- 8
AFDMsummary$Time_Point[AFDMsummary$Time_Point == 2] <- 41
AFDMsummary$Time_Point[AFDMsummary$Time_Point == 3] <- 68
AFDMsummary$Time_Point[AFDMsummary$Time_Point == 4] <- 98
AFDMsummary$Reach<-revalue(AFDMsummary$Reach, c("US"="Upstream"))
AFDMsummary$Reach<-revalue(AFDMsummary$Reach, c("DS"="Downstream"))

#plot on y axis % remaining and x axis days exposure
ggplot(AFDMsummary, aes(x=Time_Point, y=percAFDMremain+1, color=Leaf_Type)) +
  geom_smooth(aes(group=Leaf_Type), method=lm, se=F, fullrange=TRUE, size=2)+
  geom_errorbar(aes(ymin=percAFDMremain+1-se, ymax=percAFDMremain+1+se), width=2) +
  geom_point(size=3) +
  xlab("Days of Exposure") +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  scale_y_log10(name="%AFDM Remaining +1 (± SEM)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)
ggplot(AFDMsummary, aes(x=Time_Point, y=percAFDMremain, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=percAFDMremain-se, ymax=percAFDMremain+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean %AFDM Remaining (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#Upload kvalue means and se's 
KMeans<-read.csv("kvalues.csv", sep = ",", header = T)
KMeans$negk<--1*(KMeans$k)
KMeans$Reach<-factor(KMeans$Reach, levels=c("Upstream","Gap","Downstream"))
ggplot(KMeans, aes(x=Reach, y=negk, fill=Leaf_Type)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=negk-SE, ymax=negk+SE), width=.2,
                position=position_dodge(.9))+
  xlab("Gap Location") +
  ylab("Mean Decomposition Rate (-k ± SE)") +
  scale_fill_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=16),
        legend.title=element_text(size=22),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 20))

#Upload historical leaf processing data
HistK<-read.csv("Processing_coefficients.csv", sep = ",", header = T)
HistK$Year<-as.factor(HistK$Year)
HistK[is.na(HistK)] <- 0
HistK$Location<-revalue(HistK$Location, c("US"="Upstream"))
HistK$Location<-revalue(HistK$Location, c("DS"="Downstream"))
HistK$Location<- factor(HistK$Location, levels=rev(levels(HistK$Location)))
ggplot(HistK, aes(x=Year, y=k, color=Taxa)) +
  geom_pointrange(size=1, aes(ymin=k-SE, ymax=k+SE)) +
  xlab("Year of study") +
  scale_y_continuous(name="Leaf Processing Coefficient (-k ± SEM)") +
  scale_color_manual(values=twoleaftaxacolvec,name = "Leaf Taxa") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=24),axis.title.y=element_text(size=22),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=18),
        legend.title=element_text(size=24),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  facet_wrap(.~Location, scales="free_x")

#Upload macroinvertebrate data
LPMacros<-read.csv("Leaf_Pack_Macroinvertebrates.csv", sep = ",", header = T)
LPMacrosM<-merge(LPMacros, LPMetadata, by="Pack_ID")
LPMacrosM$rowsum<-rowSums(LPMacrosM[,2:22])

LPMacrosMno<-subset(LPMacrosM, rowsum!=0)
LPMacroCom<-LPMacrosMno[,2:22]
sum(LPMacroCom)
#644 individuals
colSums(LPMacroCom)
#Taeniopterygidae most abundant with 167 individuals

LPMacroEnv<-LPMacrosMno[,c(1,23:ncol(LPMacrosMno))]
LPMacroComtot<-LPMacrosM[,2:22]
LPMacroComtotbray0<-as.matrix(bray0(LPMacroComtot))
LPMacroComRA<-data.frame(make_relative(as.matrix(LPMacrosM[,2:22])))
LPMacroComRA[is.nan(LPMacroComRA)] <- 0
LPMacroEnvtot<-LPMacrosM[,c(1,23:ncol(LPMacrosM))]
LPMacroRA<-merge(LPMacroComRA, LPMacroEnvtot, by=0)
Reachtot<-droplevels(as.factor(LPMacroEnvtot$Reach))

#delete time point 0 days
LPMacrosM_no<-subset(LPMacrosM, Time_Point!=0)

LPMacro_noComtot<-LPMacrosM_no[,2:22]
LPMacro_noComtotbray0<-as.matrix(bray0(LPMacro_noComtot))
LPMacro_noComRA<-data.frame(make_relative(as.matrix(LPMacrosM_no[,2:22])))
LPMacro_noComRA[is.nan(LPMacro_noComRA)] <- 0
LPMacro_noEnvtot<-LPMacrosM_no[,c(1,23:ncol(LPMacrosM_no))]
LPMacro_noRA<-merge(LPMacro_noComRA, LPMacro_noEnvtot, by=0)

stat.desc(LPMacro_noRA$Taeniopterygidae)

#permanova
macperm<-adonis2(as.dist(LPMacro_noComtotbray0)~Reach*Leaf_Type*Days_Exposure, data=LPMacro_noEnvtot, permutations=999)
#days, reach, leaf type, and reach:leaf type signficant

Reachtot<-droplevels(as.factor(LPMacroEnvtot$Reach))
Leaf_Taxatot<-as.factor(LPMacroEnvtot$Leaf_Type)
ReachLeaftot<-as.factor(paste(LPMacro_noEnvtot$Reach,LPMacroEnvtot$Leaf_Type))
ReachDaystot<-as.factor(paste(LPMacro_noEnvtot$Reach,LPMacroEnvtot$Time_Point))

Reachtot_no<-droplevels(as.factor(LPMacro_noEnvtot$Reach))
Leaf_Taxatot_no<-as.factor(LPMacro_noEnvtot$Leaf_Type)
ReachLeaftot_no<-as.factor(paste(LPMacro_noEnvtot$Reach,LPMacro_noEnvtot$Leaf_Type))
ReachDaystot_no<-as.factor(paste(LPMacro_noEnvtot$Reach,LPMacro_noEnvtot$Time_Point))

#NMDS analysis
AC_Macroinvertebrate_NMDS<-metaMDS(as.dist(LPMacroComtotbray0))
#stress=0.14
stressplot(AC_Macroinvertebrate_NMDS)
ordiplot(AC_Macroinvertebrate_NMDS, type="n")
with(AC_Macroinvertebrate_NMDS, points(AC_Macroinvertebrate_NMDS, display="sites", col=leaftaxacolvec[Leaf_Taxatot], pch=19, pt.bg=leaftaxacolvec))
with(AC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(Leaf_Taxatot), bty="n", col=leaftaxacolvec, pch=19, pt.bg=leaftaxacolvec))
with(AC_Macroinvertebrate_NMDS, ordiellipse(AC_Macroinvertebrate_NMDS, Leaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_Macroinvertebrate_NMDS, ordiellipse(AC_Macroinvertebrate_NMDS, Leaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_Macroinvertebrate_NMDS, ordiellipse(AC_Macroinvertebrate_NMDS, Leaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Cotton"))
with(AC_Macroinvertebrate_NMDS, ordiellipse(AC_Macroinvertebrate_NMDS, Leaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

Reachtot<-revalue(Reachtot, c("US"="Upstream"))
Reachtot<-revalue(Reachtot, c("DS"="Downstream"))
ordiplot(AC_Macroinvertebrate_NMDS, type="n")
with(AC_Macroinvertebrate_NMDS, points(AC_Macroinvertebrate_NMDS, display="sites", col=reach_col_vec[Reachtot], pch=19))
with(AC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(Reachtot), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_Macroinvertebrate_NMDS, ordiellipse(AC_Macroinvertebrate_NMDS, Reachtot, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "Upstream"))
with(AC_Macroinvertebrate_NMDS, ordiellipse(AC_Macroinvertebrate_NMDS, Reachtot, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_Macroinvertebrate_NMDS, ordiellipse(AC_Macroinvertebrate_NMDS, Reachtot, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "Downstream"))
#############

#Diversity metrics
#Calculate richess
LPMacro_noEnvtot$Richness<-rowSums(LPMacro_noComtot > 0)
LPMacro_noEnvtot$Time_Point_cat<-as.factor(LPMacro_noEnvtot$Time_Point)
#summary stats
LPMacro_noEnvtot %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(Richness, type = "mean_sd")
#check for outliers
LPMacro_noEnvtot %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Richness)
#no outliers
#visualize
ggboxplot(LPMacro_noEnvtot, x = "Reach", y = "Richness",
  color = "Time_Point_cat", palette = "jco",
  facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LPMacro_noEnvtot %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(Richness)
#not normal
ggqqplot(LPMacro_noEnvtot, "Richness", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#Transform Richness using boxcox
MR.lm<- lm((Richness+1) ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot)
ggqqplot(residuals(MR.lm))
shapiro_test(residuals(MR.lm))
#0.000000556 not normal, see boxcox
boxcox(MR.lm)
#labmda=0. log transform
LPMacro_noEnvtot$log10Richness<-log10(LPMacro_noEnvtot$Richness+1)
LPMacro_noEnvtot %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(log10Richness)
#not normal
ggqqplot(LPMacro_noEnvtot, "log10Richness", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#Use this transformation
logRichmac.lm<- lm(log10Richness ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot)
ggqqplot(residuals(logRichmac.lm))
shapiro_test(residuals(logRichmac.lm))
#residuals fall aprox in qqplot, but shapiro test still <0.05
#Note that, if your sample size is greater than 50, the normal QQ plot is preferred 
#because at larger sample sizes the Shapiro-Wilk test becomes very sensitive even to a
#minor deviation from normality.
#looks good, test homogeniety of variance
leveneTest(log10Richness ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot)
#not significant, therefore assume homogeniety of variance
res.aov.logRich <- LPMacro_noEnvtot %>% anova_test(log10Richness ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.logRich
#reachxdays and leaftypexdays exposure significant, so group the data by reach and time, and fit  anova
LPMacro_noEnvtot %>% group_by(Reach) %>%
  anova_test(log10Richness ~ Days_Exposure, error = logRichmac.lm)
LPMacro_noEnvtot %>% group_by(Time_Point_cat) %>%
  anova_test(log10Richness ~ Reach, error = logRichmac.lm)
LPMacro_noEnvtot %>% group_by(Leaf_Type) %>%
  anova_test(log10Richness ~ Days_Exposure, error = logRichmac.lm)
LPMacro_noEnvtot %>% group_by(Time_Point_cat) %>%
  anova_test(log10Richness ~ Leaf_Type, error = logRichmac.lm)
# Pairwise comparisons
LPMacro_noEnvtot %>% group_by(Reach) %>% emmeans_test(log10Richness ~ Time_Point_cat, 
                                                      p.adjust.method = "bonferroni",
                                                      model=logRichmac.lm)
LPMacro_noEnvtot %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Richness ~ Reach, p.adjust.method = "bonferroni", 
               model=logRichmac.lm)
LPMacro_noEnvtot %>% group_by(Leaf_Type) %>%
  emmeans_test(log10Richness ~ Time_Point_cat, p.adjust.method = "bonferroni", 
               model=logRichmac.lm)
MRTPxLT<-LPMacro_noEnvtot %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Richness ~ Leaf_Type, p.adjust.method = "bonferroni",
               model=logRichmac.lm)
#emm for non interactions
LPMacro_noEnvtot %>% emmeans_test(log10Richness ~ Leaf_Type, p.adjust.method = "bonferroni",
                                  model = logRichmac.lm)
#ash greather than cotton and buckthorn, oak greather than buckthorn
LPMacro_noEnvtot %>% emmeans_test(log10Richness ~ Reach, p.adjust.method = "bonferroni",
                                  model = logRichmac.lm)
#Gap greater than us or ds
LPMacro_noEnvtot %>% emmeans_test(log10Richness ~ Time_Point_cat, p.adjust.method = "bonferroni",
                                  model = logRichmac.lm)
#See summary statistics for significant groups
LPMacro_noEnvtot_US<-subset(LPMacro_noEnvtot, Reach=="US")
LPMacro_noEnvtot_G<-subset(LPMacro_noEnvtot, Reach=="Gap")
LPMacro_noEnvtot_DS<-subset(LPMacro_noEnvtot, Reach=="DS")
stat.desc(LPMacro_noEnvtot_US$Richness)
stat.desc(LPMacro_noEnvtot_G$Richness)
stat.desc(LPMacro_noEnvtot_DS$Richness)
LPMacro_noEnvtot_A<-subset(LPMacro_noEnvtot, Leaf_Type=="Ash")
LPMacro_noEnvtot_B<-subset(LPMacro_noEnvtot, Leaf_Type=="Buckthorn")
LPMacro_noEnvtot_C<-subset(LPMacro_noEnvtot, Leaf_Type=="Cotton")
LPMacro_noEnvtot_O<-subset(LPMacro_noEnvtot, Leaf_Type=="Oak")
stat.desc(LPMacro_noEnvtot_A$Richness)
stat.desc(LPMacro_noEnvtot_B$Richness)
stat.desc(LPMacro_noEnvtot_C$Richness)
stat.desc(LPMacro_noEnvtot_O$Richness)

#visualize reach and leaf type
LPMacro_noEnvtot$Reach<-revalue(LPMacro_noEnvtot$Reach, c("US"="Upstream"))
LPMacro_noEnvtot$Reach<-revalue(LPMacro_noEnvtot$Reach, c("DS"="Downstream"))
ggplot(LPMacro_noEnvtot, aes(x=Reach, y=Richness, fill=Leaf_Type)) +
  geom_boxplot() +
  xlab("Gap Location") +
  scale_y_continuous(name="Genus Richness") +
  scale_fill_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))
#visualize reach, leaf type and time
LPMacroEnvtot$Richness<-rowSums(LPMacroComtot > 0)
MRsummary<-summarySE(LPMacroEnvtot, measurevar=c("Richness"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
MRsummary$Time_Point[MRsummary$Time_Point == 1] <- 8
MRsummary$Time_Point[MRsummary$Time_Point == 2] <- 41
MRsummary$Time_Point[MRsummary$Time_Point == 3] <- 68
MRsummary$Time_Point[MRsummary$Time_Point == 4] <- 98
MRsummary$Reach<-revalue(MRsummary$Reach, c("US"="Upstream"))
MRsummary$Reach<-revalue(MRsummary$Reach, c("DS"="Downstream"))

#plot on y richness and x axis days exposure
ggplot(MRsummary, aes(x=Time_Point, y=Richness, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Genus Richness (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#Calculate diversity
LPMacro_noEnvtot$Simp<-diversity(LPMacro_noComtot, index="invsimpson")
LPMacro_noEnvtot$Simp[!is.finite(LPMacro_noEnvtot$Simp)] <- 0
#Check outliers
LPMacro_noEnvtot %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Simp)
#no outliers
#Check model assumptions
Simpmac.lm<- lm(Simp+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot)
ggqqplot(residuals(Simpmac.lm))
shapiro_test(residuals(Simpmac.lm))
#qqplot has points in middle that fall out of zone, significant, not normal
#Transform simp using boxcox
boxcox(Simpmac.lm)
#labmda=0. log transform
LPMacro_noEnvtot$log10Simp<-log10(LPMacro_noEnvtot$Simp+1)
logSimpmac.lm<- lm(log10Simp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot)
ggqqplot(residuals(logSimpmac.lm))
shapiro_test(residuals(logSimpmac.lm))
#qqplot looks better, shaprio test is significant
LPMacro_noEnvtot %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(log10Simp)
#shaprio test significant
ggqqplot(LPMacro_noEnvtot, "log10Simp", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(log10Simp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot)
#not significant, therefore assume homogeniety of variance
LPMacro_noEnvtot %>% anova_test(log10Simp ~ Reach*Leaf_Type*Time_Point_cat)
anova(Simpmac.lm)
#Reach, leaf type, reach x time and leaf type x time significant
LPMacro_noEnvtot %>% emmeans_test(log10Simp ~ Reach, p.adjust.method = "bonferroni",
                              model = logSimpmac.lm)
#US and DS different than G
LPMacro_noEnvtot %>% emmeans_test(log10Simp ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = logSimpmac.lm)
#buckthorn and cotton less than ash
LPMacro_noEnvtot %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Simp ~ Reach, p.adjust.method = "bonferroni", 
               model=logSimpmac.lm)
LRSTPLT<-LPMacro_noEnvtot %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Simp ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=logSimpmac.lm)

#Descriptive stats
LPMacro_noEnvtot_US<-subset(LPMacro_noEnvtot, Reach=="Upstream")
LPMacro_noEnvtot_G<-subset(LPMacro_noEnvtot, Reach=="Gap")
LPMacro_noEnvtot_DS<-subset(LPMacro_noEnvtot, Reach=="Downstream")
stat.desc(LPMacro_noEnvtot_US$Simp)
stat.desc(LPMacro_noEnvtot_G$Simp)
stat.desc(LPMacro_noEnvtot_DS$Simp)
LPMacro_noEnvtot_A<-subset(LPMacro_noEnvtot, Leaf_Type=="Ash")
LPMacro_noEnvtot_B<-subset(LPMacro_noEnvtot, Leaf_Type=="Buckthorn")
LPMacro_noEnvtot_C<-subset(LPMacro_noEnvtot, Leaf_Type=="Cotton")
LPMacro_noEnvtot_O<-subset(LPMacro_noEnvtot, Leaf_Type=="Oak")
stat.desc(LPMacro_noEnvtot_A$Simp)
stat.desc(LPMacro_noEnvtot_B$Simp)
stat.desc(LPMacro_noEnvtot_C$Simp)
stat.desc(LPMacro_noEnvtot_O$Simp)

#visualize reach, leaf type and time
LPMacroEnvtot$Simp<-diversity(LPMacroComtot, index="invsimpson")
LPMacroEnvtot$Simp[!is.finite(LPMacroEnvtot$Simp)] <- 0
Divsummary<-summarySE(LPMacroEnvtot, measurevar=c("Simp"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
Divsummary$Time_Point[Divsummary$Time_Point == 1] <- 8
Divsummary$Time_Point[Divsummary$Time_Point == 2] <- 41
Divsummary$Time_Point[Divsummary$Time_Point == 3] <- 68
Divsummary$Time_Point[Divsummary$Time_Point == 4] <- 98
Divsummary$Reach<-revalue(Divsummary$Reach, c("US"="Upstream"))
Divsummary$Reach<-revalue(Divsummary$Reach, c("DS"="Downstream"))

#plot on y axis diversity and x axis days exposure
ggplot(Divsummary, aes(x=Time_Point, y=Simp, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=Simp-se, ymax=Simp+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Inverse Simpson's Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

AC_Inverts_Com_indic_leaf<-data.frame(signassoc(LPMacro_noComtot, cluster=Leaf_Taxatot_no,  mode=0, 
                                                alternative = "two.sided",control = how(nperm=999)))
AC_Inverts_Com_indic_leaf_sig<-subset(AC_Inverts_Com_indic_leaf, psidak<=0.05)
#Chironomidae (cotton), Gammaridae (ash), Nemouridae (ash), Simuliidae(Oak), Taeniopteridae (buckthorn)
AC_Inverts_Com_reach_indic<-data.frame(signassoc(LPMacro_noComtot, cluster=Reachtot_no,  mode=0, alternative = "two.sided",
                                                 control = how(nperm=999)))
AC_Inverts_Com_reach_indic_sig<-subset(AC_Inverts_Com_reach_indic, psidak<=0.05)
#Chironomid (US), ephemerellidae, Nemouridae, simuliidae Taeniopteridae- gap

#Run models for these 6 populations
#Because there are so many zeros, use dataset that deletes communities with no individuals
#Try using relative abundances
LPMacroComRAno<-data.frame(make_relative(as.matrix(LPMacroCom)))

#need to have at least 30 observations of presence to model with most limited dataset

LPMacrosM_no$Time_Point_cat<-as.factor(LPMacrosM_no$Time_Point)

#Chironomidae
##60 zeros out of 91
#Check outliers
LPMacrosM_no %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Chironomidae)
#no outliers
#Check model assumptions
Ch.lm<- lm(Chironomidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(Ch.lm))
shapiro_test(residuals(Ch.lm))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Ch.lm)
#lambda =-2
LPMacrosM_no$boxChi<-(LPMacrosM_no$Chironomidae+1)^-2
boxCh.lm<- lm(boxChi ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(boxCh.lm))
shapiro_test(residuals(boxCh.lm))
#shaprio test is significant
LPMacrosM_no %>% group_by(Reach) %>% shapiro_test(boxChi)
#shaprio test significant
ggqqplot(LPMacrosM_no, "boxChi", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxChi ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no %>% anova_test(boxChi ~ Reach*Leaf_Type*Time_Point_cat)
#reach, leaf type, time, reach x leaf type, reach x time significant
LPMacrosM_no %>% group_by(Reach) %>% anova_test(boxChi ~ Leaf_Type, error = boxCh.lm)
#gap significant
LPMacrosM_no %>% group_by(Leaf_Type) %>%
  emmeans_test(boxChi ~ Reach, p.adjust.method = "bonferroni", model=boxCh.lm)
#In gap, buckthorn and oak less than cotton, buckthorn less than ash
LPMacrosM_no %>% emmeans_test(boxChi ~ Reach, p.adjust.method = "bonferroni",
    model = boxCh.lm)
#gap different than us and ds
LPMacrosM_no %>% emmeans_test(boxChi ~ Leaf_Type, p.adjust.method = "bonferroni",
    model = boxCh.lm)
#Cotton greater than buckthorn 
LPMacrosM_no %>% emmeans_test(boxChi ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxCh.lm)
#1 < 2,3,4
LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxChi ~ Reach, p.adjust.method = "bonferroni", model=boxCh.lm)
#times 2 - 4- Gap > DS and US

LPMacrosM_no_US<-subset(LPMacrosM_no, Reach=="US")
LPMacrosM_no_G<-subset(LPMacrosM_no, Reach=="Gap")
LPMacrosM_no_DS<-subset(LPMacrosM_no, Reach=="DS")
stat.desc(LPMacrosM_no_US$Chironomidae)
stat.desc(LPMacrosM_no_G$Chironomidae)
stat.desc(LPMacrosM_no_DS$Chironomidae)
LPMacrosM_no_G_A<-subset(LPMacrosM_no_G, Leaf_Type=="Ash")
LPMacrosM_no_G_B<-subset(LPMacrosM_no_G, Leaf_Type=="Buckthorn")
LPMacrosM_no_G_C<-subset(LPMacrosM_no_G, Leaf_Type=="Cotton")
LPMacrosM_no_G_O<-subset(LPMacrosM_no_G, Leaf_Type=="Oak")
stat.desc(LPMacrosM_no_G_A$Chironomidae)
stat.desc(LPMacrosM_no_G_B$Chironomidae)
stat.desc(LPMacrosM_no_G_C$Chironomidae)
stat.desc(LPMacrosM_no_G_O$Chironomidae)

#Chironomids gap cotton
ChiSumm<-summarySE(LPMacrosM, measurevar=c("Chironomidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
ChiSumm$Time_Point[ChiSumm$Time_Point == 1] <- 8
ChiSumm$Time_Point[ChiSumm$Time_Point == 2] <- 41
ChiSumm$Time_Point[ChiSumm$Time_Point == 3] <- 68
ChiSumm$Time_Point[ChiSumm$Time_Point == 4] <- 98
ChiSumm$Reach<-revalue(ChiSumm$Reach, c("US"="Upstream"))
ChiSumm$Reach<-revalue(ChiSumm$Reach, c("DS"="Downstream"))
ggplot(ChiSumm, aes(x=Time_Point, y=Chironomidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Chironomidae-se, ymax=Chironomidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name="Mean Chironomidae Abundance (± SE)") +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#ephemerellidae
#78 zeros (out of 91)
#Check outliers
LPMacrosM_no %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Ephemerellidae)
#no outliers
#Check model assumptions
Ep.lm<- lm(Ephemerellidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(Ep.lm))
shapiro_test(residuals(Ep.lm))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Ep.lm)
#lambda =-2
LPMacrosM_no$boxEp<-(LPMacrosM_no$Ephemerellidae+1)^-2
boxEp.lm<- lm(boxEp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(boxEp.lm))
shapiro_test(residuals(boxEp.lm))
#shaprio test is significant
LPMacrosM_no %>% group_by(Leaf_Type) %>% shapiro_test(boxEp)
#shaprio test significant
ggqqplot(LPMacrosM_no, "boxEp", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxEp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no %>% anova_test(boxEp ~ Reach*Leaf_Type*Time_Point_cat)
#Reach, time and reach x time significant
LPMacrosM_no %>% emmeans_test(boxEp ~ Reach, p.adjust.method = "bonferroni",
                              model = boxEp.lm)
#Gap greather than US and DS
LPMacrosM_no %>% emmeans_test(boxEp ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxEp.lm)
#3 > 1
LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxEp ~ Reach, p.adjust.method = "bonferroni", model = boxEp.lm)
#time 3 and 4: Gap > us and ds

stat.desc(LPMacrosM_no_US$Ephemerellidae)
stat.desc(LPMacrosM_no_G$Ephemerellidae)
stat.desc(LPMacrosM_no_DS$Ephemerellidae)

#ephemerellidae gap 
EpSumm<-summarySE(LPMacrosM, measurevar=c("Ephemerellidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
EpSumm$Time_Point[EpSumm$Time_Point == 1] <- 8
EpSumm$Time_Point[EpSumm$Time_Point == 2] <- 41
EpSumm$Time_Point[EpSumm$Time_Point == 3] <- 68
EpSumm$Time_Point[EpSumm$Time_Point == 4] <- 98
EpSumm$Reach<-revalue(EpSumm$Reach, c("US"="Upstream"))
EpSumm$Reach<-revalue(EpSumm$Reach, c("DS"="Downstream"))
Eplab<-expression(paste("Mean ",italic("Ephemerella"), " Abundance (± SE)"))
ggplot(EpSumm, aes(x=Time_Point, y=Ephemerellidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Ephemerellidae-se, ymax=Ephemerellidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Eplab)) +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Gammaridae 81 zeros (out of 144)
#First absolute abundances
#Check outliers
LPMacrosM_no %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Gammaridae)
#no outliers
#Check model assumptions
Gam.lm<- lm(Gammaridae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(Gam.lm))
shapiro_test(residuals(Gam.lm))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Gam.lm)
#lambda =-2
LPMacrosM_no$boxGam<-(LPMacrosM_no$Gammaridae+1)^-2
boxGam.lm<- lm(boxGam ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(boxGam.lm))
shapiro_test(residuals(boxGam.lm))
#shaprio test is significant
LPMacrosM_no %>% group_by(Reach) %>% shapiro_test(boxGam)
#shaprio test significant
ggqqplot(LPMacrosM_no, "boxGam", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxGam ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no %>% anova_test(boxGam ~ Reach*Leaf_Type*Time_Point_cat)
#leaf, time, reach x time, leaf type x time and 3 way all significant
#just do overal leaf and reach
LPMacrosM_no %>% emmeans_test(boxGam ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxGam.lm)
#Ash greater than buckthorn and oak
LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxGam ~ Reach, p.adjust.method = "bonferroni", model=boxGam.lm)
#US>DS at time point 2, US<gap at time point 3
GTPLT<-LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxGam ~ Leaf_Type, p.adjust.method = "bonferroni", model=boxGam.lm)
#2 a > b,c,o
GTPLTR<-LPMacrosM_no %>% group_by(Time_Point_cat, Leaf_Type) %>%
  emmeans_test(boxGam ~ Reach, p.adjust.method = "bonferroni", model=boxGam.lm)
LPMacrosM_no %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(Gammaridae, type = "mean_se")


#Gammardiae indicates ash leaves
GamSumm<-summarySE(LPMacrosM, measurevar=c("Gammaridae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
GamSumm$Time_Point[GamSumm$Time_Point == 1] <- 8
GamSumm$Time_Point[GamSumm$Time_Point == 2] <- 41
GamSumm$Time_Point[GamSumm$Time_Point == 3] <- 68
GamSumm$Time_Point[GamSumm$Time_Point == 4] <- 98
GamSumm$Reach<-revalue(GamSumm$Reach, c("US"="Upstream"))
GamSumm$Reach<-revalue(GamSumm$Reach, c("DS"="Downstream"))
Gamlab<-expression(paste("Mean ",italic("Gammarus"), " Abundance (± SE)"))
ggplot(GamSumm, aes(x=Time_Point, y=Gammaridae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Gammaridae-se, ymax=Gammaridae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab(Gamlab) +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#nemouridae 64 zeros (out of 91)
#Check outliers
LPMacrosM_no %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Nemouridae)
#no outliers
#Check model assumptions
Ne.lm<- lm(Nemouridae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(Ne.lm))
shapiro_test(residuals(Ne.lm))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Ne.lm)
#lambda =-2
LPMacrosM_no$boxNe<-(LPMacrosM_no$Nemouridae+1)^-2
boxNe.lm<- lm(boxNe ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(boxNe.lm))
shapiro_test(residuals(boxNe.lm))
#shaprio test is significant
LPMacrosM_no %>% group_by(Reach) %>% shapiro_test(boxNe)
#shaprio test significant
ggqqplot(LPMacrosM_no, "boxNe", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxNe ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no %>% anova_test(boxNe ~ Reach*Leaf_Type*Time_Point_cat)
#Reach, leaftype, time, reach x leaftype, reach x time, leaftype x time significant
LPMacrosM_no %>% emmeans_test(boxNe ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxNe.lm)
#ash greather than buckthorn and cotton, oak greater than cotton
LPMacrosM_no %>% emmeans_test(boxNe ~ Reach, p.adjust.method = "bonferroni",
                              model = boxNe.lm)
#gap greater than us and ds
LPMacrosM_no %>% emmeans_test(boxNe ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxNe.lm)
#2 > 1,3
LPMacrosM_no %>% group_by(Leaf_Type) %>%
  emmeans_test(boxNe ~ Reach, p.adjust.method = "bonferroni", model = boxNe.lm)
LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxNe ~ Reach, p.adjust.method = "bonferroni", model = boxNe.lm)
NTPLT<-LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxNe ~ Leaf_Type, p.adjust.method = "bonferroni", model = boxNe.lm)
LPMacrosM_no %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(Nemouridae, type = "mean_se")
stat.desc(LPMacrosM_no_US$Nemouridae)
stat.desc(LPMacrosM_no_G$Nemouridae)
stat.desc(LPMacrosM_no_DS$Nemouridae)

#Nemouridae gap ash
NeSumm<-summarySE(LPMacrosM, measurevar=c("Nemouridae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
NeSumm$Time_Point[NeSumm$Time_Point == 1] <- 8
NeSumm$Time_Point[NeSumm$Time_Point == 2] <- 41
NeSumm$Time_Point[NeSumm$Time_Point == 3] <- 68
NeSumm$Time_Point[NeSumm$Time_Point == 4] <- 98
NeSumm$Reach<-revalue(NeSumm$Reach, c("US"="Upstream"))
NeSumm$Reach<-revalue(NeSumm$Reach, c("DS"="Downstream"))
Nelab<-expression(paste("Mean ",italic("Nemoura"), " Abundance (± SE)"))
ggplot(NeSumm, aes(x=Time_Point, y=Nemouridae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Nemouridae-se, ymax=Nemouridae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Nelab)) +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#simuliidae 73 zeros (out of 91)
#First absolute abundances
#Check outliers
LPMacrosM_no %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Simuliidae)
#no outliers
#Check model assumptions
Sim.lm<- lm(Simuliidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(Sim.lm))
shapiro_test(residuals(Sim.lm))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Sim.lm)
#lambda =-2
LPMacrosM_no$boxSim<-(LPMacrosM_no$Simuliidae+1)^-2
boxSim.lm<- lm(boxSim ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(boxSim.lm))
shapiro_test(residuals(boxSim.lm))
#shaprio test is significant
LPMacrosM_no %>% group_by(Leaf_Type) %>% shapiro_test(boxSim)
#shaprio test significant
ggqqplot(LPMacrosM_no, "boxSim", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxSim ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no %>% anova_test(boxSim ~ Reach*Leaf_Type*Time_Point_cat)
#reach, leaf type, time, reach x time, and 3 way interaction significant
LPMacrosM_no %>% emmeans_test(boxSim ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxSim.lm)
#oak greater than ash and buckthorn
LPMacrosM_no %>% emmeans_test(boxSim ~ Reach, p.adjust.method = "bonferroni",
                              model = boxSim.lm)
#gap greater than us and ds
LPMacrosM_no %>% emmeans_test(boxSim ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxSim.lm)
#2 > 1,4
LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxSim ~ Reach, p.adjust.method = "bonferroni", model = boxSim.lm)
STPLYR<-LPMacrosM_no %>% group_by(Time_Point_cat, Leaf_Type) %>%
  emmeans_test(boxSim ~ Reach, p.adjust.method = "bonferroni", model = boxSim.lm)

stat.desc(LPMacrosM_no_US$Simuliidae)
stat.desc(LPMacrosM_no_G$Simuliidae)
stat.desc(LPMacrosM_no_DS$Simuliidae)

#gap oak
SiSumm<-summarySE(LPMacrosM, measurevar=c("Simuliidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
SiSumm$Time_Point[SiSumm$Time_Point == 1] <- 8
SiSumm$Time_Point[SiSumm$Time_Point == 2] <- 41
SiSumm$Time_Point[SiSumm$Time_Point == 3] <- 68
SiSumm$Time_Point[SiSumm$Time_Point == 4] <- 98
SiSumm$Reach<-revalue(SiSumm$Reach, c("US"="Upstream"))
SiSumm$Reach<-revalue(SiSumm$Reach, c("DS"="Downstream"))
Silab<-expression(paste("Mean ",italic("Prosimulium"), " Abundance (± SE)"))
ggplot(SiSumm, aes(x=Time_Point, y=Simuliidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Simuliidae-se, ymax=Simuliidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Silab)) +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Taeniopteridae indicates buckthorn and gap
#First absolute abundances
#Check outliers
LPMacrosM_no %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Taeniopterygidae)
#no outliers
#Check model assumptions
Ta.lm<- lm(Taeniopterygidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(Ta.lm))
shapiro_test(residuals(Ta.lm))
#qqplot has many points that fall outside zone, significant, not normal 
boxcox(Ta.lm)
#lambda =-1
LPMacrosM_no$boxTa<-(LPMacrosM_no$Taeniopterygidae+1)^-1
boxTa.lm<- lm(boxTa ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
ggqqplot(residuals(boxTa.lm))
shapiro_test(residuals(boxTa.lm))
#shaprio test is significant
LPMacrosM_no %>% group_by(Leaf_Type) %>% shapiro_test(boxTa)
#shaprio test significant
ggqqplot(LPMacrosM_no, "boxTa", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxTa ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no %>% anova_test(boxTa ~ Reach*Leaf_Type*Time_Point_cat)
#reach, leaf type, time, reach x time and leaf type x time significant
LPMacrosM_no %>% emmeans_test(boxTa ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxTa.lm)
#ash greater than buckthorn and cotton, oak greater than buckthorn
LPMacrosM_no %>% emmeans_test(boxTa ~ Reach, p.adjust.method = "bonferroni",
                              model = boxTa.lm)
#gap greater than us and ds
LPMacrosM_no %>% emmeans_test(boxTa ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxTa.lm)
#(1,2) > 4, 2>3>4
LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxTa ~ Reach, p.adjust.method = "bonferroni", model = boxTa.lm)
# times 2 and 3 Gap > US and DS
TTPLT<-LPMacrosM_no %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxTa ~ Leaf_Type, p.adjust.method = "bonferroni", model = boxTa.lm)

stat.desc(LPMacrosM_no_US$Taeniopterygidae)
stat.desc(LPMacrosM_no_G$Taeniopterygidae)
stat.desc(LPMacrosM_no_DS$Taeniopterygidae)

#visualize taen gap
#absolute abundances
TaenSumm<-summarySE(LPMacrosM, measurevar=c("Taeniopterygidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
TaenSumm$Time_Point[TaenSumm$Time_Point == 1] <- 8
TaenSumm$Time_Point[TaenSumm$Time_Point == 2] <- 41
TaenSumm$Time_Point[TaenSumm$Time_Point == 3] <- 68
TaenSumm$Time_Point[TaenSumm$Time_Point == 4] <- 98
TaenSumm$Reach<-revalue(TaenSumm$Reach, c("US"="Upstream"))
TaenSumm$Reach<-revalue(TaenSumm$Reach, c("DS"="Downstream"))
Taenlab<-expression(paste("Mean ",italic("Taeniopteryx"), " Abundance (± SE)"))
ggplot(TaenSumm, aes(x=Time_Point, y=Taeniopterygidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Taeniopterygidae-se, ymax=Taeniopterygidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Taenlab)) +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Functional feeding groups
#Melt dataset
LPMacrosmelt<- melt(LPMacros, id=c("Pack_ID")) 
names(LPMacrosmelt)[names(LPMacrosmelt) == "variable"] <- "Taxa"
#Combine FFG with dataset
MacroFFGs<-read.csv("AC_Macro_FFGs.csv", sep = ",", header = T)
LeafPackMacrosMFFG<-merge(LPMacrosmelt, MacroFFGs, by="Taxa")
str(LeafPackMacrosMFFG)
LPMacrosFFG<-cast(LeafPackMacrosMFFG, Pack_ID~FFG, sum)
LeafPackFFG<-merge(LPMacrosFFG, LPMetadata, by="Pack_ID")
str(LeafPackFFG)

LeafPackFFG$rowsum<-rowSums(LeafPackFFG[,2:6])

LeafPackFFGno<-subset(LeafPackFFG, rowsum!=0)
LeafPackFFGCom<-LeafPackFFGno[,2:6]

LPFFGEnv<-LeafPackFFGno[,c(1,7:ncol(LeafPackFFGno))]
LPFFGComtot<-LeafPackFFG[,2:6]
LPFFGComtotbray0<-as.matrix(bray0(LPFFGComtot))
LPFFGComRA<-data.frame(make_relative(as.matrix(LeafPackFFG[,2:6])))
LPFFGComRA[is.nan(LPFFGComRA)] <- 0
LPFFGEnvtot<-LeafPackFFG[,c(1,7:ncol(LeafPackFFG))]
LPFFGRA<-merge(LPFFGComRA, LPFFGEnvtot, by=0)
FFGReachtot<-droplevels(as.factor(LPFFGEnvtot$Reach))

#delete time point 0 days
LeafPackFFG_no<-subset(LeafPackFFG, Time_Point!=0)

LeafPackFFG_noComtot<-LeafPackFFG_no[,2:6]
LeafPackFFG_noComtotbray0<-as.matrix(bray0(LeafPackFFG_noComtot))
LeafPackFFG_noComRA<-data.frame(make_relative(as.matrix(LeafPackFFG_no[,2:6])))
LeafPackFFG_noComRA[is.nan(LeafPackFFG_noComRA)] <- 0
LeafPackFFG_noEnvtot<-LeafPackFFG_no[,c(1,7:ncol(LeafPackFFG_no))]
LeafPackFFG_noRA<-merge(LeafPackFFG_noComRA, LeafPackFFG_noEnvtot, by=0)

#permanova
adonis2(as.dist(LeafPackFFG_noComtotbray0)~Reach*Leaf_Type*Days_Exposure, data=LeafPackFFG_noEnvtot, permutations=999)
#days, reach, leaf type signficant

FFGReachtot<-droplevels(as.factor(LPFFGEnvtot$Reach))
FFGLeaf_Taxatot<-as.factor(LPFFGEnvtot$Leaf_Type)

FFGReachtot_no<-droplevels(as.factor(LeafPackFFG_noEnvtot$Reach))
FFGLeaf_Taxatot_no<-as.factor(LeafPackFFG_noEnvtot$Leaf_Type)

#NMDS analysis
AC_FFG_NMDS<-metaMDS(as.dist(LPFFGComtotbray0))
#STRESS 0.11
stressplot(AC_FFG_NMDS)
ordiplot(AC_FFG_NMDS, type="n")
with(AC_FFG_NMDS, points(AC_FFG_NMDS, display="sites", col=leaftaxacolvec[FFGLeaf_Taxatot], pch=19))
with(AC_FFG_NMDS, legend("topleft", legend=levels(FFGLeaf_Taxatot), bty="n", col=leaftaxacolvec, pch=19, pt.bg=leaftaxacolvec))
with(AC_FFG_NMDS, ordiellipse(AC_FFG_NMDS, FFGLeaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_FFG_NMDS, ordiellipse(AC_FFG_NMDS, FFGLeaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_FFG_NMDS, ordiellipse(AC_FFG_NMDS, FFGLeaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Cotton"))
with(AC_FFG_NMDS, ordiellipse(AC_FFG_NMDS, FFGLeaf_Taxatot, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

FFGReachtot<-revalue(FFGReachtot, c("US"="Upstream"))
FFGReachtot<-revalue(FFGReachtot, c("DS"="Downstream"))
ordiplot(AC_FFG_NMDS, type="n")
with(AC_FFG_NMDS, points(AC_FFG_NMDS, display="sites", col=reach_col_vec[FFGReachtot], pch=19))
with(AC_FFG_NMDS, legend("topleft", legend=levels(FFGReachtot), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_FFG_NMDS, ordiellipse(AC_FFG_NMDS, Reachtot, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "Upstream"))
with(AC_FFG_NMDS, ordiellipse(AC_FFG_NMDS, Reachtot, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_FFG_NMDS, ordiellipse(AC_FFG_NMDS, Reachtot, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "Downstream"))

AC_FFG_Com_indic_leaf<-data.frame(signassoc(LeafPackFFG_noComtot, cluster=FFGLeaf_Taxatot_no,  mode=0, 
                                            alternative = "two.sided",control = how(nperm=999)))
AC_FFG_Com_indic_leaf_sig<-subset(AC_FFG_Com_indic_leaf, psidak<=0.05)
#shredder ash
AC_FFG_Com_reach_indic<-data.frame(signassoc(LeafPackFFG_noComtot, cluster=FFGReachtot_no,  mode=0, 
                                             alternative = "two.sided",control = how(nperm=999)))
AC_FFG_Com_reach_indic_sig<-subset(AC_FFG_Com_reach_indic, psidak<=0.05)
#collector filt, shredder gap, gatherer us

#model collector filt, collector gath, shredder because of indicator species analy
#also grazer because addressed in hypotheses

#Collector filterer
LeafPackFFG_noRA$Time_Point_cat<-as.factor(LeafPackFFG_noRA$Time_Point)
LeafPackFFGnoRA<-subset(LeafPackFFG_noRA, rowsum>0)
#summary stats
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(CollectorFilterer, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(CollectorFilterer)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA, x = "Reach", y = "CollectorFilterer",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA %>%
  group_by(Reach) %>%
  shapiro_test(CollectorFilterer)
#not normal
ggqqplot(LeafPackFFGnoRA, "CollectorFilterer", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
CF.lm<- lm((CollectorFilterer+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA)
ggqqplot(residuals(CF.lm))
shapiro_test(residuals(CF.lm))
#0.0000000423 not normal, see boxcox
boxcox(CF.lm)
#labmda=-2. transform
LeafPackFFGnoRA$boxCF<-(LeafPackFFGnoRA$CollectorFilterer+1)^-2
#check normality on transformed CF
LeafPackFFGnoRA %>%
  group_by(Reach) %>%
  shapiro_test(boxCF)
#not normal
LeafPackFFGnoRA %>% group_by(Leaf_Type) %>% shapiro_test(boxCF)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA, "boxCF", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxCF ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA %>% anova_test(boxCF ~ Reach*Leaf_Type*Time_Point_cat)
#Nothing significant

#Visualize
CFSumm<-summarySE(LPFFGRA, measurevar=c("CollectorFilterer"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
CFSumm$Time_Point[CFSumm$Time_Point == 1] <- 8
CFSumm$Time_Point[CFSumm$Time_Point == 2] <- 41
CFSumm$Time_Point[CFSumm$Time_Point == 3] <- 68
CFSumm$Time_Point[CFSumm$Time_Point == 4] <- 98
CFSumm$Reach<-revalue(CFSumm$Reach, c("US"="Upstream"))
CFSumm$Reach<-revalue(CFSumm$Reach, c("DS"="Downstream"))
ggplot(CFSumm, aes(x=Time_Point, y=CollectorFilterer, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=CollectorFilterer-se, ymax=CollectorFilterer+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name="Collector-Filterer Relative Abundance (± SEM)") +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Collector gatherer
#summary stats
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(CollectorGatherer, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(CollectorGatherer)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA, x = "Reach", y = "CollectorGatherer",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA %>%
  group_by(Reach) %>%
  shapiro_test(CollectorGatherer)
#not normal
ggqqplot(LeafPackFFGnoRA, "CollectorGatherer", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
CG.lm<- lm((CollectorGatherer+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA)
ggqqplot(residuals(CG.lm))
shapiro_test(residuals(CG.lm))
#0.00000000237 not normal, see boxcox
boxcox(CG.lm)
#labmda=-2. transform
LeafPackFFGnoRA$boxCG<-(LeafPackFFGnoRA$CollectorGatherer+1)^-2
#check normality on transformed CF
LeafPackFFGnoRA %>%
  group_by(Reach) %>%
  shapiro_test(boxCG)
#not normal
LeafPackFFGnoRA %>% group_by(Leaf_Type) %>% shapiro_test(boxCG)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA, "boxCG", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxCG ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA %>% anova_test(boxCG ~ Reach*Leaf_Type*Time_Point_cat)
#Leaf type, time and 3 way interaction significant
boxCG.lm<-lm(boxCG ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA)
LeafPackFFGnoRA %>%
  group_by(Time_Point_cat) %>%
  anova_test(boxCG ~ Reach*Leaf_Type, error = boxCG.lm)
#leaf type significant on time 2, 3, and 4
#Reach x leaf type significant on time 2, 3 and 4
LeafPackFFGnoRA %>% group_by(Time_Point_cat) %>%
  anova_test(boxCG ~ Reach, error = boxCG.lm)
#Time point 3, reach significant
threewayboxCG<-LeafPackFFGnoRA  %>%
  group_by(Time_Point_cat, Leaf_Type) %>%
  emmeans_test(boxCG ~ Reach, p.adjust.method = "bonferroni", model=boxCG.lm)
#time point 2 buckthorn, gap less than ds
#time point 4 cotton, US and DS less than gap
LeafPackFFGnoRA %>% emmeans_test(boxCG ~ Leaf_Type, p.adjust.method = "bonferroni")
#cotton greather than ash, buckthorn and oak
LeafPackFFGnoRA %>% emmeans_test(boxCG ~ Time_Point_cat, p.adjust.method = "bonferroni")
#3 > 1

threewayCG<-LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(CollectorGatherer, type = "mean_se")
LeafPackFFGnoRA %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(CollectorGatherer, type = "mean_se")

#Visualize
CGSumm<-summarySE(LPFFGRA, measurevar=c("CollectorGatherer"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
CGSumm$Time_Point[CGSumm$Time_Point == 1] <- 8
CGSumm$Time_Point[CGSumm$Time_Point == 2] <- 41
CGSumm$Time_Point[CGSumm$Time_Point == 3] <- 68
CGSumm$Time_Point[CGSumm$Time_Point == 4] <- 98
CGSumm$Reach<-revalue(CGSumm$Reach, c("US"="Upstream"))
CGSumm$Reach<-revalue(CGSumm$Reach, c("DS"="Downstream"))
ggplot(CGSumm, aes(x=Time_Point, y=CollectorGatherer, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=CollectorGatherer-se, ymax=CollectorGatherer+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name="Mean Collector-Gatherer Relative Abundance (± SE)") +
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Shredder
#summary stats
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(Shredder, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Shredder)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA, x = "Reach", y = "Shredder",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA %>% group_by(Reach) %>% shapiro_test(Shredder)
#not normal
ggqqplot(LeafPackFFGnoRA, "Shredder", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
Sh.lm<- lm((Shredder+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA)
ggqqplot(residuals(Sh.lm))
shapiro_test(residuals(Sh.lm))
#0.0129 not normal, see boxcox
boxcox(Sh.lm)
#labmda=0.5. transform
LeafPackFFGnoRA$boxSh<-(LeafPackFFGnoRA$Shredder+1)^0.5
#check normality on transformed CF
LeafPackFFGnoRA %>% group_by(Reach) %>% shapiro_test(boxSh)
#not normal
LeafPackFFGnoRA %>% group_by(Leaf_Type) %>% shapiro_test(boxSh)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA, "boxSh", ggtheme = theme_bw()) +
  facet_grid(Reach  ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxSh ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA %>% anova_test(boxSh ~ Reach*Leaf_Type*Time_Point_cat)
#time significant
LeafPackFFGnoRA %>% emmeans_test(boxSh ~ Time_Point_cat, p.adjust.method = "bonferroni")

#Grazer
#summary stats
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(Grazer, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Grazer)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA, x = "Reach", y = "Grazer",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA %>% group_by(Reach) %>% shapiro_test(Grazer)
#not normal
ggqqplot(LeafPackFFGnoRA, "Grazer", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
Gr.lm<- lm((Grazer+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA)
ggqqplot(residuals(Gr.lm))
shapiro_test(residuals(Gr.lm))
#4.07e-11 not normal, see boxcox
boxcox(Gr.lm)
#labmda=-2. transform
LeafPackFFGnoRA$boxGr<-(LeafPackFFGnoRA$Grazer+1)^-2
#check normality on transformed CF
LeafPackFFGnoRA %>% group_by(Reach) %>% shapiro_test(boxGr)
#not normal
LeafPackFFGnoRA %>% group_by(Leaf_Type) %>% shapiro_test(boxGr)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA, "boxGr", ggtheme = theme_bw()) +
  facet_grid(Reach  ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxGr ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA %>% anova_test(boxGr ~ Reach*Leaf_Type*Time_Point_cat)
#Nothing significant

#Make stacked bar graphs
LPFFGmelt<-merge(LeafPackMacrosMFFG, LPMetadata, by="Pack_ID")
str(LPFFGmelt)
LPFFGmelt$Reach<-revalue(LPFFGmelt$Reach, c("US"="Upstream"))
LPFFGmelt$Reach<-revalue(LPFFGmelt$Reach, c("DS"="Downstream"))
LPFFGmelt$FFG<-revalue(LPFFGmelt$FFG, c("CollectorFilterer"="Collector-Filterer"))
LPFFGmelt$FFG<-revalue(LPFFGmelt$FFG, c("CollectorGatherer"="Collector-Gatherer"))
ggplot(LPFFGmelt, aes(fill=FFG, y=value, x=Leaf_Type)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=ffgcolvec,name = "Functional\nFeeding Group") +
  facet_grid( ~ Reach) +
  xlab("Leaf Type") +
  scale_y_continuous(name="Relative Abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))
ggplot(LPFFGmelt, aes(fill=FFG, y=value, x=Reach)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=ffgcolvec,name = "Functional\nFeeding Group") +
  facet_grid(Leaf_Type ~.) +
  scale_y_continuous(name="Relative Abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.y = element_text(size = 16))

#########################################################
#Now run without cotton
###################################################

#subset cotton out of leaf litter decay dataset
LeafPackMi_nc<-subset(LeafPackMi, Leaf_Type!="Cotton")

afdm.lm.ni.nc<-lm(lnpercAFDMremain~Days_Exposure+Days_Exposure:Leaf_Type+Days_Exposure:Reach+
                 Days_Exposure:Reach:Leaf_Type, data=LeafPackMi_nc)
summary(afdm.lm.ni.nc)
#intercept, time, time:buckthorn and time:cotton significant
#output to excel file
#use variance covariance matrix to calculate SE for -k estimates
vcovafdm.nc<-vcov(afdm.lm.ni.nc)
write.csv(vcovafdm.nc,'vcovafdm.nc.csv')
#create k value and ses .csv based on calculations
#or use this code
anova(afdm.lm.ni.nc)
#only leaf type significant
afdm.lst.ni.nc<- lstrends(afdm.lm.ni.nc, "Leaf_Type", var="Days_Exposure")
pairs(afdm.lst.ni.nc)

#Calculate mean and sd for each leaf type and sampling day
AFDMsummary_nc<-summarySE(LeafPackMi, measurevar=c("percAFDMremain"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
AFDMsummary_nc<-subset(AFDMsummary, percAFDMremain!="NaN")
AFDMsummary_nc$Time_Point[AFDMsummary$Time_Point == 1] <- 8
AFDMsummary_nc$Time_Point[AFDMsummary$Time_Point == 2] <- 41
AFDMsummary_nc$Time_Point[AFDMsummary$Time_Point == 3] <- 68
AFDMsummary_nc$Time_Point[AFDMsummary$Time_Point == 4] <- 98
AFDMsummary_nc$Reach<-revalue(AFDMsummary$Reach, c("US"="Upstream"))
AFDMsummary_nc$Reach<-revalue(AFDMsummary$Reach, c("DS"="Downstream"))

AFDMsummary_nc<-subset(AFDMsummary, Leaf_Type!="Cotton")

#plot on y axis % remaining and x axis days exposure
ggplot(AFDMsummary_nc, aes(x=Time_Point, y=percAFDMremain+1, color=Leaf_Type)) +
  geom_smooth(aes(group=Leaf_Type), method=lm, se=F, fullrange=TRUE, size=2)+
  geom_errorbar(aes(ymin=percAFDMremain+1-se, ymax=percAFDMremain+1+se), width=2) +
  geom_point(size=3) +
  xlab("Days of Exposure") +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  scale_y_log10(name="%AFDM Remaining +1 (± SEM)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)
ggplot(AFDMsummary_nc, aes(x=Time_Point, y=percAFDMremain, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=percAFDMremain-se, ymax=percAFDMremain+se), width=7) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean %AFDM Remaining (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=18),
        legend.title=element_text(size=24),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 22)) +
  facet_wrap(.~Reach)
ggplot(AFDMsummary_nc, aes(x=Time_Point, y=percAFDMremain, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=percAFDMremain-se, ymax=percAFDMremain+se), width=7) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean %AFDM Remaining (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=18),
        legend.title=element_text(size=24),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 22),legend.position="bottom") +
  facet_wrap(.~Reach)

#Upload kvalue means and se's 
KMeans<-read.csv("kvalues.csv", sep = ",", header = T)
KMeans$negk<--1*(KMeans$k)
KMeans$Reach<-factor(KMeans$Reach, levels=c("Upstream","Gap","Downstream"))
KMeans_nc<-subset(KMeans, Leaf_Type!="Cotton")
ggplot(KMeans_nc, aes(x=Reach, y=negk, fill=Leaf_Type)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=negk-SE, ymax=negk+SE), width=.2,
                position=position_dodge(.9))+
  xlab("Gap Location") +
  ylab("Decomposition Rate (-k ± SE)") +
  scale_fill_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=16),axis.text.y = element_text(size=16),
        legend.title=element_text(size=20),legend.text = element_text(size=20))

LPMacrosM_nc<-subset(LPMacrosM, Leaf_Type!="Cotton")


LPMacrosMno_nc<-subset(LPMacrosM_nc, rowsum!=0)
LPMacroCom_nc<-LPMacrosMno_nc[,2:22]
sum(LPMacroCom_nc)
#548 individuals
colSums(LPMacroCom_nc)
#Taeniopterygidae most abundant with 147 individuals

LPMacroEnv_nc<-LPMacrosMno_nc[,c(1,23:ncol(LPMacrosMno_nc))]
LPMacroComtot_nc<-LPMacrosM_nc[,2:22]
LPMacroComtotbray0_nc<-as.matrix(bray0(LPMacroComtot_nc))
LPMacroComRA_nc<-data.frame(make_relative(as.matrix(LPMacrosM_nc[,2:22])))
LPMacroComRA_nc[is.nan(LPMacroComRA_nc)] <- 0
LPMacroEnvtot_nc<-LPMacrosM_nc[,c(1,23:ncol(LPMacrosM_nc))]
LPMacroRA_nc<-merge(LPMacroComRA_nc, LPMacroEnvtot_nc, by=0)
Reachtot_nc<-droplevels(as.factor(LPMacroEnvtot_nc$Reach))

#delete time point 0 days
LPMacrosM_no_nc<-subset(LPMacrosM_nc, Time_Point!=0)

LPMacro_noComtot_nc<-LPMacrosM_no_nc[,2:22]
LPMacro_noComtotbray0_nc<-as.matrix(bray0(LPMacro_noComtot_nc))
LPMacro_noComRA_nc<-data.frame(make_relative(as.matrix(LPMacrosM_no_nc[,2:22])))
LPMacro_noComRA_nc[is.nan(LPMacro_noComRA_nc)] <- 0
LPMacro_noEnvtot_nc<-LPMacrosM_no_nc[,c(1,23:ncol(LPMacrosM_no_nc))]
LPMacro_noRA_nc<-merge(LPMacro_noComRA_nc, LPMacro_noEnvtot_nc, by=0)

stat.desc(LPMacro_noRA_nc$Taeniopterygidae)
#0.21 +/- 0.03

#permanova
adonis2(as.dist(LPMacro_noComtotbray0_nc)~Reach*Leaf_Type*Days_Exposure, data=LPMacro_noEnvtot_nc, permutations=999)
#days, reach, and leaf type signficant

Reachtot_nc<-droplevels(as.factor(LPMacroEnvtot_nc$Reach))
Leaf_Taxatot_nc<-factor(LPMacroEnvtot_nc$Leaf_Type)
ReachLeaftot_nc<-as.factor(paste(LPMacro_noEnvtot_nc$Reach,LPMacroEnvtot_nc$Leaf_Type))
ReachDaystot_nc<-as.factor(paste(LPMacro_noEnvtot_nc$Reach,LPMacroEnvtot_nc$Time_Point))

Reachtot_no_nc<-droplevels(as.factor(LPMacro_noEnvtot_nc$Reach))
Leaf_Taxatot_no_nc<-droplevels(as.factor(LPMacro_noEnvtot_nc$Leaf_Type))
ReachLeaftot_no_nc<-factor(paste(LPMacro_noEnvtot_nc$Reach,LPMacro_noEnvtot_nc$Leaf_Type))
ReachDaystot_no_nc<-factor(paste(LPMacro_noEnvtot_nc$Reach,LPMacro_noEnvtot_nc$Time_Point))

#NMDS analysis
AC_Macroinvertebrate_NMDS_nc<-metaMDS(as.dist(LPMacroComtotbray0_nc))
stressplot(AC_Macroinvertebrate_NMDS_nc)
ordiplot(AC_Macroinvertebrate_NMDS_nc, type="n")
with(AC_Macroinvertebrate_NMDS_nc, points(AC_Macroinvertebrate_NMDS_nc, display="sites", col=leaftaxacolvec_nc[Leaf_Taxatot_nc], pch=19, pt.bg=leaftaxacolvec_nc))
with(AC_Macroinvertebrate_NMDS_nc, legend("topleft", legend=levels(Leaf_Taxatot_nc), bty="n", col=leaftaxacolvec_nc, pch=19, pt.bg=leaftaxacolvec_nc))
with(AC_Macroinvertebrate_NMDS_nc, ordiellipse(AC_Macroinvertebrate_NMDS_nc, Leaf_Taxatot_nc, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_Macroinvertebrate_NMDS_nc, ordiellipse(AC_Macroinvertebrate_NMDS_nc, Leaf_Taxatot_nc, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_Macroinvertebrate_NMDS_nc, ordiellipse(AC_Macroinvertebrate_NMDS_nc, Leaf_Taxatot_nc, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

Reachtot_nc<-revalue(Reachtot_nc, c("US"="Upstream"))
Reachtot_nc<-revalue(Reachtot_nc, c("DS"="Downstream"))
ordiplot(AC_Macroinvertebrate_NMDS_nc, type="n")
with(AC_Macroinvertebrate_NMDS_nc, points(AC_Macroinvertebrate_NMDS_nc, display="sites", col=reach_col_vec[Reachtot_nc], pch=19))
with(AC_Macroinvertebrate_NMDS_nc, legend("topleft", legend=levels(Reachtot_nc), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_Macroinvertebrate_NMDS_nc, ordiellipse(AC_Macroinvertebrate_NMDS_nc, Reachtot_nc, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "Upstream"))
with(AC_Macroinvertebrate_NMDS_nc, ordiellipse(AC_Macroinvertebrate_NMDS_nc, Reachtot_nc, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_Macroinvertebrate_NMDS_nc, ordiellipse(AC_Macroinvertebrate_NMDS_nc, Reachtot_nc, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "Downstream"))
#############

#Diversity metrics
LPMacro_noEnvtot_nc$Leaf_Type<-factor(LPMacro_noEnvtot_nc$Leaf_Type)
#Calculate richess
LPMacro_noEnvtot_nc$Richness<-rowSums(LPMacro_noComtot_nc > 0)
LPMacro_noEnvtot_nc$Time_Point_cat<-as.factor(LPMacro_noEnvtot_nc$Time_Point)
#summary stats
LPMacro_noEnvtot_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(Richness, type = "mean_sd")
#check for outliers
LPMacro_noEnvtot_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Richness)
#no outliers
#visualize
ggboxplot(LPMacro_noEnvtot_nc, x = "Reach", y = "Richness",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LPMacro_noEnvtot_nc %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(Richness)
#not normal
ggqqplot(LPMacro_noEnvtot_nc, "Richness", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#Transform Richness using boxcox
MR.lm_nc<- lm((Richness+1) ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot_nc)
ggqqplot(residuals(MR.lm_nc))
shapiro_test(residuals(MR.lm_nc))
#0.000000556 not normal, see boxcox
boxcox(MR.lm_nc)
#labmda=0. log transform
LPMacro_noEnvtot_nc$log10Richness<-log10(LPMacro_noEnvtot_nc$Richness+1)
LPMacro_noEnvtot_nc %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(log10Richness)
#not normal
ggqqplot(LPMacro_noEnvtot_nc, "log10Richness", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#Use this transformation
logRichmac.lm_nc<- lm(log10Richness ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot_nc)
ggqqplot(residuals(logRichmac.lm_nc))
shapiro_test(residuals(logRichmac.lm_nc))
#residuals fall aprox in qqplot, but shapiro test still <0.05
#Note that, if your sample size is greater than 50, the normal QQ plot is preferred 
#because at larger sample sizes the Shapiro-Wilk test becomes very sensitive even to a
#minor deviation from normality.
#looks good, test homogeniety of variance
leveneTest(log10Richness ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot_nc)
#not significant, therefore assume homogeniety of variance
res.aov.logRich_nc<- LPMacro_noEnvtot_nc %>% anova_test(log10Richness ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.logRich_nc
#reach, leaf type, time point, reachxdays and leaftypexdays exposure significant, so group the data by reach and time, and fit  anova
LPMacro_noEnvtot_nc %>% group_by(Reach) %>%
  anova_test(log10Richness ~ Days_Exposure, error = logRichmac.lm_nc)
LPMacro_noEnvtot_nc %>% group_by(Time_Point_cat) %>%
  anova_test(log10Richness ~ Reach, error = logRichmac.lm_nc)
LPMacro_noEnvtot_nc %>% group_by(Leaf_Type) %>%
  anova_test(log10Richness ~ Days_Exposure, error = logRichmac.lm_nc)
LPMacro_noEnvtot_nc %>% group_by(Time_Point_cat) %>%
  anova_test(log10Richness ~ Leaf_Type, error = logRichmac.lm_nc)
# Pairwise comparisons
LPMacro_noEnvtot_nc %>% group_by(Reach) %>% emmeans_test(log10Richness ~ Time_Point_cat, 
                                                      p.adjust.method = "bonferroni",
                                                      model=logRichmac.lm_nc)
LPMacro_noEnvtot_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Richness ~ Reach, p.adjust.method = "bonferroni", 
               model=logRichmac.lm_nc)
LPMacro_noEnvtot_nc %>% group_by(Leaf_Type) %>%
  emmeans_test(log10Richness ~ Time_Point_cat, p.adjust.method = "bonferroni", 
               model=logRichmac.lm_nc)
MRTPxLT_nc<-LPMacro_noEnvtot_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Richness ~ Leaf_Type, p.adjust.method = "bonferroni",
               model=logRichmac.lm_nc)
#emm for non interactions
LPMacro_noEnvtot_nc %>% emmeans_test(log10Richness ~ Leaf_Type, p.adjust.method = "bonferroni",
                                  model = logRichmac.lm_nc)
#ash greather oak greater than buckthorn
LPMacro_noEnvtot_nc %>% emmeans_test(log10Richness ~ Reach, p.adjust.method = "bonferroni",
                                  model = logRichmac.lm_nc)
#Gap greater than us or ds
LPMacro_noEnvtot_nc %>% emmeans_test(log10Richness ~ Time_Point_cat, p.adjust.method = "bonferroni",
                                  model = logRichmac.lm_nc)
#1&3<2,
#See summary statistics for significant groups
LPMacro_noEnvtot_nc %>%
  group_by(Reach) %>%
  get_summary_stats(log10Richness, type = "mean_se")
LPMacro_noEnvtot_nc %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(log10Richness, type = "mean_se")

#visualize reach and leaf type
LPMacro_noEnvtot_nc$Reach<-revalue(LPMacro_noEnvtot_nc$Reach, c("US"="Upstream"))
LPMacro_noEnvtot_nc$Reach<-revalue(LPMacro_noEnvtot_nc$Reach, c("DS"="Downstream"))
ggplot(LPMacro_noEnvtot_nc, aes(x=Reach, y=Richness, fill=Leaf_Type)) +
  geom_boxplot() +
  xlab("Gap Location") +
  scale_y_continuous(name="Genus Richness") +
  scale_fill_manual(values=leaftaxacolvec_nc,name = "Leaf Taxa") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))
#visualize reach, leaf type and time
LPMacroEnvtot_nc$Richness<-rowSums(LPMacroComtot_nc > 0)
MRsummary_nc<-summarySE(LPMacroEnvtot_nc, measurevar=c("Richness"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
MRsummary_nc$Time_Point[MRsummary_nc$Time_Point == 1] <- 8
MRsummary_nc$Time_Point[MRsummary_nc$Time_Point == 2] <- 41
MRsummary_nc$Time_Point[MRsummary_nc$Time_Point == 3] <- 68
MRsummary_nc$Time_Point[MRsummary_nc$Time_Point == 4] <- 98
MRsummary_nc$Reach<-revalue(MRsummary_nc$Reach, c("US"="Upstream"))
MRsummary_nc$Reach<-revalue(MRsummary_nc$Reach, c("DS"="Downstream"))

#plot on y richness and x axis days exposure
ggplot(MRsummary_nc, aes(x=Time_Point, y=Richness, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Macroinvertebrate Genus Richness (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=24),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=18),
        legend.title=element_text(size=24),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 22),legend.position = "bottom") +
  facet_wrap(.~Reach)

#Calculate diversity
LPMacro_noEnvtot_nc$Simp<-diversity(LPMacro_noComtot_nc, index="invsimpson")
LPMacro_noEnvtot_nc$Simp[!is.finite(LPMacro_noEnvtot_nc$Simp)] <- 0
#Check outliers
LPMacro_noEnvtot_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Simp)
#no outliers
#Check model assumptions
Simpmac.lm_nc<- lm(Simp+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot_nc)
ggqqplot(residuals(Simpmac.lm_nc))
shapiro_test(residuals(Simpmac.lm_nc))
#qqplot has points in middle that fall out of zone, significant, not normal
#Transform simp using boxcox
boxcox(Simpmac.lm_nc)
#labmda=0. log transform
LPMacro_noEnvtot_nc$log10Simp<-log10(LPMacro_noEnvtot_nc$Simp+1)
logSimpmac.lm_nc<- lm(log10Simp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot_nc)
ggqqplot(residuals(logSimpmac.lm_nc))
shapiro_test(residuals(logSimpmac.lm_nc))
#qqplot looks better, shaprio test is significant
LPMacro_noEnvtot_nc %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(log10Simp)
#shaprio test significant
ggqqplot(LPMacro_noEnvtot_nc, "log10Simp", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(log10Simp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacro_noEnvtot_nc)
#not significant, therefore assume homogeniety of variance
LPMacro_noEnvtot_nc %>% anova_test(log10Simp ~ Reach*Leaf_Type*Time_Point_cat)
anova(Simpmac.lm_nc)
#Reach, leaf type, reach x time and leaf type x time significant
LPMacro_noEnvtot_nc %>% emmeans_test(log10Simp ~ Reach, p.adjust.method = "bonferroni",
                                  model = logSimpmac.lm_nc)
#US and DS different than G
LPMacro_noEnvtot_nc %>% emmeans_test(log10Simp ~ Leaf_Type, p.adjust.method = "bonferroni",
                                  model = logSimpmac.lm_nc)
#ash > oak > buckthorn
LPMacro_noEnvtot_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Simp ~ Reach, p.adjust.method = "bonferroni", 
               model=logSimpmac.lm_nc)
LRSTPLT_nc<-LPMacro_noEnvtot_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(log10Simp ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=logSimpmac.lm_nc)

#Descriptive stats
LPMacro_noEnvtot_nc %>%
  group_by(Reach) %>%
  get_summary_stats(log10Simp, type = "mean_se")
LPMacro_noEnvtot_nc %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(log10Simp, type = "mean_se")

#visualize reach, leaf type and time
LPMacroEnvtot_nc$Simp<-diversity(LPMacroComtot_nc, index="invsimpson")
LPMacroEnvtot_nc$Simp[!is.finite(LPMacroEnvtot_nc$Simp)] <- 0
Divsummary_nc<-summarySE(LPMacroEnvtot_nc, measurevar=c("Simp"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
Divsummary_nc$Time_Point[Divsummary_nc$Time_Point == 1] <- 8
Divsummary_nc$Time_Point[Divsummary_nc$Time_Point == 2] <- 41
Divsummary_nc$Time_Point[Divsummary_nc$Time_Point == 3] <- 68
Divsummary_nc$Time_Point[Divsummary_nc$Time_Point == 4] <- 98
Divsummary_nc$Reach<-revalue(Divsummary_nc$Reach, c("US"="Upstream"))
Divsummary_nc$Reach<-revalue(Divsummary_nc$Reach, c("DS"="Downstream"))

#plot on y axis diversity and x axis days exposure
ggplot(Divsummary_nc, aes(x=Time_Point, y=Simp, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=Simp-se, ymax=Simp+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Inverse Simpson's Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

AC_Inverts_Com_indic_leaf_nc<-data.frame(signassoc(LPMacro_noComtot_nc, cluster=Leaf_Taxatot_no_nc,  
                                                   mode=0, alternative = "two.sided",control = how(nperm=999)))
AC_Inverts_Com_indic_leaf_sig_nc<-subset(AC_Inverts_Com_indic_leaf_nc, psidak<=0.05)
#Gammaridae (ash), Nemouridae (ash), Taeniopteridae (buckthorn)
AC_Inverts_Com_reach_indic_nc<-data.frame(signassoc(LPMacro_noComtot_nc, cluster=Reachtot_no_nc,  mode=0, 
                                                    alternative = "two.sided",control = how(nperm=999)))
AC_Inverts_Com_reach_indic_sig_nc<-subset(AC_Inverts_Com_reach_indic_nc, psidak<=0.05)
#Chironomid (US), ephemerellidae, Nemouridae, simuliidae Taeniopteridae- gap

#Run models for these 6 populations
#Because there are so many zeros, use dataset that deletes communities with no individuals
#Try using relative abundances
LPMacroComRAno_nc<-data.frame(make_relative(as.matrix(LPMacroCom_nc)))

#need to have at least 30 observations of presence to model with most limited dataset

LPMacrosM_no_nc$Time_Point_cat<-as.factor(LPMacrosM_no_nc$Time_Point)

#Chironomidae
##60 zeros out of 91
#Check outliers
LPMacrosM_no_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Chironomidae)
#no outliers
#Check model assumptions
Ch.lm_nc<- lm(Chironomidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(Ch.lm_nc))
shapiro_test(residuals(Ch.lm_nc))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Ch.lm_nc)
#lambda =-2
LPMacrosM_no_nc$boxChi<-(LPMacrosM_no_nc$Chironomidae+1)^-2
boxCh.lm_nc<- lm(boxChi ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(boxCh.lm_nc))
shapiro_test(residuals(boxCh.lm_nc))
#shaprio test is significant
LPMacrosM_no_nc %>% group_by(Reach) %>% shapiro_test(boxChi)
#shaprio test significant
ggqqplot(LPMacrosM_no_nc, "boxChi", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxChi ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no_nc %>% anova_test(boxChi ~ Reach*Leaf_Type*Time_Point_cat)
#reach significant
LPMacrosM_no_nc %>% emmeans_test(boxChi ~ Reach, p.adjust.method = "bonferroni",
                              model = boxCh.lm_nc)
#gap different than us and ds
LPMacrosM_no_nc %>%
  group_by(Reach) %>%
  get_summary_stats(Chironomidae, type = "mean_se")
#2 times higher in gap

#Chironomids gap
ChiSumm_nc<-summarySE(LPMacrosM_nc, measurevar=c("Chironomidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
ChiSumm_nc$Time_Point[ChiSumm_nc$Time_Point == 1] <- 8
ChiSumm_nc$Time_Point[ChiSumm_nc$Time_Point == 2] <- 41
ChiSumm_nc$Time_Point[ChiSumm_nc$Time_Point == 3] <- 68
ChiSumm_nc$Time_Point[ChiSumm_nc$Time_Point == 4] <- 98
ChiSumm_nc$Reach<-revalue(ChiSumm_nc$Reach, c("US"="Upstream"))
ChiSumm_nc$Reach<-revalue(ChiSumm_nc$Reach, c("DS"="Downstream"))
ggplot(ChiSumm_nc, aes(x=Time_Point, y=Chironomidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Chironomidae-se, ymax=Chironomidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name="Chironomidae Abundance (± SE)") +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#ephemerellidae
#78 zeros (out of 91)
#Check outliers
LPMacrosM_no_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Ephemerellidae)
#no outliers
#Check model assumptions
Ep.lm_nc<- lm(Ephemerellidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(Ep.lm_nc))
shapiro_test(residuals(Ep.lm_nc))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Ep.lm_nc)
#lambda =-2
LPMacrosM_no_nc$boxEp<-(LPMacrosM_no_nc$Ephemerellidae+1)^-2
boxEp.lm_nc<- lm(boxEp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(boxEp.lm_nc))
shapiro_test(residuals(boxEp.lm_nc))
#shaprio test is significant
LPMacrosM_no_nc %>% group_by(Leaf_Type) %>% shapiro_test(boxEp)
#shaprio test significant
ggqqplot(LPMacrosM_no_nc, "boxEp", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxEp ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no_nc %>% anova_test(boxEp ~ Reach*Leaf_Type*Time_Point_cat)
#Reach, time and reach x time significant
LPMacrosM_no_nc %>% emmeans_test(boxEp ~ Reach, p.adjust.method = "bonferroni",
                              model = boxEp.lm_nc)
#Gap greather than US and DS
LPMacrosM_no_nc %>% emmeans_test(boxEp ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxEp.lm_nc)
#not sig
LPMacrosM_no_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxEp ~ Reach, p.adjust.method = "bonferroni", model = boxEp.lm_nc)
#time 3: Gap > us and ds, time 4: gap>us
LPMacrosM_no_nc %>%
  group_by(Reach) %>%
  get_summary_stats(Ephemerellidae, type = "mean_se")

#ephemerellidae gap 
EpSumm_nc<-summarySE(LPMacrosM_nc, measurevar=c("Ephemerellidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
EpSumm_nc$Time_Point[EpSumm_nc$Time_Point == 1] <- 8
EpSumm_nc$Time_Point[EpSumm_nc$Time_Point == 2] <- 41
EpSumm_nc$Time_Point[EpSumm_nc$Time_Point == 3] <- 68
EpSumm_nc$Time_Point[EpSumm_nc$Time_Point == 4] <- 98
EpSumm_nc$Reach<-revalue(EpSumm_nc$Reach, c("US"="Upstream"))
EpSumm_nc$Reach<-revalue(EpSumm_nc$Reach, c("DS"="Downstream"))
Eplab_nc<-expression(paste(italic("Ephemerella"), " Abundance (± SE)"))
ggplot(EpSumm_nc, aes(x=Time_Point, y=Ephemerellidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Ephemerellidae-se, ymax=Ephemerellidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Eplab_nc)) +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Gammaridae 81 zeros (out of 144)
LPMacrosM_no_nc$Leaf_Type<-factor(LPMacrosM_no_nc$Leaf_Type)
#First absolute abundances
#Check outliers
LPMacrosM_no_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Gammaridae)
#no outliers
#Check model assumptions
Gam.lm_nc<- lm(Gammaridae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(Gam.lm_nc))
shapiro_test(residuals(Gam.lm_nc))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Gam.lm_nc)
#lambda =-2
LPMacrosM_no_nc$boxGam<-(LPMacrosM_no_nc$Gammaridae+1)^-2
boxGam.lm_nc<- lm(boxGam ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(boxGam.lm_nc))
shapiro_test(residuals(boxGam.lm_nc))
#shaprio test is significant
LPMacrosM_no_nc %>% group_by(Reach) %>% shapiro_test(boxGam)
#shaprio test significant
ggqqplot(LPMacrosM_no_nc, "boxGam", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxGam ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no_nc %>% anova_test(boxGam ~ Reach*Leaf_Type*Time_Point_cat)
#leaf, time, reach x time, leaf type x time and 3 way all significant
#just do overal leaf and reach
LPMacrosM_no_nc %>% emmeans_test(boxGam ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxGam.lm_nc)
#buckthorn and oak less than ash
LPMacrosM_no_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxGam ~ Reach, p.adjust.method = "bonferroni", model=boxGam.lm_nc)
#US<gap and ds time 2
GTPLT_nc<-LPMacrosM_no_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxGam ~ Leaf_Type, p.adjust.method = "bonferroni", model=boxGam.lm_nc)
#2 a > b,c,o
GTPLTR_nc<-LPMacrosM_no_nc %>% group_by(Time_Point_cat, Leaf_Type) %>%
  emmeans_test(boxGam ~ Reach, p.adjust.method = "bonferroni", model=boxGam.lm_nc)
LPMacrosM_no_nc %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(Gammaridae, type = "mean_se")

#Gammardiae indicates ash leaves
GamSumm_nc<-summarySE(LPMacrosM_nc, measurevar=c("Gammaridae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
GamSumm_nc$Time_Point[GamSumm_nc$Time_Point == 1] <- 8
GamSumm_nc$Time_Point[GamSumm_nc$Time_Point == 2] <- 41
GamSumm_nc$Time_Point[GamSumm_nc$Time_Point == 3] <- 68
GamSumm_nc$Time_Point[GamSumm_nc$Time_Point == 4] <- 98
GamSumm_nc$Reach<-revalue(GamSumm_nc$Reach, c("US"="Upstream"))
GamSumm_nc$Reach<-revalue(GamSumm_nc$Reach, c("DS"="Downstream"))
Gamlab_nc<-expression(paste("Mean ",italic("Gammarus"), " Abundance (± SE)"))
ggplot(GamSumm_nc, aes(x=Time_Point, y=Gammaridae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Gammaridae-se, ymax=Gammaridae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Gamlab_nc)) +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=18),
        legend.title=element_text(size=24),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 22))
  
#nemouridae 64 zeros (out of 91)
#Check outliers
LPMacrosM_no_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Nemouridae)
#no outliers
#Check model assumptions
Ne.lm_nc<- lm(Nemouridae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(Ne.lm_nc))
shapiro_test(residuals(Ne.lm_nc))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Ne.lm_nc)
#lambda =-2
LPMacrosM_no_nc$boxNe<-(LPMacrosM_no_nc$Nemouridae+1)^-2
boxNe.lm_nc<- lm(boxNe ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(boxNe.lm_nc))
shapiro_test(residuals(boxNe.lm_nc))
#shaprio test is significant
LPMacrosM_no_nc %>% group_by(Reach) %>% shapiro_test(boxNe)
#shaprio test significant
ggqqplot(LPMacrosM_no_nc, "boxNe", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxNe ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no_nc %>% anova_test(boxNe ~ Reach*Leaf_Type*Time_Point_cat)
#Reach, leaftype, time, and reach x time
LPMacrosM_no_nc %>% emmeans_test(boxNe ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxNe.lm_nc)
#ash greather than buckthorn
LPMacrosM_no_nc %>% emmeans_test(boxNe ~ Reach, p.adjust.method = "bonferroni",
                              model = boxNe.lm_nc)
#gap greater than us and ds
LPMacrosM_no_nc%>% emmeans_test(boxNe ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxNe.lm_nc)
#2 > 1,3
LPMacrosM_no_nc%>% group_by(Time_Point_cat) %>%
  emmeans_test(boxNe ~ Reach, p.adjust.method = "bonferroni", model = boxNe.lm_nc)
#time 2 gap>DS, time 4 gap > us or ds
LPMacrosM_no_nc%>%
  group_by(Leaf_Type) %>%
  get_summary_stats(Nemouridae, type = "mean_se")

#Nemouridae gap ash
NeSumm_nc<-summarySE(LPMacrosM_nc, measurevar=c("Nemouridae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
NeSumm_nc$Time_Point[NeSumm_nc$Time_Point == 1] <- 8
NeSumm_nc$Time_Point[NeSumm_nc$Time_Point == 2] <- 41
NeSumm_nc$Time_Point[NeSumm_nc$Time_Point == 3] <- 68
NeSumm_nc$Time_Point[NeSumm_nc$Time_Point == 4] <- 98
NeSumm_nc$Reach<-revalue(NeSumm_nc$Reach, c("US"="Upstream"))
NeSumm_nc$Reach<-revalue(NeSumm_nc$Reach, c("DS"="Downstream"))
Nelab_nc<-expression(paste("Mean ",italic("Nemoura"), " Abundance (± SE)"))
ggplot(NeSumm_nc, aes(x=Time_Point, y=Nemouridae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Nemouridae-se, ymax=Nemouridae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Nelab_nc)) +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=18),
        legend.title=element_text(size=24),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 22))

#simuliidae 73 zeros (out of 91)
#First absolute abundances
#Check outliers
LPMacrosM_no_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Simuliidae)
#no outliers
#Check model assumptions
Sim.lm_nc<- lm(Simuliidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(Sim.lm_nc))
shapiro_test(residuals(Sim.lm_nc))
#qqplot has many points that fall outside zone, significant, not normal
boxcox(Sim.lm_nc)
#lambda =-2
LPMacrosM_no_nc$boxSim<-(LPMacrosM_no_nc$Simuliidae+1)^-2
boxSim.lm_nc<- lm(boxSim ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(boxSim.lm_nc))
shapiro_test(residuals(boxSim.lm_nc))
#shaprio test is significant
LPMacrosM_no_nc %>% group_by(Leaf_Type) %>% shapiro_test(boxSim)
#shaprio test significant
ggqqplot(LPMacrosM_no_nc, "boxSim", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxSim ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no_nc %>% anova_test(boxSim ~ Reach*Leaf_Type*Time_Point_cat)
#reach, leaf type, time, reach x time, leaf type x time, and 3 way interaction significant
LPMacrosM_no_nc %>% emmeans_test(boxSim ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxSim.lm_nc)
#oak greater than ash and buckthorn
LPMacrosM_no_nc %>% emmeans_test(boxSim ~ Reach, p.adjust.method = "bonferroni",
                              model = boxSim.lm_nc)
#gap greater than us and ds
LPMacrosM_no_nc %>% emmeans_test(boxSim ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxSim.lm_nc)
#2 > 1,4
LPMacrosM_no_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxSim ~ Reach, p.adjust.method = "bonferroni", model = boxSim.lm_nc)
#US=gap>DS in 2, gap>USandDS 3
LPMacrosM_no_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxSim ~ Leaf_Type, p.adjust.method = "bonferroni", model = boxSim.lm_nc)
#oak > ash and buckthorn in 2 and 3
STPLYR_nc<-LPMacrosM_no_nc%>% group_by(Time_Point_cat, Leaf_Type) %>%
  emmeans_test(boxSim ~ Reach, p.adjust.method = "bonferroni", model = boxSim.lm_nc)
LPMacrosM_no_nc%>%
  group_by(Leaf_Type) %>%
  get_summary_stats(Simuliidae, type = "mean_se")

#gap oak
SiSumm_nc<-summarySE(LPMacrosM_nc, measurevar=c("Simuliidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
SiSumm_nc$Time_Point[SiSumm_nc$Time_Point == 1] <- 8
SiSumm_nc$Time_Point[SiSumm_nc$Time_Point == 2] <- 41
SiSumm_nc$Time_Point[SiSumm_nc$Time_Point == 3] <- 68
SiSumm_nc$Time_Point[SiSumm_nc$Time_Point == 4] <- 98
SiSumm_nc$Reach<-revalue(SiSumm_nc$Reach, c("US"="Upstream"))
SiSumm_nc$Reach<-revalue(SiSumm_nc$Reach, c("DS"="Downstream"))
Silab_nc<-expression(paste(italic("Prosimulium"), " Abundance (± SE)"))
ggplot(SiSumm_nc, aes(x=Time_Point, y=Simuliidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Simuliidae-se, ymax=Simuliidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Silab_nc)) +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Taeniopteridae indicates buckthorn and gap
#First absolute abundances
#Check outliers
LPMacrosM_no_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Taeniopterygidae)
#no outliers
#Check model assumptions
Ta.lm_nc<- lm(Taeniopterygidae+1 ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(Ta.lm_nc))
shapiro_test(residuals(Ta.lm_nc))
#qqplot has many points that fall outside zone, significant, not normal 
boxcox(Ta.lm_nc)
#lambda =-1
LPMacrosM_no_nc$boxTa<-(LPMacrosM_no_nc$Taeniopterygidae+1)^-1
boxTa.lm_nc<- lm(boxTa ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
ggqqplot(residuals(boxTa.lm_nc))
shapiro_test(residuals(boxTa.lm_nc))
#shaprio test is significant
LPMacrosM_no_nc %>% group_by(Leaf_Type) %>% shapiro_test(boxTa)
#shaprio test significant
ggqqplot(LPMacrosM_no_nc, "boxTa", ggtheme = theme_bw()) +
  facet_grid(Reach + Leaf_Type ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxTa ~ Reach*Leaf_Type*Time_Point_cat, data = LPMacrosM_no_nc)
#not significant, therefore assume homogeniety of variance
LPMacrosM_no_nc %>% anova_test(boxTa ~ Reach*Leaf_Type*Time_Point_cat)
#reach, leaf type, time, and leaf type x time significant
LPMacrosM_no_nc %>% emmeans_test(boxTa ~ Leaf_Type, p.adjust.method = "bonferroni",
                              model = boxTa.lm_nc)
#ash greater than buckthorn and cotton, oak greater than buckthorn
LPMacrosM_no_nc %>% emmeans_test(boxTa ~ Reach, p.adjust.method = "bonferroni",
                              model = boxTa.lm_nc)
#gap greater than us and ds
LPMacrosM_no_nc %>% emmeans_test(boxTa ~ Time_Point_cat, p.adjust.method = "bonferroni",
                              model = boxTa.lm_nc)
#(1,2) > 4, 2>3>4
TTPLT_nc<-LPMacrosM_no_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(boxTa ~ Leaf_Type, p.adjust.method = "bonferroni", model = boxTa.lm_nc)
LPMacrosM_no_nc%>%
  group_by(Leaf_Type) %>%
  get_summary_stats(Taeniopterygidae, type = "mean_se")

#visualize taen gap
#absolute abundances
TaenSumm_nc<-summarySE(LPMacrosM_nc, measurevar=c("Taeniopterygidae"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
TaenSumm_nc$Time_Point[TaenSumm_nc$Time_Point == 1] <- 8
TaenSumm_nc$Time_Point[TaenSumm_nc$Time_Point == 2] <- 41
TaenSumm_nc$Time_Point[TaenSumm_nc$Time_Point == 3] <- 68
TaenSumm_nc$Time_Point[TaenSumm_nc$Time_Point == 4] <- 98
TaenSumm_nc$Reach<-revalue(TaenSumm_nc$Reach, c("US"="Upstream"))
TaenSumm_nc$Reach<-revalue(TaenSumm_nc$Reach, c("DS"="Downstream"))
Taenlab_nc<-expression(paste(italic("Taeniopteryx"), " Abundance (± SE)"))
ggplot(TaenSumm_nc, aes(x=Time_Point, y=Taeniopterygidae, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=Taeniopterygidae-se, ymax=Taeniopterygidae+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name=c(Taenlab_nc)) +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Functional feeding groups
str(LeafPackFFG)
LeafPackFFG_nc<-subset(LeafPackFFG, Leaf_Type!="Cotton")
LeafPackFFG_nc$rowsum<-rowSums(LeafPackFFG_nc[,2:6])

LeafPackFFGno_nc<-subset(LeafPackFFG_nc, rowsum!=0)
LeafPackFFGCom_nc<-LeafPackFFGno_nc[,2:6]

LPFFGEnv_nc<-LeafPackFFGno_nc[,c(1,7:ncol(LeafPackFFGno_nc))]
LPFFGComtot_nc<-LeafPackFFG_nc[,2:6]
LPFFGComtotbray0_nc<-as.matrix(bray0(LPFFGComtot_nc))
LPFFGComRA_nc<-data.frame(make_relative(as.matrix(LeafPackFFG_nc[,2:6])))
LPFFGComRA_nc[is.nan(LPFFGComRA_nc)] <- 0
LPFFGEnvtot_nc<-LeafPackFFG_nc[,c(1,7:ncol(LeafPackFFG_nc))]
LPFFGRA_nc<-merge(LPFFGComRA_nc, LPFFGEnvtot_nc, by=0)
FFGReachtot_nc<-droplevels(as.factor(LPFFGEnvtot_nc$Reach))

#delete time point 0 days
LeafPackFFG_no_nc<-subset(LeafPackFFG_nc, Time_Point!=0)

LeafPackFFG_noComtot_nc<-LeafPackFFG_no_nc[,2:6]
LeafPackFFG_noComtotbray0_nc<-as.matrix(bray0(LeafPackFFG_noComtot_nc))
LeafPackFFG_noComRA_nc<-data.frame(make_relative(as.matrix(LeafPackFFG_no_nc[,2:6])))
LeafPackFFG_noComRA_nc[is.nan(LeafPackFFG_noComRA_nc)] <- 0
LeafPackFFG_noEnvtot_nc<-LeafPackFFG_no_nc[,c(1,7:ncol(LeafPackFFG_no_nc))]
LeafPackFFG_noRA_nc<-merge(LeafPackFFG_noComRA_nc, LeafPackFFG_noEnvtot_nc, by=0)

#permanova
adonis2(as.dist(LeafPackFFG_noComtotbray0_nc)~Reach*Leaf_Type*Days_Exposure, data=LeafPackFFG_noEnvtot_nc, permutations=999)
#reach, leaf type signficant

FFGReachtot_nc<-droplevels(as.factor(LPFFGEnvtot_nc$Reach))
FFGLeaf_Taxatot_nc<-droplevels(as.factor(LPFFGEnvtot_nc$Leaf_Type))

FFGReachtot_no_nc<-droplevels(as.factor(LeafPackFFG_noEnvtot_nc$Reach))
FFGLeaf_Taxatot_no_nc<-droplevels(as.factor(LeafPackFFG_noEnvtot_nc$Leaf_Type))

#NMDS analysis
AC_FFG_NMDS_nc<-metaMDS(as.dist(LPFFGComtotbray0_nc))
stressplot(AC_FFG_NMDS_nc)
ordiplot(AC_FFG_NMDS_nc, type="n")
with(AC_FFG_NMDS_nc, points(AC_FFG_NMDS_nc, display="sites", col=leaftaxacolvec_nc[FFGLeaf_Taxatot_nc], pch=19))
with(AC_FFG_NMDS_nc, legend("topleft", legend=levels(FFGLeaf_Taxatot_nc), bty="n", col=leaftaxacolvec_nc, pch=19, pt.bg=leaftaxacolvec_nc))
with(AC_FFG_NMDS_nc, ordiellipse(AC_FFG_NMDS_nc, FFGLeaf_Taxatot_nc, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_FFG_NMDS_nc, ordiellipse(AC_FFG_NMDS_nc, FFGLeaf_Taxatot_nc, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_FFG_NMDS_nc, ordiellipse(AC_FFG_NMDS_nc, FFGLeaf_Taxatot_nc, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

FFGReachtot_nc<-revalue(FFGReachtot_nc, c("US"="Upstream"))
FFGReachtot_nc<-revalue(FFGReachtot_nc, c("DS"="Downstream"))
ordiplot(AC_FFG_NMDS_nc, type="n")
with(AC_FFG_NMDS_nc, points(AC_FFG_NMDS_nc, display="sites", col=reach_col_vec[FFGReachtot_nc], pch=19))
with(AC_FFG_NMDS_nc, legend("topleft", legend=levels(FFGReachtot_nc), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_FFG_NMDS_nc, ordiellipse(AC_FFG_NMDS_nc, Reachtot_nc, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "Upstream"))
with(AC_FFG_NMDS_nc, ordiellipse(AC_FFG_NMDS_nc, Reachtot_nc, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_FFG_NMDS_nc, ordiellipse(AC_FFG_NMDS_nc, Reachtot_nc, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "Downstream"))

AC_FFG_Com_indic_leaf_nc<-data.frame(signassoc(LeafPackFFG_noComtot_nc, cluster=FFGLeaf_Taxatot_no_nc,  mode=0, 
                                               alternative = "two.sided",control = how(nperm=999)))
AC_FFG_Com_indic_leaf_sig_nc<-subset(AC_FFG_Com_indic_leaf_nc, psidak<=0.05)
#shredder buckthorn
AC_FFG_Com_reach_indic_nc<-data.frame(signassoc(LeafPackFFG_noComtot_nc, cluster=FFGReachtot_no_nc,  mode=0, 
                                                alternative = "two.sided",control = how(nperm=999)))
AC_FFG_Com_reach_indic_sig_nc<-subset(AC_FFG_Com_reach_indic_nc, psidak<=0.05)
#collector filt, shredder gap, gatherer us

#model collector filt, collector gath, shredder because of indicator species analy
#also grazer because addressed in hypotheses

#Collector filterer
LeafPackFFG_noRA_nc$Time_Point_cat<-as.factor(LeafPackFFG_noRA_nc$Time_Point)
LeafPackFFGnoRA_nc<-subset(LeafPackFFG_noRA_nc, rowsum>0)
#summary stats
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(CollectorFilterer, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(CollectorFilterer)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA_nc, x = "Reach", y = "CollectorFilterer",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA_nc %>%
  group_by(Reach) %>%
  shapiro_test(CollectorFilterer)
#not normal
ggqqplot(LeafPackFFGnoRA_nc, "CollectorFilterer", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
CF.lm_nc<- lm((CollectorFilterer+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA_nc)
ggqqplot(residuals(CF.lm))
shapiro_test(residuals(CF.lm))
#0.0000000423 not normal, see boxcox
boxcox(CF.lm_nc)
#labmda=-2. transform
LeafPackFFGnoRA_nc$boxCF<-(LeafPackFFGnoRA_nc$CollectorFilterer+1)^-2
#check normality on transformed CF
LeafPackFFGnoRA_nc %>%
  group_by(Reach) %>%
  shapiro_test(boxCF)
#not normal
LeafPackFFGnoRA_nc %>% group_by(Leaf_Type) %>% shapiro_test(boxCF)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA_nc, "boxCF", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxCF ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA_nc)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA_nc %>% anova_test(boxCF ~ Reach*Leaf_Type*Time_Point_cat)
#Nothing significant

#Collector gatherer
LeafPackFFGnoRA_nc$Leaf_Type<-factor(LeafPackFFGnoRA_nc$Leaf_Type)
#summary stats
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(CollectorGatherer, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(CollectorGatherer)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA_nc, x = "Reach", y = "CollectorGatherer",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA_nc %>%
  group_by(Reach) %>%
  shapiro_test(CollectorGatherer)
#not normal
ggqqplot(LeafPackFFGnoRA_nc, "CollectorGatherer", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
CG.lm_nc<- lm((CollectorGatherer+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA_nc)
ggqqplot(residuals(CG.lm_nc))
shapiro_test(residuals(CG.lm_nc))
#0.00000000237 not normal, see boxcox
boxcox(CG.lm_nc)
#labmda=-2. transform
LeafPackFFGnoRA_nc$boxCG<-(LeafPackFFGnoRA_nc$CollectorGatherer+1)^-2
#check normality on transformed CF
LeafPackFFGnoRA_nc %>%
  group_by(Reach) %>%
  shapiro_test(boxCG)
#not normal
LeafPackFFGnoRA_nc %>% group_by(Leaf_Type) %>% shapiro_test(boxCG)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA_nc, "boxCG", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxCG ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA_nc)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA_nc %>% anova_test(boxCG ~ Reach*Leaf_Type*Time_Point_cat)
#3 way interaction significant
boxCG.lm_nc<-lm(boxCG ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA_nc)
LeafPackFFGnoRA_nc %>%
  group_by(Time_Point_cat) %>%
  anova_test(boxCG ~ Reach*Leaf_Type, error = boxCG.lm_nc)
#leaf type significant on time 2
#Reach x leaf type significant on time 2, 3
LeafPackFFGnoRA_nc %>% group_by(Time_Point_cat) %>%
  anova_test(boxCG ~ Reach, error = boxCG.lm_nc)
#nothing significant
threewayboxCG_nc<-LeafPackFFGnoRA_nc  %>%
  group_by(Time_Point_cat, Leaf_Type) %>%
  emmeans_test(boxCG ~ Reach, p.adjust.method = "bonferroni", model=boxCG.lm_nc)
#time point 2 buckthorn, gap less than ds
threewayboxCG2_nc<-LeafPackFFGnoRA_nc  %>%
  group_by(Time_Point_cat, Reach) %>%
  emmeans_test(boxCG ~ Leaf_Type, p.adjust.method = "bonferroni", model=boxCG.lm_nc)

threewayCG_nc<-LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(CollectorGatherer, type = "mean_se")

#Visualize
CGSumm_nc<-summarySE(LPFFGRA_nc, measurevar=c("CollectorGatherer"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
CGSumm_nc$Time_Point[CGSumm_nc$Time_Point == 1] <- 8
CGSumm_nc$Time_Point[CGSumm_nc$Time_Point == 2] <- 41
CGSumm_nc$Time_Point[CGSumm_nc$Time_Point == 3] <- 68
CGSumm_nc$Time_Point[CGSumm_nc$Time_Point == 4] <- 98
CGSumm_nc$Reach<-revalue(CGSumm_nc$Reach, c("US"="Upstream"))
CGSumm_nc$Reach<-revalue(CGSumm_nc$Reach, c("DS"="Downstream"))
ggplot(CGSumm_nc, aes(x=Time_Point, y=CollectorGatherer, color=Leaf_Type)) +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=CollectorGatherer-se, ymax=CollectorGatherer+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_y_continuous(name="Mean Collector-Gatherer Rel. Abund. (± SE)") +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Taxa") +
  facet_wrap(.~Reach) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=24),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=18),axis.text.y = element_text(size=18),
        legend.title=element_text(size=24),legend.text = element_text(size=20),
        strip.text.x = element_text(size = 22))

#Shredder
#summary stats
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(Shredder, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Shredder)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA_nc, x = "Reach", y = "Shredder",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA_nc %>% group_by(Reach) %>% shapiro_test(Shredder)
#not normal
ggqqplot(LeafPackFFGnoRA_nc, "Shredder", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
Sh.lm_nc<- lm((Shredder+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA_nc)
ggqqplot(residuals(Sh.lm_nc))
shapiro_test(residuals(Sh.lm_nc))
#0.0129 not normal, see boxcox
boxcox(Sh.lm_nc)
#labmda=0.5. transform
LeafPackFFGnoRA_nc$boxSh<-(LeafPackFFGnoRA_nc$Shredder+1)^0.5
#check normality on transformed CF
LeafPackFFGnoRA_nc %>% group_by(Reach) %>% shapiro_test(boxSh)
#not normal
LeafPackFFGnoRA_nc %>% group_by(Leaf_Type) %>% shapiro_test(boxSh)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA_nc, "boxSh", ggtheme = theme_bw()) +
  facet_grid(Reach  ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxSh ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA_nc)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA_nc %>% anova_test(boxSh ~ Reach*Leaf_Type*Time_Point_cat)
#time significant
LeafPackFFGnoRA_nc %>% emmeans_test(boxSh ~ Time_Point_cat, p.adjust.method = "bonferroni")

#Grazer
#summary stats
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(Grazer, type = "mean_sd")
#check for outliers
LeafPackFFGnoRA_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(Grazer)
#no outliers
#visualize
ggboxplot(LeafPackFFGnoRA_nc, x = "Reach", y = "Grazer",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
LeafPackFFGnoRA_nc %>% group_by(Reach) %>% shapiro_test(Grazer)
#not normal
ggqqplot(LeafPackFFGnoRA_nc, "Grazer", ggtheme = theme_bw()) +
  facet_grid(Reach ~ Time_Point_cat, labeller = "label_both")
Gr.lm_nc<- lm((Grazer+1) ~ Reach*Leaf_Type*Time_Point_cat, 
           data = LeafPackFFGnoRA_nc)
ggqqplot(residuals(Gr.lm_nc))
shapiro_test(residuals(Gr.lm_nc))
#4.07e-11 not normal, see boxcox
boxcox(Gr.lm_nc)
#labmda=-2. transform
LeafPackFFGnoRA_nc$boxGr<-(LeafPackFFGnoRA_nc$Grazer+1)^-2
#check normality on transformed CF
LeafPackFFGnoRA_nc%>% group_by(Reach) %>% shapiro_test(boxGr)
#not normal
LeafPackFFGnoRA_nc%>% group_by(Leaf_Type) %>% shapiro_test(boxGr)
#shaprio test significant
ggqqplot(LeafPackFFGnoRA_nc, "boxGr", ggtheme = theme_bw()) +
  facet_grid(Reach  ~ Time_Point_cat, labeller = "label_both")
#looks good, test homogeniety of variance
leveneTest(boxGr ~ Reach*Leaf_Type*Time_Point_cat, data = LeafPackFFGnoRA_nc)
#not significant, therefore assume homogeniety of variance
LeafPackFFGnoRA_nc%>% anova_test(boxGr ~ Reach*Leaf_Type*Time_Point_cat)
#Nothing significant

#Make stacked bar graphs
LPFFGmelt_nc<-subset(LPFFGmelt, Leaf_Type!="Cotton")
LPFFGmelt_nc$Leaf_Type<-factor(LPFFGmelt_nc$Leaf_Type)
ggplot(LPFFGmelt_nc, aes(fill=FFG, y=value, x=Leaf_Type)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=ffgcolvec,name = "Functional\nFeeding Group") +
  facet_grid( ~ Reach) +
  xlab("Leaf Type") +
  scale_y_continuous(name="Relative Abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))
ggplot(LPFFGmelt_nc, aes(fill=FFG, y=value, x=Reach)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=ffgcolvec,name = "Functional\nFeeding Group") +
  facet_grid(Leaf_Type ~.) +
  scale_y_continuous(name="Relative Abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.y = element_text(size = 16))

###############################
#Microbes
############################

#Use QIIME2 for sequence processing

#Upload unifrac distance table
AC_16S_uni<-read.csv("AC_WUnifrac.csv", header = T, check.names=F)

#Combine distance matrix with environmental variables
LPMicM<-merge(AC_16S_uni, LPMetadata, by="Pack_ID")
row.names(LPMicM)<-LPMicM[,1]
LPMicM<-LPMicM[,-c(1)]
names(LPMicM)
AC_16S_uni<-as.matrix(LPMicM[,c(1:138)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
LPMicM_env<-LPMicM[,c(139:144)]
LPMicM_env_R<-LPMicM_env$Reach
LPMicM_env_LT<-as.factor(LPMicM_env$Leaf_Type)

#UNI-Overall permanova with unifrac distances
adonis2(as.dist(AC_16S_uni) ~ Leaf_Type*Reach*Days_Exposure, data=LPMicM_env,
       permutations=999)
#Leaf type, reach, time, leafxtime and reach x time significant

#Visualize via nmds
AC_Mic_NMDS<-metaMDS(as.dist(AC_16S_uni))
#stress 0.1

#Stressplot macroinvertebrate Nmds
stressplot(AC_Mic_NMDS)

#NMDS plot for Reach
ordiplot(AC_Mic_NMDS, type="n")
with(AC_Mic_NMDS, points(AC_Mic_NMDS, display="sites", col=reach_col_vec[LPMicM_env_R], pch=19))
with(AC_Mic_NMDS, legend("topleft", legend=levels(LPMicM_env_R), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_R, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "US"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_R, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_R, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "DS"))

#NMDS plot for leaf type
ordiplot(AC_Mic_NMDS, type="n")
with(AC_Mic_NMDS, points(AC_Mic_NMDS, display="sites", col=leaftaxacolvec[LPMicM_env_LT], pch=19))
with(AC_Mic_NMDS, legend("topleft", legend=levels(LPMicM_env_LT), bty="n", col=leaftaxacolvec, pch=19, pt.bg=leaftaxacolvec))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Cotton"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

#Upload family level taxonomy
AC_16S_f<-read.csv("AC_f.csv", header=T, check.names=F)
#Format data frame so the family is row name
row.names(AC_16S_f)<-AC_16S_f[,1]
#Delete otu id column, now that otu id is row name
AC_16S_f$Pack_ID<-NULL
AC_16S_f_t<-t(AC_16S_f)
AC_16S_f_t<-data.frame(AC_16S_f_t, check.names = FALSE)
AC_16S_f_t$Pack_ID<-row.names(AC_16S_f_t)
#Combine with environmental variables
AC_16S_f_map <-merge(LPMetadata,AC_16S_f_t, by="Pack_ID")
names(AC_16S_f_map)
sort(colSums(AC_16S_f_t[,1:536]))
#most common family is CComamonadaceae 
stat.desc(AC_16S_f_t$`k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae`)
#269 +/- 14
AC_16S_f_map_R<-AC_16S_f_map$Reach
AC_16S_f_map_L<-AC_16S_f_map$Leaf_Type

#indicator species analysis for reach
AC_Mic_Com_R_indic<-data.frame(signassoc(AC_16S_f_t[,1:536], cluster=AC_16S_f_map_R,  mode=0, alternative = "two.sided",
                                         control = how(nperm=999)))
AC_Mic_Com_R_indic_sig<-subset(AC_Mic_Com_R_indic, psidak<=0.05)
#27 indicator families for watershed

#indicator species analysis for leaf type
AC_Mic_Com_L_indic<-data.frame(signassoc(AC_16S_f_t[,1:536], cluster=AC_16S_f_map_L,  mode=0, alternative = "two.sided",control = how(nperm=999)))
AC_Mic_Com_L_indic_sig<-subset(AC_Mic_Com_L_indic, psidak<=0.05)
#122 indicator families for leaf type

#Upload phylogenetic diversity
AC_16S_fpd<- read.csv("AC_fpd.csv", header=T)
#Combine with environmental variables
AC_16S_fpd_map<-merge(LPMetadata, AC_16S_fpd, by="Pack_ID")
AC_16S_fpd_map$Time_Point_cat<-as.factor(AC_16S_fpd_map$Time_Point)

#Faith's PD
range(AC_16S_fpd_map$faith_pd)
#8.76-54.38
#summary stats
AC_16S_fpd_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(faith_pd, type = "mean_sd")
#check for outliers
AC_16S_fpd_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(faith_pd)
#no outliers
#visualize
ggboxplot(AC_16S_fpd_map, x = "Reach", y = "faith_pd",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_fpd_map %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(faith_pd)
#normal
fpd.lm<- lm(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map)
ggqqplot(residuals(fpd.lm))
shapiro_test(residuals(fpd.lm))
#residuals fall aprox in qqplot, and shapiro test not significant
#looks good, test homogeniety of variance
leveneTest(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map)
#not significant, therefore assume homogeniety of variance
res.aov.fpd<- AC_16S_fpd_map %>% anova_test(faith_pd ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.fpd
#leaf type, time point, and leaf typextime point significant, so group the data by leaf type and time, and fit  anova
AC_16S_fpd_map %>% group_by(Leaf_Type) %>%
  anova_test(faith_pd~ Days_Exposure, error = fpd.lm)
AC_16S_fpd_map%>% group_by(Time_Point_cat) %>%
  anova_test(faith_pd ~ Leaf_Type, error = fpd.lm)
# Pairwise comparisons
AC_16S_fpd_leafxtimecat<-AC_16S_fpd_map %>% group_by(Leaf_Type) %>% emmeans_test(faith_pd ~ Time_Point_cat, 
                                                      p.adjust.method = "bonferroni",
                                                      model=fpd.lm)
AC_16S_fpd_timecatxleaf<-AC_16S_fpd_map %>% group_by(Time_Point_cat) %>%
  emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=fpd.lm)
#emm for non interactions
AC_16S_fpd_map%>% emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni")
#cotton greater than ash and buckthorn, ash and cotton greater than oak
AC_16S_fpd_map%>% emmeans_test(faith_pd ~ Time_Point_cat, p.adjust.method = "bonferroni")
#0<1=2=3=4
#See summary statistics for significant groups
AC_16S_fpd_map_A<-subset(AC_16S_fpd_map, Leaf_Type=="Ash")
AC_16S_fpd_map_B<-subset(AC_16S_fpd_map, Leaf_Type=="Buckthorn")
AC_16S_fpd_map_C<-subset(AC_16S_fpd_map, Leaf_Type=="Cotton")
AC_16S_fpd_map_O<-subset(AC_16S_fpd_map, Leaf_Type=="Oak")
stat.desc(AC_16S_fpd_map_A$faith_pd)
stat.desc(AC_16S_fpd_map_B$faith_pd)
stat.desc(AC_16S_fpd_map_C$faith_pd)
stat.desc(AC_16S_fpd_map_O$faith_pd)
AC_16S_fpd_map_T0<-subset(AC_16S_fpd_map, Time_Point_cat=="0")
AC_16S_fpd_map_T1<-subset(AC_16S_fpd_map, Time_Point_cat=="1")
AC_16S_fpd_map_T2<-subset(AC_16S_fpd_map, Time_Point_cat=="2")
AC_16S_fpd_map_T3<-subset(AC_16S_fpd_map, Time_Point_cat=="3")
AC_16S_fpd_map_T4<-subset(AC_16S_fpd_map, Time_Point_cat=="4")
stat.desc(AC_16S_fpd_map_T0$faith_pd)
stat.desc(AC_16S_fpd_map_T1$faith_pd)
stat.desc(AC_16S_fpd_map_T2$faith_pd)
stat.desc(AC_16S_fpd_map_T3$faith_pd)
stat.desc(AC_16S_fpd_map_T4$faith_pd)

#visualize reach, leaf type and time
fpdsummary<-summarySE(AC_16S_fpd_map, measurevar=c("faith_pd"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
fpdsummary$Time_Point[fpdsummary$Time_Point == 1] <- 8
fpdsummary$Time_Point[fpdsummary$Time_Point == 2] <- 41
fpdsummary$Time_Point[fpdsummary$Time_Point == 3] <- 68
fpdsummary$Time_Point[fpdsummary$Time_Point == 4] <- 98
fpdsummary$Reach<-revalue(fpdsummary$Reach, c("US"="Upstream"))
fpdsummary$Reach<-revalue(fpdsummary$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(fpdsummary, aes(x=Time_Point, y=faith_pd, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=faith_pd-se, ymax=faith_pd+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Phylogenetic Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

fpdsummarynr<-summarySE(AC_16S_fpd_map, measurevar=c("faith_pd"), groupvars=c("Leaf_Type","Time_Point"), na.rm=TRUE)
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 1] <- 8
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 2] <- 41
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 3] <- 68
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 4] <- 98
# plot on y fpd and x days, no facet
ggplot(fpdsummarynr, aes(x=Time_Point, y=faith_pd, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=faith_pd-se, ymax=faith_pd+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Phylogenetic Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Chao1
#Upload chao1
AC_16S_ch<- read.csv("AC_chao.csv", header=T)
#Combine with environmental variables
AC_16S_ch_map<-merge(LPMetadata, AC_16S_ch, by="Pack_ID")
AC_16S_ch_map$Time_Point_cat<-as.factor(AC_16S_ch_map$Time_Point)

#Chao1
range(AC_16S_ch_map$chao1)
#40-903.9
#summary stats
AC_16S_ch_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(chao1, type = "mean_sd")
#check for outliers
AC_16S_ch_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(chao1)
#no outliers
#visualize
ggboxplot(AC_16S_ch_map, x = "Reach", y = "chao1",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_ch_map %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(chao1)
#normal
ch.lm<- lm(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map)
ggqqplot(residuals(ch.lm))
shapiro_test(residuals(ch.lm))
#residuals fall aprox in qqplot, and shapiro test significant
#test homogeniety of variance
leveneTest(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map)
#not significant, therefore assume homogeniety of variance
res.aov.ch<- AC_16S_ch_map %>% anova_test(chao1 ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.ch
#leaf type, time point, and leaf typextime point significant, so group the data by leaf type and time, and fit  anova
AC_16S_ch_map %>% group_by(Leaf_Type) %>%
  anova_test(chao1~ Days_Exposure, error = fpd.lm)
AC_16S_ch_map%>% group_by(Time_Point_cat) %>%
  anova_test(chao1 ~ Leaf_Type, error = fpd.lm)
# Pairwise comparisons
AC_16S_ch_map %>% group_by(Leaf_Type) %>% emmeans_test(chao1 ~ Time_Point_cat, 
                                                        p.adjust.method = "bonferroni",
                                                        model=ch.lm)
AC_16S_ch_timexleaf<-AC_16S_ch_map %>% group_by(Time_Point_cat) %>%
  emmeans_test(chao1 ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=ch.lm)
#emm for non interactions
AC_16S_ch_map%>% emmeans_test(chao1 ~ Leaf_Type, p.adjust.method = "bonferroni")
#Cotton>Ash=buckthorn>Oak
AC_16S_ch_map%>% emmeans_test(chao1 ~ Time_Point_cat, p.adjust.method = "bonferroni")
#0<1=2=3=4
#See summary statistics for significant groups
AC_16S_ch_map_A<-subset(AC_16S_ch_map, Leaf_Type=="Ash")
AC_16S_ch_map_B<-subset(AC_16S_ch_map, Leaf_Type=="Buckthorn")
AC_16S_ch_map_C<-subset(AC_16S_ch_map, Leaf_Type=="Cotton")
AC_16S_ch_map_O<-subset(AC_16S_ch_map, Leaf_Type=="Oak")
stat.desc(AC_16S_ch_map_A$chao1)
stat.desc(AC_16S_ch_map_B$chao1)
stat.desc(AC_16S_ch_map_C$chao1)
stat.desc(AC_16S_ch_map_O$chao1)

#visualize reach, leaf type and time
chsummary<-summarySE(AC_16S_ch_map, measurevar=c("chao1"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
chsummary$Time_Point[chsummary$Time_Point == 1] <- 8
chsummary$Time_Point[chsummary$Time_Point == 2] <- 41
chsummary$Time_Point[chsummary$Time_Point == 3] <- 68
chsummary$Time_Point[chsummary$Time_Point == 4] <- 98
chsummary$Reach<-revalue(chsummary$Reach, c("US"="Upstream"))
chsummary$Reach<-revalue(chsummary$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(chsummary, aes(x=Time_Point, y=chao1, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=chao1-se, ymax=chao1+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Richness (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#now move on to venn diagrams to find unique OTUs introduced in time 0
#create venn diagram for shared OTUs among leaf type on day 0
#upload OTU table
AC_OTU<- read.csv("AC_otu_table.csv", header=T,check.names=F)
names(AC_OTU)
#Format data frame so the family is row name
row.names(AC_OTU)<-AC_OTU[,1]
#Delete otu id column, now that otu id is row name
AC_OTU$Pack_ID<-NULL
AC_OTU_t<-t(AC_OTU)
AC_OTU_t<-data.frame(AC_OTU_t, check.names = FALSE)
AC_OTU_t$Pack_ID<-row.names(AC_OTU_t)
#Combine with environmental variables
AC_OTU_map<-merge(LPMetadata, AC_OTU_t, by="Pack_ID")
#subset into only time point 0
AC_OTU_map_0<-subset(AC_OTU_map, Time_Point==0)
#first list OTU names found in each group
names(AC_OTU_map_0)
AC_OTU_map_0_ag_venn<-aggregate(AC_OTU_map_0[8:ncol(AC_OTU_map_0)], 
                                   by=list(Leaf_Type=AC_OTU_map_0$Leaf_Type),
                                   FUN=sum)
str(AC_OTU_map_0_ag_venn)
row.names(AC_OTU_map_0_ag_venn)<-AC_OTU_map_0_ag_venn$Leaf_Type
AC_OTU_map_0_ag_venn$Leaf_Type<-NULL
row.names(AC_OTU_map_0_ag_venn)
Ash<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Ash",] > 0]
Buckthorn<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Buckthorn",] > 0]
Cotton<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Cotton",] > 0]
Oak<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Oak",] > 0]
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")
LTVen<-venn_diagram4(Ash, Buckthorn, Cotton, Oak,
                      "Ash", "Buckthorn", "Cotton", "Oak",
                      colors=leaftaxacolvec)
#Each leaf type has unique OTUs
#visualize without cotton
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
quartz()
LTVen_nc<-venn_diagram3(Ash,Buckthorn,Oak, "Ash",
                       "Buckthorn","Oak",colors=leaftaxacolvec_nc)
#each leaf type has unique OTUs

##################
#now do other microbial analyses without cotton
#################
LPMIcM_nc<-subset(LPMicM, Leaf_Type!="Cotton")
names(LPMIcM_nc)
LPMIcM_nc<- select(LPMIcM_nc, -contains("-C-"))
names(LPMIcM_nc)
AC_16S_uni_nc<-as.matrix(LPMIcM_nc[,c(1:106)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
LPMicM_env_nc<-LPMIcM_nc[,c(107:112)]
LPMicM_env_nc$Reach<-revalue(LPMicM_env_nc$Reach, c("US"="Upstream", "DS"="Downstream"))
LPMicM_env_R_nc<-LPMicM_env_nc$Reach
LPMicM_env_LT_nc<-droplevels(as.factor(LPMicM_env_nc$Leaf_Type))

#UNI-Overall permanova with unifrac distances
adonis(as.dist(AC_16S_uni_nc) ~ Leaf_Type*Reach*Days_Exposure, data=LPMicM_env_nc,
       permutations=999)
#Leaf type, reach, time, leafxtime and reach x time significant

#Visualize via nmds
AC_Mic_NMDS_nc<-metaMDS(as.dist(AC_16S_uni_nc))
#stress 0.1

#Stressplot microbial Nmds
stressplot(AC_Mic_NMDS_nc)

#NMDS plot for Reach
ordiplot(AC_Mic_NMDS_nc, type="n")
with(AC_Mic_NMDS_nc, points(AC_Mic_NMDS_nc, display="sites", col=reach_col_vec[LPMicM_env_R_nc], pch=19))
with(AC_Mic_NMDS_nc, legend("topleft", legend=levels(LPMicM_env_R_nc), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_R_nc, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "Upstream"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_R_nc, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_R_nc, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "Downstream"))

#NMDS plot for leaf type
ordiplot(AC_Mic_NMDS_nc, type="n")
with(AC_Mic_NMDS_nc, points(AC_Mic_NMDS_nc, display="sites", col=leaftaxacolvec_nc[LPMicM_env_LT_nc], pch=19))
with(AC_Mic_NMDS_nc, legend("topleft", legend=levels(LPMicM_env_LT_nc), bty="n", col=leaftaxacolvec_nc, pch=19, pt.bg=leaftaxacolvec_nc))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_LT_nc, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_LT_nc, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_LT_nc, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

#family level taxonomy
AC_16S_f_map_nc<-subset(AC_16S_f_map, Leaf_Type!="Cotton")
names(AC_16S_f_map_nc)
sort(colSums(AC_16S_f_map_nc[,8:543]))
#most common family is Comamonadaceae 
stat.desc(AC_16S_f_map_nc$`k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae`)
#261 +/- 18
AC_16S_f_map_R_nc<-AC_16S_f_map_nc$Reach
AC_16S_f_map_L_nc<-droplevels(AC_16S_f_map_nc$Leaf_Type)

#indicator species analysis for reach
AC_Mic_Com_R_indic_nc<-signassoc(AC_16S_f_map_nc[,8:543], cluster=AC_16S_f_map_R_nc,  mode=0, alternative = "two.sided",control = how(nperm=999))
AC_Mic_Com_R_indic_sig_nc<-subset(AC_Mic_Com_R_indic_nc, psidak<=0.05)
#23 indicator families for reach, 2 for gap
#k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f_
#k__Bacteria;p__Bacteroidetes;c__Sphingobacteriia;o__Sphingobacteriales;f__
AC_16S_f_map_nc%>%
  group_by(Reach) %>%
  get_summary_stats("k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__", type = "mean_se")
#lowered in gaps
AC_16S_f_map_nc%>%
  group_by(Reach) %>%
  get_summary_stats("k__Bacteria;p__Bacteroidetes;c__Sphingobacteriia;o__Sphingobacteriales;f__", type = "mean_se")
#lowered in gaps

#indicator species analysis for leaf type
AC_Mic_Com_L_indic_nc<-signassoc(AC_16S_f_map_nc[,8:543], cluster=AC_16S_f_map_L_nc,  mode=0, alternative = "two.sided",control = how(nperm=999))
AC_Mic_Com_L_indic_sig_nc<-subset(AC_Mic_Com_L_indic_nc, psidak<=0.05)
#48 indicator families for leaf type

#Upload phylogenetic diversity
AC_16S_fpd_map_nc<-subset(AC_16S_fpd_map, Leaf_Type!="Cotton")
AC_16S_fpd_map_nc$Leaf_Type<-droplevels(AC_16S_fpd_map_nc$Leaf_Type)
#Faith's PD
range(AC_16S_fpd_map_nc$faith_pd)
#8.76-51.66
#summary stats
AC_16S_fpd_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(faith_pd, type = "mean_sd")
#check for outliers
AC_16S_fpd_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(faith_pd)
#no outliers
#visualize
ggboxplot(AC_16S_fpd_map_nc, x = "Reach", y = "faith_pd",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_fpd_map_nc %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(faith_pd)
#normal
fpd.lm_nc<- lm(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map_nc)
ggqqplot(residuals(fpd.lm_nc))
shapiro_test(residuals(fpd.lm_nc))
#residuals fall aprox in qqplot, and shapiro test not significant
#looks good, test homogeniety of variance
leveneTest(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map_nc)
#not significant, therefore assume homogeniety of variance
res.aov.fpd_nc<- AC_16S_fpd_map_nc %>% anova_test(faith_pd ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.fpd_nc
#leaf type, time point, and leaf typextime point significant, so group the data by leaf type and time, and fit  anova
AC_16S_fpd_map_nc %>% group_by(Leaf_Type) %>%
  anova_test(faith_pd~ Days_Exposure, error = fpd.lm_nc)
AC_16S_fpd_map_nc%>% group_by(Time_Point_cat) %>%
  anova_test(faith_pd ~ Leaf_Type, error = fpd.lm_nc)
# Pairwise comparisons
AC_16S_fpd_map_nc %>% group_by(Leaf_Type) %>% emmeans_test(faith_pd ~ Time_Point_cat, 
                                                        p.adjust.method = "bonferroni",
                                                        model=fpd.lm_nc)
AC_16S_fpd_map_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=fpd.lm_nc)
#emm for non interactions
AC_16S_fpd_map_nc%>% emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni")
#(Ash=Buckthorn)>Oak
#See summary statistics for significant groups
AC_16S_fpd_map_nc %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(faith_pd, type = "mean_se")

#visualize reach, leaf type and time
fpdsummary_nc<-summarySE(AC_16S_fpd_map_nc, measurevar=c("faith_pd"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 1] <- 8
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 2] <- 41
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 3] <- 68
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 4] <- 98
fpdsummary_nc$Reach<-revalue(fpdsummary_nc$Reach, c("US"="Upstream"))
fpdsummary_nc$Reach<-revalue(fpdsummary_nc$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(fpdsummary_nc, aes(x=Time_Point, y=faith_pd, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=faith_pd-se, ymax=faith_pd+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Phylogenetic Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#Chao1
#Upload chao1
AC_16S_ch_map_nc<-subset(AC_16S_ch_map, Leaf_Type!="Cotton")
AC_16S_ch_map_nc$Leaf_Type<-droplevels(AC_16S_ch_map_nc$Leaf_Type)
#Chao1
range(AC_16S_ch_map_nc$chao1)
#40-903.9
#summary stats
AC_16S_ch_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(chao1, type = "mean_sd")
#check for outliers
AC_16S_ch_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(chao1)
#no outliers
#visualize
ggboxplot(AC_16S_ch_map_nc, x = "Reach", y = "chao1",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_ch_map_nc %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(chao1)
#normal
ch.lm_nc<- lm(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map_nc)
ggqqplot(residuals(ch.lm_nc))
shapiro_test(residuals(ch.lm_nc))
#residuals fall aprox in qqplot, and shapiro test significant
#test homogeniety of variance
leveneTest(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map_nc)
#not significant, therefore assume homogeniety of variance
res.aov.ch_nc<- AC_16S_ch_map_nc %>% anova_test(chao1 ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.ch_nc
#leaf type and time point significant, so group the data by leaf type and time, and fit  anova
#emm for non interactions
AC_16S_ch_map_nc%>% emmeans_test(chao1 ~ Leaf_Type, p.adjust.method = "bonferroni")
#Ash=buckthorn>Oak
#See summary statistics for significant groups
AC_16S_ch_map_nc %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(chao1, type = "mean_se")

#visualize reach, leaf type and time
chsummary_nc<-summarySE(AC_16S_ch_map_nc, measurevar=c("chao1"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 1] <- 8
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 2] <- 41
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 3] <- 68
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 4] <- 98
chsummary_nc$Reach<-revalue(chsummary_nc$Reach, c("US"="Upstream"))
chsummary_nc$Reach<-revalue(chsummary_nc$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(chsummary_nc, aes(x=Time_Point, y=chao1, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=chao1-se, ymax=chao1+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Richness (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#Facetted graph for SFS 2021 poster
chsummary_nc$facet<-rep("Mean Microbial Richness (± SE)", 42)
names(chsummary_nc)[names(chsummary_nc) == "chao1"] <- "Measurement"
MRsummary_nc$facet<-rep("Mean Macroinvertebrate Richness (± SE)", 45)
names(MRsummary_nc)[names(MRsummary_nc) == "Richness"] <- "Measurement"
AFDMsummary_nc$facet<-rep("Mean %AFDM Remaining (± SE)", 45)
names(AFDMsummary_nc)[names(AFDMsummary_nc) == "percAFDMremain"] <- "Measurement"
poster<-rbind(chsummary_nc,MRsummary_nc,AFDMsummary_nc)
poster$facet<- factor(poster$facet, levels = c("Mean %AFDM Remaining (± SE)", "Mean Microbial Richness (± SE)", "Mean Macroinvertebrate Richness (± SE)"))
ggplot(poster, aes(x=Time_Point, y=Measurement, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=Measurement-se, ymax=Measurement+se), width=7) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=14),axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=12),
        legend.title=element_text(size=14),legend.text = element_text(size=12),
        strip.text.x = element_text(size = 14),strip.text.y=element_text(size=10),legend.position="bottom") +
  facet_grid(facet~Reach,scales="free_y", labeller = label_wrap_gen(width=17))


###############################
#Microbes
############################

#Use QIIME2 for sequence processing

#Upload unifrac distance table
AC_16S_uni<-read.csv("~/Documents/Research/Field_Experiment/Microbes/AC_WUnifrac.csv", header = T, check.names=F)

#Combine distance matrix with environmental variables
LPMicM<-merge(AC_16S_uni, LPMetadata, by="Pack_ID")
row.names(LPMicM)<-LPMicM[,1]
LPMicM<-LPMicM[,-c(1)]
names(LPMicM)
AC_16S_uni<-as.matrix(LPMicM[,c(1:138)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
LPMicM_env<-LPMicM[,c(139:144)]
LPMicM_env_R<-LPMicM_env$Reach
LPMicM_env_LT<-as.factor(LPMicM_env$Leaf_Type)

#UNI-Overall permanova with unifrac distances
adonis2(as.dist(AC_16S_uni) ~ Leaf_Type*Reach*Days_Exposure, data=LPMicM_env,
        permutations=999)
#Leaf type, reach, time, leafxtime and reach x time significant

#Visualize via nmds
AC_Mic_NMDS<-metaMDS(as.dist(AC_16S_uni))
#stress 0.1

#Stressplot macroinvertebrate Nmds
stressplot(AC_Mic_NMDS)

#NMDS plot for Reach
ordiplot(AC_Mic_NMDS, type="n")
with(AC_Mic_NMDS, points(AC_Mic_NMDS, display="sites", col=reach_col_vec[LPMicM_env_R], pch=19))
with(AC_Mic_NMDS, legend("topleft", legend=levels(LPMicM_env_R), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_R, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "US"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_R, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_R, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "DS"))

#NMDS plot for leaf type
ordiplot(AC_Mic_NMDS, type="n")
with(AC_Mic_NMDS, points(AC_Mic_NMDS, display="sites", col=leaftaxacolvec[LPMicM_env_LT], pch=19))
with(AC_Mic_NMDS, legend("topleft", legend=levels(LPMicM_env_LT), bty="n", col=leaftaxacolvec, pch=19, pt.bg=leaftaxacolvec))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Cotton"))
with(AC_Mic_NMDS, ordiellipse(AC_Mic_NMDS, LPMicM_env_LT, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

#Upload family level taxonomy
AC_16S_f<-read.csv("~/Documents/Research/Field_Experiment/Microbes/AC_f.csv", header=T, check.names=F)
#Format data frame so the family is row name
row.names(AC_16S_f)<-AC_16S_f[,1]
#Delete otu id column, now that otu id is row name
AC_16S_f$Pack_ID<-NULL
AC_16S_f_t<-t(AC_16S_f)
AC_16S_f_t<-data.frame(AC_16S_f_t, check.names = FALSE)
AC_16S_f_t$Pack_ID<-row.names(AC_16S_f_t)
#Combine with environmental variables
AC_16S_f_map <-merge(LPMetadata,AC_16S_f_t, by="Pack_ID")
names(AC_16S_f_map)
sort(colSums(AC_16S_f_t[,1:536]))
#most common family is CComamonadaceae 
stat.desc(AC_16S_f_t$`k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae`)
#269 +/- 14
AC_16S_f_map_R<-AC_16S_f_map$Reach
AC_16S_f_map_L<-AC_16S_f_map$Leaf_Type

#indicator species analysis for reach
AC_Mic_Com_R_indic<-data.frame(signassoc(AC_16S_f_t[,1:536], cluster=AC_16S_f_map_R,  mode=0, alternative = "two.sided",control = how(nperm=999)))
AC_Mic_Com_R_indic_sig<-subset(AC_Mic_Com_R_indic, psidak<=0.05)
#27 indicator families for watershed

#indicator species analysis for leaf type
AC_Mic_Com_L_indic<-data.frame(signassoc(AC_16S_f_t[,1:536], cluster=AC_16S_f_map_L,  mode=0, alternative = "two.sided",control = how(nperm=999)))
AC_Mic_Com_L_indic_sig<-subset(AC_Mic_Com_L_indic, psidak<=0.05)
#122 indicator families for leaf type

#Upload phylogenetic diversity
AC_16S_fpd<- read.csv("~/Documents/Research/Field_Experiment/Microbes/AC_fpd.csv", header=T)
#Combine with environmental variables
AC_16S_fpd_map<-merge(LPMetadata, AC_16S_fpd, by="Pack_ID")
AC_16S_fpd_map$Time_Point_cat<-as.factor(AC_16S_fpd_map$Time_Point)

#Faith's PD
range(AC_16S_fpd_map$faith_pd)
#8.76-54.38
#summary stats
AC_16S_fpd_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(faith_pd, type = "mean_sd")
#check for outliers
AC_16S_fpd_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(faith_pd)
#no outliers
#visualize
ggboxplot(AC_16S_fpd_map, x = "Reach", y = "faith_pd",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_fpd_map %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(faith_pd)
#normal
fpd.lm<- lm(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map)
ggqqplot(residuals(fpd.lm))
shapiro_test(residuals(fpd.lm))
#residuals fall aprox in qqplot, and shapiro test not significant
#looks good, test homogeniety of variance
leveneTest(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map)
#not significant, therefore assume homogeniety of variance
res.aov.fpd<- AC_16S_fpd_map %>% anova_test(faith_pd ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.fpd
#leaf type, time point, and leaf typextime point significant, so group the data by leaf type and time, and fit  anova
AC_16S_fpd_map %>% group_by(Leaf_Type) %>%
  anova_test(faith_pd~ Days_Exposure, error = fpd.lm)
AC_16S_fpd_map%>% group_by(Time_Point_cat) %>%
  anova_test(faith_pd ~ Leaf_Type, error = fpd.lm)
# Pairwise comparisons
AC_16S_fpd_map %>% group_by(Leaf_Type) %>% emmeans_test(faith_pd ~ Time_Point_cat, 
                                                        p.adjust.method = "bonferroni",
                                                        model=fpd.lm)
AC_16S_fpd_map %>% group_by(Time_Point_cat) %>%
  emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=fpd.lm)
#emm for non interactions
AC_16S_fpd_map%>% emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni")
#cotton greater than ash and buckthorn, ash and cotton greater than oak
#See summary statistics for significant groups
AC_16S_fpd_map_A<-subset(AC_16S_fpd_map, Leaf_Type=="Ash")
AC_16S_fpd_map_B<-subset(AC_16S_fpd_map, Leaf_Type=="Buckthorn")
AC_16S_fpd_map_C<-subset(AC_16S_fpd_map, Leaf_Type=="Cotton")
AC_16S_fpd_map_O<-subset(AC_16S_fpd_map, Leaf_Type=="Oak")
stat.desc(AC_16S_fpd_map_A$faith_pd)
stat.desc(AC_16S_fpd_map_B$faith_pd)
stat.desc(AC_16S_fpd_map_C$faith_pd)
stat.desc(AC_16S_fpd_map_O$faith_pd)

#visualize reach, leaf type and time
fpdsummary<-summarySE(AC_16S_fpd_map, measurevar=c("faith_pd"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
fpdsummary$Time_Point[fpdsummary$Time_Point == 1] <- 8
fpdsummary$Time_Point[fpdsummary$Time_Point == 2] <- 41
fpdsummary$Time_Point[fpdsummary$Time_Point == 3] <- 68
fpdsummary$Time_Point[fpdsummary$Time_Point == 4] <- 98
fpdsummary$Reach<-revalue(fpdsummary$Reach, c("US"="Upstream"))
fpdsummary$Reach<-revalue(fpdsummary$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(fpdsummary, aes(x=Time_Point, y=faith_pd, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=faith_pd-se, ymax=faith_pd+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Phylogenetic Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

fpdsummarynr<-summarySE(AC_16S_fpd_map, measurevar=c("faith_pd"), groupvars=c("Leaf_Type","Time_Point"), na.rm=TRUE)
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 1] <- 8
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 2] <- 41
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 3] <- 68
fpdsummarynr$Time_Point[fpdsummarynr$Time_Point == 4] <- 98
# plot on y fpd and x days, no facet
ggplot(fpdsummarynr, aes(x=Time_Point, y=faith_pd, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=faith_pd-se, ymax=faith_pd+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Phylogenetic Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18))

#Chao1
#Upload chao1
AC_16S_ch<- read.csv("~/Documents/Research/Field_Experiment/Microbes/AC_chao.csv", header=T)
#Combine with environmental variables
AC_16S_ch_map<-merge(LPMetadata, AC_16S_ch, by="Pack_ID")
AC_16S_ch_map$Time_Point_cat<-as.factor(AC_16S_ch_map$Time_Point)

#Chao1
range(AC_16S_ch_map$chao1)
#40-903.9
#summary stats
AC_16S_ch_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(chao1, type = "mean_sd")
#check for outliers
AC_16S_ch_map %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(chao1)
#no outliers
#visualize
ggboxplot(AC_16S_ch_map, x = "Reach", y = "chao1",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_ch_map %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(chao1)
#normal
ch.lm<- lm(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map)
ggqqplot(residuals(ch.lm))
shapiro_test(residuals(ch.lm))
#residuals fall aprox in qqplot, and shapiro test significant
#test homogeniety of variance
leveneTest(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map)
#not significant, therefore assume homogeniety of variance
res.aov.ch<- AC_16S_ch_map %>% anova_test(chao1 ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.ch
#leaf type, time point, and leaf typextime point significant, so group the data by leaf type and time, and fit  anova
AC_16S_ch_map %>% group_by(Leaf_Type) %>%
  anova_test(chao1~ Days_Exposure, error = fpd.lm)
AC_16S_ch_map%>% group_by(Time_Point_cat) %>%
  anova_test(chao1 ~ Leaf_Type, error = fpd.lm)
# Pairwise comparisons
AC_16S_ch_map %>% group_by(Leaf_Type) %>% emmeans_test(chao1 ~ Time_Point_cat, 
                                                       p.adjust.method = "bonferroni",
                                                       model=ch.lm)
AC_16S_ch_map %>% group_by(Time_Point_cat) %>%
  emmeans_test(chao1 ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=ch.lm)
#emm for non interactions
AC_16S_ch_map%>% emmeans_test(chao1 ~ Leaf_Type, p.adjust.method = "bonferroni")
#Cotton>Ash=buckthorn>Oak
#See summary statistics for significant groups
AC_16S_ch_map_A<-subset(AC_16S_ch_map, Leaf_Type=="Ash")
AC_16S_ch_map_B<-subset(AC_16S_ch_map, Leaf_Type=="Buckthorn")
AC_16S_ch_map_C<-subset(AC_16S_ch_map, Leaf_Type=="Cotton")
AC_16S_ch_map_O<-subset(AC_16S_ch_map, Leaf_Type=="Oak")
stat.desc(AC_16S_ch_map_A$chao1)
stat.desc(AC_16S_ch_map_B$chao1)
stat.desc(AC_16S_ch_map_C$chao1)
stat.desc(AC_16S_ch_map_O$chao1)

#visualize reach, leaf type and time
chsummary<-summarySE(AC_16S_ch_map, measurevar=c("chao1"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
chsummary$Time_Point[chsummary$Time_Point == 1] <- 8
chsummary$Time_Point[chsummary$Time_Point == 2] <- 41
chsummary$Time_Point[chsummary$Time_Point == 3] <- 68
chsummary$Time_Point[chsummary$Time_Point == 4] <- 98
chsummary$Reach<-revalue(chsummary$Reach, c("US"="Upstream"))
chsummary$Reach<-revalue(chsummary$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(chsummary, aes(x=Time_Point, y=chao1, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=chao1-se, ymax=chao1+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Richness (± SE)")+
  scale_color_manual(values=leaftaxacolvec,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#now move on to venn diagrams to find unique OTUs introduced in time 0
#create venn diagram for shared OTUs among leaf type on day 0
#upload OTU table
AC_OTU<- read.csv("~/Documents/MSU/Research/Field_Experiment/Microbes/AC_otu_table.csv", header=T,check.names=F)
names(AC_OTU)
#Format data frame so the family is row name
row.names(AC_OTU)<-AC_OTU[,1]
#Delete otu id column, now that otu id is row name
AC_OTU$Pack_ID<-NULL
AC_OTU_t<-t(AC_OTU)
AC_OTU_t<-data.frame(AC_OTU_t, check.names = FALSE)
AC_OTU_t$Pack_ID<-row.names(AC_OTU_t)
#Combine with environmental variables
AC_OTU_map<-merge(LPMetadata, AC_OTU_t, by="Pack_ID")
#subset into only time point 0
AC_OTU_map_0<-subset(AC_OTU_map, Time_Point==0)
#first list OTU names found in each group
names(AC_OTU_map_0)
AC_OTU_map_0_ag_venn<-aggregate(AC_OTU_map_0[8:ncol(AC_OTU_map_0)], 
                                by=list(Leaf_Type=AC_OTU_map_0$Leaf_Type),
                                FUN=sum)
str(AC_OTU_map_0_ag_venn)
row.names(AC_OTU_map_0_ag_venn)<-AC_OTU_map_0_ag_venn$Leaf_Type
AC_OTU_map_0_ag_venn$Leaf_Type<-NULL
row.names(AC_OTU_map_0_ag_venn)
Ash<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Ash",] > 0]
Buckthorn<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Buckthorn",] > 0]
Cotton<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Cotton",] > 0]
Oak<-colnames(AC_OTU_map_0_ag_venn)[AC_OTU_map_0_ag_venn["Oak",] > 0]
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")
quartz()
LTVen<-venn_diagram4(Ash, Buckthorn, Cotton, Oak,
                     "Ash", "Buckthorn", "Cotton", "Oak",
                     colors=leaftaxacolvec)
#Each leaf type has unique OTUs
#visualize without cotton
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
quartz()
LTVen_nc<-venn_diagram3(Ash,Buckthorn,Oak, "Ash",
                        "Buckthorn","Oak",colors=leaftaxacolvec_nc)
#each leaf type has unique OTUs

##################
#now do other microbial analyses without cotton
#################
LPMIcM_nc<-subset(LPMicM, Leaf_Type!="Cotton")
names(LPMIcM_nc)
LPMIcM_nc<- select(LPMIcM_nc, -contains("-C-"))
names(LPMIcM_nc)
AC_16S_uni_nc<-as.matrix(LPMIcM_nc[,c(1:106)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
LPMicM_env_nc<-LPMIcM_nc[,c(107:112)]
LPMicM_env_nc$Reach<-revalue(LPMicM_env_nc$Reach, c("US"="Upstream", "DS"="Downstream"))
LPMicM_env_R_nc<-LPMicM_env_nc$Reach
LPMicM_env_LT_nc<-droplevels(as.factor(LPMicM_env_nc$Leaf_Type))

#UNI-Overall permanova with unifrac distances
adonis(as.dist(AC_16S_uni_nc) ~ Leaf_Type*Reach*Days_Exposure, data=LPMicM_env_nc,
       permutations=999)
#Leaf type, reach, time, leafxtime and reach x time significant

#Visualize via nmds
AC_Mic_NMDS_nc<-metaMDS(as.dist(AC_16S_uni_nc))
#stress 0.1

#Stressplot microbial Nmds
stressplot(AC_Mic_NMDS_nc)

#NMDS plot for Reach
ordiplot(AC_Mic_NMDS_nc, type="n")
with(AC_Mic_NMDS_nc, points(AC_Mic_NMDS_nc, display="sites", col=reach_col_vec[LPMicM_env_R_nc], pch=19))
with(AC_Mic_NMDS_nc, legend("topleft", legend=levels(LPMicM_env_R_nc), bty="n", col=reach_col_vec, pch=19, pt.bg=reach_col_vec))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_R_nc, kind="se", conf=0.95, lwd=2, col="#1b9e77", show.groups = "Upstream"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_R_nc, kind="se", conf=0.95, lwd=2, col="#d95f02", show.groups = "Gap"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_R_nc, kind="se", conf=0.95, lwd=2, col="#7570b3", show.groups = "Downstream"))

#NMDS plot for leaf type
ordiplot(AC_Mic_NMDS_nc, type="n")
with(AC_Mic_NMDS_nc, points(AC_Mic_NMDS_nc, display="sites", col=leaftaxacolvec_nc[LPMicM_env_LT_nc], pch=19))
with(AC_Mic_NMDS_nc, legend("topleft", legend=levels(LPMicM_env_LT_nc), bty="n", col=leaftaxacolvec_nc, pch=19, pt.bg=leaftaxacolvec_nc))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_LT_nc, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Ash"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_LT_nc, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Buckthorn"))
with(AC_Mic_NMDS_nc, ordiellipse(AC_Mic_NMDS_nc, LPMicM_env_LT_nc, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Oak"))

#family level taxonomy
AC_16S_f_map_nc<-subset(AC_16S_f_map, Leaf_Type!="Cotton")
names(AC_16S_f_map_nc)
sort(colSums(AC_16S_f_map_nc[,8:543]))
#most common family is Comamonadaceae 
stat.desc(AC_16S_f_map_nc$`k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae`)
#261 +/- 18
AC_16S_f_map_R_nc<-AC_16S_f_map_nc$Reach
AC_16S_f_map_L_nc<-droplevels(AC_16S_f_map_nc$Leaf_Type)

#indicator species analysis for reach
AC_Mic_Com_R_indic_nc<-signassoc(AC_16S_f_map_nc[,8:543], cluster=AC_16S_f_map_R_nc,  mode=0, alternative = "two.sided",control = how(nperm=999))
AC_Mic_Com_R_indic_sig_nc<-subset(AC_Mic_Com_R_indic_nc, psidak<=0.05)
#23 indicator families for reach, 2 for gap
#k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f_
#k__Bacteria;p__Bacteroidetes;c__Sphingobacteriia;o__Sphingobacteriales;f__
AC_16S_f_map_nc%>%
  group_by(Reach) %>%
  get_summary_stats("k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__", type = "mean_se")
#lowered in gaps
AC_16S_f_map_nc%>%
  group_by(Reach) %>%
  get_summary_stats("k__Bacteria;p__Bacteroidetes;c__Sphingobacteriia;o__Sphingobacteriales;f__", type = "mean_se")
#lowered in gaps

#indicator species analysis for leaf type
AC_Mic_Com_L_indic_nc<-signassoc(AC_16S_f_map_nc[,8:543], cluster=AC_16S_f_map_L_nc,  mode=0, alternative = "two.sided",control = how(nperm=999))
AC_Mic_Com_L_indic_sig_nc<-subset(AC_Mic_Com_L_indic_nc, psidak<=0.05)
#48 indicator families for leaf type

#Upload phylogenetic diversity
AC_16S_fpd_map_nc<-subset(AC_16S_fpd_map, Leaf_Type!="Cotton")
AC_16S_fpd_map_nc$Leaf_Type<-droplevels(AC_16S_fpd_map_nc$Leaf_Type)
#Faith's PD
range(AC_16S_fpd_map_nc$faith_pd)
#8.76-51.66
#summary stats
AC_16S_fpd_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(faith_pd, type = "mean_sd")
#check for outliers
AC_16S_fpd_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(faith_pd)
#no outliers
#visualize
ggboxplot(AC_16S_fpd_map_nc, x = "Reach", y = "faith_pd",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_fpd_map_nc %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(faith_pd)
#normal
fpd.lm_nc<- lm(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map_nc)
ggqqplot(residuals(fpd.lm_nc))
shapiro_test(residuals(fpd.lm_nc))
#residuals fall aprox in qqplot, and shapiro test not significant
#looks good, test homogeniety of variance
leveneTest(faith_pd ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_fpd_map_nc)
#not significant, therefore assume homogeniety of variance
res.aov.fpd_nc<- AC_16S_fpd_map_nc %>% anova_test(faith_pd ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.fpd_nc
#leaf type, time point, and leaf typextime point significant, so group the data by leaf type and time, and fit  anova
AC_16S_fpd_map_nc %>% group_by(Leaf_Type) %>%
  anova_test(faith_pd~ Days_Exposure, error = fpd.lm_nc)
AC_16S_fpd_map_nc%>% group_by(Time_Point_cat) %>%
  anova_test(faith_pd ~ Leaf_Type, error = fpd.lm_nc)
# Pairwise comparisons
AC_16S_fpd_map_nc %>% group_by(Leaf_Type) %>% emmeans_test(faith_pd ~ Time_Point_cat, 
                                                           p.adjust.method = "bonferroni",
                                                           model=fpd.lm_nc)
AC_16S_fpd_map_nc %>% group_by(Time_Point_cat) %>%
  emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni", 
               model=fpd.lm_nc)
#emm for non interactions
AC_16S_fpd_map_nc%>% emmeans_test(faith_pd ~ Leaf_Type, p.adjust.method = "bonferroni")
#(Ash=Buckthorn)>Oak
#See summary statistics for significant groups
AC_16S_fpd_map_nc %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(faith_pd, type = "mean_se")

#visualize reach, leaf type and time
fpdsummary_nc<-summarySE(AC_16S_fpd_map_nc, measurevar=c("faith_pd"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 1] <- 8
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 2] <- 41
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 3] <- 68
fpdsummary_nc$Time_Point[fpdsummary_nc$Time_Point == 4] <- 98
fpdsummary_nc$Reach<-revalue(fpdsummary_nc$Reach, c("US"="Upstream"))
fpdsummary_nc$Reach<-revalue(fpdsummary_nc$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(fpdsummary_nc, aes(x=Time_Point, y=faith_pd, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=faith_pd-se, ymax=faith_pd+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Phylogenetic Diversity (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#Chao1
#Upload chao1
AC_16S_ch_map_nc<-subset(AC_16S_ch_map, Leaf_Type!="Cotton")
AC_16S_ch_map_nc$Leaf_Type<-droplevels(AC_16S_ch_map_nc$Leaf_Type)
#Chao1
range(AC_16S_ch_map_nc$chao1)
#40-903.9
#summary stats
AC_16S_ch_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  get_summary_stats(chao1, type = "mean_sd")
#check for outliers
AC_16S_ch_map_nc %>%
  group_by(Reach, Leaf_Type, Time_Point_cat) %>%
  identify_outliers(chao1)
#no outliers
#visualize
ggboxplot(AC_16S_ch_map_nc, x = "Reach", y = "chao1",
          color = "Time_Point_cat", palette = "jco",
          facet.by = "Leaf_Type", short.panel.labs = FALSE)
#Check model assumptions
#normality, not enough replicates to do all three at same time
AC_16S_ch_map_nc %>%
  group_by(Reach, Leaf_Type) %>%
  shapiro_test(chao1)
#normal
ch.lm_nc<- lm(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map_nc)
ggqqplot(residuals(ch.lm_nc))
shapiro_test(residuals(ch.lm_nc))
#residuals fall aprox in qqplot, and shapiro test significant
#test homogeniety of variance
leveneTest(chao1 ~ Reach*Leaf_Type*Time_Point_cat, data = AC_16S_ch_map_nc)
#not significant, therefore assume homogeniety of variance
res.aov.ch_nc<- AC_16S_ch_map_nc %>% anova_test(chao1 ~ Reach*Leaf_Type*Time_Point_cat)
res.aov.ch_nc
#leaf type and time point significant, so group the data by leaf type and time, and fit  anova
#emm for non interactions
AC_16S_ch_map_nc%>% emmeans_test(chao1 ~ Leaf_Type, p.adjust.method = "bonferroni")
#Ash=buckthorn>Oak
#See summary statistics for significant groups
AC_16S_ch_map_nc %>%
  group_by(Leaf_Type) %>%
  get_summary_stats(chao1, type = "mean_se")

#visualize reach, leaf type and time
chsummary_nc<-summarySE(AC_16S_ch_map_nc, measurevar=c("chao1"), groupvars=c("Leaf_Type","Time_Point","Reach"), na.rm=TRUE)
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 1] <- 8
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 2] <- 41
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 3] <- 68
chsummary_nc$Time_Point[chsummary_nc$Time_Point == 4] <- 98
chsummary_nc$Reach<-revalue(chsummary_nc$Reach, c("US"="Upstream"))
chsummary_nc$Reach<-revalue(chsummary_nc$Reach, c("DS"="Downstream"))

#plot on y fpd and x axis days exposure
ggplot(chsummary_nc, aes(x=Time_Point, y=chao1, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=chao1-se, ymax=chao1+se), width=1) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  ylab("Mean Microbial Richness (± SE)")+
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18)) +
  facet_wrap(.~Reach)

#Facetted graph for SFS 2021 poster
chsummary_nc$facet<-rep("Mean Microbial Richness (± SE)", 42)
names(chsummary_nc)[names(chsummary_nc) == "chao1"] <- "Measurement"
MRsummary_nc$facet<-rep("Mean Macroinvertebrate Richness (± SE)", 45)
names(MRsummary_nc)[names(MRsummary_nc) == "Richness"] <- "Measurement"
AFDMsummary_nc$facet<-rep("Mean %AFDM Remaining (± SE)", 45)
names(AFDMsummary_nc)[names(AFDMsummary_nc) == "percAFDMremain"] <- "Measurement"
poster<-rbind(chsummary_nc,MRsummary_nc,AFDMsummary_nc)
poster$facet<- factor(poster$facet, levels = c("Mean %AFDM Remaining (± SE)", "Mean Microbial Richness (± SE)", "Mean Macroinvertebrate Richness (± SE)"))
ggplot(poster, aes(x=Time_Point, y=Measurement, color=Leaf_Type)) +
  geom_line(aes(group=Leaf_Type), size=1.5)+
  geom_errorbar(aes(ymin=Measurement-se, ymax=Measurement+se), width=7) +
  geom_point(size=1.5) +
  xlab("Days of Exposure") +
  scale_color_manual(values=leaftaxacolvec_nc,name = "Leaf Type") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=14),axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=12),
        legend.title=element_text(size=14),legend.text = element_text(size=12),
        strip.text.x = element_text(size = 14),strip.text.y=element_text(size=10),legend.position="bottom") +
  facet_grid(facet~Reach,scales="free_y", labeller = label_wrap_gen(width=17))

