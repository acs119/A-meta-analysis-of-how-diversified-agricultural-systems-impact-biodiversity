##Library###
library("magrittr") # for piping (%>%)
library(dplyr) # for data manipulation
library(stringr) # string manipulation
library(tidyr)
library(reshape2)
library(Matrix)
library(data.table)
library(metafor) #meta-analysis##
library(Formula)
library(ggplot2)
library(histogram)
library(funModeling)
library(frequency)
library(forestplot)
library(rlang)
library(ggpubr)


###############################################################################################################################
#Import database used for Sánchez et al "A meta-analysis of how diversified agricultural systems impact biodiversity"
# Access link:  https://doi.org/10.7910/DVN/NBPIPW
data-> read.csv("database.csv", header = TRUE,  sep = ",")

##SUBSET DATA ABUNDANCE 
abundance_logRR <- data %>% subset(B_measure == "Abundance")
abundance_logRR$ES_ID <- as.numeric(1:nrow(abundance_logRR)) #add a new column with the effect size ID number

hist(abundance_logRR$LRR) ##Frequency of ES##

##SUBSET DATA SPECIES RICHNESS
richness_logRR <- data %>% subset(B_measure == "Species Richness")
richness_logRR$ES_ID <- as.numeric(1:nrow(richness_logRR)) #add a new column with the effect size ID number
hist(richness_logRR$LRR) ##Frequency of ES##

# Number of effect sizes after filtering LRRvar = 0
nrow(abundance_logRR) #number of effect sizes after the LRRvar filter
nrow(richness_logRR) #number of effect sizes after the LRRvar filter

length(sort(unique(abundance_logRR$ID))) #Number of articles
length(sort(unique(richness_logRR$ID))) #Number of articles


####################################################################################################################

#################################### META-ANALYSIS RESULTS #################################################################################

##########################################################################################################################3
#####Equation: Cheung (2014) Formula to calculate the estimate sampling variance (formula 14)####
#b= LRR_var
estimated.sampling.variance.func <- function (b) {  
  result<- ((length(b)-1) * sum(1/b))/ (((sum(1/b))^2)-(sum(1/(b^2))))
  return(result)
}

######################################################################################################################
############################## A B U N D A N C E ######################################################################
#######################################################################################################################

########---- INTERCEPT ONLY MODEL -------###########
# Estimate the overall mean effect size by fitting an intercept-only model.
abun.overall <- rma.mv(y= LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                       tdist= TRUE, data=abundance_logRR, method="REML")
summary(abun.overall, digits=3)

# Heterogeneity of within-study variance (level 2)###
abun.modelnovar2 <- rma.mv(y=LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           sigma2=c(0,NA), tdist=TRUE, data=abundance_logRR, method="REML", verbose=TRUE)
summary(abun.modelnovar2)

# Perform a likelihood-ratio-test to determine the significance of the within-study variance (level2).
anova(abun.overall,abun.modelnovar2) #ABUNDANCE

# Heterogeneity of between-study variance (level 3)
abun.modelnovar3 <- rma.mv(y=LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID),
                           sigma2=c(NA,0), tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.modelnovar3)

# Perform a likelihood-ratio-test to determine the significance of the between-study variance.
anova(abun.overall,abun.modelnovar3) #ABUNDANCE

# The distribution of the variance over the three levels of the meta-analytic model
abun.estimated.sampling.variance<- estimated.sampling.variance.func(abundance_logRR$LRR_var)

## Each of the three variance components (I2_1, I2_2, I2_3) is divided by the total amount of variance
# Sampling variance (Amount of variance at level 1)
((abun.estimated.sampling.variance)/(abun.overall$sigma2[1]+abun.overall$sigma2[2]+abun.estimated.sampling.variance))*100

# Within-study variance (Amount of variance at level 2)
((abun.overall$sigma2[1]) / (abun.overall$sigma2[1] + abun.overall$sigma2[2] + abun.estimated.sampling.variance))*100

# Between-study variance (Amount of variance at level 3)
((abun.overall$sigma2[2]) / (abun.overall$sigma2[1] + abun.overall$sigma2[2] + abun.estimated.sampling.variance))*100

getMKLthreads()

# SENSITIVITY ANALYSIS
## COOK´S DISTANCE 
system.time(abun.overall_cook <- cooks.distance(abun.overall, reestimate=FALSE))
abun.overall_cook_df<- as.data.frame(abun.overall_cook)
abun.overall_cook_df$ES_ID <- as.numeric(1:nrow(abun.overall_cook_df))

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1)
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(abun.overall))))
plot(x=abun.overall_cook_df$ES_ID, y= abun.overall_cook_df$abun.overall_cook, xlim=c(0,2000),
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="A) Database: Abundance. Intercept only model")
lines(x=abun.overall_cook_df$ES_ID,y= abun.overall_cook_df$abun.overall_cook, col = "black")
abline(h = qchisq(0.5, df=1), lty=2, lwd=2, col="red") #0.4549364;

#Identify possible effect size outliers and exclude outliers from the abundance database
abundance_logRR_sensitivity.overall_cook<- left_join(abun.overall_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.overall_cook<0.01)

abundance_logRR_sensitivity.overall_cook_outliers<- left_join(abun.overall_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.overall_cook>0.01)
#plot outliers points in different color
abline(h = 0.01, lty=2, lwd=2, col="red") #0.4549364;
points(x=abundance_logRR_sensitivity.overall_cook_outliers$ES_ID, y=abundance_logRR_sensitivity.overall_cook_outliers$abun.overall_cook,
       col="red",pch=19)
text(abundance_logRR_sensitivity.overall_cook_outliers$ES_ID, abundance_logRR_sensitivity.overall_cook_outliers$abun.overall_cook, 
     labels = abundance_logRR_sensitivity.overall_cook_outliers$ES_ID, cex= 1, pos = c(2,4), col="red")

##Meta-analysis model without effect size outliers
abun.overall.sensitivity_cook <- rma.mv(y= LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                        tdist= TRUE, data=abundance_logRR_sensitivity.overall_cook, method="REML")
summary(abun.overall.sensitivity_cook, digits=3)
summary(abun.overall, digits=3)

########---- META-REGRESSION MODEL (MODERATOR: FUNCTIONAL GROUPS) -------###########
abun.FG <- rma.mv(y=LRR, V=LRR_var, mods = ~ Functional_group_recla, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.FG, digits=3)

## SENSITIVITY ANALYSIS 
# COOKS DISTANCE
system.time(abun.FG_cook <- cooks.distance(abun.FG, reestimate=FALSE))
abun.FG_cook_df<- as.data.frame(abun.FG_cook)
abun.FG_cook_df$ES_ID <- as.numeric(1:nrow(abun.FG_cook_df)) #add a new column with the effect size ID number

#Recommended cutoff for Cooks distance chi-squared distribution with df=(k+1)  
#or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(abun.FG))))
plot(x=abun.FG_cook_df$ES_ID, y= abun.FG_cook_df$abun.FG_cook, xlim=c(0,2000),
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="B) Database: Abundance. Meta-regression model (Moderator=Functional groups)")
lines(x=abun.FG_cook_df$ES_ID,y= abun.FG_cook_df$abun.FG_cook, col = "black")
abline(h = (qchisq(0.5, df=(length(coef(abun.FG))))), lty=2, lwd=2, col="red") #5.348121

#Identify possible effect size outliers and exclude outliers from the abundance database
abundance_logRR_sensitivity.FG_cook<- left_join(abun.FG_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.FG_cook<0.3)

abundance_logRR_sensitivity.FG_cook_outliers<- left_join(abun.FG_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.FG_cook>0.3)
#plot outliers points in different color
abline(h = 0.3, lty=2, lwd=2, col="red") #0.3;
points(x=abundance_logRR_sensitivity.FG_cook_outliers$ES_ID, y=abundance_logRR_sensitivity.FG_cook_outliers$abun.FG_cook, col="red",pch=19)
text(abundance_logRR_sensitivity.FG_cook_outliers$ES_ID, abundance_logRR_sensitivity.FG_cook_outliers$abun.FG_cook, 
     labels = abundance_logRR_sensitivity.FG_cook_outliers$ES_ID, cex= 1, pos = c(2,4), col="red")

# Meta-analysis model without effect size outliers
abun.FG.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~ Functional_group_recla, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                   tdist=TRUE, data=abundance_logRR_sensitivity.FG_cook, method="REML")
summary(abun.FG.sensitivity_cook, digits=3)

#------- PLOT RESULTS INTERCEPT ONLY MODEL AND META-REGRESSION (MODERATORS = Functional groups)
abun.FG.out.int<-rma.mv(y=LRR, V=LRR_var, mods = ~ (Functional_group_recla)-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.FG.out.int, digits=3)

# PLOT EFFECT SIZES ABUNDANCE AND FUNCTIONAL GROUPS
abun.FG_studies_ES<- abundance_logRR %>%
  group_by(Functional_group_recla)%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies", "Others", "Pests", "Pollinators"))

abun.FG.comb<- coef(summary(abun.FG.out.int))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies", "Others", "Pests", "Pollinators"),
         FG = as.factor(FG))%>%
  right_join(y=abun.FG_studies_ES, by="FG")

abun.overall_studies_ES <- abundance_logRR%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Overall"))

abun.overall.comb<- coef(summary(abun.overall))%>%
  mutate(Functional_group_recla = c("Summary"),
         FG = c("Overall"),
         FG= as.factor(FG))%>%
  right_join(y=abun.overall_studies_ES, by="FG")

abun.comb<- rbind(abun.FG.comb, abun.overall.comb)%>%
  rename(ES_LRR = estimate, 
         SE = se)%>%
  mutate(ES_percent = (100*(exp(ES_LRR)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         ES_LRR = round(ES_LRR, digits = 2),
         SE = round(SE, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =0),
         ci.lb_percent = round(ci.lb_percent, digits = 0),
         ci.ub_percent = round(ci.ub_percent, digits = 0))%>%
  select("FG", "ES_LRR", "SE","ci.lb","ci.ub","tval" ,"pval", "n_studies","n_effectsizes","ES_percent",
         "ci.lb_percent", "ci.ub_percent")

abun.comb.graph<- abun.comb%>%
  add_row(ES_LRR = c(NA), SE = c(NA), tval= c(NA), pval= c(NA), ci.lb= c(NA), ci.ub= c(NA),
          FG= c(NA), ES_percent= c(NA), ci.lb_percent = c(NA), ci.ub_percent= c(NA), .before = 1)%>%
  mutate(ID = row_number())

abun.tabletext<-abun.comb%>%
  mutate(ES_LRR = paste(ES_LRR, " ","(", ES_percent, ")", sep=""),
         ci.lb = paste(ci.lb, " ","(", ci.lb_percent, ")", sep=""),
         ci.ub = paste(ci.ub, " ","(", ci.ub_percent, ")", sep=""),
         SE = as.character(SE),
         tval = as.character(tval),
         pval = as.character(pval),
         n_studies = as.character(n_studies),
         n_effectsizes = as.character(n_effectsizes),
         pval= if_else(pval == "0", "<0.001", pval),
         pval = if_else(FG == "Pests"| FG == "Decomposers"|
                          FG =="Pollinators", paste(pval, "**", sep= ""),pval))%>%
  select(FG,ES_LRR, SE ,ci.lb, ci.ub, tval, pval,n_studies, n_effectsizes)%>%
  add_row(FG= c("Functional groups"),ES_LRR = c("Mean ES"), SE = c("SE"), tval= c("t-value"), 
          pval= c("p-value"), ci.lb= c("LC"), ci.ub= c("UC"), n_studies= c("#Articles"), 
          n_effectsizes =c("#ES"), .before = 1)%>%
  mutate(FG = c("Functional groups","Autotrophs","Decomposers", "Natural enemies",
                "Others", "Pests", "Pollinators",  "Overall"))
abun.tabletext

# Forest plot figure_2
figure_2<-forestplot(abun.tabletext, mean= abun.comb.graph$ES_percent, lower= abun.comb.graph$ci.lb_percent, 
                     upper= abun.comb.graph$ci.ub_percent, graph.pos=2,
                     new_page = TRUE, is.summary = c(TRUE,rep(FALSE,6), TRUE),  
                     xlab = "% Effect size (±95% CI)", clip = c(-30,200),
                     ci.vertices = TRUE, boxsize= c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4),
                     col= fpColors(box="#22211d", line="#22211d", summary="#686D35", zero ="gray50"),
                     txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", col= "black", cex = 1.3),
                                      ticks= gpar(fontfamily = "sans", cex = 1, col= "black"),
                                      xlab= gpar(fontfamily = "sans", cex = 1.2, col= "black")),
                     hrzl_lines = list("2" = gpar(lwd=2, col = "black"),
                                       "8" = gpar(lwd=2, col = "black")),
                     vertices=TRUE, align = c("l",rep("c",8)),
                     colgap = unit(3.5,"mm"),
                     xticks = c(-50, -25, 0, 25, 50,  100, 150, 200),
                     lineheight=unit(1.4,'cm'),
                     graphwidth = unit(4.2, "cm"))
figure_2

#------- PLOT RESULTS FROM SENSITIVITY ANALYSIS -> INTERCEPT ONLY MODEL AND META-REGRESSION (MODERATORS = Functional groups)
abun.FG.out.int.sensitivity_cook<-rma.mv(y=LRR, V=LRR_var, mods = ~ (Functional_group_recla)-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                         tdist=TRUE, data=abundance_logRR_sensitivity.FG_cook, method="REML")
summary(abun.FG.out.int.sensitivity_cook, digits=3)

#Plot data
abun.FG_studies_ES.sensitivity<- abundance_logRR_sensitivity.FG_cook %>%
  group_by(Functional_group_recla)%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies","Other", "Pests", "Pollinators"))

abun.FG.comb.sensitivity<- coef(summary(abun.FG.out.int.sensitivity_cook))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies", "Other", "Pests", "Pollinators"),
         FG = as.factor(FG))%>%
  right_join(y=abun.FG_studies_ES.sensitivity, by="FG")
abun.FG.comb.sensitivity

abun.overall_studies_ES.sensitivity <- abundance_logRR_sensitivity.overall_cook%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Overall"))

abun.overall.comb.sensitivity<- coef(summary(abun.overall.sensitivity_cook))%>%
  mutate(Functional_group_recla = c("Summary"),
         FG = c("Overall"),
         FG= as.factor(FG))%>%
  right_join(y=abun.overall_studies_ES.sensitivity, by="FG")
abun.overall.comb.sensitivity

abun.comb.sensitivity<- rbind(abun.FG.comb.sensitivity,abun.overall.comb.sensitivity)%>%
  rename(ES_LRR = estimate, 
         SE = se)%>%
  mutate(ES_percent = (100*(exp(ES_LRR)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         ES_LRR = round(ES_LRR, digits = 2),
         SE = round(SE, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =0),
         ci.lb_percent = round(ci.lb_percent, digits = 0),
         ci.ub_percent = round(ci.ub_percent, digits = 0))%>%
  select("FG", "ES_LRR", "SE","ci.lb","ci.ub","tval" ,"pval", "n_studies","n_effectsizes","ES_percent",
         "ci.lb_percent", "ci.ub_percent")

abun.comb.graph.sensitivity<- abun.comb.sensitivity%>%
  add_row(ES_LRR = c(NA), SE = c(NA), tval= c(NA), pval= c(NA), ci.lb= c(NA), ci.ub= c(NA),
          FG= c(NA), ES_percent= c(NA), ci.lb_percent = c(NA), ci.ub_percent= c(NA), .before = 1)%>%
  mutate(ID = row_number())

abun.tabletext.sensitivity<-abun.comb.sensitivity%>%
  mutate(ES_LRR = paste(ES_LRR, " ","(", ES_percent, ")", sep=""),
         ci.lb = paste(ci.lb, " ","(", ci.lb_percent, ")", sep=""),
         ci.ub = paste(ci.ub, " ","(", ci.ub_percent, ")", sep=""),
         SE = as.character(SE),
         tval = as.character(tval),
         pval = as.character(pval),
         n_studies = as.character(n_studies),
         n_effectsizes = as.character(n_effectsizes),
         pval= if_else(pval == "0", "<0.001", pval),
         pval = if_else(FG =="Autotrophs"|FG == "Pests"| FG == "Decomposers"| FG=="Pollinators", paste(pval, "**", sep= ""),pval))%>%
  select(FG,ES_LRR, SE ,ci.lb, ci.ub, tval, pval,n_studies, n_effectsizes)%>%
  add_row(FG= c("Functional groups"),ES_LRR = c("Mean ES"), SE = c("SE"), tval= c("t-value"), 
          pval= c("p-value"), ci.lb= c("LC"), ci.ub= c("UC"), n_studies= c("#Articles"), 
          n_effectsizes =c("#ES"), .before = 1)%>%
  mutate(FG = c("Functional groups","Autotrophs","Decomposers", "Natural enemies", "Other",
                "Pests", "Pollinators",  "Overall"))

##Forest plot
Figure_S5<- forestplot(abun.tabletext.sensitivity, mean= abun.comb.graph.sensitivity$ES_percent, lower= abun.comb.graph.sensitivity$ci.lb_percent, 
                       upper= abun.comb.graph.sensitivity$ci.ub_percent, graph.pos=2,
                       new_page = TRUE, is.summary = c(TRUE,rep(FALSE,6), TRUE),  
                       xlab = "% Effect size (±95% CI)", clip = c(-30,250),
                       ci.vertices = TRUE, boxsize= c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4),
                       col= fpColors(box="#22211d", line="#22211d", summary="#686D35", zero ="gray50"),
                       txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", col= "black", cex = 1.3),
                                        ticks= gpar(fontfamily = "sans", cex = 1, col= "black"),
                                        xlab= gpar(fontfamily = "sans", cex = 1.2, col= "black")),
                       hrzl_lines = list("2" = gpar(lwd=2, col = "black"),
                                         "8" = gpar(lwd=2, col = "black")),
                       vertices=TRUE, align = c("l",rep("c",8)),
                       colgap = unit(3.5,"mm"),
                       xticks = c(-60, -30, 0, 25,50, 150, 250),
                       lineheight=unit(1.4,'cm'),
                       graphwidth = unit(4.3, "cm"))

#---------------- META-REGRESSION (MODERATORS = LANDSCAPE METRICS)
# Determine the potential moderating effect of the assessed landscape metrics 

#------- MODERATOR = % OF NATURAL OR SEMI-NATURAL HABITATS
abun.natural_percentage <- rma.mv(y=LRR, V=LRR_var, mods = ~ natural_percentage_mean, random = 
                                    list(~ 1 | ES_ID, ~ 1 | ID),tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.natural_percentage, digits=3)

#Results meta-regression
coef(summary(abun.natural_percentage))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =1),ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

## SENSITIVITY ANALYSIS
# COOK´S DISTANCE 
system.time(abun.natural_percentage_cook <- cooks.distance(abun.natural_percentage, reestimate=FALSE))
abun.natural_percentage_cook_df<- as.data.frame(abun.natural_percentage_cook)
abun.natural_percentage_cook_df$ES_ID <- as.numeric(1:nrow(abun.natural_percentage_cook_df))

#Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1)
#or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(abun.natural_percentage))))
plot(x=abun.natural_percentage_cook_df$ES_ID, y= abun.natural_percentage_cook_df$abun.natural_percentage_cook, 
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="C) Database: Abundance. Meta-regression model (Moderator= % Natural habitats)")
lines(x=abun.natural_percentage_cook_df$ES_ID,y= abun.natural_percentage_cook_df$abun.natural_percentage_cook, col = "black")
abline(h = qchisq(0.5, df=1), lty=2, lwd=2, col="red") #1.386294;

#Identify possible effect size outliers and exclude outliers from the abundance database
abundance_logRR_sensitivity.natural_percentage_cook<- left_join(abun.natural_percentage_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.natural_percentage_cook<0.01)

abundance_logRR_sensitivity.natural_percentage_cook_outliers<- left_join(abun.natural_percentage_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.natural_percentage_cook>0.01)
#plot outliers points in different color
abline(h = 0.01, lty=2, lwd=2, col="red") #0.01
points(x=abundance_logRR_sensitivity.natural_percentage_cook_outliers$ES_ID, y=abundance_logRR_sensitivity.natural_percentage_cook_outliers$abun.natural_percentage_cook,
       col="red",pch=19)
text(abundance_logRR_sensitivity.natural_percentage_cook_outliers$ES_ID, abundance_logRR_sensitivity.natural_percentage_cook_outliers$abun.natural_percentage_cook, 
     labels = abundance_logRR_sensitivity.natural_percentage_cook_outliers$ES_ID, cex= 1, pos = c(2,4), col="red")

##Meta-analysis model without effect size outliers
abun.natural_percentage.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~ natural_percentage_mean, random = 
                                                     list(~ 1 | ES_ID, ~ 1 | ID),tdist=TRUE, data=abundance_logRR_sensitivity.natural_percentage_cook, method="REML")
summary(abun.natural_percentage.sensitivity_cook, digits=3)
summary(abun.natural_percentage, digits=3)

#Results meta-regression (after sensitivity analysis)
coef(summary(abun.natural_percentage.sensitivity))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =4),ci.lb_percent = round(ci.lb_percent, digits = 4),
         ci.ub_percent = round(ci.ub_percent, digits = 4))

#------- MODERATOR = % OF AGRICULTURAL AREAS
abun.arable_percentage <- rma.mv(y=LRR, V=LRR_var, mods = ~ arable_percentage_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                 tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.arable_percentage, digits=5)

# Results meta-regression
coef(summary(abun.arable_percentage))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =3),ci.lb_percent = round(ci.lb_percent, digits = 3),
         ci.ub_percent = round(ci.ub_percent, digits = 3))

## SENSITIVITY ANALYSIS
# COOK´S DISTANCE 
system.time(abun.arable_percentage_cook <- cooks.distance(abun.arable_percentage, reestimate=FALSE))
abun.arable_percentage_cook_df<- as.data.frame(abun.arable_percentage_cook)
abun.arable_percentage_cook_df$ES_ID <- as.numeric(1:nrow(abun.arable_percentage_cook_df))

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1) and alpha=0.5, 
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(abun.arable_percentage))))
plot(x=abun.arable_percentage_cook_df$ES_ID, y= abun.arable_percentage_cook_df$abun.arable_percentage_cook, 
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="D) Database: Abundance. Meta-regression model (Moderator= % Agricultural land)")
lines(x=abun.arable_percentage_cook_df$ES_ID,y= abun.arable_percentage_cook_df$abun.arable_percentage_cook, col = "black")
abline(h = qchisq(0.5, df=1), lty=2, lwd=2, col="red") #1.386294;

# Identify possible effect size outliers and exclude outliers from the abundance database
abundance_logRR_sensitivity.arable_percentage_cook<- left_join(abun.arable_percentage_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.arable_percentage_cook<0.02)

abundance_logRR_sensitivity.arable_percentage_cook_outliers<- left_join(abun.arable_percentage_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.arable_percentage_cook>0.02)
# plot outliers points in different color
abline(h = 0.02, lty=2, lwd=2, col="red") #0.02;
points(x=abundance_logRR_sensitivity.arable_percentage_cook_outliers$ES_ID, y=abundance_logRR_sensitivity.arable_percentage_cook_outliers$abun.arable_percentage_cook,
       col="red",pch=19)
text(abundance_logRR_sensitivity.arable_percentage_cook_outliers$ES_ID, abundance_logRR_sensitivity.arable_percentage_cook_outliers$abun.arable_percentage_cook, 
     labels = abundance_logRR_sensitivity.arable_percentage_cook_outliers$ES_ID, cex= 1, pos = c(2,4), col="red")

## Meta-analysis model without effect size outliers
abun.arable_percentage.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~ arable_percentage_mean, random = 
                                                    list(~ 1 | ES_ID, ~ 1 | ID),tdist=TRUE, data=abundance_logRR_sensitivity.arable_percentage_cook, method="REML")
summary(abun.arable_percentage.sensitivity_cook, digits=3)

# Results meta-regression after sensitivity analysis
coef(summary(abun.arable_percentage.sensitivity_cook))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =1),ci.lb_percent = round(ci.lb_percent, digits = 3),
         ci.ub_percent = round(ci.ub_percent, digits = 3))

#------- MODERATOR = DISTANCE TO NATURAL HABITATS
## META-REGRESSION MODEL: using min_distance_mean in meters as moderator:
abun.distance <- rma.mv(y=LRR, V=LRR_var, mods = ~min_distance_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.distance, digits=5)

# Calculate the residuals
rs.abun.distance<-rstandard.rma.mv(abun.distance, type="rstandard")

# Calculate predictor values for the model
preds.abun.distance <-predict(abun.distance,newmods=c(0:1001), addx=TRUE)
preds.abun.distance<-as.data.frame(preds.abun.distance)

## META-REGRESSION MODEL: using min_log_distance_mean (ln(distance+1)) as moderator :
abun.log_distance <- rma.mv(y=LRR, V=LRR_var, mods = ~min_log_distance_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.log_distance, digits=5)

# Calulate the residuals from the metaregression model using distance in log
rs.abun.log_distance<-rstandard.rma.mv(abun.log_distance, type="rstandard")

# Calculate predictor values for the model
preds.abun.log_distance <-predict(abun.log_distance,newmods=c(0:7), addx=TRUE)
preds.abun.log_distance<-as.data.frame(preds.abun.log_distance)
preds.abun.log_distance

# Plot residuals vs. distance in meters
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))#two plots in the same graph
plot(x=abundance_logRR$min_distance_mean, y=rs.abun.distance$resid,
     ylab="Standardized residuals", xlab= "Distance (m)", main = "A")
lines(x=preds.abun.distance$X.min_distance_mean ,y=preds.abun.distance$pred, col="#686D35", lwd=2)

# Plot residuals vs. distance (ln(distance+1))
plot(abundance_logRR$min_log_distance_mean, rs.abun.log_distance$resid,
     ylab="Standardized residuals", xlab= "ln(Distance + 1)", main = "B")
lines(x=preds.abun.log_distance$X.min_log_distance_mean ,y=preds.abun.log_distance$pred, col="#686D35", lwd=2)
mtext("Database: Abundance", outer = TRUE, cex = 1.5)

# SENSITIVITY ANALYSIS
## COOK´S DISTANCE 
system.time(abun.log_distance_cook <- cooks.distance(abun.log_distance, reestimate=FALSE))
abun.log_distance_cook_df<- as.data.frame(abun.log_distance_cook)
abun.log_distance_cook_df$ES_ID <- as.numeric(1:nrow(abun.log_distance_cook_df))
View(abun.log_distance_cook_df)

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1) and alpha=0.5, 
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(abun.log_distance))))
plot(x=abun.log_distance_cook_df$ES_ID, y= abun.log_distance_cook_df$abun.log_distance_cook, 
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="E) Database: Abundance. Meta-regression model (Moderator= Log Distance to Natural habitats)")
lines(x=abun.log_distance_cook_df$ES_ID,y= abun.log_distance_cook_df$abun.log_distance_cook, col = "black")
abline(h = qchisq(0.5, df=1), lty=2, lwd=2, col="red") #1.386294;

# Identify possible effect size outliers and exclude outliers from the abundance database
abundance_logRR_sensitivity.log_distance_cook<- left_join(abun.log_distance_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.log_distance_cook<0.02)

abundance_logRR_sensitivity.log_distance_cook_outliers<- left_join(abun.log_distance_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(abun.log_distance_cook>0.02)

# plot outliers points in different color
abline(h = 0.02, lty=2, lwd=2, col="red") #0.02;
points(x=abundance_logRR_sensitivity.log_distance_cook_outliers$ES_ID, y=abundance_logRR_sensitivity.log_distance_cook_outliers$abun.log_distance_cook,
       col="red",pch=19)
text(abundance_logRR_sensitivity.log_distance_cook_outliers$ES_ID, abundance_logRR_sensitivity.log_distance_cook_outliers$abun.log_distance_cook, 
     labels = abundance_logRR_sensitivity.log_distance_cook_outliers$ES_ID, cex= 1, pos = c(2,4), col="red")

## Meta-analysis model without effect size outliers
abun.log_distance.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~min_log_distance_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                             tdist=TRUE, data=abundance_logRR_sensitivity.log_distance_cook, method="REML")

summary(abun.log_distance.sensitivity_cook, digits=3)

# Results meta-regression after sensitivity analysis
coef(summary(abun.log_distance.sensitivity_cook))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =1),ci.lb_percent = round(ci.lb_percent, digits = 3),
         ci.ub_percent = round(ci.ub_percent, digits = 3))

##------ Plot RESULTS META-REGRESSION MODEL (MODERATORS = LANDSCAPE METRICS)
#A) %Natural habitats
preds.abun.natural.percentage_plot <-predict(abun.natural_percentage,newmods=c(0:100), addx=TRUE)
preds.abun.natural.percentage_plot<-as.data.frame(preds.abun.natural.percentage_plot)%>%
  mutate(pred_percentage = (100*(exp(pred)-1)))
class(preds.abun.natural.percentage_plot)
preds.abun.natural.percentage_plot
View(preds.abun.natural.percentage_plot)

summary(abun.natural_percentage)
abun.natural.percentage.summary<- coef(summary(abun.natural_percentage))%>%
  mutate(estimate= round(estimate, digits = 4),
         tval= round(tval, digits =2),
         pval= round(pval, digits =5),
         ci.lb = round(ci.lb, digits =4),
         ci.ub = round(ci.ub, digits =4),
         name = c("intercept", "slope"),
         x = c(1,4),
         y= c(2.5,3.7),
         label = paste("slope (%) = ", estimate, " (95% CI: ", ci.lb, ", ", ci.ub, ")",sep=""))%>%
  filter(name == "slope")
abun.natural.percentage.summary

figure_3_A<- ggplot(preds.abun.natural.percentage_plot,aes(x=X.natural_percentage_mean,y=pred))+
  geom_point(data=abundance_logRR,aes(x=natural_percentage_mean,y=LRR),colour="grey70",fill="white",
             position="jitter",alpha=0.4,shape= 20,size=2)+
  scale_y_discrete(limits = c(-8,-6,-4,-2,0,2,4,6))+
  geom_hline(yintercept = 0, colour = "black",)+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1.5, colour="#686D35")+
  geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), alpha=0.3)+
  labs(x=expression(bold(paste("Percentage of natural/semi-natural habitats"))),
       y= expression(bold(paste("Effect size (LRR)"))))+
  ggtitle("A)")+
  theme(
    axis.text.x = element_text(color="#22211d",size=8,  family = "sans"),
    axis.text.y = element_text(color="#22211d",size=8, family = "sans"),
    text = element_text(color = "#22211d", size =8, face = "bold", family = "sans"),
    plot.background = element_rect(fill = "White", color = "White"), 
    panel.background = element_rect(fill = "White", color = "White"), 
    axis.line = element_line(colour = "black"),
    plot.title = element_text(color="#22211d",size=9, family = "sans", vjust = 0.5, hjust=-0.1))
figure_3_A

# %Arable land
preds.abun.arable.percentage_plot <-predict(abun.arable_percentage,newmods=c(0:100), addx=TRUE)
preds.abun.arable.percentage_plot<-as.data.frame(preds.abun.arable.percentage_plot)%>%
  mutate(pred_percentage = (100*(exp(pred)-1)))
class(preds.abun.arable.percentage_plot)
preds.abun.arable.percentage_plot
View(preds.abun.arable.percentage_plot)

summary(abun.arable_percentage)
abun.arable.percentage.summary<- coef(summary(abun.arable_percentage))%>%
  mutate(estimate= round(estimate, digits = 4),
         tval= round(tval, digits =2),
         pval= round(pval, digits =5),
         ci.lb = round(ci.lb, digits =4),
         ci.ub = round(ci.ub, digits =4),
         name = c("intercept", "slope"),
         x = c(1,4),
         y= c(2.5,3.7),
         label = paste("slope (%) = ", estimate, " (95% CI: ", ci.lb, ", ", ci.ub, ")",sep=""))%>%
  filter(name == "slope")
abun.arable.percentage.summary

figure_3_B<- ggplot(preds.abun.arable.percentage_plot,aes(x=X.arable_percentage_mean,y=pred))+
  geom_point(data=abundance_logRR,aes(x=arable_percentage_mean,y=LRR),colour="grey70",fill="white",
             position="jitter",alpha=0.4,shape= 20,size=2)+
  scale_y_discrete(limits = c(-8,-6,-4,-2,0,2,4,6))+
  geom_hline(yintercept = 0, colour = "black",)+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1.5, colour="#686D35")+
  geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), alpha=0.3)+
  labs(x=expression(bold(paste("Percentage of agricultural land"))),
       y= expression(bold(paste("Effect size (LRR)"))))+
  ggtitle("B)")+
  theme(
    axis.text.x = element_text(color="#22211d",size=8,  family = "sans"),
    axis.text.y = element_text(color="#22211d",size=8, family = "sans"),
    text = element_text(color = "#22211d", size =8, face = "bold", family = "sans"),
    plot.background = element_rect(fill = "White", color = "White"), 
    panel.background = element_rect(fill = "White", color = "White"), 
    axis.line = element_line(colour = "black"),
    plot.title = element_text(color="#22211d",size=9, family = "sans", vjust = 0.5, hjust=-0.1))
figure_3_B

# log distance to natural or semi-natural habitats
preds.abun.log.distance_plot <-predict(abun.log_distance,newmods=c(0:8), addx=TRUE)
preds.abun.log.distance_plot<-as.data.frame(preds.abun.log.distance_plot)%>%
  mutate(pred_percentage = (100*(exp(pred)-1)))
class(preds.abun.log.distance_plot)
preds.abun.log.distance_plot
View(preds.abun.log.distance_plot)

summary(abun.log_distance)
abun.log.distance.summary<- coef(summary(abun.log_distance))%>%
  mutate(estimate= round(estimate, digits = 4),
         tval= round(tval, digits =2),
         pval= round(pval, digits =5),
         ci.lb = round(ci.lb, digits =4),
         ci.ub = round(ci.ub, digits =4),
         name = c("intercept", "slope"),
         x = c(1,4),
         y= c(2.5,3.7),
         label = paste("slope (%) = ", estimate, " (95% CI: ", ci.lb, ", ", ci.ub, ")",sep=""))%>%
  filter(name == "slope")
abun.log.distance.summary

figure_3_C<- ggplot(preds.abun.log.distance_plot,aes(x=X.min_log_distance_mean,y=pred))+
  geom_point(data=abundance_logRR,aes(x=min_log_distance_mean,y=LRR),colour="grey70",fill="white",
             position="jitter",alpha=0.4,shape= 20,size=2)+
  scale_y_discrete(limits = c(-8,-6,-4,-2,0,2,4,6))+
  scale_x_discrete(limits=c(0, 2,4,6,8),labels=c("0", "2", "4", "6", "8", ""))+
  geom_hline(yintercept = 0, colour = "black",)+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1.5, colour="#686D35")+
  geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), alpha=0.3)+
  labs(x=expression(bold(paste("Distance to natural or semi-natural habitats [ln(distance + 1)]"))),
       y= expression(bold(paste("Effect size (LRR)"))))+
  ggtitle("C)")+
  theme(
    axis.text.x = element_text(color="#22211d",size=8,  family = "sans"),
    axis.text.y = element_text(color="#22211d",size=8, family = "sans"),
    text = element_text(color = "#22211d", size =8, face = "bold", family = "sans"),
    plot.background = element_rect(fill = "White", color = "White"), 
    panel.background = element_rect(fill = "White", color = "White"), 
    axis.line = element_line(colour = "black"),
    plot.title = element_text(color="#22211d",size=9, family = "sans", vjust = 0.5, hjust=-0.1))
figure_3_C

ggpubr::ggarrange(figure_3_A,figure_3_B, figure_3_C, nrow = 3)

######################################################################################################################
############################## S P E C I E S  R I C H N E S S ######################################################################
#######################################################################################################################

########---- INTERCEPT ONLY MODEL -------###########
richness.overall <- rma.mv(y= LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist= TRUE, data=richness_logRR, method="REML")
summary(richness.overall, digits=3)

# Heterogeneity of within-study variance (level 2)
richness.modelnovar2 <- rma.mv(y=LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID), sigma2=c(0,NA), tdist=TRUE, 
                               data=richness_logRR, method="REML")
summary(richness.modelnovar2)

# Perform a likelihood-ratio-test to determine the significance of the within-study variance (level2).
anova(richness.overall,richness.modelnovar2)

# Heterogeneity of between-study variance (level 3)
richness.modelnovar3 <- rma.mv(y=LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID),
                               sigma2=c(NA,0), tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.modelnovar3)

# Perform a likelihood-ratio-test to determine the significance of the between-study variance.
anova(richness.overall,richness.modelnovar3) #SPECIES RICHNESS

# The distribution of the variance over the three levels of the meta-analytic model
richness.estimated.sampling.variance<- estimated.sampling.variance.func(richness_logRR$LRR_var)

###Each of the three variance components (I2_1, I2_2, I2_3) is divided by the total amount of variance:
#Sampling variance (Amount of variance at level 1)
((richness.estimated.sampling.variance)/(richness.overall$sigma2[1]+richness.overall$sigma2[2]+richness.estimated.sampling.variance))*100
#Within-study variance (Amount of variance at level 2)
((richness.overall$sigma2[1]) / (richness.overall$sigma2[1] + richness.overall$sigma2[2] + richness.estimated.sampling.variance))*100
#Between-study variance (Amount of variance at level 3)
((richness.overall$sigma2[2]) / (richness.overall$sigma2[1] + richness.overall$sigma2[2] + richness.estimated.sampling.variance))*100

# SENSITIVITY ANALYSIS 
## COOK DISTANCE
system.time(richness.overall_cook <- cooks.distance(richness.overall, reestimate=FALSE))
richness.overall_cook_df<- as.data.frame(richness.overall_cook)
richness.overall_cook_df$ES_ID <- as.numeric(1:nrow(richness.overall_cook_df))

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1) 
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(richness.overall))))
plot(x=richness.overall_cook_df$ES_ID, y= richness.overall_cook_df$richness.overall_cook, xlim=c(0,400),
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="F) Database: Species Richness. Intercept only model")
lines(x=richness.overall_cook_df$ES_ID,y= richness.overall_cook_df$richness.overall_cook, col = "black")
abline(h = qchisq(0.5, 1), lty=2, lwd=2, col="red") #0.4549364 

#Identify possible effect size outliers and exclude outliers from the abundance database
richness_logRR_sensitivity.overall_cook<- left_join(richness.overall_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(richness.overall_cook<0.02)

richness_logRR_sensitivity.overall_cook_outliers<- left_join(richness.overall_cook_df,abundance_logRR, by="ES_ID")%>%
  filter(richness.overall_cook>0.02)

#plot outliers points in different color
abline(h = 0.02, lty=2, lwd=2, col="red") #0.02;
points(x=richness_logRR_sensitivity.overall_cook_outliers$ES_ID, y=richness_logRR_sensitivity.overall_cook_outliers$richness.overall_cook,
       col="red",pch=19)
text(richness_logRR_sensitivity.overall_cook_outliers$ES_ID, richness_logRR_sensitivity.overall_cook_outliers$richness.overall_cook, 
     labels = richness_logRR_sensitivity.overall_cook_outliers$ES_ID, cex= 1, pos = c(2,4), col="red")

##Meta-analysis model without effect size outliers
richness.overall.sensitivity_cook <- rma.mv(y= LRR, V=LRR_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                            tdist= TRUE, data=richness_logRR_sensitivity.overall_cook, method="REML")
summary(richness.overall.sensitivity_cook, digits=3)

########---- META-REGRESSION MODEL (MODERATOR: FUNCTIONAL GROUPS) -------###########
richness.FG <- rma.mv(y=LRR, V=LRR_var, mods = ~ Functional_group_recla, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.FG, digits=3)

## SENSITIVITY ANALYSIS
# COOK DISTANCE
system.time(richness.FG_cook <- cooks.distance(richness.FG, reestimate=FALSE))
names(richness.FG_cook)
richness.FG_cook_df<- as.data.frame(richness.FG_cook)
richness.FG_cook_df$ES_ID <- as.numeric(1:nrow(richness.FG_cook_df)) #add a new column with the effect size ID number

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1)
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(richness.FG))))
plot(x=richness.FG_cook_df$ES_ID, y= richness.FG_cook_df$richness.FG_cook, xlim=c(0,400),
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="G) Database: Species Richness. Meta-regression model (Moderator=Functional groups)")
lines(x=richness.FG_cook_df$ES_ID,y= richness.FG_cook_df$richness.FG_cook, col = "black")
abline(h = qchisq(0.5, 6), lty=2, lwd=2, col="red") #5.348121

#Identify possible effect size outliers and exclude outliers from the abundance database
richness_logRR_sensitivity.FG_cook<- left_join(richness.FG_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.FG_cook<0.3)

richness_logRR_sensitivity.FG_cook_outlier<- left_join(richness.FG_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.FG_cook>0.3)

# plot outliers points in different color
abline(h = 0.3, lty=2, lwd=2, col="red") #0.3;
points(x=richness_logRR_sensitivity.FG_cook_outlier$ES_ID, y=richness_logRR_sensitivity.FG_cook_outlier$richness.FG_cook,
       col="red",pch=19)
text(richness_logRR_sensitivity.FG_cook_outlier$ES_ID, richness_logRR_sensitivity.FG_cook_outlier$richness.FG_cook, 
     labels = richness_logRR_sensitivity.FG_cook_outlier$ES_ID, cex= 1, pos = c(2,4), col="red")

## Meta-analysis model without effect size outliers
richness.FG.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~ Functional_group_recla, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                       tdist=TRUE, data=richness_logRR_sensitivity.FG_cook, method="REML")
summary(richness.FG.sensitivity_cook, digits=3)

#------- PLOT RESULTS INTERCEPT ONLY MODEL AND META-REGRESSION MODEL (MODERATORS = Functional groups)
richness.FG.out.int<-rma.mv(y=LRR, V=LRR_var, mods = ~ (Functional_group_recla)-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.FG.out.int, digits=3)

# Plot effect sizes Species richness and Functional groups
rich.FG_studies_ES<- richness_logRR %>%
  group_by(Functional_group_recla)%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies", "Others", "Pests", "Pollinators"))

rich.FG.comb<- coef(summary(richness.FG.out.int))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies", "Others", "Pests", "Pollinators"),
         FG = as.factor(FG))%>%
  right_join(y=rich.FG_studies_ES, by="FG")

rich.overall_studies_ES <- richness_logRR%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Overall"))

rich.overall.comb<- coef(summary(richness.overall))%>%
  mutate(Functional_group_recla = c("Summary"),
         FG = c("Overall"),
         FG= as.factor(FG))%>%
  right_join(y=rich.overall_studies_ES, by="FG")

rich.comb<- rbind(rich.FG.comb,rich.overall.comb)%>%
  rename(ES_LRR = estimate, 
         SE = se)%>%
  mutate(ES_percent = (100*(exp(ES_LRR)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         ES_LRR = round(ES_LRR, digits = 2),
         SE = round(SE, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =0),
         ci.lb_percent = round(ci.lb_percent, digits = 0),
         ci.ub_percent = round(ci.ub_percent, digits = 0))%>%
  select("FG", "ES_LRR", "SE","ci.lb","ci.ub","tval" ,"pval", "n_studies","n_effectsizes","ES_percent",
         "ci.lb_percent", "ci.ub_percent")

rich.comb.graph<- rich.comb%>%
  add_row(ES_LRR = c(NA), SE = c(NA), tval= c(NA), pval= c(NA), ci.lb= c(NA), ci.ub= c(NA),
          FG= c(NA), ES_percent= c(NA), ci.lb_percent = c(NA), ci.ub_percent= c(NA), .before = 1)%>%
  mutate(ID = row_number())
rich.comb.graph

rich.tabletext<-rich.comb%>%
  mutate(ES_LRR = paste(ES_LRR, " ","(", ES_percent, ")", sep=""),
         ci.lb = paste(ci.lb, " ","(", ci.lb_percent, ")", sep=""),
         ci.ub = paste(ci.ub, " ","(", ci.ub_percent, ")", sep=""),
         SE = as.character(SE),
         tval = as.character(tval),
         pval = as.character(pval),
         n_studies = as.character(n_studies),
         n_effectsizes = as.character(n_effectsizes),
         pval= if_else(pval == "0", "<0.001", pval),
         pval = if_else(FG == "Autotrophs"| FG == "Natural enemies"| FG == "Decomposers"| 
                          FG =="Pollinators"| FG== "Overall", paste(pval, "**", sep= ""),pval))%>%
  select(FG,ES_LRR, SE ,ci.lb, ci.ub, tval, pval,n_studies, n_effectsizes)%>%
  add_row(FG= c("Functional groups"),ES_LRR = c("Mean ES"), SE = c("SE"), tval= c("t-value"), 
          pval= c("p-value"), ci.lb= c("LC"), ci.ub= c("UC"), n_studies= c("#Articles"), 
          n_effectsizes =c("#ES"), .before = 1)%>%
  mutate(FG = c("Functional groups","Autotrophs","Decomposers", "Natural enemies",
                "Others", "Pests", "Pollinators",  "Overall"))
rich.tabletext

## Figure: Forest plot 
figure_4<-forestplot(rich.tabletext, mean= rich.comb.graph$ES_percent, lower= rich.comb.graph$ci.lb_percent, 
                     upper= rich.comb.graph$ci.ub_percent, graph.pos=2,
                     new_page = TRUE, is.summary = c(TRUE,rep(FALSE,6), TRUE),  
                     xlab = "% Effect size (±95% CI)", clip = c(-10,100),
                     ci.vertices = TRUE, boxsize= c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4),
                     col= fpColors(box="#22211d", line="#22211d", summary="#A63117", zero ="gray50"),
                     txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", col= "black", cex = 1.3),
                                      ticks= gpar(fontfamily = "sans", cex = 1, col= "black"),
                                      xlab= gpar(fontfamily = "sans", cex = 1.2, col= "black")),
                     hrzl_lines = list("2" = gpar(lwd=2, col = "black"),
                                       "8" = gpar(lwd=2, col = "black")),
                     vertices=TRUE, align = c("l",rep("c",8)),
                     colgap = unit(3.5,"mm"),
                     xticks = c(-25, 0, 25, 50,  75),
                     lineheight=unit(1.4,'cm'),
                     graphwidth = unit(4.3, "cm"))

#------- PLOT RESULTS FROM SENSITIVITY ANALYSIS -> INTERCEPT ONLY MODEL AND META-REGRESSION (MODERATORS = Functional groups)
richness.FG.out.int_sensitivity<-rma.mv(y=LRR, V=LRR_var, mods = ~ (Functional_group_recla)-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                        tdist=TRUE, data=richness_logRR_sensitivity.FG_cook, method="REML")
summary(richness.FG.out.int_sensitivity, digits=3)

#Plot effect sizes Species richness and Functional groups
rich.FG_studies_ES_sensitivity<- richness_logRR_sensitivity.FG_cook %>%
  group_by(Functional_group_recla)%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies", "Other", "Pests", "Pollinators"))

rich.FG.comb_sensitivity<- coef(summary(richness.FG.out.int_sensitivity))%>%
  mutate(FG = c("Autotrophs","Decomposers", "Natural enemies", "Other","Pests", "Pollinators"),
         FG = as.factor(FG))%>%
  right_join(y=rich.FG_studies_ES_sensitivity, by="FG")

rich.overall_studies_ES_sensitivity <- richness_logRR_sensitivity.overall_cook%>%
  summarise(n_studies = n_distinct(ID),
            n_effectsizes = n_distinct(ES_ID))%>%
  mutate(FG = c("Overall"))

rich.overall.comb_sensitivity<- coef(summary(richness.overall.sensitivity_cook))%>%
  mutate(Functional_group_recla = c("Summary"),
         FG = c("Overall"),
         FG= as.factor(FG))%>%
  right_join(y=rich.overall_studies_ES_sensitivity, by="FG")

rich.comb_sensitivity<- rbind(rich.FG.comb_sensitivity,rich.overall.comb_sensitivity)%>%
  rename(ES_LRR = estimate, 
         SE = se)%>%
  mutate(ES_percent = (100*(exp(ES_LRR)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         ES_LRR = round(ES_LRR, digits = 2),
         SE = round(SE, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =0),
         ci.lb_percent = round(ci.lb_percent, digits = 0),
         ci.ub_percent = round(ci.ub_percent, digits = 0))%>%
  select("FG", "ES_LRR", "SE","ci.lb","ci.ub","tval" ,"pval", "n_studies","n_effectsizes","ES_percent",
         "ci.lb_percent", "ci.ub_percent")

rich.comb.graph_sensitivity<- rich.comb_sensitivity%>%
  add_row(ES_LRR = c(NA), SE = c(NA), tval= c(NA), pval= c(NA), ci.lb= c(NA), ci.ub= c(NA),
          FG= c(NA), ES_percent= c(NA), ci.lb_percent = c(NA), ci.ub_percent= c(NA), .before = 1)%>%
  mutate(ID = row_number())
rich.comb.graph_sensitivity

rich.tabletext_sensitivity<-rich.comb_sensitivity%>%
  mutate(ES_LRR = paste(ES_LRR, " ","(", ES_percent, ")", sep=""),
         ci.lb = paste(ci.lb, " ","(", ci.lb_percent, ")", sep=""),
         ci.ub = paste(ci.ub, " ","(", ci.ub_percent, ")", sep=""),
         SE = as.character(SE),
         tval = as.character(tval),
         pval = as.character(pval),
         n_studies = as.character(n_studies),
         n_effectsizes = as.character(n_effectsizes),
         pval= if_else(pval == "0", "<0.001", pval),
         pval = if_else(FG == "Autotrophs"| FG == "Natural enemies"| FG=="Decomposers"|
                          FG =="Pollinators"| FG== "Overall", paste(pval, "**", sep= ""),pval))%>%
  select(FG,ES_LRR, SE ,ci.lb, ci.ub, tval, pval,n_studies, n_effectsizes)%>%
  add_row(FG= c("Functional groups"),ES_LRR = c("ES"), SE = c("SE"), tval= c("t-value"), 
          pval= c("p-value"), ci.lb= c("LC"), ci.ub= c("UC"), n_studies= c("#Articles"), 
          n_effectsizes =c("# ES"), .before = 1)%>%
  mutate(FG = c("Functional groups","Autotrophs","Decomposers", "Natural enemies", "Other",
                "Pests", "Pollinators",  "Overall"))
rich.tabletext_sensitivity

# Figure_S: forest plot 
Figure_S6<- forestplot(rich.tabletext_sensitivity, mean= rich.comb.graph_sensitivity$ES_percent, lower= rich.comb.graph_sensitivity$ci.lb_percent, 
                       upper= rich.comb.graph_sensitivity$ci.ub_percent, graph.pos=2,
                       new_page = TRUE, is.summary = c(TRUE,rep(FALSE,6), TRUE),  
                       xlab = "% Effect size (±95% CI)", clip = c(-10,100),
                       ci.vertices = TRUE, boxsize= c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4),
                       col= fpColors(box="#22211d", line="#22211d", summary="#A63117", zero ="gray50"),
                       txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", col= "black", cex = 1.3),
                                        ticks= gpar(fontfamily = "sans", cex = 1, col= "black"),
                                        xlab= gpar(fontfamily = "sans", cex = 1.2, col= "black")),
                       hrzl_lines = list("2" = gpar(lwd=2, col = "black"),
                                         "8" = gpar(lwd=2, col = "black")),
                       vertices=TRUE, align = c("l",rep("c",8)),
                       colgap = unit(3.5,"mm"),
                       xticks = c(-25, 0, 25, 50,  75),
                       lineheight=unit(1.4,'cm'),
                       graphwidth = unit(4.3, "cm"))

#---------------- META-REGRESSION (MODERATORS = LANDSCAPE METRICS)
# Determine the potential moderating effect of the assessed landscape metrics 

#------- MODERATOR = % OF NATURAL OR SEMI-NATURAL HABITATS
richness.natural_percentage <- rma.mv(y=LRR, V=LRR_var, mods = ~ natural_percentage_mean,random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                      tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.natural_percentage, digits=5)

# Results meta-regression
coef(summary(richness.natural_percentage))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =1),ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

## SENSITIVITY ANALYSIS
# COOK DISTANCE
system.time(richness.natural_percentage_cook <- cooks.distance(richness.natural_percentage, reestimate=FALSE))
richness.natural_percentage_cook_df<- as.data.frame(richness.natural_percentage_cook)
richness.natural_percentage_cook_df$ES_ID <- as.numeric(1:nrow(richness.natural_percentage_cook_df)) #add a new column with the effect size ID number

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1) and ???=0.5, 
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(richness.natural_percentage))))
plot(x=richness.natural_percentage_cook_df$ES_ID, y= richness.natural_percentage_cook_df$richness.natural_percentage_cook, 
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="H) Database: Species richness. Meta-regression model (Moderator= %Natural habitats)")
lines(x=richness.natural_percentage_cook_df$ES_ID,y= richness.natural_percentage_cook_df$richness.natural_percentage_cook, col = "black")
abline(h = qchisq(0.5, df=(length(coef(richness.natural_percentage)))), lty=2, lwd=2, col="red") #1.386294

# Identify possible effect size outliers and exclude outliers from the abundance database
richness_logRR_sensitivity.natural_percentage_cook<- left_join(richness.natural_percentage_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.natural_percentage_cook<0.02)

richness_logRR_sensitivity.natural_percentage_cook_outlier<- left_join(richness.natural_percentage_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.natural_percentage_cook>0.02)

# plot outliers points in different color
abline(h = 0.02, lty=2, lwd=2, col="red") #0.3;
points(x=richness_logRR_sensitivity.natural_percentage_cook_outlier$ES_ID, y=richness_logRR_sensitivity.natural_percentage_cook_outlier$richness.natural_percentage_cook,
       col="red",pch=19)
text(richness_logRR_sensitivity.natural_percentage_cook_outlier$ES_ID, richness_logRR_sensitivity.natural_percentage_cook_outlier$richness.natural_percentage_cook, 
     labels = richness_logRR_sensitivity.natural_percentage_cook_outlier$ES_ID, cex= 1, pos = c(2,4), col="red")

## Meta-analysis model without effect size outliers
richness.natural_percentage.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~ natural_percentage_mean, 
                                                       random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                                       tdist=TRUE, data=richness_logRR_sensitivity.natural_percentage_cook, method="REML")
summary(richness.natural_percentage.sensitivity_cook, digits=3)

# Results meta-regression
coef(summary(richness.natural_percentage.sensitivity_cook))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =1),ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

#------- MODERATOR = % OF ARABLE HABITATS
richness.arable_percentage <- rma.mv(y=LRR, V=LRR_var, mods = ~ arable_percentage_mean, random = list(~ 1 | ES_ID, ~ 1 | ID),
                                     tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.arable_percentage, digits=3)

# Results meta-regression
coef(summary(richness.arable_percentage))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =1),ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

## SENSITIVITY ANALYSIS
# COOK DISTANCE
system.time(richness.arable_percentage_cook <- cooks.distance(richness.arable_percentage, reestimate=FALSE))
names(richness.arable_percentage_cook)
richness.arable_percentage_cook_df<- as.data.frame(richness.arable_percentage_cook)
richness.arable_percentage_cook_df$ES_ID <- as.numeric(1:nrow(richness.arable_percentage_cook_df)) #add a new column with the effect size ID number

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1) and ???=0.5, 
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(richness.arable_percentage))))
plot(x=richness.arable_percentage_cook_df$ES_ID, y= richness.arable_percentage_cook_df$richness.arable_percentage_cook, 
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="I) Database: Species richness. Meta-regression model (Moderator= %Agricultural land)")
lines(x=richness.arable_percentage_cook_df$ES_ID,y= richness.arable_percentage_cook_df$richness.arable_percentage_cook, col = "black")
abline(h = qchisq(0.5, df=(length(coef(richness.arable_percentage)))), lty=2, lwd=2, col="red") #1.386294

# Identify possible effect size outliers and exclude outliers from the abundance database
richness_logRR_sensitivity.arable_percentage_cook<- left_join(richness.arable_percentage_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.arable_percentage_cook<0.024)

richness_logRR_sensitivity.arable_percentage_cook_outlier<- left_join(richness.arable_percentage_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.arable_percentage_cook>0.024)

# plot outliers points in different color
abline(h = 0.024, lty=2, lwd=2, col="red") #0.3;
points(x=richness_logRR_sensitivity.arable_percentage_cook_outlier$ES_ID, y=richness_logRR_sensitivity.arable_percentage_cook_outlier$richness.arable_percentage_cook,
       col="red",pch=19)
text(richness_logRR_sensitivity.arable_percentage_cook_outlier$ES_ID, richness_logRR_sensitivity.arable_percentage_cook_outlier$richness.arable_percentage_cook, 
     labels = richness_logRR_sensitivity.arable_percentage_cook_outlier$ES_ID, cex= 1, pos = c(2,2), col="red")

## Meta-analysis model without effect size outliers
richness.arable_percentage.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~ arable_percentage_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                                      tdist=TRUE, data=richness_logRR_sensitivity.arable_percentage_cook, method="REML")
summary(richness.arable_percentage.sensitivity_cook, digits=3)

# Results meta-regression (after sensitivity analysis)
coef(summary(richness.arable_percentage.sensitivity_cook))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3),
         ES_percent =round(ES_percent, digits =1),ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

#------- MODERATOR = DISTANCE TO NATURAL HABITATS
## Meta-regression model: using min_distance_mean (meters) as moderator:
richness.distance <- rma.mv(y=LRR, V=LRR_var, mods = ~min_distance_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.distance, digits=4)

# Calulate the residuals from the metaregression model using distance in meters
rs.richness.distance<-rstandard.rma.mv(richness.distance, type="rstandard")

# Calculate predictor values for the model
preds.richness.distance <-predict(richness.distance,newmods=c(0:1001), addx=TRUE)
preds.richness.distance<-as.data.frame(preds.richness.distance)

## Meta-regression model: using min_log_distance_mean (ln(distance + 1)) as moderator
richness.log_distance <- rma.mv(y=LRR, V=LRR_var, mods = ~min_log_distance_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.log_distance, digits=4)

# Calulate the residuals from the metaregression model using distance in log10
rs.richness.log_distance<-rstandard.rma.mv(richness.log_distance, type="rstandard")

# Calculate predictor values for the model
preds.richness.log_distance <-predict(richness.log_distance,newmods=c(0:7), addx=TRUE)
preds.richness.log_distance<-as.data.frame(preds.richness.log_distance)

par(mfrow=c(1,2),oma = c(0, 0, 2, 0)) #two plots in the same graph
# Plot residuals vs. distance in meters
plot(richness_logRR$min_distance_mean, rs.richness.distance$resid,
     ylab="Standardized residuals", xlab= "Distance (m)", main="A")
lines(x=preds.richness.distance$X.min_distance_mean ,y=preds.richness.distance$pred, col="#A63117", lwd=2)
# Plot residuals vs. distance in (ln(Distance + 1)
plot(richness_logRR$min_log_distance_mean, rs.richness.log_distance$resid,
     ylab="Standardized residuals", xlab= "ln(Distance + 1)", main="B")
lines(x=preds.richness.log_distance$X.min_log_distance_mean ,y=preds.richness.log_distance$pred, col="#A63117", lwd=2)
mtext("Database: Species richness", outer = TRUE, cex = 1.5)

## SENSITIVITY ANALYSIS
# COOK DISTANCE
system.time(richness.log_distance_cook <- cooks.distance(richness.log_distance, reestimate=FALSE))
richness.log_distance_cook_df<- as.data.frame(richness.log_distance_cook)
richness.log_distance_cook_df$ES_ID <- as.numeric(1:nrow(richness.log_distance_cook_df)) #add a new column with the effect size ID number

# Recomended cutoff for Cooks distance chi-squared distribution with df=(k+1) 
# or comparatively larger when compared to the other effect sizes (Viechtbauera & Cheung, 2010)
qchisq(0.5, df=(length(coef(richness.log_distance))))
plot(x=richness.log_distance_cook_df$ES_ID, y= richness.log_distance_cook_df$richness.log_distance_cook, 
     ylab= "Cook´s distance", xlab= "Effect size ID", cex.main=0.9, main="J) Database: Species richness. Meta-regression model (Moderator= Log Distance to natural habitats)")
lines(x=richness.log_distance_cook_df$ES_ID,y= richness.log_distance_cook_df$richness.log_distance_cook, col = "black")
abline(h = qchisq(0.5, df=(length(coef(richness.log_distance)))), lty=2, lwd=2, col="red") #1.386294

# Identify possible effect size outliers and exclude outliers from the abundance database
richness_logRR_sensitivity.log_distance_cook<- left_join(richness.log_distance_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.log_distance_cook<0.06)

richness_logRR_sensitivity.log_distance_cook_outlier<- left_join(richness.log_distance_cook_df,richness_logRR, by="ES_ID")%>%
  filter(richness.log_distance_cook>0.06)

# plot outliers points in different color
abline(h = 0.06, lty=2, lwd=2, col="red") #0.3;
points(x=richness_logRR_sensitivity.log_distance_cook_outlier$ES_ID, y=richness_logRR_sensitivity.log_distance_cook_outlier$richness.log_distance_cook,
       col="red",pch=19)
text(richness_logRR_sensitivity.log_distance_cook_outlier$ES_ID, richness_logRR_sensitivity.log_distance_cook_outlier$richness.log_distance_cook, 
     labels = richness_logRR_sensitivity.log_distance_cook_outlier$ES_ID, cex= 1, pos = c(2,4), col="red")

## Meta-analysis model without effect size outliers
richness.log_distance.sensitivity_cook <- rma.mv(y=LRR, V=LRR_var, mods = ~ min_log_distance_mean, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                                 tdist=TRUE, data=richness_logRR_sensitivity.log_distance_cook, method="REML")
summary(richness.log_distance.sensitivity_cook, digits=3)

# Results meta-regression (after sensitivity analysis)
coef(summary(richness.log_distance.sensitivity_cook))%>%
  mutate(estimate = round(estimate, digits = 3),
         se = round(se, digits =3),tval= round(tval, digits =3),pval= round(pval, digits =3),
         ci.lb = round(ci.lb, digits =3),ci.ub = round(ci.ub, digits =3))

##------ Plot RESULTS META-REGRESSION MODEL (MODERATORS = LANDSCAPE METRICS)
#A) %Natural habitats
preds.rich.natural.percentage_plot <-predict(richness.natural_percentage,newmods=c(0:100), addx=TRUE)
preds.rich.natural.percentage_plot<-as.data.frame(preds.rich.natural.percentage_plot)%>%
  mutate(pred_percentage = (100*(exp(pred)-1)))
class(preds.rich.natural.percentage_plot)
preds.rich.natural.percentage_plot
View(preds.rich.natural.percentage_plot)

summary(richness.natural_percentage)
rich.natural.percentage.summary<- coef(summary(richness.natural_percentage))%>%
  mutate(estimate= round(estimate, digits = 4),
         tval= round(tval, digits =2),
         pval= round(pval, digits =5),
         ci.lb = round(ci.lb, digits =4),
         ci.ub = round(ci.ub, digits =4),
         name = c("intercept", "slope"),
         x = c(1,4),
         y= c(2.5,3.7),
         label = paste("slope (%) = ", estimate, " (95% CI: ", ci.lb, ", ", ci.ub, ")",sep=""))%>%
  filter(name == "slope")
rich.natural.percentage.summary

figure_5_A<- ggplot(preds.rich.natural.percentage_plot,aes(x=X.natural_percentage_mean,y=pred))+
  geom_point(data=richness_logRR,aes(x=natural_percentage_mean,y=LRR),colour="grey70",fill="white",
             position="jitter",alpha=0.4,shape= 20,size=2)+
  #scale_y_continuous(limits = c(-100, 500))+
  geom_hline(yintercept = 0, colour = "black",)+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1.5, colour="#A63117")+
  geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), alpha=0.3)+
  labs(x=expression(bold(paste("Percentage of natural/semi-natural habitats"))),
       y= expression(bold(paste("Effect size (LRR)"))))+
  ggtitle("A)")+
  theme(
    axis.text.x = element_text(color="#22211d",size=8,  family = "sans"),
    axis.text.y = element_text(color="#22211d",size=8, family = "sans"),
    text = element_text(color = "#22211d", size =8, face = "bold", family = "sans"),
    plot.background = element_rect(fill = "White", color = "White"), 
    panel.background = element_rect(fill = "White", color = "White"), 
    axis.line = element_line(colour = "black"),
    plot.title = element_text(color="#22211d",size=9, family = "sans", vjust = 0.5, hjust=-0.1))
figure_5_A

#B) %Agricultural land
preds.rich.arable.percentage_plot <-predict(richness.arable_percentage,newmods=c(0:100), addx=TRUE)
preds.rich.arable.percentage_plot<-as.data.frame(preds.rich.arable.percentage_plot)
class(preds.rich.arable.percentage_plot)
preds.rich.arable.percentage_plot

summary(richness.arable_percentage)
rich.arable.percentage.summary<- coef(summary(richness.arable_percentage))%>%
  mutate(estimate= round(estimate, digits = 4),
         tval= round(tval, digits =2),
         pval= round(pval, digits =5),
         ci.lb = round(ci.lb, digits =4),
         ci.ub = round(ci.ub, digits =4),
         name = c("intercept", "slope"),
         x = c(1,4),
         y= c(2.5,3.7),
         label = paste("slope (%) = ", estimate, " (95% CI: ", ci.lb, ", ", ci.ub, ")",sep=""))%>%
  filter(name == "slope")
rich.arable.percentage.summary

figure_5_B<- ggplot2::ggplot(preds.rich.arable.percentage_plot,aes(x=X.arable_percentage_mean,y=pred))+
  geom_point(data=richness_logRR,aes(x=arable_percentage_mean,y=LRR),colour="grey70",fill="white",
             position="jitter",alpha=0.4,shape= 20,size=2)+
  geom_hline(yintercept = 0, colour = "black",)+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1.5, colour="#A63117")+
  geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), alpha=0.3)+
  labs(x=expression(bold(paste("Percentage of agricultural land"))),
       y= expression(bold(paste("Effect size (LRR)"))))+
  ggtitle("B)")+
  theme(
    axis.text.x = element_text(color="#22211d",size=8,  family = "sans"),
    axis.text.y = element_text(color="#22211d",size=8, family = "sans"),
    text = element_text(color = "#22211d", size =8, face = "bold", family = "sans"),
    plot.background = element_rect(fill = "White", color = "White"), 
    panel.background = element_rect(fill = "White", color = "White"), 
    axis.line = element_line(colour = "black"),
    plot.title = element_text(color="#22211d",size=9, family = "sans", vjust = 0.5, hjust=-0.1))
figure_5_B

#C) Distance to natural habitats
preds.rich.log.distance_plot <-predict(richness.log_distance,newmods=c(0,1,1.15, 1.2,1.3,1.4,1.5,2,8), addx=TRUE)
preds.rich.log.distance_plot<-as.data.frame(preds.rich.log.distance_plot)
class(preds.rich.log.distance_plot)
View(preds.rich.log.distance_plot)

summary(richness.log_distance)
rich.log.distance.summary<- coef(summary(richness.log_distance))%>%
  mutate(estimate= round(estimate, digits = 2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =3),
         ci.ub = round(ci.ub, digits =3),
         name = c("intercept", "slope"),
         x = c(1,4),
         y= c(2.5,3.7),
         label = paste("slope (%) = ", estimate, " (95% CI: ", ci.lb, ", ", ci.ub, ")",sep=""))%>%
  filter(name == "slope")
rich.log.distance.summary

figure_5_C<- ggplot(preds.rich.log.distance_plot,aes(x=X.min_log_distance_mean,y=pred))+
  geom_point(data=richness_logRR,aes(x=min_log_distance_mean,y=LRR),colour="grey70",fill="white",
             position="jitter",alpha=0.4,shape= 20,size=2)+
  geom_hline(yintercept = 0, colour = "black",)+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1.5, colour="#A63117")+
  geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), alpha=0.3)+
  scale_x_discrete(limits=c(0, 2,4,6,8),labels=c("0", "2", "4", "6", "8", ""))+
  labs(x=expression(bold(paste("Distance to natural or semi-natural habitats [ln(distance + 1)]"))),
       y= expression(bold(paste("Effect size (LRR)"))))+
  ggtitle("C)")+
  theme(
    axis.text.x = element_text(color="#22211d",size=8,  family = "sans"),
    axis.text.y = element_text(color="#22211d",size=8, family = "sans"),
    text = element_text(color = "#22211d", size =8, face = "bold", family = "sans"),
    plot.background = element_rect(fill = "White", color = "White"), 
    panel.background = element_rect(fill = "White", color = "White"), 
    axis.line = element_line(colour = "black"),
    plot.title = element_text(color="#22211d",size=9, family = "sans", vjust = 0.5, hjust=-0.1))
figure_5_C
summary(richness.log_distance)

ggpubr::ggarrange(figure_5_A,figure_5_B, figure_5_C, nrow = 3)

##############################################################################################################################333
################################  PUBLICATION BIAS ANALYSIS ###################################################################3----------------------------------------------------#
#####################################################################################################################################33

############################ A B U N D A N C E ######################################
##Egger regression test
##Replace the moderator by SE (see Habeck & Schultz 2015)
abun.bias.egger <-rma.mv(y=LRR, V=LRR_var, mods = ~ LRR_se, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                         tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.bias.egger, digits=3)
coef(summary(abun.bias.egger, digits=3))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         estimate = round(estimate, digits = 2),
         se= round(se, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =0),
         ci.lb_percent = round(ci.lb_percent, digits = 0),
         ci.ub_percent = round(ci.ub_percent, digits = 0))

resid.1<-rstandard.rma.mv(abun.bias.egger, type="rstandard")
resid.1
### Meta-regression
## MODERATOR: YEAR OF PUBLICATION
abun.bias.year <- rma.mv(y=LRR, V=LRR_var, mods = ~ Year, method="REML",
                         random = list(~ 1 | ES_ID, ~1 | ID), tdist=TRUE, data=abundance_logRR)
summary(abun.bias.year, digits=3)
coef(summary(abun.bias.year, digits=3))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         estimate = round(estimate, digits = 2),
         se= round(se, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =1),
         ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

## MODERATOR: LATITUDE
names(abundance_logRR)
abun.latitude <- rma.mv(y=LRR, V=LRR_var, mods = ~ Lat_C_2, 
                        random = list(~ 1 | ES_ID, ~1 | ID), tdist=TRUE, data=abundance_logRR, method="REML")
summary(abun.latitude, digits=3)
coef(summary(abun.latitude, digits=3))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         estimate = round(estimate, digits = 2),
         se= round(se, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =1),
         ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

############################ S P E C I E S  R I C H N E S S ######################################
##Egger regression test
##Replace the moderator by SE (see Habeck & Schultz 2015)
richness.bias.egger <-rma.mv(y=LRR, V=LRR_var, mods = ~ LRR_se, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                             tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.bias.egger, digits=3)
coef(summary(richness.bias.egger, digits=3))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         estimate = round(estimate, digits = 2),
         se= round(se, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =0),
         ci.lb_percent = round(ci.lb_percent, digits = 0),
         ci.ub_percent = round(ci.ub_percent, digits = 0))

##-- Meta-regression
#### MODERATOR: YEAR OF PUBLICATION
richness.year <- rma.mv(y=LRR, V=LRR_var, mods = ~ Year, method="REML",
                        random = list(~ 1 | ES_ID, ~1 | ID), tdist=TRUE, data= richness_logRR)
summary(richness.year, digits=3)
coef(summary(richness.year, digits=3))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         estimate = round(estimate, digits = 2),
         se= round(se, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =1),
         ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))

#### MODERATOR: LATITUDE
richness.latitude <- rma.mv(y=LRR, V=LRR_var, mods = ~ Lat_C_2, 
                            random = list(~ 1 | ES_ID, ~1 | ID), tdist=TRUE, data=richness_logRR, method="REML")
summary(richness.latitude, digits=3)
coef(summary(richness.latitude, digits=3))%>%
  mutate(ES_percent = (100*(exp(estimate)-1)),
         ci.lb_percent = (100*(exp(ci.lb)-1)),
         ci.ub_percent = (100*(exp(ci.ub)-1)),
         estimate = round(estimate, digits = 2),
         se= round(se, digits =2),
         tval= round(tval, digits =2),
         pval= round(pval, digits =2),
         ci.lb = round(ci.lb, digits =2),
         ci.ub = round(ci.ub, digits =2),
         ES_percent =round(ES_percent, digits =1),
         ci.lb_percent = round(ci.lb_percent, digits = 1),
         ci.ub_percent = round(ci.ub_percent, digits = 1))