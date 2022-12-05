#library(dplyr)

# Load datasets
clinical <- read.csv("input-datasets/Clinical Data miRNA Cohort.csv")
mirna <- read.csv("input-datasets/miRNA expression raw data.csv")
mri <- read.csv("input-datasets/MRI Imaging  results_2022-04-20.csv")

# RNA DATA PREPROCESSING
# Remove control samples
mirna = subset(mirna, select = -c(fb1, fb.2, FB.RNA3, FB.RNA4))
# Make miRNA IDs as index
rownames(mirna) = mirna$Sample.Name
# Remove miRNA IDs column
mirna = mirna[,-1]

# Separate plasma and CSF miRNA samples
mirna.p = mirna[ , grepl( "p" , names(mirna))]
mirna.c = mirna[ , grepl( "c" , names(mirna))]

# Pre-process plasma samples
mirna.p = as.matrix(mirna.p)
mirna.p = t(mirna.p)
mirna.p = as.data.frame(mirna.p)
colnames(mirna.p) <- paste(colnames(mirna.p), "p", sep = "_")
SDG_ID = rownames(mirna.p); mirna.p = cbind(SDG_ID, mirna.p)
library(stringr)
mirna.p$SDG_ID = str_sub(mirna.p$SDG_ID, 2, -2)
mirna.p$SDG_ID = gsub("\\_", ".", mirna.p$SDG_ID)
mirna.p$SDG_ID = gsub("\\.", "-", mirna.p$SDG_ID)

# Pre-process CSF samples
mirna.c = as.matrix(mirna.c)
mirna.c = t(mirna.c)
mirna.c = as.data.frame(mirna.c)
colnames(mirna.c) <- paste(colnames(mirna.c), "c", sep = "_")
SDG_ID = rownames(mirna.c); mirna.c = cbind(SDG_ID, mirna.c)
library(stringr)
mirna.c$SDG_ID = str_sub(mirna.c$SDG_ID, 2, -2)
mirna.c$SDG_ID = gsub("\\.", "-", mirna.c$SDG_ID)

#exclude sample 15635-21p, 15635-221p, 15635-246p, 15635-132c, 15635-225c but keep 15635-132c2
mirna.p = mirna.p[-which(mirna.p$SDG_ID == "15635-21"),]
mirna.p = mirna.p[-which(mirna.p$SDG_ID == "15635-221"),]
mirna.p = mirna.p[-which(mirna.p$SDG_ID == "15635-246"),]
mirna.c = mirna.c[-which(mirna.c$SDG_ID == "15635-132"),]
mirna.c = mirna.c[-which(mirna.c$SDG_ID == "15635-225"),]
#exclude two plasma samples that have “no evidence of disease” at time of collection based on MRI findings: 15635-254p, 15635-251p
mirna.p = mirna.p[-which(mirna.p$SDG_ID == "15635-254"),]
mirna.p = mirna.p[-which(mirna.p$SDG_ID == "15635-251"),]
# Rename Sample ID "15635-132c" to "15635-132"
mirna.c[which(mirna.c$SDG_ID == "15635-132c"),"SDG_ID"] = "15635-132"

# Remove control miRNAs and housekeeping miRNAs
mirna.p <- subset(mirna.p, select = -c(CTRL_ANT1_p, CTRL_ANT2_p, CTRL_ANT3_p, CTRL_ANT4_p, CTRL_ANT5_p, 
                            CTRL_miR_POS_p, HK_ACTB_p, HK_B2M_p, HK_GAPDH_p, HK_PPIA_p,
                            HK_RNU47_p, HK_RNU75_p, HK_RNY3_p, HK_RPL19_p, HK_RPL27_p,
                            HK_RPS12_p, HK_RPS20_p, HK_SNORA66_p, HK_YWHAZ_p))
mirna.c <- subset(mirna.c, select = -c(CTRL_ANT1_c, CTRL_ANT2_c, CTRL_ANT3_c, CTRL_ANT4_c, CTRL_ANT5_c, 
                                       CTRL_miR_POS_c, HK_ACTB_c, HK_B2M_c, HK_GAPDH_c, HK_PPIA_c,
                                       HK_RNU47_c, HK_RNU75_c, HK_RNY3_c, HK_RPL19_c, HK_RPL27_c,
                                       HK_RPS12_c, HK_RPS20_c, HK_SNORA66_c, HK_YWHAZ_c))

# CLINICAL DATA PREPROCESSING
clinical = clinical[,c("StudySubject_ID", "SDG_ID", "Specimen_Type", "Short_histology", "Age_at_LBcollection", 
                   "Age_at_.Initial_Diagnosis", "Age_LastKnownClinicalStatus", "OverallSurvival_Time_From_Initial_Tumor.Diagnosis", 
                   "PFS_from_InitialDiagnosis", "OverallSurvival_LBcollection")]

# MRI DATA PREPROCESSING
# Remove Patient LB00070 rows
mri = mri[-which(mri$StudySubject_ID == "LB00070"),]

mri = mri[,c("StudySubject_ID", "SDG_ID", "Short_histology", "Imaging.timepoint..age.in.days.",
             "Total.Tumor.Volume..mm.3.", "Enhancing..mm.3.", "Non.enhancing..mm.3.", "Cystic..core...mm.3.",
             "Cystic..reactive...mm.3.", "AmountTumor_Enhancing", "AmountTumor_Nonenhancing",
             "AmountTumor_CysticCore", "Adjacent.ventricular.system..y.n.", "Adjacent.surface..y.n.",
             "Leptomeningeal.disease..y.n.")]

################################################################################
# MRI MSGL
library(msgl)
mri.fit.cv <- msgl::cv(as.matrix(mri[,4:15]), mri$Short_histology, fold = 4, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
mri.fit.cv

mri.fit <- msgl::fit(as.matrix(mri[,4:15]), mri$Short_histology, alpha = 0.5, lambda = 0.1, standardize = FALSE)
mri.fit

#find the best model index from cv
features(mri.fit)[[best_model(mri.fit.cv)]]
parameters(mri.fit)[[best_model(mri.fit.cv)]]
coef(mri.fit, best_model(mri.fit.cv))
mri.coef = coef(mri.fit, best_model(mri.fit.cv))
mri.coef = as.matrix(mri.coef)
mri.coef[which(mri.coef == 0)] = NA
mri.coef = data.frame(mri.coef)
mri.coef = mri.coef[,-1]
mri.coef$tumor_type = rownames(mri.coef)

library(reshape2)
mri.coef = melt(mri.coef)
mri.coef$variable = as.character(mri.coef$variable)
mri.coef$variable[which(mri.coef$variable == "Total.Tumor.Volume..mm.3.")] = "Total Tumor Volume (mm^3)"
mri.coef$variable[which(mri.coef$variable == "Enhancing..mm.3.")] = "Enhancing (mm^3)"
mri.coef$variable[which(mri.coef$variable == "Non.enhancing..mm.3.")] = "Non-enhancing (mm^3)"
mri.coef$variable[which(mri.coef$variable == "Cystic..core...mm.3.")] = "Cystic (core) (mm^3)"


options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(mri.coef, aes(variable, tumor_type, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coefficient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(mri, aes(x = Short_histology, y = `Total.Tumor.Volume..mm.3.`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Total Tumor Volume (mm^3)", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p3 = ggplot(mri, aes(x = Short_histology, y = `Enhancing..mm.3.`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Enhancing (mm^3)", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p4 = ggplot(mri, aes(x = Short_histology, y = `Non.enhancing..mm.3.`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Non-enhancing (mm^3)", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p5 = ggplot(mri, aes(x = Short_histology, y = `Cystic..core...mm.3.`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Cystic (core) (mm^3)", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

library(ggpubr)
library(grid)
ggarrange(p1,
          labels = c("A"))
ggarrange(p2, p3, p4, p5,
          labels = c("B","C","D","E"))

require(nnet)
mri.multinom.fit <- multinom(Short_histology ~ `Total.Tumor.Volume..mm.3.` + `Enhancing..mm.3.` + `Non.enhancing..mm.3.` + `Cystic..core...mm.3.`, data = mri)
summary(mri.multinom.fit)
round((1 - pnorm(abs(summary(mri.multinom.fit)$coefficients/summary(mri.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(mri.multinom.fit))

head(mri.pp <- fitted(mri.multinom.fit))
library(sjPlot)
tab_model(mri.multinom.fit, digits = 6)

mri.prediction = predict(mri.multinom.fit, newdata = mri, "class")
sum(mri.prediction == mri$Short_histology)/nrow(mri)


mri.result1 = predict(mri.multinom.fit, newdata = mri, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(mri$Short_histology),
  predicted = as.matrix(mri.result1)
))
################################################################################
# Clinical MSGL
library(msgl)
clinical.fit.cv <- msgl::cv(as.matrix(clinical[,5:10]), clinical$Short_histology, fold = 4, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
clinical.fit.cv

clinical.fit <- msgl::fit(as.matrix(clinical[,5:10]), clinical$Short_histology, alpha = 0.5, lambda = 0.1, standardize = FALSE)
clinical.fit

#find the best model index from cv
features(clinical.fit)[[best_model(clinical.fit.cv)]]
parameters(clinical.fit)[[best_model(clinical.fit.cv)]]
coef(clinical.fit, best_model(clinical.fit.cv))
clinical.coef = coef(clinical.fit, best_model(clinical.fit.cv))
clinical.coef = as.matrix(clinical.coef)
clinical.coef[which(clinical.coef == 0)] = NA
clinical.coef = data.frame(clinical.coef)
clinical.coef = clinical.coef[,-1]
clinical.coef$tumor_type = rownames(clinical.coef)

library(reshape2)
clinical.coef = melt(clinical.coef)
clinical.coef$variable = as.character(clinical.coef$variable)
clinical.coef$variable[which(clinical.coef$variable == "Age_at_.Initial_Diagnosis")] = "Age at Initial Diagnosis"
clinical.coef$variable[which(clinical.coef$variable == "Age_LastKnownClinicalStatus")] = "Age at Last Known Clinical Status"
clinical.coef$variable[which(clinical.coef$variable == "OverallSurvival_Time_From_Initial_Tumor.Diagnosis")] = "Overall Survival Time From Initial Tumor Diagnosis"
clinical.coef$variable[which(clinical.coef$variable == "PFS_from_InitialDiagnosis")] = "PFS from Initial Diagnosis"


options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(clinical.coef, aes(variable, tumor_type, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coefficient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(clinical, aes(x = Short_histology, y = `Age_at_.Initial_Diagnosis`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Age at Initial Diagnosis", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p3 = ggplot(clinical, aes(x = Short_histology, y = `Age_LastKnownClinicalStatus`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Age at Last Known Clinical Status", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p4 = ggplot(clinical, aes(x = Short_histology, y = `OverallSurvival_Time_From_Initial_Tumor.Diagnosis`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Overall Survival Time From Initial Tumor Diagnosis", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p5 = ggplot(clinical, aes(x = Short_histology, y = `PFS_from_InitialDiagnosis`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "PFS from Initial Diagnosis", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

library(ggpubr)
library(grid)
ggarrange(p1,
          labels = c("A"))
ggarrange(p2, p3, p4, p5,
          labels = c("B","C","D","E"))

require(nnet)
clinical.multinom.fit <- multinom(Short_histology ~ `Age_at_.Initial_Diagnosis` + `Age_LastKnownClinicalStatus` + `OverallSurvival_Time_From_Initial_Tumor.Diagnosis` + `PFS_from_InitialDiagnosis`, data = clinical)
summary(clinical.multinom.fit)
round((1 - pnorm(abs(summary(clinical.multinom.fit)$coefficients/summary(clinical.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(clinical.multinom.fit))

head(clinical.pp <- fitted(clinical.multinom.fit))
library(sjPlot)
tab_model(clinical.multinom.fit, digits = 6)

clinical.prediction = predict(clinical.multinom.fit, newdata = clinical, "class")
sum(clinical.prediction == clinical$Short_histology)/nrow(clinical)


clinical.result1 = predict(clinical.multinom.fit, newdata = clinical, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(clinical$Short_histology),
  predicted = as.matrix(clinical.result1)
))

################################################################################
# Plasma miRNA MSGL

# Add Short Histology to plasma miRNA dataset
df1 = clinical[,c("StudySubject_ID", "SDG_ID", "Short_histology")]
pmirna = merge(df1, mirna.p, by = "SDG_ID", all = T)

#exclude two plasma samples that have “no evidence of disease” at time of collection based on MRI findings: 15635-254p, 15635-251p
pmirna = pmirna[-which(pmirna$SDG_ID == "15635-254"),]
pmirna = pmirna[-which(pmirna$SDG_ID == "15635-251"),]



library(msgl)
pmirna.fit.cv <- msgl::cv(as.matrix(pmirna[,4:2086]), pmirna$Short_histology, fold = 4, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
pmirna.fit.cv

pmirna.fit <- msgl::fit(as.matrix(pmirna[,4:2086]), pmirna$Short_histology, alpha = 0.5, lambda = 0.1, standardize = FALSE)
pmirna.fit

#find the best model index from cv
features(pmirna.fit)[[best_model(pmirna.fit.cv)]]
parameters(pmirna.fit)[[best_model(pmirna.fit.cv)]]
coef(pmirna.fit, best_model(pmirna.fit.cv))
pmirna.coef = coef(pmirna.fit, best_model(pmirna.fit.cv))
pmirna.coef = as.matrix(pmirna.coef)
#pmirna.coef[which(pmirna.coef == 0)] = NA
pmirna.coef = data.frame(pmirna.coef)
pmirna.coef$tumor_type = rownames(pmirna.coef)
# Drop first column containing intercepts...
pmirna.coef = pmirna.coef[,-1]

library(reshape2)
pmirna.coef = melt(pmirna.coef)
pmirna.coef$variable = as.character(pmirna.coef$variable)
pmirna.coef$variable[which(pmirna.coef$variable == "miR.451a_p")] = "miR-451a"


options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(pmirna.coef, aes(variable, tumor_type, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coefficient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(pmirna, aes(x = Short_histology, y = `miR-451a_p`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-451a", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

library(ggpubr)
library(grid)
ggarrange(p1,
          labels = c("A"))
ggarrange(p2,
          labels = c("B"))

require(nnet)
pmirna.multinom.fit <- multinom(Short_histology ~ `miR-451a_p`, data = pmirna)
summary(pmirna.multinom.fit)
round((1 - pnorm(abs(summary(pmirna.multinom.fit)$coefficients/summary(pmirna.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(pmirna.multinom.fit))

head(pmirna.pp <- fitted(pmirna.multinom.fit))
library(sjPlot)
tab_model(pmirna.multinom.fit, digits = 6)

pmirna.prediction = predict(pmirna.multinom.fit, newdata = pmirna, "class")
sum(pmirna.prediction == pmirna$Short_histology)/nrow(pmirna)


pmirna.result1 = predict(pmirna.multinom.fit, newdata = pmirna, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(pmirna$Short_histology),
  predicted = as.matrix(pmirna.result1)
))

################################################################################
# CSF miRNA MSGL

# Add Short Histology to CSF miRNA dataset
df1 = clinical[,c("StudySubject_ID", "SDG_ID", "Short_histology")]
cmirna = merge(df1, mirna.c, by = "SDG_ID", all = T)

cmirna = cmirna[-which(cmirna$SDG_ID == "15635-108"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-127"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-129"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-134"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-148"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-154"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-158"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-161"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-182"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-187"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-215"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-225"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-228"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-234"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-251"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-254"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-31"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-38"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-42"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-53"),]
cmirna = cmirna[-which(cmirna$SDG_ID == "15635-90"),]

# Only one ATRT sample, so it will be dropped
#table(cmirna$Short_histology)
cmirna = cmirna[-which(cmirna$Short_histology == "ATRT"),]

library(msgl)
cmirna.fit.cv <- msgl::cv(as.matrix(cmirna[,4:2086]), cmirna$Short_histology, fold = 3, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
cmirna.fit.cv


cmirna.fit <- msgl::fit(as.matrix(cmirna[,4:2086]), cmirna$Short_histology, alpha = 0.5, lambda = 0.1, standardize = FALSE)
cmirna.fit

#find the best model index from cv
features(cmirna.fit)[[best_model(cmirna.fit.cv)]]
parameters(cmirna.fit)[[best_model(cmirna.fit.cv)]]
coef(cmirna.fit, best_model(cmirna.fit.cv))
cmirna.coef = coef(cmirna.fit, best_model(cmirna.fit.cv))
cmirna.coef = as.matrix(cmirna.coef)
cmirna.coef[which(cmirna.coef == 0)] = NA
cmirna.coef = data.frame(cmirna.coef)
cmirna.coef$tumor_type = rownames(cmirna.coef)
# Drop first column containing intercepts
cmirna.coef = cmirna.coef[,-1]

library(reshape2)
cmirna.coef = melt(cmirna.coef)
cmirna.coef$variable = as.character(cmirna.coef$variable)
cmirna.coef$variable[which(cmirna.coef$variable == "miR.124.3p_c")] = "miR-124-3p"
cmirna.coef$variable[which(cmirna.coef$variable == "miR.3197_c")] = "miR-3197"
cmirna.coef$variable[which(cmirna.coef$variable == "miR.451a_c")] = "miR-451a"
cmirna.coef$variable[which(cmirna.coef$variable == "miR.6126_c")] = "miR-6126"
cmirna.coef$variable[which(cmirna.coef$variable == "miR.6131_c")] = "miR-6131"
cmirna.coef$variable[which(cmirna.coef$variable == "miR.6870.3p_c")] = "miR-6870-3p"

options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(cmirna.coef, aes(variable, tumor_type, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coefficient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(cmirna, aes(x = Short_histology, y = `miR-124-3p_c`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-124-3p", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p3 = ggplot(cmirna, aes(x = Short_histology, y = `miR-3197_c`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-3197", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p4 = ggplot(cmirna, aes(x = Short_histology, y = `miR-451a_c`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-451a", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p5 = ggplot(cmirna, aes(x = Short_histology, y = `miR-6126_c`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-6126", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p6 = ggplot(cmirna, aes(x = Short_histology, y = `miR-6131_c`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-6131", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p7 = ggplot(cmirna, aes(x = Short_histology, y = `miR-6870-3p_c`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-6870-3p", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

library(ggpubr)
library(grid)
ggarrange(p1,
          labels = c("A"))
# ggarrange(p2,p3,
#           labels = c("B","C"))
# ggarrange(p4,p5,
#           labels = c("D","E"))
# ggarrange(p6,p7,
#           labels = c("F","G"))
ggarrange(p2,p3,p4,p5,p6,p7,
          labels = c("B","C","D","E","F","G"))

require(nnet)
cmirna.multinom.fit <- multinom(Short_histology ~ `miR-124-3p_c`+`miR-3197_c`+`miR-451a_c`+`miR-6126_c`+`miR-6131_c`+`miR-6870-3p_c`, data = cmirna)
summary(cmirna.multinom.fit)
round((1 - pnorm(abs(summary(cmirna.multinom.fit)$coefficients/summary(cmirna.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(cmirna.multinom.fit))

head(cmirna.pp <- fitted(cmirna.multinom.fit))
library(sjPlot)
tab_model(cmirna.multinom.fit, digits = 6)

cmirna.prediction = predict(cmirna.multinom.fit, newdata = cmirna, "class")
sum(cmirna.prediction == cmirna$Short_histology)/nrow(cmirna)


cmirna.result1 = predict(cmirna.multinom.fit, newdata = cmirna, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(cmirna$Short_histology),
  predicted = as.matrix(cmirna.result1)
))

################################################################################
# Plasma miRNA + MRI
mri.pmirna = merge(mri, mirna.p, by = "SDG_ID", all = T)
mri.pmirna = mri.pmirna[-which(mri.pmirna$SDG_ID == "15635-127"),]



library(msgl)
mri.pmirna.fit.cv <- msgl::cv(as.matrix(mri.pmirna[,4:2098]), mri.pmirna$Short_histology, fold = 4, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
mri.pmirna.fit.cv

mri.pmirna.fit <- msgl::fit(as.matrix(mri.pmirna[,4:2098]), mri.pmirna$Short_histology, alpha = 0.5, lambda = 0.1, standardize = FALSE)
mri.pmirna.fit

#find the best model index from cv
features(mri.pmirna.fit)[[best_model(mri.pmirna.fit.cv)]]
parameters(mri.pmirna.fit)[[best_model(mri.pmirna.fit.cv)]]
coef(mri.pmirna.fit, best_model(mri.pmirna.fit.cv))
mri.pmirna.coef = coef(mri.pmirna.fit, best_model(mri.pmirna.fit.cv))
mri.pmirna.coef = as.matrix(mri.pmirna.coef)
mri.pmirna.coef[which(mri.pmirna.coef == 0)] = NA
mri.pmirna.coef = data.frame(mri.pmirna.coef)
mri.pmirna.coef$tumor_type = rownames(mri.pmirna.coef)
# Drop first column containing intercepts
mri.pmirna.coef = mri.pmirna.coef[,-1]

library(reshape2)
mri.pmirna.coef = melt(mri.pmirna.coef)
mri.pmirna.coef$variable = as.character(mri.pmirna.coef$variable)
mri.pmirna.coef$variable[which(mri.pmirna.coef$variable == "miR.451a_p")] = "miR-451a"


options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(mri.pmirna.coef, aes(variable, tumor_type, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coefficient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(mri.pmirna, aes(x = Short_histology, y = `miR-451a_p`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-451a", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

library(ggpubr)
library(grid)
ggarrange(p1,
          labels = c("A"))
ggarrange(p2,
          labels = c("B"))

require(nnet)
mri.pmirna.multinom.fit <- multinom(Short_histology ~ `Cystic..core...mm.3.` + `Enhancing..mm.3.` + `Non.enhancing..mm.3.` + `Total.Tumor.Volume..mm.3.` + `miR-451a_p`, data = mri.pmirna)
summary(mri.pmirna.multinom.fit)
round((1 - pnorm(abs(summary(mri.pmirna.multinom.fit)$coefficients/summary(mri.pmirna.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(mri.pmirna.multinom.fit))

head(mri.pmirna.pp <- fitted(mri.pmirna.multinom.fit))
library(sjPlot)
tab_model(mri.pmirna.multinom.fit, digits = 6)

mri.pmirna.prediction = predict(mri.pmirna.multinom.fit, newdata = mri.pmirna, "class")
sum(mri.pmirna.prediction == mri.pmirna$Short_histology)/nrow(mri.pmirna)


mri.pmirna.result1 = predict(mri.pmirna.multinom.fit, newdata = mri.pmirna, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(mri.pmirna$Short_histology),
  predicted = as.matrix(mri.pmirna.result1)
))

################################################################################
# CSF miRNA + MRI
mri.cmirna = merge(mri, mirna.c, by = "SDG_ID", all = T)

mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-108"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-129"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-134"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-148"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-154"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-158"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-161"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-182"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-187"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-215"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-225"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-228"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-234"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-31"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-38"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-42"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-53"),]
mri.cmirna = mri.cmirna[-which(mri.cmirna$SDG_ID == "15635-90"),]


# Drop ATRT single sample
# table(mri.cmirna$Short_histology)
mri.cmirna = mri.cmirna[-which(mri.cmirna$Short_histology == "ATRT"),]

library(msgl)
mri.cmirna.fit.cv <- msgl::cv(as.matrix(mri.cmirna[,4:2098]), mri.cmirna$Short_histology, fold = 3, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
mri.cmirna.fit.cv

mri.cmirna.fit <- msgl::fit(as.matrix(mri.cmirna[,4:2098]), mri.cmirna$Short_histology, alpha = 0.5, lambda = 0.1, standardize = FALSE)
mri.cmirna.fit

#find the best model index from cv
features(mri.cmirna.fit)[[best_model(mri.cmirna.fit.cv)]]
parameters(mri.cmirna.fit)[[best_model(mri.cmirna.fit.cv)]]
coef(mri.cmirna.fit, best_model(mri.cmirna.fit.cv))
mri.cmirna.coef = coef(mri.cmirna.fit, best_model(mri.cmirna.fit.cv))
mri.cmirna.coef = as.matrix(mri.cmirna.coef)
mri.cmirna.coef[which(mri.cmirna.coef == 0)] = NA
mri.cmirna.coef = data.frame(mri.cmirna.coef)
mri.cmirna.coef$tumor_type = rownames(mri.cmirna.coef)
# Drop first column containing intercepts
mri.cmirna.coef = mri.cmirna.coef[,-1]

library(reshape2)
mri.cmirna.coef = melt(mri.cmirna.coef)
mri.cmirna.coef$variable = as.character(mri.cmirna.coef$variable)
mri.cmirna.coef$variable[which(mri.cmirna.coef$variable == "miR.3197_c")] = "miR-3197"


options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(mri.cmirna.coef, aes(variable, tumor_type, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coefficient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(mri.cmirna, aes(x = Short_histology, y = `miR-3197_c`, group = Short_histology, color = Short_histology)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-3197", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

library(ggpubr)
library(grid)
ggarrange(p1,
          labels = c("A"))
ggarrange(p2,
          labels = c("B"))

require(nnet)
mri.cmirna.multinom.fit <- multinom(Short_histology ~ `Cystic..core...mm.3.` + `Enhancing..mm.3.` + `Non.enhancing..mm.3.` + `Total.Tumor.Volume..mm.3.` + `miR-124-3p_c` + `miR-3197_c` + `miR-451a_c` + `miR-6126_c` + `miR-6131_c` + `miR-6870-3p_c`, data = mri.cmirna)
summary(mri.cmirna.multinom.fit)
round((1 - pnorm(abs(summary(mri.cmirna.multinom.fit)$coefficients/summary(mri.cmirna.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(mri.cmirna.multinom.fit))

head(mri.cmirna.pp <- fitted(mri.cmirna.multinom.fit))
library(sjPlot)
tab_model(mri.cmirna.multinom.fit, digits = 6)

mri.cmirna.prediction = predict(mri.cmirna.multinom.fit, newdata = mri.cmirna, "class")
sum(mri.cmirna.prediction == mri.cmirna$Short_histology)/nrow(mri.cmirna)


mri.cmirna.result1 = predict(mri.cmirna.multinom.fit, newdata = mri.cmirna, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(mri.cmirna$Short_histology),
  predicted = as.matrix(mri.cmirna.result1)
))


###################
# Export Plasma miRNA data
#write.csv(pmirna,"preprocessed-datasets/plasma_mirna.csv", row.names = FALSE)
# Export CSF miRNA data
#write.csv(cmirna,"preprocessed-datasets/csf_mirna.csv", row.names = FALSE)
# Export clinical data
#write.csv(clinical,"preprocessed-datasets/clinical_data.csv", row.names = FALSE)
# Export MRI data
#write.csv(mri,"preprocessed-datasets/mri_data.csv", row.names = FALSE)
# Export MRI + Plasma miRNA data
#write.csv(mri.pmirna,"preprocessed-datasets/mri_plasma_mirna.csv", row.names = FALSE)
# Export MRI + CSF miRNA data
#write.csv(mri.cmirna,"preprocessed-datasets/mri_csf_mirna.csv", row.names = FALSE)

###############################################################################################
###############################################################################################
# Add ATRT back into CSF Dataset
# CSF miRNA MSGL

# Add Short Histology to CSF miRNA dataset
df1 = clinical[,c("StudySubject_ID", "SDG_ID", "Short_histology")]
cmirna.v2 = merge(df1, mirna.c, by = "SDG_ID", all = T)

cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-108"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-127"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-129"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-134"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-148"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-154"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-158"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-161"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-182"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-187"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-215"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-225"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-228"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-234"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-251"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-254"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-31"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-38"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-42"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-53"),]
cmirna.v2 = cmirna.v2[-which(cmirna.v2$SDG_ID == "15635-90"),]

require(nnet)
cmirna.v2.multinom.fit <- multinom(Short_histology ~ `miR-124-3p_c`+`miR-3197_c`+`miR-451a_c`+`miR-6126_c`+`miR-6131_c`+`miR-6870-3p_c`, data = cmirna.v2)
summary(cmirna.v2.multinom.fit)
round((1 - pnorm(abs(summary(cmirna.v2.multinom.fit)$coefficients/summary(cmirna.v2.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(cmirna.v2.multinom.fit))

head(cmirna.v2.pp <- fitted(cmirna.v2.multinom.fit))
library(sjPlot)
tab_model(cmirna.v2.multinom.fit, digits = 6)

cmirna.v2.prediction = predict(cmirna.v2.multinom.fit, newdata = cmirna.v2, "class")
sum(cmirna.v2.prediction == cmirna.v2$Short_histology)/nrow(cmirna.v2)


cmirna.v2.result1 = predict(cmirna.v2.multinom.fit, newdata = cmirna.v2, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(cmirna.v2$Short_histology),
  predicted = as.matrix(cmirna.v2.result1)
))

###############################################################################################
###############################################################################################
# Manually selecting more miRNAs for Plasma miRNA Dataset
require(nnet)
pmirna.v2.multinom.fit <- multinom(Short_histology ~ `miR-16-5p_p` + `miR-3197_p` + `miR-451a_p` + `miR-4745-3p_p` + `miR-6126_p` + `miR-6870-3p_p`, data = pmirna)
summary(pmirna.v2.multinom.fit)
round((1 - pnorm(abs(summary(pmirna.v2.multinom.fit)$coefficients/summary(pmirna.v2.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(pmirna.v2.multinom.fit))

head(pmirna.v2.pp <- fitted(pmirna.v2.multinom.fit))
library(sjPlot)
tab_model(pmirna.v2.multinom.fit, digits = 6)

pmirna.v2.prediction = predict(pmirna.v2.multinom.fit, newdata = pmirna, "class")
sum(pmirna.v2.prediction == pmirna$Short_histology)/nrow(pmirna)


pmirna.v2.result1 = predict(pmirna.v2.multinom.fit, newdata = pmirna, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(pmirna$Short_histology),
  predicted = as.matrix(pmirna.v2.result1)
))

###############################################################################################
###############################################################################################
# Re-Doing MRI + Plasma miRNA
require(nnet)
mri.pmirna.v2.multinom.fit <- multinom(Short_histology ~ `Cystic..core...mm.3.` + `Enhancing..mm.3.` + `Non.enhancing..mm.3.` + `Total.Tumor.Volume..mm.3.` + `miR-16-5p_p` + `miR-3197_p` + `miR-451a_p` + `miR-4745-3p_p` + `miR-6126_p` + `miR-6870-3p_p`, data = mri.pmirna)
summary(mri.pmirna.v2.multinom.fit)
round((1 - pnorm(abs(summary(mri.pmirna.v2.multinom.fit)$coefficients/summary(mri.pmirna.v2.multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(mri.pmirna.v2.multinom.fit))

head(mri.pmirna.v2.pp <- fitted(mri.pmirna.v2.multinom.fit))
library(sjPlot)
tab_model(mri.pmirna.v2.multinom.fit, digits = 6)

mri.pmirna.v2.prediction = predict(mri.pmirna.v2.multinom.fit, newdata = mri.pmirna, "class")
sum(mri.pmirna.v2.prediction == mri.pmirna$Short_histology)/nrow(mri.pmirna)


mri.pmirna.v2.result1 = predict(mri.pmirna.v2.multinom.fit, newdata = mri.pmirna, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(mri.pmirna$Short_histology),
  predicted = as.matrix(mri.pmirna.v2.result1)
))


################################################################################
# Bar Plot of Results
#rm(list = ls())
df = read.csv("Accuracy.csv")
df$Features = factor(df$Features,
                     levels = c("MRI", "Plasma miRNA", "CSF miRNA", "MRI and Plasma miRNA", "MRI and CSF miRNA"))
df$Value = round(df$Value, digits = 2)

ggplot(df, aes(x = Features, y = Value, fill = Variable)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label = Value), vjust = -1.6, position = position_dodge(0.9)) + theme_bw() + 
  labs(x = "", y = "Prediction accuracy", fill = "") + ylim(0,1.05) + 
  # scale_x_discrete(labels=c("MRI" = "4 MRI vars", 
  #                           "Plasma miRNA" = "6 miRNA vars",
  #                           "CSF miRNA" = "6 miRNA vars",
  #                           "Plasma miRNA + MRI" = "4 MRI vars + 6 Plasma miRNA",
  #                           "CSF miRNA + MRI" = "4 MRI vars + 6 Plasma miRNA"
  #                           )) + 
  scale_x_discrete(labels=c("MRI" = "MRI (4)", 
                            "Plasma miRNA" = "Plasma miRNA (6)",
                            "CSF miRNA" = "CSF miRNA (6)",
                            "MRI and Plasma miRNA" = "MRI (4) + Plasma miRNA (6)",
                            "MRI and CSF miRNA" = "MRI (4) + CSF miRNA (6)")) + 
  theme(legend.position = "right",
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, face = "plain")) + 
  scale_fill_manual(values = c("black","grey"))
