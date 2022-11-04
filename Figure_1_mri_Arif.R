rm(list = ls())
setwd("/Users/shehbeelarif/Desktop/chao_repeat/Datasets")

##MRI###
mri = read.csv("mri_data.csv")

length(unique(mri$StudySubject_ID));unique(mri$StudySubject_ID)

## Use the variables recommended by Mateusz
mri = mri[,c("StudySubject_ID","SDG_ID","Short_histology","Total.Tumor.Volume..mm.3.","Enhancing..mm.3.",
             "Non.enhancing..mm.3.","Cystic..core...mm.3.","Cystic..reactive...mm.3.",
             "Edema..mm.3.","Adjacent.surface..y.n.","Leptomeningeal.disease..y.n.",
             "Tumor.location..general.")] 

library(msgl)
fit.cv <- msgl::cv(as.matrix(mri[,4:12]), mri$Short_histology, fold = 4, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
fit.cv

fit <- msgl::fit(as.matrix(mri[,4:12]), mri$Short_histology, alpha = 0.5, lambda = 0.1, standardize = FALSE)
fit
#find the best model index from cv
features(fit)[[best_model(fit.cv)]]
parameters(fit)[[best_model(fit.cv)]]
coef(fit, best_model(fit.cv))
coef = coef(fit, best_model(fit.cv))
coef = as.matrix(coef)
coef[which(coef == 0)] = NA
coef = data.frame(coef)
coef = coef[,-1]
coef$cancer_type = rownames(coef)

library(reshape2)
coef = melt(coef)
coef$variable = as.character(coef$variable)
coef$variable[which(coef$variable == "Total.Tumor.Volume..mm.3.")] = "Total Tumor Volume (mm^3)"
coef$variable[which(coef$variable == "Enhancing..mm.3.")] = "Enhancing (mm^3)"
coef$variable[which(coef$variable == "Non.enhancing..mm.3.")] = "Non-enhancing (mm^3)"
coef$variable[which(coef$variable == "Cystic..core...mm.3.")] = "Cystic (core) (mm^3)"
coef$variable[which(coef$variable == "Edema..mm.3.")] = "Edema (mm^3)"

options(scipen=10000)
library(corrplot)
library(viridis)
library(ggplot2)
library(RColorBrewer)
p1 = ggplot(coef, aes(variable, cancer_type, fill = value > 0)) + 
  geom_tile() + 
  geom_text(aes(label = formatC(value, format = "e", digits = 2))) + 
  scale_fill_discrete(labels = c("Coeffcient < 0", "Coefficient > 0", "Not selected")) + 
  labs(x = "", y = "", fill = "") + 
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5, color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain"),
        strip.text.x = element_text(size = 12))

p2 = ggplot(mri, aes(x = reorder(Short_histology, `Total.Tumor.Volume..mm.3.`, mean), y = `Total.Tumor.Volume..mm.3.`, group = Short_histology, color = reorder(Short_histology, -`Total.Tumor.Volume..mm.3.`, mean))) + 
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

p3 = ggplot(mri, aes(x = reorder(Short_histology, `Enhancing..mm.3.`, mean), y = `Enhancing..mm.3.`, group = Short_histology, color = reorder(Short_histology, -`Enhancing..mm.3.`, mean))) + 
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

p4 = ggplot(mri, aes(x = reorder(Short_histology, `Non.enhancing..mm.3.`, mean), y = `Non.enhancing..mm.3.`, group = Short_histology, color = reorder(Short_histology, -`Non.enhancing..mm.3.`, mean))) + 
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

p5 = ggplot(mri, aes(x = reorder(Short_histology, `Cystic..core...mm.3.`, mean), y = `Cystic..core...mm.3.`, group = Short_histology, color = reorder(Short_histology, -`Cystic..core...mm.3.`, mean))) + 
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

p5 = ggplot(mri, aes(x = reorder(Short_histology, `Edema..mm.3.`, mean), y = `Edema..mm.3.`, group = Short_histology, color = reorder(Short_histology, -`Edema..mm.3.`, mean))) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "Edema (mm^3)", color = "") + theme_bw() +
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
multinom.fit <- multinom(Short_histology ~ `Total.Tumor.Volume..mm.3.`+`Enhancing..mm.3.` + `Non.enhancing..mm.3.` + `Cystic..core...mm.3.` + `Edema..mm.3.`, data = mri)
summary(multinom.fit)
round((1 - pnorm(abs(summary(multinom.fit)$coefficients/summary(multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(multinom.fit))

head(pp <- fitted(multinom.fit))
library(sjPlot)
tab_model(multinom.fit, digits = 6)

prediction = predict(multinom.fit, newdata = mri, "class")
sum(prediction == mri$Short_histology)/nrow(mri)


result1 = predict(multinom.fit, newdata = mri, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(mri$Short_histology),
  predicted = as.matrix(result1)
))




