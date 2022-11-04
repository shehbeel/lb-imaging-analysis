rm(list = ls())
setwd("/Users/shehbeelarif/Desktop/chao_repeat/Datasets")

##Plasma miRNA###
pmirna = read.csv("pmirna_data.csv")

length(unique(pmirna$SDG_ID));unique(pmirna$SDG_ID)

library(msgl)
fit.cv <- msgl::cv(as.matrix(pmirna[,3:2084]), pmirna$Short_histology, fold = 4, alpha = 0.5, lambda = 0.1, use_parallel = TRUE, standardize = FALSE)
fit.cv

fit <- msgl::fit(as.matrix(pmirna[,3:2084]), pmirna$Short_histology, alpha = 0.1, lambda = 0.1, standardize = FALSE)
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
coef$variable[which(coef$variable == "miR.16.5p_p")] = "miR-16-5p"
coef$variable[which(coef$variable == "miR.3197_p")] = "miR-3197"
coef$variable[which(coef$variable == "miR.451a_p")] = "miR-451a"
coef$variable[which(coef$variable == "miR.4745.3p_p")] = "miR-4745-3p"
coef$variable[which(coef$variable == "miR.6126_p")] = "miR-6126"
coef$variable[which(coef$variable == "miR.6870.3p_p")] = "miR-6870-3p"

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

p2 = ggplot(pmirna, aes(x = reorder(Short_histology, `miR.16.5p_p`, mean), y = `miR.16.5p_p`, group = Short_histology, color = reorder(Short_histology, -`miR.16.5p_p`, mean))) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-16-5p", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p3 = ggplot(pmirna, aes(x = reorder(Short_histology, `miR.3197_p`, mean), y = `miR.3197_p`, group = Short_histology, color = reorder(Short_histology, -`miR.3197_p`, mean))) + 
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

p4 = ggplot(pmirna, aes(x = reorder(Short_histology, `miR.451a_p`, mean), y = `miR.451a_p`, group = Short_histology, color = reorder(Short_histology, -`miR.451a_p`, mean))) + 
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

p5 = ggplot(pmirna, aes(x = reorder(Short_histology, `miR.4745.3p_p`, mean), y = `miR.4745.3p_p`, group = Short_histology, color = reorder(Short_histology, -`miR.4745.3p_p`, mean))) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, size = 6,position=position_jitter(0.2)) + 
  labs(x = "", y = "miR-4745-3p", color = "") + theme_bw() +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "plain")) + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="red", fill="red") + 
  scale_color_manual(values = viridis(6))

p6 = ggplot(pmirna, aes(x = reorder(Short_histology, `miR.6126_p`, mean), y = `miR.6126_p`, group = Short_histology, color = reorder(Short_histology, -`miR.6126_p`, mean))) + 
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

p7 = ggplot(pmirna, aes(x = reorder(Short_histology, `miR.6870.3p_p`, mean), y = `miR.6870.3p_p`, group = Short_histology, color = reorder(Short_histology, -`miR.6870.3p_p`, mean))) + 
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
ggarrange(p2, p3,
          labels = c("B","C"))
ggarrange(p4, p5,
          labels = c("D","E"))
ggarrange(p6, p7,
          labels = c("F", "G"))


coef$variable[which(coef$variable == "miR.16.5p_p")] = "miR-16-5p"
coef$variable[which(coef$variable == "miR.3197_p")] = "miR-3197"
coef$variable[which(coef$variable == "miR.451a_p")] = "miR-451a"
coef$variable[which(coef$variable == "miR.4745.3p_p")] = "miR-4745-3p"
coef$variable[which(coef$variable == "miR.6126_p")] = "miR-6126"
coef$variable[which(coef$variable == "miR.6870.3p_p")] = "miR-6870-3p"





require(nnet)
multinom.fit <- multinom(Short_histology ~ `miR.16.5p_p`+`miR.3197_p` + `miR.451a_p` + `miR.4745.3p_p` + `miR.6126_p`, `miR.6870.3p_p`, data = pmirna)
summary(multinom.fit)
round((1 - pnorm(abs(summary(multinom.fit)$coefficients/summary(multinom.fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(multinom.fit))

head(pp <- fitted(multinom.fit))
library(sjPlot)
tab_model(multinom.fit, digits = 6)

prediction = predict(multinom.fit, newdata = pmirna, "class")
sum(prediction == pmirna$Short_histology)/nrow(pmirna)


result1 = predict(multinom.fit, newdata = pmirna, type='probs')
library(HandTill2001)
auc(multcap(
  response = factor(pmirna$Short_histology),
  predicted = as.matrix(result1)
))

