#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Corr_rois_cv_1.1.R                                                       #
# Author: Yang Hu                                                          #
# Established Date: May. 17, 2018                                          #
# Last Change: Feb. 5, 2020                                                #
# Project: Corruption (Corr)                                               #
# PI: Yang Hu                                                              #
# Collaborators: Chen Qu, Jean-Claude Dreher                               #
# This is a custom R script to draw plots of contrast values (cv) of rois  #
# from GLM analyses (fMRI study).                                          #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#**************************************************************************#
#                          Step 1: Preparation                             #
#**************************************************************************#
library(readxl)
library(ggplot2)
# library(ggpmisc)
require(GGally)
require(R.matlab)
require(scales)
# require(robust)
require(robustbase)
require(complmrob)
require(betas) # for standard correlation coef of lmrob
require(MASS) # for stat_method('rlm')
require(ggExtra)
library(lme4)
library(psych)
library(reshape2)
library(plyr)

sbjs <- c('sbj01','sbj02','sbj03','sbj05','sbj06','sbj07','sbj08','sbj09','sbj10',
          'sbj11','sbj12','sbj13','sbj14','sbj15','sbj16','sbj17','sbj18','sbj19','sbj20',
          'sbj21','sbj22','sbj23','sbj24','sbj25','sbj26','sbj27','sbj28','sbj29','sbj30',
          'sbj31','sbj32','sbj33','sbj34','sbj35','sbj36','sbj37','sbj38','sbj39','sbj40'); # 39 sbjs
# % Excessive Headmotion (>3mm): sbj04

sbjs_reduced <- c('sbj01','sbj02','sbj03','sbj05','sbj06','sbj07','sbj09','sbj10',
                  'sbj11','sbj13','sbj14','sbj16','sbj17','sbj18','sbj19','sbj20',
                  'sbj21','sbj22','sbj23','sbj24','sbj25','sbj26','sbj27','sbj28','sbj29',
                  'sbj31','sbj32','sbj33','sbj35','sbj36','sbj37','sbj38','sbj40'); # 33 pps for this analyses
# % Missing pps:
#   % Excessive headmotion: sbj04;
# % No accept in SL and DL: sbj12;
# % No reject in SH: sbj15, sbj30;
# % No reject in SH and DH: sbj08, sbj39;
# % No reject in SH, SL, DH and DL: sbj34

rootdir <- '/gpfs/share/home/1906388343/ISC_local/Studies_CNRS_Lyon/2017_Corruption/fMRI/'
# df_ind_pm <- read.table("./outputs/estimation/7models_201906/winmodel_pm_spm_39sbjs.txt", header = F)
df_ind_pm <- read.table(paste0(rootdir,"outputs/estimation/7models_201906/winmodel_pm_spm_39sbjs.txt"), header = F)

colnames(df_ind_pm) <- c('beta_s','beta_p','theta','omega','gamma')
df_ind_pm$SbjID <- sbjs

myconfig <- theme(#axis.line.x = element_line(size = 2),
                  #axis.line.y = element_line(size = 2),
                  axis.text = element_text(size = 20),
                  axis.title.x = element_text(size = 23,vjust = -1),
                  axis.title.y = element_text(size = 23,vjust = 1),
                  legend.text = element_text(size = 20),
                  legend.title = element_text(size = 23),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 23),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.spacing = unit("0.5", "cm"),
                  plot.title = element_text(size = 25, vjust = 3),
                  plot.margin = unit(c(1, 1, 1, 1), "cm"))

# plot_savedir <- "./outputs/figures/roiplots_39sbjs/"
plot_savedir <- paste0(rootdir,"outputs/figures/roiplots_39sbjs/")

# Built-in function: mylmrob #
mylmrob <- function(x,df_input) {
  if(is.data.frame(df_input) == F) {
    df_input <- as.data.frame(df_input)
  }
  
  p_rob <- vector()
  beta_rob <- vector()
  ci_rob <- matrix(0, nrow = length(names(df_input)), ncol = 2) 
  idx <- 0
  
  for (i in 1:length(names(df_input))) {
      idx <- idx + 1
      y_tmp <- df_input[,i]
      fit1 <- lmrob(scale(y_tmp) ~ scale(x)) # data = df_input
      p_rob[idx] <- summary(fit1)$coef[2,4]
      beta_rob[idx] <- coef(fit1)[2]
      ci_rob[idx,1] <- confint(fit1)[[2,1]]
      ci_rob[idx,2] <- confint(fit1)[[2,2]]
  }
  # return(list(p_rob = p_rob,
  #             std_r = beta_rob,
  #             std_ci = ci_rob))
  
  df_lmrob_summary <- data.frame(beta_rob, p_rob, ci_rob)
  colnames(df_lmrob_summary) <- c("std_r","p_rob","ci_2.5","ci_97.5")
  #df_lmrob_summary$bene <- rep(c("S-GH","C-GH"), time = length(names(df_input))/2)
  
  return(df_lmrob_summary)
}

myggline <- function(df_input,GLMs) {
  
  if (GLMs == "GLM7e") {
    colnum <- 2
    width_tmp <- 8  # 4 rois: 10
    height_tmp <- 5 # 4 rois: 7
  } else {
    colnum <- 1
    width_tmp <- 8
    height_tmp <- 6
  }
  
  g1 <- ggplot(data = df_input, 
               aes(x = f1, y = mean, color = Condition, group = f2)) + 
        facet_wrap(~ roi_name, ncol = colnum) +
        # geom_bar(position = position_dodge(), stat = "identity", size = 1.5, color = "black") +
        geom_point(position = position_dodge(0.1), size = 1.5) +
        geom_line(size = 0.5, color = "black") +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                      width = .1, size = 1, position = position_dodge(0.1)) +
        scale_color_manual(name = "", values=c("#66FF99", "#FF3399", 
                                               "#00FF00", "#FF00FF")) + # light green/magenta
        # scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
        labs(x = "", y = "Contrast Value\n", size = "") +
        theme_bw() + myconfig +
        theme(strip.text.x = element_text(margin = margin(b = 10, t = 0)))
  print(g1)
  ggsave(paste0(plot_savedir,"cv_", GLMs, "_line.jpeg"), 
         width = width_tmp, height = height_tmp, dpi = 600)
}


############################################################################
#                           Step 2: Bar graphs                             #
############################################################################

#%%%%%%%%%%%%%%%%    ROI in GLM7e: bilateral dlPFC     %%%%%%%%%%%%%%%%#
# df_cvs <- readMat("./data_roi/cv/Corr_ROI_cv_GLM7e.mat")
df_cvs <- readMat(paste0(rootdir,"data_roi/cv/Corr_ROI_cv_GLM7e.mat"))
rois <- c("L_dlPFC","R_dlPFC")# "L_IPL","R_IPL"
roi_mat_idx <- c(1,2)# 3,4
cv_modi <- data.frame()

for (i in 1:length(rois)) {
  
  cv_tmp <- as.data.frame(unlist(df_cvs$betas[[roi_mat_idx[i]]][[1]]))
  cv_modi <- rbind(cv_modi, cv_tmp)
  
}

cv_modi$SbjID <- rep(sbjs_reduced, time = 2)
cv_modi$roi_name <- rep(rois, each = length(sbjs_reduced))
colnames(cv_modi)[1:8] <- c("dec_SH_accept","dec_SH_reject","dec_SL_accept","dec_SL_reject",
                             "dec_DH_accept","dec_DH_reject","dec_DL_accept","dec_DL_reject")
cv_modi$anti_SH <- cv_modi$dec_SH_reject - cv_modi$dec_SH_accept
cv_modi$anti_SL <- cv_modi$dec_SL_reject - cv_modi$dec_SL_accept
cv_modi$anti_DH <- cv_modi$dec_DH_reject - cv_modi$dec_DH_accept
cv_modi$anti_DL <- cv_modi$dec_DL_reject - cv_modi$dec_DL_accept

cv_modi$anti_solo <- cv_modi$anti_SL - cv_modi$anti_SH
cv_modi$anti_dyad <- cv_modi$anti_DL - cv_modi$anti_DH

cv_modi_long <- melt(cv_modi[,c(9:14)], id.vars = c("SbjID", "roi_name"))

# cv_modi_long$roi_name <- factor(cv_modi_long$roi_name, levels = rois[1:4])
cv_modi_long$f1 <- ifelse(cv_modi_long$variable == "anti_SH" | 
                          cv_modi_long$variable == "anti_SL","solo","dyad")
cv_modi_long$f2 <- ifelse(cv_modi_long$variable == "anti_SH" | 
                          cv_modi_long$variable == "anti_DH","honesty","lie")

cv_modi_long$f1 <-factor(cv_modi_long$f1, levels = c("solo","dyad"))
cv_modi_long$Condition <- substr(cv_modi_long$variable, 6, 7)
cv_modi_long$Condition <- factor(cv_modi_long$Condition, levels = c("SH","SL","DH","DL"))

cv_modi_summary <- ddply(cv_modi_long, c("roi_name","f1","f2","Condition"), summarise,
                         N = length(SbjID), mean = mean(value), 
                         sd = sd(value), se = sd/sqrt(N), .drop = FALSE)

cv_modi_summary <- cv_modi_summary[cv_modi_summary$N != 0, ] # in this way we still keep Condition
# cv_modi_summary$Condition <- factor(cv_modi_summary$Condition, levels = c("SH","SL","DH","DL"))
# myggline(cv_modi_summary[cv_modi_summary$roi_name == "L_dlPFC" |
#                            cv_modi_summary$roi_name == "R_dlPFC",],"GLM7e")
myggline(cv_modi_summary,"GLM7e")

# Make the same bar plot
g1 <- ggplot(data = cv_modi_summary, 
             aes(x = f1, y = mean, color = Condition, fill = Condition, group = f2)) + 
      facet_wrap(~ roi_name, ncol = 2) +
      geom_bar(position = position_dodge(), stat = "identity", size = 1.5) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                    width = .1, size = 1, position = position_dodge(.9)) +
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", 
                                             "#00FF00", "#FF00FF")) + # light green/magenta
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", 
                                             "#00FF00", "#FF00FF")) + # light green/magenta
      # scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
      labs(x = "", y = "Contrast Value\n", size = "") +
      theme_bw() + myconfig +
      theme(strip.text.x = element_text(margin = margin(b = 10, t = 0)))
print(g1)
ggsave(paste0(plot_savedir,"cv_GLM7e_bar.jpeg"), 
       width = 9, height = 6, dpi = 600)


# A tidy way to run post-hoc analyses
library(tidyverse)   # 1.00
post_hoc <- cv_modi_long %>%
            group_by(roi_name,f1) %>%
            nest() %>%
            mutate(fit = map(data, ~t.test(value ~ f2, paired = T, data = .x)),
                   t_val = map_dbl(fit, "statistic"),
                   p_val = map_dbl(fit, "p.value"))
# see https://drsimonj.svbtle.com/running-a-model-on-separate-groups for examples
write.table(post_hoc[,c(1,2,5,6)], paste0(plot_savedir,"post-hoc_Ttest_GLM7e.txt"))


# New figures of dlPFC: anti-corruption signals
cv_modi_long1 <- melt(cv_modi[,c(9,10,15,16)], id.vars = c("SbjID", "roi_name"))
cv_modi_long1$f1 <- ifelse(cv_modi_long1$variable == "anti_solo","solo","dyad")
cv_modi_long1$f1 <-factor(cv_modi_long1$f1, levels = c("solo","dyad"))
cv_modi_summary <- ddply(cv_modi_long1, c("roi_name","f1"), summarise,
                         N = length(SbjID), mean = mean(value), 
                         sd = sd(value), se = sd/sqrt(N), .drop = FALSE)

g1 <- ggplot(data = cv_modi_summary, 
             aes(x = roi_name, y = mean, color = f1, fill = f1, group = f1)) + 
      geom_bar(position = position_dodge(), stat = "identity", size = 1.5, alpha = 0.75) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                    width = .1, size = 1, position = position_dodge(.9)) +
      scale_color_manual(name = "", values=c("white", "white")) + # light green/magenta
      scale_fill_manual(name = "", values=c("skyblue", "blue")) + # light green/magenta
      # scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
      labs(x = "", y = "Anti-Corruption Signal\n", size = "") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"cv_GLM7e_bar_2conds.jpeg"), 
       width = 8, height = 5, dpi = 600)

post_hoc <- cv_modi_long1 %>%
            group_by(roi_name,f1) %>%
            nest() %>%
            mutate(fit = map(data, ~t.test(.x$value, mu = 0)),
                   t_val = map_dbl(fit, "statistic"),
                   p_val = map_dbl(fit, "p.value"))


# New figures of dlPFC: raw signals in 8 conditions
cv_modi_long2 <- melt(cv_modi[,c(1:10)], id.vars = c("SbjID", "roi_name"))
cv_modi_long2$roi_name <- factor(cv_modi_long2$roi_name, levels = c("L_dlPFC","R_dlPFC"))
cv_modi_long2$Context <- ifelse(cv_modi_long2$variable == "dec_SH_accept" | 
                                  cv_modi_long2$variable == "dec_SH_reject" |
                                  cv_modi_long2$variable == "dec_SL_accept" | 
                                  cv_modi_long2$variable == "dec_SL_reject","solo","dyad")
cv_modi_long2$Sender_Act <- ifelse(cv_modi_long2$variable == "dec_SH_accept" | 
                                     cv_modi_long2$variable == "dec_SH_reject" |
                                     cv_modi_long2$variable == "dec_DH_accept" | 
                                     cv_modi_long2$variable == "dec_DH_reject","honesty","lie")

cv_modi_long2$Decision <- ifelse(cv_modi_long2$variable == "dec_SH_accept" | 
                                   cv_modi_long2$variable == "dec_SL_accept" |
                                   cv_modi_long2$variable == "dec_DH_accept" | 
                                   cv_modi_long2$variable == "dec_DL_accept","accept","reject")

cv_modi_long2$Context <-factor(cv_modi_long2$Context, levels = c("solo","dyad"))
cv_modi_long2$Condition <- substr(cv_modi_long2$variable, 5,6)

cv_modi_summary <- ddply(cv_modi_long2, c("roi_name","Decision","Context","Sender_Act"), summarise,
                          N = length(SbjID), mean = mean(value), 
                          sd = sd(value), se = sd/sqrt(N), .drop = FALSE)

cv_modi_summary$variable <- paste0(substr(cv_modi_summary$Context,1,1),
                                   substr(cv_modi_summary$Sender_Act,1,1),"_",
                                   cv_modi_summary$Decision)
cv_modi_summary$variable <- factor(cv_modi_summary$variable,
                                    levels = c("sh_accept","sh_reject",
                                               "sl_accept","sl_reject",
                                               "dh_accept","dh_reject",
                                               "dl_accept","dl_reject"))

cv_modi_summary$Condition <- paste0(substr(cv_modi_summary$Context,1,1),
                                    substr(cv_modi_summary$Sender_Act,1,1))
cv_modi_summary$Condition <- factor(cv_modi_summary$Condition,
                                    levels = c("sh","sl","dh","dl"))

g1 <- ggplot(data = cv_modi_summary, 
             aes(x = Condition, y = mean, fill = variable, color = variable)) + 
  facet_wrap(~ roi_name, ncol = 2) +
  geom_bar(position = position_dodge(), stat = "identity", size = 1.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                width = .1, size = 1, position = position_dodge(.9)) +
  scale_fill_manual(name = "", values=c("#66FF99", "white", "#FF3399", "white", 
                                        "#00FF00", "white", "#FF00FF", "white")) + # light green/magenta
  scale_color_manual(name = "", values=c("#66FF99", "#66FF99", "#FF3399", "#FF3399",
                                         "#00FF00", "#00FF00", "#FF00FF", "#FF00FF")) + # dark green/magenta
  # scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
  labs(x = "", y = "Contrast Value\n", size = "") +
  theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"cv_GLM7e_bar_3ways.jpeg"), 
       width = 12, height = 5, dpi = 600)

cv_modi_summary$variable <- factor(cv_modi_summary$variable,
                                   levels = c("sh_accept","sl_accept",
                                              "sh_reject","sl_reject",
                                              "dh_accept", "dl_accept",
                                              "dh_reject","dl_reject"))
g1 <- ggplot(data = cv_modi_summary, 
             aes(x = Decision, y = mean, fill = variable, color = variable)) + 
  facet_wrap(roi_name ~ Context, ncol = 2) +
  geom_bar(position = position_dodge(), stat = "identity", size = 1.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                width = .1, size = 1, position = position_dodge(.9)) +
  scale_fill_manual(name = "", values=c("#66FF99", "#FF3399","white","white",
                                        "#00FF00", "#FF00FF","white","white")) + # light green/magenta
  scale_color_manual(name = "", values=c("#66FF99", "#FF3399", "#66FF99", "#FF3399", 
                                         "#00FF00", "#FF00FF","#00FF00", "#FF00FF")) + # dark green/magenta
  # scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
  labs(x = "", y = "Contrast Value\n", size = "") +
  theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"cv_GLM7e_bar_3ways_4panels.jpeg"), 
       width = 11, height = 9, dpi = 600)

# detach("package:plyr", unload=TRUE) # if you load the plyr, it is better to detach it first
# library(tidyverse)
cv_modi_long2 <- cv_modi_long2 %>% mutate(SbjID = as.factor(SbjID))
cv_solo <- cv_modi_long2 %>% filter(Context == "solo")
cv_dyad <- cv_modi_long2 %>% filter(Context == "dyad")
aov_fit_solo_LdlPFC <- aov(value ~ Decision*Sender_Act + Error(SbjID/Decision*Sender_Act),
                          data = cv_solo[cv_solo$roi_name == "L_dlPFC",]) 
aov_fit_dyad_LdlPFC <- aov(value ~ Decision*Sender_Act + Error(SbjID/Decision*Sender_Act),
                           data = cv_dyad[cv_dyad$roi_name == "L_dlPFC",]) 
aov_fit_solo_RdlPFC <- aov(value ~ Decision*Sender_Act + Error(SbjID/Decision*Sender_Act),
                           data = cv_solo[cv_solo$roi_name == "R_dlPFC",]) 
aov_fit_dyad_RdlPFC <- aov(value ~ Decision*Sender_Act + Error(SbjID/Decision*Sender_Act),
                           data = cv_dyad[cv_dyad$roi_name == "R_dlPFC",])

post_hoc <- cv_modi_long2 %>%
  group_by(roi_name, Context, Decision) %>%
  nest() %>%
  mutate(fit = map(data, ~t.test(.x$value ~ .x$Sender_Act, paired = T)),
         t_val = map_dbl(fit, "statistic"),
         p_val = map_dbl(fit, "p.value"))

#%%%%%%%%%%%%%%%%    ROI in GLM8d: dACC     %%%%%%%%%%%%%%%%#
df_cvs <- readMat(paste0(rootdir,"data_roi/cv/Corr_ROI_cv_GLM8d.mat"))
rois <- c("dACC")
roi_mat_idx <- c(1)

#%%%%%%%%%%%%%%%%    ROI in GLM813: bilateral vAI     %%%%%%%%%%%%%%%%#
df_cvs <- readMat(paste0(rootdir,"data_roi/cv/Corr_ROI_cv_GLM13e.mat"))
rois <- c("lvAI","rvAI")
roi_mat_idx <- c(1,2)
cv_modi <- data.frame()
for (i in 1:length(rois)) {
  cv_tmp <- as.data.frame(unlist(df_cvs$betas[[roi_mat_idx[i]]][[1]]))
  cv_modi <- rbind(cv_modi, cv_tmp)
}

cv_modi$SbjID <- rep(sbjs, time = length(rois))
cv_modi$roi_name <- rep(rois, each = length(sbjs))
colnames(cv_modi)[1:4] <- c("dec_SH","dec_SL","dec_DH","dec_DL")

cv_modi_long <- melt(cv_modi, id.vars = c("SbjID", "roi_name"))
cv_modi_long$f1 <- ifelse(cv_modi_long$variable == "dec_SH" | 
                            cv_modi_long$variable == "dec_SL","solo","dyad")
cv_modi_long$f2 <- ifelse(cv_modi_long$variable == "dec_SH" | 
                            cv_modi_long$variable == "dec_DH","honesty","lie")
cv_modi_long$f1 <- factor(cv_modi_long$f1, levels = c("solo","dyad"))

cv_modi_long$Condition <- substr(cv_modi_long$variable, 5, 6)
cv_modi_long$Condition <- factor(cv_modi_long$Condition, levels = c("SH","SL","DH","DL"))

cv_modi_summary <- ddply(cv_modi_long, c("f1","f2","roi_name"), summarise,
                         N = length(SbjID), mean = mean(value), 
                         sd = sd(value), se = sd/sqrt(N), .drop = FALSE)
cv_modi_summary$Condition <- ifelse(cv_modi_summary$f1 == "solo",ifelse(cv_modi_summary$f2 == "honesty","SH","SL"),
                                    ifelse(cv_modi_summary$f2 == "honesty","DH","DL"))
cv_modi_summary$Condition <- factor(cv_modi_summary$Condition, levels = c("SH","SL","DH","DL"))

g1 <- ggplot(data = cv_modi_summary, 
             aes(x = f1, y = mean, color = Condition, fill = Condition, group = f2)) + 
      facet_wrap(~ roi_name, ncol = 2) +
      geom_bar(position = position_dodge(), stat = "identity", size = 1.5) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                    width = .1, size = 1, position = position_dodge(.9)) +
      # geom_point(data = cv_modi_long[cv_modi_long$SbjID != "sbj34",], 
      #            aes(x = f1, y = value, group = f2, size = 1),
      #            position = position_dodge(0.5), color = "black", alpha = 0.2) +
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", 
                                             "#00FF00", "#FF00FF")) + # light green/magenta
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", 
                                            "#00FF00", "#FF00FF")) + # light green/magenta
      # scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
      labs(x = "", y = "Contrast Value\n", size = "") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"cv_GLM8d_bar.jpeg"), 
       width = 9, height = 6, dpi = 600)

ggsave(paste0(plot_savedir,"cv_GLM13e_bar.jpeg"), 
       width = 12, height = 6, dpi = 600)

############################################################################
#                         Step 2: Scatter Plots                            #
############################################################################

#%%%%%%%%%%%%%%%%    ROI in GLM8d: gPPI: vmPFC-dlPFC ~ theta  %%%%%%%%%%%%%%%%#
rm(list =  c("df_cvs","cv_modi"))
df_cvs <- readMat("./data_roi/cv/Corr_gPPI_theta_GLM8d_ROI.mat")
rois <- c("L_dlPFC","R_dlPFC")
roi_mat_idx <- c(1,2)

#%%%%%%%%%%%%%%%%    ROI in GLM8d: gPPI: control analyses  %%%%%%%%%%%%%%%%#
df_cvs <- readMat("./data_roi/cv/Corr_gPPI_thetaomega_GLM8d_ROI.mat")
rois <- c("R_dlPFC")
roi_mat_idx <- c(1)

cv_modi <- data.frame()

for (i in 1:length(rois)) {
  
  cv_tmp <- as.data.frame(unlist(df_cvs$betas[[roi_mat_idx[i]]][[1]]))
  cv_modi <- rbind(cv_modi, cv_tmp)
  
}

cv_modi$SbjID <- rep(sbjs, time = length(rois))
cv_modi$roi_name <- rep(rois, each = length(sbjs))
colnames(cv_modi)[1:4] <- c("dec_SH","dec_SL","dec_DH","dec_DL")
cv_modi$honesty <- (cv_modi$dec_SH + cv_modi$dec_DH)*0.5
cv_modi$lie <- (cv_modi$dec_SL + cv_modi$dec_DL)*0.5
cv_modi$LH <- cv_modi$lie - cv_modi$honesty
cv_modi <- merge(cv_modi, df_ind_pm[,c(3,6)],by = c("SbjID"))

# Difference
cv_modi_wide <- dcast(cv_modi[,c(1,6,9,10)],SbjID + theta ~ roi_name, value.var="LH")

lmrob_summary <-mylmrob(cv_modi_wide$theta, cv_modi_wide[,c(3)]) #cv_modi_wide[,c(3,4)]
print(lmrob_summary)
lmrob_summary$roi_name <- rep(rois, each = 1)
# lmrob_summary$x_pos <- rep(c(8), time = length(rois))
# lmrob_summary$y_pos <- rep(c(1.7), time = length(rois))
lmrob_summary$x_pos <- rep(c(8,8), time = 1)
lmrob_summary$y_pos <- rep(c(1.7,1.9), time = 1)
lmrob_summary$condition <- rep(c("lie vs. honesty"), time = length(rois))

for (i in 1:length(rois)) {
  
  g1 <- ggplot(data = cv_modi[cv_modi$roi_name == rois[i],], 
               aes(x = theta, y = LH)) +
        # facet_wrap(~ roi_name, ncol = length(rois)) +
        geom_point(shape = 16, size = 5, alpha = 0.3) +
        geom_smooth(method = "rlm", se = T) +
        # scale_x_continuous(limits=c(-2,12), breaks = seq(-2,12,2)) +
        labs(x = "\ntheta", y = "PPI: Lie vs. Honesty\n(Contrast Value)") +
        theme_bw() + myconfig +
        theme(plot.title = element_text(size = 15, hjust = 0.5, vjust = 3))
  
  # fig_name <- paste0("scatter_PPI_cv_theta_vmPFC-",rois[i],"_LH.jpeg")
  fig_name <- paste0("scatter_PPI_cv_thetaomega_vmPFC-",rois[i],"_LH.jpeg")

  ggsave(paste0(plot_savedir,fig_name),
         ggMarginal(g1, xparams = list(fill = "grey", alpha = 0.3), 
                        yparams = list(fill = "grey", alpha = 0.3)), 
         width = 6, height = 5, dpi = 600)
}


# Seperate conditions 
cv_modi_wide_H <- dcast(cv_modi[,c(1,6,7,10)],SbjID + theta ~ roi_name,
                        value.var = "honesty")
colnames(cv_modi_wide_H)[3:4] <- c("L_dlPFC_H","R_dlPFC_H")
cv_modi_wide_L <- dcast(cv_modi[,c(1,6,8,10)],SbjID + theta ~ roi_name,
                        value.var = "lie")
colnames(cv_modi_wide_L)[3:4] <- c("L_dlPFC_L","R_dlPFC_L")
cv_modi_wide1 <- merge(cv_modi_wide_H, cv_modi_wide_L, by = c("SbjID","theta"))

lmrob_summary <-mylmrob(cv_modi_wide1$theta,cv_modi_wide1[,c(3:6)])
print(lmrob_summary)

lmrob_summary$roi_name <- rep(rois, time = 2)
lmrob_summary$x_pos <- rep(c(8,8), time = length(rois))
lmrob_summary$y_pos <- rep(c(-1.8,-2.2), each = length(rois))
lmrob_summary$condition <- rep(c("honesty","lie"), each = length(rois))

cv_modi_long <- melt(cv_modi[,c(1,6:8,10)], 
                     id.vars = c("SbjID","roi_name", "theta"))
colnames(cv_modi_long)[4] <- "condition"

g1 <- ggplot(data = cv_modi_long, aes(x = theta, y = value,
                                      group = condition, color = condition, fill = condition)) +
      facet_wrap(~roi_name, ncol = 2) +
      geom_point(shape = 16, size = 5, alpha = 0.3) + 
      geom_smooth(method = "rlm", se = T) +
      geom_text(data = lmrob_summary,
                aes(label = sprintf("r = %0.3f, p = %0.3e",round(std_r, digits = 3),p_rob),
                    x = x_pos, y = y_pos), size = 5, parse = F) +
      scale_fill_manual(name = "", values=c("cyan", "magenta")) + # light green/magenta
      scale_color_manual(name = "", values=c("cyan", "magenta")) + # dark green/magenta
      labs(x = "\ntheta", y = "PPI\n") +
      # scale_y_continuous(limits=c(-3,3), breaks = seq(-3,3,1)) +
      # ggtitle("Sender's gain due to bribery: Left Insula") +
      theme_bw() + myconfig + 
      theme(plot.title = element_text(size = 15, hjust = 0.5, vjust = 3))
print(g1)

# require(ggExtra)
# ggMarginal(g1, groupColour = TRUE, groupFill = TRUE)

ggsave(paste0(plot_savedir,"scatter_PPI_cv_theta_vmPFC-dlPFC_sep.jpeg"),
       width = 10, height = 6, dpi = 600)
