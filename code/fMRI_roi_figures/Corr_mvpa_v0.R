#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Corr_mvpa_v0.R                                                     #
# Author: Yang Hu                                                    #
# Established Date: May. 17, 2018                                    #
# Last Change: Dec. 16, 2020                                         #
# Project: Corruption (Corr)                                         #
# PI: Yang Hu                                                        #
# Collaborators: Chen Qu, Jean-Claude Dreher                         #
# This script aims to manage permutation-based decoding ROI plots,   #
# and pattern similarity plots                                       #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

######################################################################
###                      Step 1: Preparation                       ###
######################################################################
# required packages
library(tidyverse)   # 1.00
library(stringr)
# library(lme4)  # lmer / glmer  1.1-16
# library(effsize) # cohens d calculator  0.7.1
# library(optimx) # additional optimizers for glmer  2013.8.7
library(psych)
library(readxl)
# library(sjPlot)
library(broom)
library(pander)
library(cowplot)
library(emmeans)
library(lmPerm) # permutation-based linear regression
library(EMAtools) # Cohen's D for each effect in an lme4 object
require(ggbeeswarm)
require(lmerTest) # post-edited by Yang Hu: to get the p-value for lmer fit
require(car) # post-edited by Yang Hu: Anova
require(scales)

# source("./code/code_R/Corr_tDCS_custom_funs.R")
# source("./code/code_R/tidy_lmer.R")

myconfig <- theme(#axis.line.x = element_line(size = 2),
                  #axis.line.y = element_line(size = 2),
                  axis.text = element_text(size = 15),
                  axis.title.x = element_text(size = 20,vjust = -1),
                  axis.title.y = element_text(size = 20,vjust = 1),
                  legend.text = element_text(size = 15),
                  legend.title = element_text(size = 15),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 15),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.spacing = unit("0.5", "cm"),
                  plot.title = element_text(size = 23, hjust = 0.5, vjust = 3),
                  plot.margin = unit(c(1, 1, 1, 1), "cm"))

# data_path <- "./data/backup/"
rootdir <- "/gpfs/share/home/1906388343/ISC_local/Studies_CNRS_Lyon/2017_Corruption/fMRI"
# plot_savedir <- "./outputs/figures/rsaplots/"
plot_savedir <- paste0(rootdir,"/outputs/figures/roiplots_39sbjs/")
# GLM_label <- c("GLM8d","GLM7e")
# rois <- c("TPJ","dlPFC")

######################################################################
###           Step 2: Plots for SVM Permutation Results            ###
######################################################################
df_svm <- data.frame()
# svm_file_name <- c("GLM8d_bribe_S_TPJ","GLM8d_bribe_D_TPJ",
#                    "GLM7e_bribe_S_ar_dlPFC","GLM7e_bribe_D_ar_dlPFC")

svm_file_name <- c("GLM8d_bribe_S_lTPJ","GLM8d_bribe_D_lTPJ",
                   "GLM8d_bribe_S_rTPJ","GLM8d_bribe_D_rTPJ") # in response to elife reviewer

for (i in 1:length(svm_file_name)) {
    df_tmp <- read.table(paste0(rootdir,"/data_roi/decoding/",svm_file_name[i],"_svm_acc_perm5000.csv"), 
                         header = T, sep = ",") %>%
              mutate(N_perm = seq(1,length(ROI)))
    df_svm <- rbind(df_svm,df_tmp)
}

# df_svm_TPJ_r <- df_svm %>%
#                 select(1,2,8:13) %>%
#                 filter(ROI == "TPJ")
# 
# df_svm_dlPFC_r <- df_svm %>%
#                   select(1,2,8:13) %>%
#                   filter(ROI == "dlPFC")
# df_svm_r <- df_svm_TPJ_r
# df_svm_r <- df_svm_dlPFC_r

df_svm <- df_svm %>% mutate(ROI_new = rep(c("lTPJ","rTPJ"), each = 5000*2)) # the old name is always TPJ; we update it
df_svm_lTPJ_r <- df_svm %>%
                 select(1,14,8:13) %>%
                 filter(ROI_new == "lTPJ")

df_svm_rTPJ_r <- df_svm %>%
                 select(1,14,8:13) %>%
                 filter(ROI_new == "rTPJ")

# df_svm_r <- df_svm_lTPJ_r
df_svm_r <- df_svm_rTPJ_r

svm_summary <- df_svm_r %>%
               group_by(Comparison, ROI_new) %>%
               summarise(N = n(), mean = mean(acc_fc), se = mean(se_fc)) %>%
               ungroup() %>%
               rename(acc_fc_perm = 4) # for the plots later

g1 <- ggplot(data = df_svm_r, aes(x = Comparison, y = acc_fc_perm, fill = Comparison, color = Comparison)) +
      # facet_wrap( ~ pars, ncol = 2) +
      geom_violin(position=position_dodge(1), alpha = 0.5) +
      geom_dotplot(binaxis='y', stackdir='center',binwidth = 1/250,
                   position = position_dodge(1), alpha = 0.2, dotsize = 0.1, fill = "black") + #
      # geom_jitter(position=position_jitter(0.4),size = 0.1) +
      geom_point(data = svm_summary, aes(x = Comparison, y = acc_fc_perm),
                 position = position_dodge(1), stat = "identity", size = 2, color = "black") +
      geom_errorbar(data = svm_summary, aes(ymin = acc_fc_perm - se, ymax = acc_fc_perm + se), color = "black",
                    width = .1, size = 1, position = position_dodge(1)) +
      geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
      scale_fill_manual(name = "", values=c("white", "white")) + # light magenta/green
      scale_color_manual(name = "", values=c("skyblue", "blue")) + # dark magenta/green
      scale_x_discrete(labels = c("Solo","Dyad"))+
      scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
      labs(x = "", y = "decoding accuracy\n") +
      theme_bw() + myconfig
print(g1)
# ggsave(paste0(plot_savedir,"GLM8d_TPJ_svm_acc.jpeg"), width = 8, height = 6, dpi = 600)
# ggsave(paste0(plot_savedir,"GLM7e_dlPFC_svm_acc.jpeg"), width = 8, height = 6, dpi = 600)
ggsave(paste0(plot_savedir,"GLM8d_lTPJ_svm_acc.jpeg"), width = 8, height = 6, dpi = 600)
ggsave(paste0(plot_savedir,"GLM8d_rTPJ_svm_acc.jpeg"), width = 8, height = 6, dpi = 600)

######################################################################
###              Step 3: Plots for Pattern Similarity              ###
######################################################################
df_patternsim <- data.frame()
patternsim_file_name <- c("GLM16_group_pattern_sim_SVdiff_HL")

for (i in 1:length(patternsim_file_name)) {
  df_tmp <- read.table(paste0(rootdir,"/data_roi/patternsim/",patternsim_file_name[i],"_perm5000.csv"), 
                       header = T, sep = ",")
  df_patternsim <- rbind(df_patternsim,df_tmp)
}

r_summary <- df_patternsim %>% group_by(ROI) %>%
             summarise(N = n(), mean = mean(PatternSim), 
                   sd = sd(PatternSim), se = sd/sqrt(N)) %>%
             ungroup()
g1 <- df_patternsim %>% ggplot(aes(x = PatternSim)) +
      geom_histogram(binwidth = 0.1, alpha = 0.5,
                     #               #position = "dodge",
                     position = "identity",fill = "grey",
                     color = "grey") +
      # geom_density(alpha = 0.2,color = "skyblue") + 
      geom_vline(data = r_summary, color = "red",
                 aes(xintercept = mean), linetype = "dashed", size = 1) +
      scale_x_continuous(limits = c(-0.8, 1),breaks = seq(-0.8, 1, 0.2)) +
      #scale_colour_manual(name = "", values = c("#FF3399", "#66FF99")) +
      #scale_fill_manual(name = "", values = c("#FF3399", "#66FF99")) +
      xlab("\nPattern Similarity") +
      ylab("Frequency\n") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir, "GLM16_vmPFC_SVdiff_patternsim.jpeg"), width = 8, height = 5, dpi = 600)


