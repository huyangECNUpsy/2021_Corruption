#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Corr_behav_v4.R                                                    #
# Author: Yang Hu                                                    #
# Date: Nov. 15, 2017                                                #
# Date of last change: Dec. 15, 2020                                 #
# Project: corruption (fMRI study)                                   #
#                                                                    #
# This script aims to make plots and conduct statistics for          #
# behavioral data (i.e., choice and decision time).                  #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

######################################################################
###                       Step 1: Preparations                     ###
######################################################################

#### Load relevant packages ###
# library(reshape2)
library(ggplot2)
# library(plyr)
library(grid)
library(psych)
library(readxl)
library(corrplot)
require(scales)
library(sjPlot)
library(lme4)
library(car)
library(EMAtools) # To compute the cohen's d for each effect in an lme4 object (lme.dscore)
require(lmerTest)
require(lsmeans)
require(optimx)
# install.packages("phia") # post-hoc interaction analysis
# library(phia)

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

####   Prepare the data   ####
rootdir <- '/gpfs/share/home/1906388343/ISC_local/Studies_CNRS_Lyon/2017_Corruption/fMRI/'
# df_merge <- read.table("./data_behav/Corr_merge_40sbjs_fMRI.csv", header = T, 
#                        sep = ",", comment.char = "$")
df_merge <- read.table(paste0(rootdir,"data_behav/Corr_merge_40sbjs_fMRI.csv"), header = T, 
                       sep = ",", comment.char = "$")
# df_merge_run1 <- df_merge[df_merge$RunNo == 1,]
# df_merge_run2 <- df_merge[df_merge$RunNo == 2,]
# df_input <- df_merge_run2
# table(df_input$SbjID,df_input$Condition,df_input$Decision)
table(df_merge$Condition,df_merge$Decision)
# drop invalid decisions and participant (sbj04; due to excessive headmotion) 
df_valid <- df_merge[df_merge$Decision != "invalid" & df_merge$SbjID != 4, ]
colnames(df_valid)[9:10]<-c("f1","f2")
df_valid$SbjID <- as.factor(df_valid$SbjID)
df_valid$f1 <- factor(df_valid$f1, levels = c("solo","dyad"))
df_valid$f2 <- as.factor(df_valid$f2)
df_valid$Condition <- factor(df_valid$Condition, levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))
df_valid$Decision <- factor(df_valid$Decision, levels = c("accept","reject"))
df_valid$IsAccept <- ifelse(df_valid$Decision == "accept", 1, 0)
df_valid$Bribeprop <- round(df_valid$Sender_Bribe*10/df_valid$Sender_Payoff)/10
sbjs  <- levels(df_valid$SbjID)
sbjno <- nlevels(df_valid$SbjID)

# plot_savedir <- "./outputs/figures/behav_39sbjs/"
plot_savedir <- paste0(rootdir,"outputs/figures/behav_39sbjs/")
options(scipen = 999) # inhibit displaying in scientific notation

######################################################################
###    Step 2: Making descriptive plots of choice proportion (%)   ###
######################################################################

#### summarize: within each pp ####
choice_ind_summary <- ddply(df_valid,c("SbjID","f1","f2","Decision"), summarise,
                            N = length(Sender), .drop = FALSE)
choice_ind_summary$prop <- choice_ind_summary$N/36 # 36 trials per condition

# # Individual acceptance proportion (for fMRI correlation analyses)
# accept_ind_summary <- choice_ind_summary %>%
#                       filter(Decision == "accept") %>%
#                       select(-Decision) %>%
#                       mutate(Condition = paste0(f1,'_',f2)) %>%
#                       select(-f1, -f2)
# accept_ind_summary_wide <- accept_ind_summary %>%
#                            select(-N) %>%
#                            group_by(SbjID) %>%
#                            spread(Condition, prop) %>%
#                            ungroup() %>%
#                            mutate(TL_solo = solo_honesty - solo_lie ,
#                                   TL_dyad = dyad_honesty - dyad_lie,
#                                   TL_all = (TL_dyad + TL_solo)/2,
#                                   inter = TL_dyad - TL_solo,
#                                   ave_all = (dyad_honesty + dyad_lie + solo_honesty + solo_lie)/4)
# write.csv(accept_ind_summary_wide, paste0(rootdir,"data_behav/Corr_acceptprop_39sbjs.csv"),
#           row.names = F)

choice_ind_summary$Condition <- paste(choice_ind_summary$f1,choice_ind_summary$f2,sep = "_")
choice_ind_summary$Condition <- factor(choice_ind_summary$Condition, 
                                       levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))
choice_ind_summary$Count <- ave(choice_ind_summary$SbjID, 
                                choice_ind_summary[,c("Condition","Decision","prop")], 
                                FUN = length)
choice_ind_summary$Count <- as.numeric(levels(choice_ind_summary$Count))[choice_ind_summary$Count] # from factor to numeric


#### summarize: across each pp ####
choice_summary <- ddply(choice_ind_summary, c("f1","f2","Decision"), summarise,
                        N = length(SbjID), mean = mean(prop), 
                        sd = sd(prop), se = sd/sqrt(N), .drop = FALSE)
choice_summary$Condition <- paste(choice_summary$f1, choice_summary$f2, sep = "_")
choice_summary$Condition <- factor(choice_summary$Condition, 
                                   levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))

#### figure 1: histogram of mean accept rate ####
g1 <- ggplot(data = choice_ind_summary[choice_ind_summary$Decision == "accept", ], 
             aes(x = prop, group = Condition, color = Condition, fill = Condition)) +
      #geom_histogram(binwidth = 0.1, alpha = 0.5, 
      #               #position = "dodge",
      #               position = "identity") +
      geom_density(alpha = 0.2) + 
      geom_vline(data = choice_summary[choice_summary$Decision == "accept", ], 
                 aes(xintercept = mean,  color = Condition),
                 linetype = "dashed", size = 1) +
      scale_x_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
      scale_colour_manual(name = "", values = c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) +
      scale_fill_manual(name = "", values = c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) +
      xlab("\nAcceptance (%)") +
      ylab("Density\n") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"acceptprop_density.jpeg"), width = 8, height = 5, dpi = 600)

#### figure 2: bar plot (with individual data on) of mean accept rate ####
g1 <- ggplot(data = choice_summary[choice_summary$Decision == "accept", ], 
             aes(x = f1, y = mean, color = Condition, fill = Condition)) + 
      # facet_wrap(~ Decision, ncol = 2) +
      geom_bar(position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                    width = .1, size = 1, position = position_dodge(.9)) +
      geom_point(data = choice_ind_summary[choice_summary$Decision == "accept", ], 
                 aes(x = f1, y = prop, group = f2, size = Count),
                 position = position_dodge(0.5), color = "black", alpha = 0.2) +
      # geom_violin(data = df_tmp, aes(x = f1, y = Score), 
      #             trim = TRUE, size = 1,
      #             position = position_dodge(.9), color = "black", alpha = 0.5) +
      # scale_size_manual(name = "", values=c(1:7)) +
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # light green/magenta
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # dark green/magenta
      scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
      labs(x = "", y = "Acceptance (%)\n", size = "") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"acceptprop_bar.jpeg"), width = 10, height = 6, dpi = 600)

#### figure 3: heat map of mean accept rate as a function of bribeprop and pool ####
# sbjno <- nlevels(df_valid$SbjID)


# # Some correlation between stimuli (payoffs or inequity)
library(tidyverse)
stim <- df_valid %>%
        filter(SbjID == 1 & Condition == "solo_honesty") %>%
        dplyr::select(13,16,17,71) %>%
        mutate(inequity = Sender_Bribe - Sender_Kept,
               inequity_abs = abs(inequity),
               potential_loss = Sender_Payoff - (100 - Sender_Payoff)) %>%
        arrange(Sender_Payoff,Sender_Bribe) # in terms of the order of RSA
# corr.test(stim)
# Correlation matrix 
#                Sender_Payoff Sender_Kept Sender_Bribe Bribeprop inequity inequity_abs potential_loss
# Sender_Payoff           1.00        0.01         0.55      0.42     0.31         0.17           1.00
# Sender_Kept             0.01        1.00        -0.83     -0.89    -0.95         0.33           0.01
# Sender_Bribe            0.55       -0.83         1.00      0.98     0.96        -0.18           0.55
# Bribeprop               0.42       -0.89         0.98      1.00     0.99        -0.28           0.42
# inequity                0.31       -0.95         0.96      0.99     1.00        -0.26           0.31
# inequity_abs            0.17        0.33        -0.18     -0.28    -0.26         1.00           0.17
# potential_loss          1.00        0.01         0.55      0.42     0.31         0.17           1.00

# print(corr.test(stim[,c(2,3,7)]), short = F) 
# Sender_Kept (expected gains for the proposer);
# Sender_Bribe (expected gains for the pp);
# potential_loss (expected losses to the third party)

df_valid$Bribeprop <- as.factor(df_valid$Bribeprop)
df_valid$Bribeint <- ifelse(df_valid$Bribeprop <= 0.5, "weak", "strong") # Bribeint: intensity of bribe
df_valid$Bribeint <- factor(df_valid$Bribeint, levels = c("weak", "strong"))
df_valid$Sender_Payoff <- as.factor(df_valid$Sender_Payoff)
choice_summary0 <- ddply(df_valid,c("f1","f2","Bribeprop","Sender_Payoff"), summarise,
                         N_valid = length(Sender), .drop = FALSE) # to make sure the valid sample size of each cell
choice_summary1 <- ddply(df_valid,c("f1","f2","Decision","Bribeprop","Sender_Payoff"), summarise,
                         N = length(Sender), .drop = FALSE)
choice_summary1 <- merge(choice_summary1, choice_summary0, by = c("f1","f2","Bribeprop","Sender_Payoff"))
choice_summary1 <- choice_summary1[choice_summary1$N_valid != 0,]
choice_summary1$prop <- choice_summary1$N/choice_summary1$N_valid
# valid_idx <- c(rep(1, time = 12), 0, rep(1, time = 5), 0, rep(1, time = 5),
#                0, 0, rep(1, time = 4), 0, 0, rep(1, time = 4),
#                rep(0, time = 3), rep(1, time = 3), rep(0, time = 4), rep(1, time = 2),
#                rep(0, time = 5), 1)
# choice_summary1$valid <- rep(valid_idx, time = 8)
g1 <- ggplot(data = choice_summary1[choice_summary1$Decision == "accept", ], 
             aes(x = Sender_Payoff, y = Bribeprop)) +
      facet_wrap(f1 ~ f2, ncol = 2) +
      geom_tile(aes(fill = prop)) +
      scale_y_discrete(labels = c("10%","20%","30%","40%","50%","60%","70%","80%","90%"))+
      scale_fill_gradient2(low = "green", high = "red", mid = "yellow", 
                           midpoint = 0.5, limit = c(0,1), space = "Lab", 
                           name = "Acceptance\n(Corruption)\n", labels = percent) +
      labs(x = "\nPool (CNY)", y = "Bribe Proportion\n") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"acceptprop_heatmap.jpeg"), width = 10, height = 6, dpi = 600)

#### figure 4: line + ribbon plot of accept rate as a function of bribeprop ####

choice_ind_summary1 <- ddply(df_valid,c("SbjID","f1","f2","Decision","Bribeprop"), summarise,
                             N = length(Sender), .drop = FALSE)
choice_ind_summary2 <- ddply(df_valid,c("SbjID","f1","f2","Bribeprop"), summarise,
                             N_valid = length(Sender), .drop = FALSE) # to make sure the valid sample size of each bribeprop
choice_ind_summary1 <- merge(choice_ind_summary1, choice_ind_summary2, by = c("SbjID","f1","f2","Bribeprop"))
choice_ind_summary1$prop <- choice_ind_summary1$N/choice_ind_summary1$N_valid

choice_summary2 <- ddply(choice_ind_summary1[choice_ind_summary1$Decision == "accept",],
                         c("f1","f2","Bribeprop"), summarise,
                         N = length(SbjID), mean = mean(prop), 
                         sd = sd(prop), se = sd/sqrt(N), .drop = FALSE)
choice_summary2$Condition <- paste(choice_summary2$f1, choice_summary2$f2, sep = "_")
choice_summary2$Condition <- factor(choice_summary2$Condition, 
                                        levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))
g1 <- ggplot(data = choice_summary2, 
             aes(x = Bribeprop, y = mean, group = Condition, color = Condition, fill = Condition)) + 
      #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
      #              width = .1, size = 1, position = position_dodge()) +
      geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.3) + 
      geom_line(size = 2) +
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # light green/magenta
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # dark green/magenta
      scale_x_discrete(labels = c("10%","20%","30%","40%","50%","60%","70%","80%","90%")) +
      scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
      labs(x = "\nBribe Proportion", y = "Acceptance(%)\n", size = "") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"acceptprop_line.jpeg"), width = 9, height = 6, dpi = 600)

# choice_summary1_wide <- dcast(choice_summary1, f1 + f2 + Bribeprop + Sender_Payoff + valid + Bribeprop1 ~ Decision, value.var="prop")

#### figure 5: sigmoid function predicting accept rate by bribeprop ####
choice_summary2$intensity <- as.numeric(levels(choice_summary2$Bribeprop))[choice_summary2$Bribeprop] # from factor to numeric
choice_summary2$IsAccept <- choice_summary2$mean 
df_valid$intensity <- as.numeric(levels(df_valid$Bribeprop))[df_valid$Bribeprop] # from factor to numeric

g1 <- ggplot(data = df_valid, aes(y = IsAccept, x = intensity, 
                                  group = Condition, color = Condition, fill = Condition)) +
      geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                  fullrange = T, se = T) + 
      geom_point(data = choice_summary2, aes(x = intensity, y = mean), 
                 shape = 1, size = 4) + # position = position_dodge(0.9)
      geom_errorbar(data = choice_summary2, aes(x = intensity, ymin = mean - se, ymax = mean + se), 
                    width = 0, size = 1) +
      scale_x_continuous(limits = c(0.1,0.9), breaks = seq(0.1,0.9,0.1),
                         labels = percent) + # oob = rescale_none
      scale_y_continuous(limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1),
                         labels = c("0", "0.2", "0.4", "0.5", "0.6", "0.8", "1")) +
      scale_color_manual(name = "", values = c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # magenta
      scale_fill_manual(name = "", values = c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # magenta
      geom_hline(aes(yintercept = 0.5),
                 colour="grey", linetype="dashed", size = 1) +
      geom_vline(aes(xintercept = 0.5),
                 colour="grey", linetype="dashed", size = 1) +
      labs(x = "Bribe Intensity", y = "P(Acceptance)\n") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"acceptprop_sigmoid.jpeg"), width = 8, height = 6, dpi = 600)

#### figure 6: bar plot of accept rate as a function of bribeint ####

choice_ind_summary3 <- ddply(df_valid,c("SbjID","f1","f2","Decision","Bribeint"), summarise,
                             N = length(Sender), .drop = FALSE)
choice_ind_summary4 <- ddply(df_valid,c("SbjID","f1","f2","Bribeint"), summarise,
                             N_valid = length(Sender), .drop = FALSE) # to make sure the valid sample size of each bribeprop
choice_ind_summary3 <- merge(choice_ind_summary3, choice_ind_summary4, by = c("SbjID","f1","f2","Bribeint"))
choice_ind_summary3$prop <- choice_ind_summary3$N/choice_ind_summary3$N_valid
choice_ind_summary3$Condition <- paste(choice_ind_summary3$f1, choice_ind_summary3$f2, sep = "_")
choice_ind_summary3$Condition <- factor(choice_ind_summary3$Condition, 
                                        levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))

choice_ind_summary3$Count <- ave(choice_ind_summary3$SbjID, 
                                 choice_ind_summary3[,c("Condition","Decision","Bribeint","prop")], 
                                 FUN = length)
choice_ind_summary3$Count <- as.numeric(levels(choice_ind_summary3$Count))[choice_ind_summary3$Count] # from factor to numeric


choice_summary3 <- ddply(choice_ind_summary3, c("f1","f2","Decision","Bribeint"), summarise,
                        N = length(SbjID), mean = mean(prop), 
                        sd = sd(prop), se = sd/sqrt(N), .drop = FALSE)
choice_summary3$Condition <- paste(choice_summary3$f1, choice_summary3$f2, sep = "_")
choice_summary3$Condition <- factor(choice_summary3$Condition, 
                                   levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))

g1 <- ggplot(data = choice_summary3[choice_summary3$Decision == "accept", ], 
             aes(x = f1, y = mean, color = Condition, fill = Condition)) + 
      facet_wrap(~ Bribeint, ncol = 2) +
      geom_bar(position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                    width = .1, size = 1, position = position_dodge(.9)) +
      geom_point(data = choice_ind_summary3[choice_ind_summary3$Decision == "accept", ], 
                 aes(x = f1, y = prop, group = f2, size = Count),
                 position = position_dodge(0.9), color = "black", alpha = 0.2) +
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # light green/magenta
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # dark green/magenta
      scale_y_continuous(labels = percent, limits = c(0,1), breaks = seq(0,1,0.2))+
      labs(x = "", y = "Acceptance\n", size = "") +
      theme_bw() + myconfig
print(g1)
ggsave(paste0(plot_savedir,"acceptprop_bribeint_bar.jpeg"), width = 10, height = 6, dpi = 600)

######################################################################
###      Step 3: Making descriptive plots of decision time (ms)    ###
######################################################################

#### summarize: within/across pps ####
DT_ind_summary <- ddply(df_valid,c("SbjID","f1","f2","Decision"),summarise,
                        N = length(DT), avgDT = mean(DT), sd = sd(DT),
                        se = sd/sqrt(N), .drop = FALSE)
DT_ind_summary$Condition <- paste(DT_ind_summary$f1, DT_ind_summary$f2, sep = "_")
DT_ind_summary$Condition <- factor(DT_ind_summary$Condition, 
                                   levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))

DT_summary <- ddply(DT_ind_summary[is.na(DT_ind_summary$avgDT) == F, ], 
                    c("f1","f2","Decision"), summarise,
                    N = length(avgDT), mean =mean(avgDT), sd = sd(avgDT),
                    se = sd/sqrt(N),.drop = F)
DT_summary$Condition <- paste(DT_summary$f1, DT_summary$f2, sep = "_")
DT_summary$Condition <- factor(DT_summary$Condition, 
                               levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))
#### figure 1: histogram of DT ####
g1 <- ggplot(data = df_valid, 
             aes(x = DT, group = Condition, color = Condition, fill = Condition)) +
      facet_wrap(~ Decision, ncol = 2) +
      geom_histogram(# aes(y = ..density..), 
                     binwidth = 500, alpha = 0.5, 
                     position = "identity") + # position = "dodge"
      # geom_density(alpha = 0.2) + 
      geom_vline(data = DT_summary, aes(xintercept = mean,  color = Condition),
                 linetype = "dashed", size = 1) +
      scale_x_continuous(limits = c(0, 8000),breaks = seq(0, 8000, 2000)) +
      scale_colour_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) +
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) +
      xlab("\nDecision Time (ms)") +
      ylab("Frequency\n") +
      theme_bw() + myconfig + 
      theme(strip.text.x = element_text(margin = margin(b = 10, t = 0)))
print(g1)
ggsave(paste0(plot_savedir,"DT_hist.jpeg"), width = 10, height = 6, dpi = 300)

#### figure 2: bar plot of DT ####
g1 <- ggplot(data = DT_summary, 
             aes(x = f1, y = mean, color = Condition, fill = Condition)) + 
      facet_wrap(~ Decision, ncol = 2) +
      geom_bar(position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                    width = .1, size = 1, position = position_dodge(.9)) +
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # light green/magenta
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # dark green/magenta
      scale_y_continuous(limits = c(1000,3000), breaks = seq(1000,3000,500), oob = rescale_none)+
      labs(x = "", y = "Decision Time (ms)\n", size = "") +
      theme_bw() + myconfig +
      theme(strip.text.x = element_text(margin = margin(b = 10, t = 0)))
print(g1)
ggsave(paste0(plot_savedir,"DT_bar.jpeg"), width = 11, height = 7, dpi = 600)


#### figure 3: violin plot of DT ####
data_summary <- function(x) {
  m <- mean(x)
  N <- length(x)
  ymin <- m-sd(x)/sqrt(N)
  ymax <- m+sd(x)/sqrt(N)
  return(c(y = m,ymin = ymin,ymax = ymax))
}

g1 <- ggplot(data = DT_ind_summary[is.na(DT_ind_summary$avgDT) == F, ], 
             aes(x = f1, y = avgDT, color = Condition, fill = Condition)) + 
      facet_wrap(~ Decision, ncol = 2) +
      geom_violin(trim = TRUE, size = 1) + 
      stat_summary(fun.data = data_summary, position = position_dodge(1),
                   geom = "pointrange", color = "black", size = 0.6) +
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # light green/magenta
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # dark green/magenta
      scale_y_continuous(limits = c(1000,7000), breaks = seq(1000,7000,1000), oob = rescale_none)+
      labs(x = "", y = "Decision Time (ms)\n", size = "") +
      theme_bw() + myconfig +
      theme(strip.text.x = element_text(margin = margin(b = 10, t = 0)))
print(g1)
ggsave(paste0(plot_savedir,"DT_violin.jpeg"), width = 11, height = 7, dpi = 300)


#### figure 4: point plot of DT as a function of bribeprop ####

DT_ind_summary1 <- ddply(df_valid,c("SbjID","Decision","f1","f2","Bribeprop"), summarise,
                         N = length(DT), avgDT = mean(DT), sd = sd(DT),
                         se = sd/sqrt(N), .drop = FALSE)
DT_ind_summary1$Condition <- paste(DT_ind_summary1$f1, DT_ind_summary1$f2, sep = "_")
DT_ind_summary1$Condition <- factor(DT_ind_summary1$Condition, 
                                   levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))

DT_summary1 <- ddply(DT_ind_summary1[is.na(DT_ind_summary1$avgDT) == F, ], 
                     c("Decision","f1","f2","Bribeprop"), summarise,
                     N = length(avgDT), mean =mean(avgDT), sd = sd(avgDT),
                     se = sd/sqrt(N),.drop = F)
DT_summary1$Condition <- paste(DT_summary1$f1, DT_summary1$f2, sep = "_")
DT_summary1$Condition <- factor(DT_summary1$Condition, 
                               levels = c("solo_honesty","solo_lie","dyad_honesty","dyad_lie"))

g1 <- ggplot(data = DT_summary1, 
             aes(x = Bribeprop, y = mean, group = Condition, color = Condition, fill = Condition)) + 
      facet_wrap(~ Decision, ncol = 2) +
      # geom_line(size = 2) +
      geom_point(size = 2, position = position_dodge(.9)) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                    width = .1, size = 1, position = position_dodge(.9)) +
      # geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.3) + 
      # geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), 
      #           color = "grey", alpha = 0.1) + 
      annotate("rect", xmin = c(1.5, 3.5, 5.5, 7.5), xmax = c(2.5, 4.5, 6.5, 8.5), 
                ymin = -Inf, ymax = Inf, color = "grey", alpha = .2) + # add grey shaded area
      scale_fill_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # light green/magenta
      scale_color_manual(name = "", values=c("#66FF99", "#FF3399", "#00FF00", "#FF00FF")) + # dark green/magenta
      scale_x_discrete(labels = c("10%","20%","30%","40%","50%","60%","70%","80%","90%")) +
      scale_y_continuous(limits = c(500,4500), breaks = seq(500, 4500, 500))+
      labs(x = "\nBribe Proportion", y = "Decition Time (ms)\n", size = "") +
      theme_bw() + myconfig +
      theme(strip.text.x = element_text(margin = margin(b = 10, t = 0)),
            axis.text.x = element_text(angle = 45, hjust = 1))
print(g1)
ggsave(paste0(plot_savedir,"DT_bribeprop_point.jpeg"), width = 15, height = 7, dpi = 600)


#####################################################################
###             Step 4: Statistical Analyses on Choices           ###
#####################################################################
library(tidyverse)
#%%%%%%%%%%%%%%%          Choice (Acceptance)        %%%%%%%%%%%%%%%%

## mean-centering the payoff covariates
# # If you run the plot first, you need to transform the factor into numberic
# df_valid$Bribeprop <- as.numeric(levels(df_valid$Bribeprop))[df_valid$Bribeprop] # from factor to numeric
# df_valid$Sender_Payoff <- as.numeric(levels(df_valid$Sender_Payoff))[df_valid$Sender_Payoff]

# Otherwise, these vars are already numeric. 
df_valid$c.Bribeprop <- as.numeric(scale(df_valid$Bribeprop, center = T, scale = F))
df_valid$c.Sender_Payoff <- as.numeric(scale(df_valid$Sender_Payoff, center = T, scale = F))

# df_valid$c.Sender_Kept <- scale(df_valid$Sender_Kept, center = T, scale = F)
# df_valid$c.Sender_Bribe <- scale(df_valid$Sender_Bribe, center = T, scale = F)
# 
# library(corrplot)
# M <- cor(df_valid[,c(73:74)])
# corrplot(M, method = "circle")

# for correlation plots and statistics of payoff-relevant variables,
# please check previous version

# Mixed-effect logistic regression on acceptance #
accept_fit1 <- glmer(IsAccept ~ 1 + f1*f2 + c.Bribeprop + c.Sender_Payoff +
                               (1|SbjID) + (1|RunNo), 
                     data = df_valid, family = "binomial",
                     control = glmerControl(optimizer = c("bobyqa")))

# This optimizer does not work
# control = glmerControl(optimizer = c("optimx"), 
#                        optCtrl = list(method = "nlminb"))

# In response to reviewer's elife:
# random slope structure: 
# s1: (1 + f1*f2|SbjID): failed to converge
# s2: (1 + f1 + f2|SbjID): converge
accept_fit1_full <- glmer(IsAccept ~ 1 + f1*f2 + c.Bribeprop + c.Sender_Payoff +
                                    (1 + f1 + f2|SbjID) + (1|RunNo), 
                     data = df_valid, family = "binomial",
                     control = glmerControl(optimizer = c("bobyqa")))

summary(accept_fit1_full)
confint(accept_fit1_full)
Anova(accept_fit1_full)

fixcoef <- fixef(accept_fit1_full)
OR <- exp(fixcoef)
print(OR)

# interaction effect: significant (N = 39: p = 0.017)
# anova(accept_fit1, accept_fit2)
set_theme(base = theme_bw() + myconfig)
g1 <- sjp.glmer(accept_fit1, type = "fe", vars = c("f1dyad","f2lie","f1dyad:f2lie",
                                                   "c.Bribeprop","c.Sender_Payoff")) 
print(g1)
ggsave(paste0(plot_savedir, "accept_sjplot_fe.jpeg"), width = 10, height = 6, dpi = 300)

# To unpack the two-way interaction
df_valid_solo <- df_valid %>% filter(f1 == "solo") %>%
                 dplyr::select(-c.Bribeprop, -c.Sender_Payoff) %>%
                 mutate(c.Bribeprop = as.numeric(scale(Bribeprop, center = T, scale = F)),
                        c.Sender_Payoff = as.numeric(scale(Sender_Payoff, center = T, scale = F)))

accept_solo_fit <- glmer(IsAccept ~ 1 + f2 + c.Bribeprop + c.Sender_Payoff + 
                                   (1 + f2|SbjID) + (1|RunNo),
                        # data = df_valid[df_valid$f1 == "solo",], 
                        data = df_valid_solo, 
                        family = "binomial",
                        control = glmerControl(optimizer = c("bobyqa")))
summary(accept_solo_fit)
confint(accept_solo_fit) # it takes quite long
Anova(accept_solo_fit)
fixcoef <- fixef(accept_solo_fit)
OR <- exp(fixcoef)
print(OR)


df_valid_dyad <- df_valid %>% filter(f1 == "dyad") %>%
  dplyr::select(-c.Bribeprop, -c.Sender_Payoff) %>%
  mutate(c.Bribeprop = as.numeric(scale(Bribeprop, center = T, scale = F)),
         c.Sender_Payoff = as.numeric(scale(Sender_Payoff, center = T, scale = F)))

accept_dyad_fit <- glmer(IsAccept ~ 1 + f2 + c.Bribeprop + c.Sender_Payoff + (1 + f2|SbjID) + (1|RunNo),
                        # data = df_valid[df_valid$f1 == "dyad",], 
                        data = df_valid_dyad, 
                        family = "binomial",
                        control = glmerControl(optimizer = c("bobyqa")))
summary(accept_dyad_fit)
confint(accept_dyad_fit) # it takes quite long
Anova(accept_dyad_fit)
fixcoef <- fixef(accept_dyad_fit)
OR <- exp(fixcoef)
print(OR)

#####################################################################
###                Step 5: Statistical Analyses on DT             ###
#####################################################################

#%%%%%%%%%%%%%%%           Decision Time (DT)         %%%%%%%%%%%%%%%%
# Testing the normality #
require(nortest)
ad.test(df_valid$DT) 
df_valid$logDT <- log(df_valid$DT) # Log-transformation

# Running regression analyses #
logDT_fit1 <- lmer(logDT ~ 1 + Decision*f1*f2 + c.Bribeprop + c.Sender_Payoff +
                          (1|SbjID) + (1|RunNo), data = df_valid)

# In response to reviewer's elife:
# random slope structure: 
# s1: (1 + Decision*f1*f2|SbjID): converge (bu singular)
# s2: (1 + Decision + f1*f2|SbjID): failed to converge
# s3: (1 + f1*f2|SbjID): failed to converge
logDT_fit1_full <- lmer(logDT ~ 1 + Decision*f1*f2 + c.Bribeprop + c.Sender_Payoff +
                               (1 + Decision*f1*f2|SbjID) + (1|RunNo), data = df_valid)

anova(logDT_fit1_full)
summary(logDT_fit1_full)
confint(logDT_fit1_full)
AIC(logDT_fit1_full)
BIC(logDT_fit1_full)

# lmerTest::anova(logDT_fit1)

g1 <- sjp.lmer(logDT_fit1, type = "fe",
               vars = c("f1dyad", "f2lie", "Decisionreject", "c.Bribeprop", "c.Sender_Payoff",
                        "Decisionreject:f1dyad", "Decisionreject:f2lie", "f1dyad:f2lie"))
print(g1)
ggsave("./outputs/figures/logDT_sjplot_fe.jpeg", width = 10, height = 6, dpi = 300)


# To unpack the two-way (Decision:f2) interaction, 
# we divided the dataset based on the decision
#%%%%%%%%%%%%%%%%%%    Accept   %%%%%%%%%%%%%%%%%%%%%%%%#
# random slope structure: 
# s1: (1 + f1*f2|SbjID): failed to converge
# s2: (1 + f1 + f2|SbjID): failed to converge
# s3: (1 + f1|SbjID): converge (AIC = 5409.9; BIC = 5479.5)
# s4: (1 + f2|SbjID): converge (AIC = 5302.4; BIC = 5372.0) # we choose this one because it has a lower AIC/BIC

df_valid_accept <- df_valid %>% filter(Decision == "accept") %>%
  dplyr::select(-c.Bribeprop, -c.Sender_Payoff) %>%
  mutate(c.Bribeprop = as.numeric(scale(Bribeprop, center = T, scale = F)),
         c.Sender_Payoff = as.numeric(scale(Sender_Payoff, center = T, scale = F)))

logDT_accept_fit1 <- lmer(logDT ~ 1 + f1*f2 + c.Bribeprop + c.Sender_Payoff +
                                 (1 + f2|SbjID) + (1|RunNo), 
                          # data = df_valid[df_valid$Decision == "accept",]
                          data = df_valid_accept)
anova(logDT_accept_fit1)
summary(logDT_accept_fit1)
confint(logDT_accept_fit1)
lme.dscore(logDT_accept_fit1, 
           df_valid[df_valid$Decision == "accept",], type = "lme4")
AIC(logDT_accept_fit1)
BIC(logDT_accept_fit1)
# main effect of f1 and f2 (no interaction)

g1 <- sjp.lmer(logDT_accept_fit1, type = "fe",
               vars = c("f1dyad", "f2lie", "c.Bribeprop", "c.Sender_Payoff"))
print(g1)
ggsave(paste0(plot_savedir, "logDT_accept_sjplot_fe.jpeg"), width = 10, height = 6, dpi = 300)

# # To be consistent with the figures
# logDT_accept_solo_fit1 <- lmer(logDT ~ 1 + f2 + c.Bribeprop + c.Sender_Payoff +
#                             (1|SbjID) + (1|RunNo), 
#                           data = df_valid[df_valid$Decision == "accept" & df_valid$f1 == "solo",])
# summary(logDT_accept_solo_fit1)
# 
# logDT_accept_dyad_fit1 <- lmer(logDT ~ 1 + f2 + c.Bribeprop + c.Sender_Payoff +
#                                  (1|SbjID) + (1|RunNo), 
#                                data = df_valid[df_valid$Decision == "accept" & df_valid$f1 == "dyad",])
# summary(logDT_accept_dyad_fit1)


#%%%%%%%%%%%%%%%%%%    Reject   %%%%%%%%%%%%%%%%%%%%%%%%#
# random slope structure: 
# s1: (1 + f1*f2|SbjID): failed to converge
# s2: (1 + f1 + f2|SbjID): converge

df_valid_reject <- df_valid %>% filter(Decision == "reject") %>%
  dplyr::select(-c.Bribeprop, -c.Sender_Payoff) %>%
  mutate(c.Bribeprop = as.numeric(scale(Bribeprop, center = T, scale = F)),
         c.Sender_Payoff = as.numeric(scale(Sender_Payoff, center = T, scale = F)))

logDT_reject_fit1 <- lmer(logDT ~ 1 + f1*f2 + c.Bribeprop + c.Sender_Payoff +
                                 (1 + f1 + f2|SbjID) + (1|RunNo), 
                          # data = df_valid[df_valid$Decision == "reject",]
                          data = df_valid_reject)
anova(logDT_reject_fit1)
summary(logDT_reject_fit1)
confint(logDT_reject_fit1)
AIC(logDT_reject_fit1)
BIC(logDT_reject_fit1)
lme.dscore(logDT_reject_fit1, 
           df_valid[df_valid$Decision == "reject",], type = "lme4")

g1 <- sjp.lmer(logDT_reject_fit1, type = "fe",
               vars = c("f1dyad", "f2lie", "c.Bribeprop", "c.Sender_Payoff"))
print(g1)
ggsave(paste0(plot_savedir, "logDT_reject_sjplot_fe.jpeg"), width = 10, height = 6, dpi = 300)
