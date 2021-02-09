#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Corr_rois_rfxplot_v0.1_GLM1b.R                                           #
# Author: Yang Hu                                                          #
# Established Date: May. 17, 2018                                          #
# Last Change: Jan. 8, 2020                                                #
# Project: Corruption (Corr)                                               #
# PI: Yang Hu                                                              #
# Collaborators: Chen Qu, Jean-Claude Dreher                               #
# This is a custom R script to draw rfx plots of rois from GLM analyses to #
# show the parametric effect (fMRI study).                                 #
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
library(lme4)
library(psych)
library(reshape2)
library(plyr)

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

rootdir <- '/gpfs/share/home/1906388343/ISC_local/Studies_CNRS_Lyon/2017_Corruption/fMRI/'
# plot_savedir <- "./outputs/figures/roiplots_39sbjs/"
plot_savedir <- paste0(rootdir,"outputs/figures/roiplots_39sbjs/")
sheet_names <- c("GLM13e_lvAI","GLM13e_rvAI","GLM13e_rvAI_solo", # 4 conditions (SH, SL, DH, DL)
                 "GLM14b_rTPJ","GLM9e_vmPFC") # 1 condition (either DL or pooled)

# Built-in function
mylm <- function (data_input) {
  
    data <- data_input
    fit <- lm(psc ~ as.numeric(as.factor(Bin)), data = data)
    
    output <- data.frame(coef(fit)[[1]],coef(fit)[[2]])
    return(output)
} 



myrfxplot <- function(GLM_label) { 
  
  # GLM_label: "GLM13e", # 4 conditions (SH, SL, DH, DL)
  # "GLM14b","GLM9e" # 1 condition (either DL or pooled)
  
  df_tmp <- df_long[df_long$GLM == GLM_label,]
  
  
  if (GLM_label == "GLM14b" | GLM_label == "GLM9e") {
    output <- ddply(df_tmp,.var = "ROI",.fun = mylm)
    colnames(output)[c(2,3)] <- c("b0","b1")
    output$bin1_hat <- output$b0 + output$b1*1
    output$bin2_hat <- output$b0 + output$b1*2
    output$bin3_hat <- output$b0 + output$b1*3
    
    output$Bin <- c("bin1")
    
    avgbetaROI <- ddply(df_tmp[is.na(df_tmp$psc) == F, ],
                        c("ROI","Bin"), summarise, 
                        N = length(psc), mean=mean(psc), sd=sd(psc), se=sd/sqrt(N),.drop = F)
    avgbetaROI$Bin <- as.factor(avgbetaROI$Bin)

    if (GLM_label == "GLM14b") {
      legend_label <- "Potential Loss"
      avgbetaROI <- avgbetaROI[avgbetaROI$ROI == "rTPJ",]
      output <- output[output$ROI == "rTPJ",]
      # output$x_start <- c(0.5,1.5)
      # output$x_end <- c(1.5,2.5)
      output$x_start <- c(0.5)
      output$x_end <- c(1.5)
    } else {
      legend_label <- "rSV"
      output$x_start <- c(0.5)
      output$x_end <- c(1.5)
    }
    
    g1 <- ggplot(avgbetaROI, aes(x = ROI, y = mean, fill = Bin))
    
  } else {
    output <- ddply(df_tmp,.var = "Cond",.fun = mylm)
    colnames(output)[c(2,3)] <- c("b0","b1")
    output$f1 <- ifelse(output$Cond == "cond1" | output$Cond == "cond2", "solo", "dyad")
    output$f1 <- factor(output$f1, levels = c("solo","dyad"))
    output$f2 <- ifelse(output$Cond == "cond1" | output$Cond == "cond3", "honesty", "lie") 
    output$bin1_hat <- output$b0 + output$b1*1
    output$bin2_hat <- output$b0 + output$b1*2
    output$bin3_hat <- output$b0 + output$b1*3
    output$x_start <- c(0.5,0.5,1.5,1.5)
    output$x_end <- c(1.5,1.5,2.5,2.5)
    output$Bin <- c("bin1") # Just need this variable
    # output_long <-melt(output[,c(1,6:8)], id.vars=c("Cond"))
    # output_long$f1 <- ifelse(output$Cond == "cond1" | output$Cond == "cond2", "solo", "dyad") 
    # output_long$f2 <- ifelse(output$Cond == "cond1" | output$Cond == "cond3", "honesty", "lie") 
    # output_long$Bin <- substr(as.character(output_long$variable),1,4)
    
    avgbetaROI <- ddply(df_tmp[is.na(df_tmp$psc) == F, ],
                        c("ROI","f1","f2","Bin"), summarise, 
                        N = length(psc), mean=mean(psc), sd=sd(psc), se=sd/sqrt(N),.drop = F)
    avgbetaROI$f1 <- factor(avgbetaROI$f1, levels = c("solo","dyad"))
    avgbetaROI$f2 <- factor(avgbetaROI$f2, levels = c("honesty","lie"))
    avgbetaROI$Bin <- as.factor(avgbetaROI$Bin)
    
    g1 <- ggplot(avgbetaROI, aes(x = f1, y = mean, fill = Bin)) + 
          facet_wrap( ~ f2, ncol = 2)
    
    if (GLM_label == "GLM13e") {
      legend_label <- "Potential Gain"
    } else {
      legend_label <- "rSV"
    }
  }
  
  g1 <- g1 + 
        geom_bar(position = position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = "black",
                      width = .1, size = 1, position = position_dodge(.9)) +
        # geom_point(data = output_long,aes(y = value), 
        #            position = position_dodge(1), color = "black") +
        geom_segment(data = output, aes(x = x_start, y = bin1_hat, 
                                        xend = x_end, yend = bin3_hat), 
                     color = "blue", size = 1.5) +
        xlab("") + ylab("% Signal Change\n") +
        # scale_x_continuous(labels = c("low","middle","high"),breaks = seq(1,4,1)) +
        # scale_y_continuous(limits = c(-0.1, 0.1), breaks = seq(-0.3,0.3,0.05)) + # GLM1
        # scale_fill_gradient(name = "",low = c("darkred"), high = c("red")) +
        # scale_fill_grey(name = legend_label,labels = c("low","middle","high"), 
        #                 start = 0.8, end = 0.2) +
        scale_fill_manual(name = legend_label,labels = c("low","middle","high"), 
                          values=c("red", "orange", "yellow")) +
        theme_bw() + myconfig 
   
        # theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(g1)
  
  rfxplot_name <- paste0("rfxplot_",GLM_label,".jpeg")
  ggsave(paste0(plot_savedir, rfxplot_name), width = 10, height = 7, dpi = 600)
  
  
}

############################################################################
#                         Step 2: Making rfx figures                       #
############################################################################
#%%%%%%%%%%%%%%%%          Extract and Prepare Data        %%%%%%%%%%%%%%%%#
df_long <- data.frame()
for (i in 1:length(sheet_names)) {
  df_tmp = read_excel(paste0(rootdir,"data_roi/rfxreg/beta_ROI_rfxplot_long.xlsx"),
                      sheet = i, col_names = F, col_types = NULL) # GLM1
  df_long <- rbind(df_long,df_tmp)
}
colnames(df_long)[1:6] <- c("GLM","ROI","SbjNo","Cond","Bin","psc")
sum(is.na(df_long$psc)) # making sure no NA 
df_long$Cond_new <- df_long$Cond
# df_long$Cond_new[df_long$GLM == "GLM14b" & 
#                  df_long$ROI == "dmPFC" & 
#                  df_long$Cond == "cond1"] <- "DL"
df_long$Cond_new[df_long$GLM == "GLM14b" & 
                   df_long$ROI == "rTPJ" & 
                   df_long$Cond == "cond1"] <- "DL"
df_long$Cond_new[df_long$GLM == "GLM9e" & 
                   df_long$ROI == "vmPFC"] <- "pooled"
attach(df_long)
df_long$f1[Cond_new == "cond1" | Cond_new == "cond2"] <- "solo"
df_long$f1[Cond_new == "cond3" | Cond_new == "cond4"] <- "dyad"
df_long$f2[Cond_new == "cond1" | Cond_new == "cond3"] <- "honesty"
df_long$f2[Cond_new == "cond2" | Cond_new == "cond4"] <- "lie"
detach(df_long)

df_long$Cond_new[df_long$Cond_new == "cond1"] <- "SH"
df_long$Cond_new[df_long$Cond_new == "cond2"] <- "SL"
df_long$Cond_new[df_long$Cond_new == "cond3"] <- "DH"
df_long$Cond_new[df_long$Cond_new == "cond4"] <- "DL"

myrfxplot("GLM14b")
myrfxplot("GLM9e")
