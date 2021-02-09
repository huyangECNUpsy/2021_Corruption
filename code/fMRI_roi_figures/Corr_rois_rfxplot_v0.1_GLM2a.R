#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Corr_rois_rfxplot_v0.1_GLM16.R                                           #
# Author: Yang Hu                                                          #
# Established Date: May. 17, 2018                                          #
# Last Change: June. 6, 2020                                               #
# Project: Corruption (Corr)                                               #
# PI: Yang Hu                                                              #
# Collaborators: Chen Qu, Jean-Claude Dreher                               #
# This is a custom R script to draw rfx plots of rois from GLM analyses to #
# show the parametric effect (fMRI study).                                 #
# GLM16 distinguished SVdiff in bribe from control condition               #
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
sheet_names <- c("GLM16_vmPFC_SVdiff_H","GLM16_vmPFC_SVdiff_L") # 1 condition (SVdiff)

# Built-in function
mylm <- function (data_input) {
  
    data <- data_input
    fit <- lm(psc ~ as.numeric(as.factor(Bin)), data = data)
    
    output <- data.frame(coef(fit)[[1]],coef(fit)[[2]])
    return(output)
} 

myrfxplot <- function(GLM_label) { 

    df_tmp <- df_long[df_long$GLM == GLM_label,]
    output <- ddply(df_tmp,.var = "Cond",.fun = mylm)
    colnames(output)[c(2,3)] <- c("b0","b1")
    output$bin1_hat <- output$b0 + output$b1*1
    output$bin2_hat <- output$b0 + output$b1*2
    output$bin3_hat <- output$b0 + output$b1*3
    output$x_start <- c(0.5,1.5)
    output$x_end <- c(1.5,2.5)
    output$Bin <- c("bin1") # Just need this variable
    
    avgbetaROI <- ddply(df_tmp[is.na(df_tmp$psc) == F, ],
                        c("ROI","Cond","Bin"), summarise, 
                        N = length(psc), mean=mean(psc), sd=sd(psc), se=sd/sqrt(N),.drop = F)
    avgbetaROI$Bin <- as.factor(avgbetaROI$Bin)
    legend_label <- "rSV"
    
    g1 <- ggplot(avgbetaROI, aes(x = Cond, y = mean, fill = Bin)) + 
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
  ggsave(paste0(plot_savedir, rfxplot_name), width = 10, height = 6, dpi = 600)

}

############################################################################
#                         Step 2: Making rfx figures                       #
############################################################################
#%%%%%%%%%%%%%%%%          Extract and Prepare Data        %%%%%%%%%%%%%%%%#
df_long <- data.frame()
for (i in 1:length(sheet_names)) {
  df_tmp = read_excel(paste0(rootdir,"data_roi/rfxreg/beta_ROI_rfxplot_long_GLM16.xlsx"),
                      sheet = i, col_names = F, col_types = NULL) # GLM1
  df_long <- rbind(df_long,df_tmp)
}
colnames(df_long)[1:6] <- c("GLM","ROI","SbjNo","Cond","Bin","psc")
sum(is.na(df_long$psc)) # making sure no NA 
# df_long$Cond_new <- df_long$Cond
df_long$Cond[df_long$ROI == "vmPFC_SVdiff_H"] <- "Control"
df_long$Cond[df_long$ROI == "vmPFC_SVdiff_L"] <- "Bribe"
df_long$Cond <- factor(df_long$Cond, levels = c("Control","Bribe"))
df_long$ROI <- "vmPFC"

myrfxplot("GLM16")
