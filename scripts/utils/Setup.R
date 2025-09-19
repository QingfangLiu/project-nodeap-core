
library(MASS) # to use rlm robust linear model, load before dplyr
library(openxlsx)
library(ggpubr) 
library(see) # to use geom_violinhalf
library(R.matlab) # to use readMat
library(ggExtra) # Add Marginal Histograms to 'ggplot2'
library(cowplot) # combine multiple panels
library(lattice) # data plotting
library(ggpattern) # to use geom_boxplot_pattern
library(ggeffects)
library(lme4)
library(plyr)
library(tidyverse)
library(ggh4x) # tweaking ggplots (e.g. facets)
library(grid)  # manually adding panel/label names with shade

project_folder <- "/Users/liuq13/project-nodeap-core"
FigPaperDir = file.path(project_folder,'figs_paper')

pd <- position_dodge(0.6) # set dodge width
#pd <- position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6, jitter.height = 0)

pd_jitter <- position_jitterdodge(jitter.width = 0.000, dodge.width = 0.2, jitter.height = 0)
pd_dodge <- position_dodge(0.2)

# some common ggplot settings
common = 
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(colour = "black"), # axis color
        axis.text.y = element_text(colour = "black"),
        text = element_text(size=16),   # text size & font
        plot.title = element_text(hjust = 0.5), # center title
        strip.background = element_blank(), # not show panel bg
        legend.title = element_blank()) # not show legend title
        
# color scales to use
use.col.conds = c("sham-sham" = "gray","cTBS-sham" = "#3bc9d9","sham-cTBS" = "#3b6fd9")
use.col.sess = c("1" = "#A1E7B8","2" = "#53BD7D","3" = "#0E600F")
use.col.conds.day1 = c("sham" = "gray","cTBS" = "#3bc9d9")
use.col.conds.day2 = c("sham" = "gray","cTBS" = "#3b6fd9")

use.col.ap.ofc = c("aOFC" = "#ED7014","pOFC" = "#8B008B")
use.col.ap.lpfc = c("aOFC" = "#FCAE1E","pOFC" = "#FF0BE8")




