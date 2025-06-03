
# ========================================================================
# Script: analyze_sphere_connectivity.R
# Description:
#   This script analyzes functional connectivity (FC) between four 
#   predefined spherical ROIs (aOFC/pOFC seeds and LPFC targets) based 
#   on preprocessed fMRI data. It includes:
#
#   1. Loading and cleaning FC matrices from .mat files, renaming spheres.
#   2. Extracting lower triangle FC values to compute within-subject 
#      averages for same-network and between-network comparisons.
#   3. Visualizing subject-level FC distributions using violin plots 
#      and paired boxplots, aligned across different connection types.
#   4. Conducting paired Wilcoxon tests:
#      - Between same-network and cross-LPFC connectivity
#      - Across different seeds targeting the same LPFC region
#   5. Creating a combined figure (p1 + p2) for publication using `ggarrange()`.
#
# Output:
#   - Summary statistics per connection type
#   - PDF figure saved as "Fig_conn_seed_stim.pdf"
#
# Requirements:
#   - FC matrices from MATLAB in "conn_matrix_4spheres.mat"
#   - Sphere naming defined by original_names → new_names map
#   - `ggplot2`, `ggpubr`, `ggh4x`, and tidyverse packages
# ========================================================================


rm(list = ls())
source('../Scripts_R_beh/Setup.R')
SubInfo = read.xlsx('../ProcessedData/SubConds.xlsx') %>%
  filter(Include == 1)
Subs = SubInfo$SubID

# Set your path to where all subject/session folders are
base_dir <- "/Volumes/X9Pro/NODEAP/FuncConn_spheres"
# Get all mat file paths
mat_files <- list.files(base_dir, pattern = "conn_matrix_4spheres.mat", 
                        recursive = TRUE, full.names = TRUE)

original_names <- c(
  "aOFC_seed_right.nii",
  "aOFC_target_right.nii",
  "pOFC_seed_right.nii",
  "pOFC_target_right.nii"
)
new_names <- c('Seed-aOFC', 'aOFC-conn-LPFC',
               'Seed-pOFC', 'pOFC-conn-LPFC')
name_map <- setNames(new_names, original_names)

# Function to load and parse FC matrix
parse_fc_file <- function(file_path, keep_lower = TRUE) {
  mat_data <- readMat(file_path)
  corr_mat <- mat_data$correlation.matrix
  spheres <- unlist(mat_data$sphere.order)
  clean_spheres <- name_map[spheres]
  
  # Create data frame and assign names
  df <- as.data.frame(as.table(corr_mat))
  colnames(df) <- c("Sphere1", "Sphere2", "FC")
  df$Sphere1 <- clean_spheres[df$Sphere1]
  df$Sphere2 <- clean_spheres[df$Sphere2]
  
  if (keep_lower) { # Keep only lower triangle
    idx1 <- match(df$Sphere1, clean_spheres)
    idx2 <- match(df$Sphere2, clean_spheres)
    df <- df[idx1 > idx2, ]
  }
  
  # Extract subject/session from file path
  path_parts <- strsplit(file_path, .Platform$file.sep)[[1]]
  df$SubID <- path_parts[6]
  df$Session <- path_parts[7]
  
  return(df)
}

# Load all FC data (lower triangle only)
fc_df <- bind_rows(lapply(mat_files, parse_fc_file)) %>%
  rowwise() %>%
  mutate(Pair = paste(sort(c(Sphere1, Sphere2)), collapse = "_")) %>%
  ungroup()

fc_subj_avg_same_network <- fc_df %>%
  subset(Pair %in% c('aOFC-conn-LPFC_Seed-aOFC','pOFC-conn-LPFC_Seed-pOFC')) %>%
  group_by(SubID) %>%
  reframe(FC = mean(FC, na.rm = TRUE)) %>%
  mutate(type = 'same network') 
fc_subj_avg_bt_LPFC <- fc_df %>%
  subset(Pair %in% 'aOFC-conn-LPFC_pOFC-conn-LPFC') %>%
  group_by(SubID) %>%
  reframe(FC = mean(FC, na.rm = TRUE)) %>%
  mutate(type = 'between LPFCs')
fc_summary_stats <- rbind(fc_subj_avg_same_network,fc_subj_avg_bt_LPFC)

# Load all data again (this time for seed-to-LPFC, no triangle filtering)
fc_df_seed_stim <- bind_rows(
  lapply(mat_files, function(f) parse_fc_file(f, keep_lower = FALSE))) %>%
  filter(grepl("Seed", Sphere1) & grepl("LPFC", Sphere2)) %>%
  group_by(SubID, Sphere1, Sphere2) %>%
  reframe(FC = mean(FC, na.rm = TRUE))

###################################

ylow = 0.82; yhigh = 0.85; ytext = 0.86

p1 = fc_summary_stats %>%
  mutate(type=factor(type,
                     levels = c('same network','between LPFCs'),
                     labels = c("Within\nnetwork", "Between\nLPFCs"))) %>%
  ggplot(aes(x = type, y = FC)) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_violinhalf(flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) + 
  scale_y_continuous(
    limits = c(-0.25, 0.88),           # Set y-axis limits
    breaks = seq(-0.25, 1, by = 0.25)  # Set specific tick marks
  ) +
  annotate("text", x=1.5, y=ytext, label='*') +
  annotate("segment",x = 1, xend = 1, y = ylow, yend = yhigh) +
  annotate("segment",x = 2, xend = 2, y = ylow, yend = yhigh) +
  annotate("segment",x = 1, xend = 2, y = yhigh, yend = yhigh) +
  labs(y=NULL,x=NULL,title=NULL) +
  common

val_same_net <- fc_summary_stats %>%
  filter(type == "same network") %>% pull(FC) 
val_bt_LPFC <- fc_summary_stats %>%
  filter(type == "between LPFCs") %>% pull(FC) 
wilcox.test(val_same_net, val_bt_LPFC, paired = TRUE, alternative = "two.sided")


strip = strip_themed(background_y = elem_list_rect(fill = use.col.ap.lpfc),
                     text_y = elem_list_text(color = 'black',face = "italic",size = 16))

p2=fc_df_seed_stim %>%
  mutate(Sphere1=factor(Sphere1,
                     levels = c('Seed-aOFC','Seed-pOFC'),
                     labels = c("Seed\naOFC", "Seed\npOFC"))) %>%
  ggplot(aes(x = Sphere1, y = FC)) +
  geom_line(aes(group = SubID), color = "gray", alpha = 0.8,
            position = position_dodge(0.2)) + 
  geom_violinhalf(flip = c(1, 3),
                  position = position_nudge(x = 0, y = 0)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(aes(group = SubID), alpha = 0.5, 
              position = position_dodge(0.2)) +
  facet_wrap2(~Sphere2,axes = 'all',strip = strip, strip.position = 'left') +
  annotate("text", x=1.5, y=ytext, label='***') +
  annotate("segment",x = 1, xend = 1, y = ylow, yend = yhigh) +
  annotate("segment",x = 2, xend = 2, y = ylow, yend = yhigh) +
  annotate("segment",x = 1, xend = 2, y = yhigh, yend = yhigh) +
  labs(y='Functional connectivity',x=NULL,title=NULL) +
  common +
  scale_y_continuous(
    limits = c(-0.25, 0.88),           # Set y-axis limits
    breaks = seq(-0.25, 1, by = 0.25)  # Set specific tick marks
  ) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16)) 

pdf(file.path(FigPaperDir,'Fig_conn_seed_stim.pdf'),9.5,4)
ggarrange(p2,p1,nrow = 1,widths = c(2.5,1),align = "h")
dev.off()


# Cross-seed test: Seed-aOFC vs. Seed-pOFC both to aOFC-conn-LPFC
aSeed_aStim <- fc_df_seed_stim %>%
  filter(Sphere1 == "Seed-aOFC", Sphere2 == "aOFC-conn-LPFC") %>% pull(FC)
pSeed_aStim <- fc_df_seed_stim %>%
  filter(Sphere1 == "Seed-pOFC", Sphere2 == "aOFC-conn-LPFC") %>% pull(FC) 
wilcox.test(aSeed_aStim, pSeed_aStim, paired = TRUE, alternative = "two.sided")

# Cross-seed test: Seed-aOFC vs. Seed-pOFC both to pOFC-conn-LPFC
aSeed_pStim <- fc_df_seed_stim %>%
  filter(Sphere1 == "Seed-aOFC", Sphere2 == "pOFC-conn-LPFC") %>% pull(FC) 
pSeed_pStim2 <- fc_df_seed_stim %>%
  filter(Sphere1 == "Seed-pOFC", Sphere2 == "pOFC-conn-LPFC") %>% pull(FC) 
wilcox.test(pSeed_pStim2, aSeed_pStim, paired = TRUE, alternative = "two.sided")

#############################################################

# use their group assignment
summary_convals_expt = dfconvals %>%
  subset(Stim=='aOFC-conn-LPFC'& StimLoc=='Anterior'|
           Stim=='pOFC-conn-LPFC'& StimLoc=='Posterior') %>%
  group_by(Seed,Stim,SubID,StimLoc) %>% 
  reframe(Conn=mean(Conn,na.rm=T)) %>% 
  arrange(SubID)

# calculate indexes
df1 = summary_convals_expt %>% subset(Seed=='Seed-aOFC')
ACI = df1$Conn
df2 = summary_convals_expt %>% subset(Seed=='Seed-pOFC')
PCI = df2$Conn

df_idx = df1 %>% 
  select(SubID,StimLoc) %>%
  mutate(ACI=ACI,PCI=PCI)

write.csv(df_idx, "../ProcessedData/ConnIdx.csv", row.names = FALSE)


  
  

