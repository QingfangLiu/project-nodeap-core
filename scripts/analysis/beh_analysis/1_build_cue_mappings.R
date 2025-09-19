
# ========================================================================
# Script: build_cue_mappings.R
#
# Purpose:
#   Construct cue–cuePair mapping tables to link cue pair indices with the
#   odor-predicting cue IDs. These mappings are needed to interpret model
#   weights (w values) that are indexed by cue pair.
#
# What it does:
#   1. Loads processed behavioral data for conditioning and post-meal choice.
#   2. Filters to trials with odor presented on the left (OdorLR == 1) and
#      Run == 1 (for conditioning).
#   3. Selects key columns (SubID, Sess, CuePair, Set, ChoiceType, CueLeft,
#      CueRight) and arranges rows for consistent indexing.
#   4. Saves two mapping tables:
#        - conditioning_cue_mapping.RData
#        - choice_cue_mapping_post.RData
#
# Output:
#   RData files in the ProcessedData folder, ready to be used to map
#   model-estimated w values back to cue IDs.
#
# Usage:
#   Edit the `project_folder` path below to match your local project location,
#   then run this script in R. The output files will be written to
#   <project_folder>/ProcessedData/.
# ========================================================================

rm(list = ls())

project_folder <- "/Users/liuq13/project-nodeap-core" 
# Minimal path helper
p <- function(...) file.path(project_folder, ...)

source(file.path(project_folder, "scripts", "utils", "Setup.R"))

# ---- Paths ----
processed_dir <- p("ProcessedData")        # keep a single, consistent folder
if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)

# ---- Conditioning: build cue mapping ----
# expects object `conditioning_dat` inside the .RData
load(p("ProcessedData", "Conditioning.RData"))

conditioning_cue_mapping = conditioning_dat %>%
  subset(OdorLR==1 & Run==1) %>%
  select(SubID,Sess,Cond,CuePair,Set,ChoiceType,CueLeft,CueRight) %>%
  arrange(SubID,Sess,CuePair)
# Note: Because we subset to OdorLR == 1 (odor on LEFT),
#       CueLeft is the odor-predicting cue ID.

save(conditioning_cue_mapping,
     file = file.path(processed_dir, "conditioning_cue_mapping.RData"))

# ---- Choice (Post): build cue mapping ----
# expects object `choice_dat` inside the .RData
load(file.path(processed_dir, "choice_dat.RData"))

choice_cue_mapping_post = choice_dat %>%
  subset(PrePost=='Post' & OdorLR==1) %>%
  select(SubID,Sess,CuePair,Set,ChoiceType,CueLeft,CueRight) %>%
  arrange(SubID,Sess,CuePair)

save(choice_cue_mapping_post,
     file = file.path(processed_dir, "choice_cue_mapping_post.RData"))

