# Behavioral Data (Cleaned)

This folder contains behavioral data that have been cleaned and organized from the raw experimental files. Sensitive participant information (e.g., identifiable demographics) has been removed.

## Files

- **`Choices_all.xlsx`**  
  Trial-level choice data from the pre- and post-meal probe tasks.

- **`Conditioning_all.xlsx`**  
  Behavioral responses during the Day 1 conditioning task.

- **`Pleasant_all.xlsx`**  
  Pleasantness ratings for each odor, collected before and after the meal.

- **`SubConds.xlsx`**  
  Cleaned subject condition table with the following:
  - Subject ID, StimLoc (Anterior/Posterior corrected for NODEAP_73)
  - Stimulation order (mapped from numeric codes to labels, e.g., `CS-SC-SS`)
  - Odor set assignment, starting odor (Sweet/Savory)
  - Age (rounded), Sex
  - Within-session and between-session intervals (in days)

- **`Survey_uncomf_strong.xlsx`**  
  Long-format data frame with one row per subject × session containing:
  - Subject ID, StimLoc, session name, session number (1–6)
  - Derived TMS type (`cTBS` vs `sham`)
  - Ratings of TMS discomfort (`Uncomfortable`) and intensity (`Strong`)
  - Z-scored versions of the ratings (within-subject standardization)


## Notes

- Metadata were sourced from `NODEAP_DataCollectionSheet_QL.xlsx` and
  `Data Organization_QL.xlsx`, then cleaned and harmonized.
- StimLoc for **NODEAP_73** was manually corrected to `Posterior`
  based on TMS coordinates.
- Calories for Day-2 sessions (S1D2, S2D2, S3D2) were extracted and aligned
  with subjects; one bad session (Subject 10, S2D2) was set to `NA` and excluded.
- TMS types were inferred from the mapped stimulation order codes, producing
  six-element sequences of `C` (cTBS) and `S` (sham) per subject.

Use these files as input for downstream statistical analysis, plots,
or modeling of behavior/TMS effects.
