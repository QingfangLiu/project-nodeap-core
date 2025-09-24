
########################################################
# Plot empirical discrimination learning overlaid with
# model-predicted ribbons from posterior w
#
# Creates
#   p1: color by Cond_day1 (cTBS/sham), facet by StimLoc
#   p2: color by Sess (1/2/3),         facet by StimLoc
#
# Model ribbon
#   • Trial prob uses pre-update active pair:
#       w[draw, sub, sess, t, CuePair[sub,sess,t]]
#   • Run summary = mean(w) over the 24 trials in each run
#   • Per draw: average across subjects to get group mean
#   • Across draws: 90% CrI for ribbon (change cred_level to adjust)
#
# Inputs (from Setup.R / results)
#   • scripts/utils/Setup.R  (paths, strip theme, color palettes:
#       use.col.conds.day1, use.col.sess, and `common` theme chunk)
#   • processed_dir/Conditioning.RData (conditioning_dat)
#   • beh_model_dir/results/posterior_learning_model.rds
#       sims.list$w [draw × sub × sess × trial_idx × cuepair]
#   • beh_model_dir/results/posterior_design_arrays.rds
#       CuePair [sub × sess × trial], run_idx [sub × sess × trial],
#       SubID_levels, meta_ss (SubID, Sess, Cond_day1, StimLoc)
########################################################

rm(list = ls())
project_folder <- "/Users/liuq13/project-nodeap-core"
source(file.path(project_folder, "scripts", "utils", "Setup.R"))

# ---- Load data ----
load(file.path(processed_dir, "Conditioning.RData"))  # conditioning_dat

post   <- readRDS(file.path(beh_model_dir, "results", "posterior_learning_model.rds"))
design <- readRDS(file.path(beh_model_dir, "results", "posterior_design_arrays.rds"))

w_sims  <- post$BUGSoutput$sims.list$w
CuePair <- design$CuePair
run_idx <- design$run_idx
meta_ss <- design$meta_ss %>% mutate(SubID = as.character(SubID), Sess = as.integer(Sess))
sub_lookup <- tibble(sub = seq_along(design$SubID_levels),
                     SubID = as.character(design$SubID_levels))

# ---- Config ----
cred_level <- 0.90
lo <- (1 - cred_level)/2; hi <- 1 - lo

# ---- Dimensions / guards ----
n_draw  <- dim(w_sims)[1]
n_sub   <- dim(w_sims)[2]
n_sess  <- dim(w_sims)[3]
ntrials <- dim(w_sims)[4] - 1L
nruns   <- max(run_idx[1, 1, ])
runs_vec <- seq_len(nruns)
stopifnot(ntrials == sum(run_idx[1,1,] == 1L) * nruns)

# ---- Build run-wise summaries per draw ----
get_run_values <- function(d, j, c) {
    # pre-update active-pair prob per trial
    wt <- vapply(seq_len(ntrials),
                 function(t) w_sims[d, j, c, t, CuePair[j, c, t]],
                 numeric(1))
    vapply(runs_vec, function(r) mean(wt[run_idx[j, c, ] == r]), numeric(1))
}

# Collect to long df: one row per (draw, sub, sess, run)
run_vals <- lapply(seq_len(n_draw), function(d) {
  do.call(rbind, lapply(seq_len(n_sub), function(j) {
    do.call(rbind, lapply(seq_len(n_sess), function(c) {
      data.frame(draw = d, sub = j, sess = c, run = runs_vec,
                 w = get_run_values(d, j, c))
    }))
  }))
}) |> bind_rows()

# Attach SubID / StimLoc / Cond_day1, and nice factors
run_vals <- run_vals %>%
  left_join(sub_lookup, by = "sub") %>%
  left_join(meta_ss, by = c("SubID", "sess" = "Sess")) %>%
  mutate(
    StimLoc   = factor(StimLoc, levels = c("aOFC","pOFC")),
    Cond_day1 = factor(Cond_day1, levels = c("cTBS","sham")),
    sess_f    = factor(sess)  # labeled 1/2/3 for plotting
  )

# ---- Posterior ribbons: subject-averaged group mean per draw → CrI ----
# by Cond_day1
ci_cond <- run_vals %>%
  group_by(draw, StimLoc, Cond_day1, run) %>%
  reframe(w = mean(w)) %>%
  group_by(StimLoc, Cond_day1, run) %>%
  reframe(w_lo = quantile(w, lo), w_hi = quantile(w, hi))

# by Sess
ci_sess <- run_vals %>%
  group_by(draw, StimLoc, sess_f, run) %>%
  reframe(w = mean(w)) %>%
  group_by(StimLoc, sess_f, run) %>%
  reframe(w_lo = quantile(w, lo), w_hi = quantile(w, hi))

# ---- Empirical summaries for overlay ----
# by Cond_day1
emp_cond <- conditioning_dat %>%
  filter(!is.na(OdorChosen)) %>%
  mutate(
    StimLoc   = factor(StimLoc, levels = c("aOFC","pOFC")),
    Cond_day1 = factor(Cond_day1, levels = c("cTBS","sham")),
    Run       = as.integer(Run)
  ) %>%
  group_by(SubID, StimLoc, Cond_day1, Run) %>%
  reframe(choice = mean(OdorChosen)) %>%
  group_by(StimLoc, Cond_day1, Run) %>%
  reframe(
    mean = mean(choice),
    se   = sd(choice)/sqrt(n()),
    lo   = mean - 1.96*se,
    hi   = mean + 1.96*se)

# by Sess
emp_sess <- conditioning_dat %>%
  filter(!is.na(OdorChosen)) %>%
  mutate(
    StimLoc = factor(StimLoc, levels = c("aOFC","pOFC")),
    Sess    = factor(Sess),
    Run     = as.integer(Run)
  ) %>%
  group_by(SubID, StimLoc, Sess, Run) %>%
  reframe(choice = mean(OdorChosen)) %>%
  group_by(StimLoc, Sess, Run) %>%
  reframe(
    mean = mean(choice),
    se   = sd(choice)/sqrt(n()),
    lo   = mean - 1.96*se,
    hi   = mean + 1.96*se)

strip <- strip_themed(background_x = elem_list_rect(fill = use.col.ap.ofc),
                      text_x = elem_list_text(color = 'white',
                                              face = "bold",
                                              size = 16))

# ---- p1: color by Cond_day1, facet by StimLoc ----
p1 <- ggplot(ci_cond) +
  facet_wrap2(~ StimLoc, nrow = 1, strip = strip, scales = 'free_y') +
  geom_ribbon(aes(x = run, ymin = w_lo, ymax = w_hi, fill = Cond_day1,
                  group = interaction(StimLoc, Cond_day1)),
              alpha = 0.5, linewidth = 0, show.legend = FALSE) +
  scale_fill_manual(values = use.col.conds.day1) +
  scale_color_manual(values = use.col.conds.day1) +
  geom_errorbar(data = emp_cond, # empirical overlay
                aes(x = Run, ymin = lo, ymax = hi, color = Cond_day1),
                width = 0.3, linewidth = 1, show.legend = FALSE) +
  stat_summary(data = emp_cond,
               aes(x = Run, y = mean, color = Cond_day1),
               fun = mean, geom = "line", linewidth = 1) +
  scale_x_continuous(breaks = runs_vec) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Runs", y = "P (Choosing rewarding stim)", color = "Day1 TMS") +
  common + theme(legend.position.inside = c(0.85, 0.3))

# ---- p2: color by Sess, facet by StimLoc ----
p2 <- ggplot(ci_sess) +
  facet_wrap2(~ StimLoc, nrow = 1, strip = strip, scales = 'free_y') +
  geom_ribbon(aes(x = run, ymin = w_lo, ymax = w_hi, fill = sess_f,
                  group = interaction(StimLoc, sess_f)),
              alpha = 0.5, linewidth = 0, show.legend = FALSE) +
  scale_fill_manual(values = use.col.sess) +
  scale_color_manual(values = use.col.sess) +
  geom_errorbar(data = emp_sess, # empirical overlay
                aes(x = Run, ymin = lo, ymax = hi, color = Sess),
                width = 0.3, linewidth = 1, show.legend = FALSE) +
  stat_summary(data = emp_sess,
               aes(x = Run, y = mean, color = Sess),
               fun = mean, geom = "line", linewidth = 1) +
  scale_x_continuous(breaks = runs_vec) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Runs", y = "P (Choosing rewarding stim)", color = "Session") +
  common + theme(legend.position.inside = c(0.85, 0.3))

pdf(file.path(FigPaperDir,'Conditioning_acc_obs_sim.pdf'),9,8)
ggarrange(p1,p2,nrow = 2)
dev.off()
