# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(insight)
library(bayestestR)
library(cowplot)
library(brms)
library(broom)
library(parameters)
library(sanon)
library(rstanarm)

# Load trial data
dat <- read.csv("../outcomes.csv")

#######################################################################
# Replicate reported results
# Primary analysis
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
mantelhaen.test(dat_primary$TreatRand, dat_primary$PrimaryEndpoint,
                dat_primary$AgeGroup)
mantelhaen.test(dat_primary$TreatRand, dat_primary$PrimaryEndpoint,
                dat_primary$AgeGroup, correct = F)

mod_primary <- glm(PrimaryEndpoint ~ TreatRand + AgeGroup,
                   data = dat_primary,
                   family = binomial(link = "log"))

# Secondary analyses
dat_secondary1 <- dat[dat$SurviveM12 != 97, ]
mantelhaen.test(dat_secondary1$TreatRand, dat_secondary1$SurviveM12,
                dat_secondary1$AgeGroup, correct = F)

dat_secondary2 <- dat[dat$SurviveM12 != 97 & !is.na(dat$DeltaVabs), ]
summary(sanon(DeltaVabs ~ grp(TreatRand) + strt(AgeGroup),
              data = dat_secondary2))

#######################################################################
# Bayesian re-analysis of primary outcome

# Flat prior analysis
# SD = 100 ensures that the prior is uninformative
dat_primary$Trt <- abs(2 - dat_primary$TreatRand)
thapca_flat <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                        data = dat_primary,
                        prior_intercept = normal(0, 100),
                        prior = normal(location = c(0, 0),
                                       scale = c(100, 100)),
                        family = binomial(link = "log"), seed = 1234)

posteriors <- insight::get_parameters(thapca_flat, iterations = 10^5)

# Some summaries used for plot legend
any_harm <- sum(exp(posteriors$Trt) < 1) / nrow(posteriors)
severe_harm <- sum(exp(posteriors$Trt) < 1/1.25) / nrow(posteriors)
rope <- sum(exp(posteriors$Trt) < 1.05 & exp(posteriors$Trt) > 1/1.05) /
        nrow(posteriors)

# Figure 2
posteriors$exptrt <- exp(posteriors$Trt)
temp_plot <- ggplot(posteriors, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 1] <- 1
plot_data$barriers[plot_data$x < 1/1.25] <- 2

fig2 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Risk", y = "") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(0.5, 4, by = 0.5),
                     breaks = seq(0.5, 4, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 4.5),
                  ylim = c(0, 1)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
                  show.legend = FALSE) +
  #geom_segment(inherit.aes = FALSE,
               #data = subset(plot_data, x > 1/1.05 & x < 1.05)[
                 #seq(1, nrow(subset(plot_data, x > 1/1.05 & x < 1.05)),
                    # 5), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  geom_segment(inherit.aes = FALSE,
               data =
                 subset(plot_data, x > 1/1.05 & x < 1.05)[
                   c(1, 3, 6, 8), ],
               aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 3.5, y = 0.9, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2),
                          "\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 3.35, xmax = 3.45, ymin = 0.943, ymax = 0.973,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 3.35, xmax = 3.45, ymin = 0.903, ymax = 0.933,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 3.25, xmax = 3.35, ymin = 0.903, ymax = 0.933,
           fill = "#3182BD") +
  annotate("rect", xmin = 3.35, xmax = 3.45, ymin = 0.863, ymax = 0.893,
           fill = "#3182BD") +
  annotate("segment", x = 3.36, xend = 3.36, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.38, xend = 3.38, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.40, xend = 3.40, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.42, xend = 3.42, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.44, xend = 3.44, y = 0.853, yend = 0.823) +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1)


# Analyses for different priors
# Helper functions to determine priors
getDiff <- function(x, priormean = 0,
                    propbelow, belowcutoff) {
  propbelow - pnorm(belowcutoff, priormean, x)
}
getPriorSD <- function(priormean = 0,
                       propbelow, belowcutoff) {
  
  root <- uniroot(getDiff, propbelow = propbelow,
                  belowcutoff = belowcutoff,
                  priormean = priormean,
                  interval = c(0, 10))$root
  return(round(root*200) / 200)
  
}

# Neutral priors
# Strong neutral prior (SD = 0.205)
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(1/1.5))
sn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, sn_sd),
                   prior = normal(location = c(0, 0),
                                  scale = c(sn_sd, sn_sd)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_sn <- insight::get_parameters(sn_mod)
any_harm_sn <- sum(exp(posteriors_sn[, 2]) < 1) / dim(posteriors_sn)[1]
severe_harm_sn <- sum(exp(posteriors_sn[, 2]) < 1/1.25) /
                  dim(posteriors_sn)[1]
quantile(exp(posteriors_sn[, 2]), c(0.025, 0.5, 0.975))

# Moderate neutral prior (SD = 0.355)
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.5))
mn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, mn_sd),
                   prior = normal(location = c(0, 0),
                                  scale = c(mn_sd, mn_sd)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mn <- insight::get_parameters(mn_mod)
quantile(exp(posteriors_mn[, 2]), c(0.025, 0.5, 0.975))

# Weak neutral prior (SD = 5)
wn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(0, 0),
                                  scale = c(5, 5)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wn <- insight::get_parameters(wn_mod)
quantile(exp(posteriors_wn[, 2]), c(0.025, 0.5, 0.975))

# Optimistic priors
# Strong optimistic prior (SD = 0.245)
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = log(1.5))
so_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(log(1.5), 0),
                                  scale = c(so_sd, 5)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_so <- insight::get_parameters(so_mod)
quantile(exp(posteriors_so[, 2]), c(0.025, 0.5, 0.975))

# Moderate optimistic prior (SD = 0.390)
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = log(1.5))
mo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(log(1.5), 0),
                                  scale = c(mo_sd, 5)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mo <- insight::get_parameters(mo_mod)
quantile(exp(posteriors_mo[, 2]), c(0.025, 0.5, 0.975))

# Weak optimistic prior (SD = 0.775)
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = log(1.5))
wo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(log(1.5), 0),
                                  scale = c(wo_sd, 5)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wo <- insight::get_parameters(wo_mod)
quantile(exp(posteriors_wo[, 2]), c(0.025, 0.5, 0.975))

# Pessimistic priors
# Strong pessimistic prior (SD = 0.245)
sp_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = -log(1/1.5))
sp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(log(1/1.5), 0),
                                  scale = c(sp_sd, 5)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_sp <- insight::get_parameters(sp_mod)
quantile(exp(posteriors_sp[, 2]), c(0.025, 0.5, 0.975))

# Moderate pessimistic prior (SD = 0.390)
mp_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = -log(1/1.5))
mp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(log(1/1.5), 0),
                                  scale = c(mp_sd, 5)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mp <- insight::get_parameters(mp_mod)
quantile(exp(posteriors_mp[, 2]), c(0.025, 0.5, 0.975))

# Weak pessimistic prior (SD = 0.775)
wp_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = -log(1/1.5))
wp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(log(1/1.5), 0),
                                  scale = c(wp_sd, 5)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wp <- insight::get_parameters(wp_mod)
quantile(exp(posteriors_wp[, 2]), c(0.025, 0.5, 0.975))

# Figure 3
# Make data set with densities for each of the 9 priors
posteriors_sn$exptrt <- exp(posteriors_sn$Trt)
temp_plot <- ggplot(posteriors_sn, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn <- p$data[[1]][, c(1, 2)]
plot_data_sn$barriers <- 0
plot_data_sn$barriers[plot_data_sn$x < 1] <- 1
plot_data_sn$strength <- " Strong"
plot_data_sn$belief <- "Neutral"
plot_data_sn$label <- paste0("RR = ",
                             round(median(posteriors_sn$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_sn$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_sn$exptrt,
                                            0.975), 2), ")")

posteriors_mn$exptrt <- exp(posteriors_mn$Trt)
temp_plot <- ggplot(posteriors_mn, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mn <- p$data[[1]][, c(1, 2)]
plot_data_mn$barriers <- 0
plot_data_mn$barriers[plot_data_mn$x < 1] <- 1
plot_data_mn$strength <- "Moderate"
plot_data_mn$belief <- "Neutral"
plot_data_mn$label <- paste0("RR = ",
                             round(median(posteriors_mn$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_mn$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_mn$exptrt,
                                            0.975), 2), ")")

posteriors_wn$exptrt <- exp(posteriors_wn$Trt)
temp_plot <- ggplot(posteriors_wn, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wn <- p$data[[1]][, c(1, 2)]
plot_data_wn$barriers <- 0
plot_data_wn$barriers[plot_data_wn$x < 1] <- 1
plot_data_wn$strength <- "Weak"
plot_data_wn$belief <- "Neutral"
plot_data_wn$label <- paste0("RR = ",
                             round(median(posteriors_wn$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_wn$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_wn$exptrt,
                                            0.975), 2), ")")

posteriors_so$exptrt <- exp(posteriors_so$Trt)
temp_plot <- ggplot(posteriors_so, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_so <- p$data[[1]][, c(1, 2)]
plot_data_so$barriers <- 0
plot_data_so$barriers[plot_data_so$x < 1] <- 1
plot_data_so$strength <- " Strong"
plot_data_so$belief <- " Optimistic"
plot_data_so$label <- paste0("RR = ",
                             round(median(posteriors_so$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_so$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_so$exptrt,
                                            0.975), 2), ")")

posteriors_mo$exptrt <- exp(posteriors_mo$Trt)
temp_plot <- ggplot(posteriors_mo, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mo <- p$data[[1]][, c(1, 2)]
plot_data_mo$barriers <- 0
plot_data_mo$barriers[plot_data_mo$x < 1] <- 1
plot_data_mo$strength <- "Moderate"
plot_data_mo$belief <- " Optimistic"
plot_data_mo$label <- paste0("RR = ",
                             round(median(posteriors_mo$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_mo$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_mo$exptrt,
                                            0.975), 2), ")")

posteriors_wo$exptrt <- exp(posteriors_wo$Trt)
temp_plot <- ggplot(posteriors_wo, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wo <- p$data[[1]][, c(1, 2)]
plot_data_wo$barriers <- 0
plot_data_wo$barriers[plot_data_wo$x < 1] <- 1
plot_data_wo$strength <- "Weak"
plot_data_wo$belief <- " Optimistic"
plot_data_wo$label <- paste0("RR = ",
                             round(median(posteriors_wo$exptrt), 2),
                             "\n", "95% CrI: (",
                             round(quantile(posteriors_wo$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_wo$exptrt,
                                            0.975), 2), ")")

posteriors_sp$exptrt <- exp(posteriors_sp$Trt)
temp_plot <- ggplot(posteriors_sp, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sp <- p$data[[1]][, c(1, 2)]
plot_data_sp$barriers <- 0
plot_data_sp$barriers[plot_data_sp$x < 1] <- 1
plot_data_sp$strength <- " Strong"
plot_data_sp$belief <- "Pessimistic"
plot_data_sp$label <- paste0("RR = ",
                             round(median(posteriors_sp$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_sp$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_sp$exptrt,
                                            0.975), 2), ")")

posteriors_mp$exptrt <- exp(posteriors_mp$Trt)
temp_plot <- ggplot(posteriors_mp, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mp <- p$data[[1]][, c(1, 2)]
plot_data_mp$barriers <- 0
plot_data_mp$barriers[plot_data_mp$x < 1] <- 1
plot_data_mp$strength <- "Moderate"
plot_data_mp$belief <- "Pessimistic"
plot_data_mp$label <- paste0("RR = ",
                             round(median(posteriors_mp$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_mp$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_mp$exptrt,
                                            0.975), 2), ")")

posteriors_wp$exptrt <- exp(posteriors_wp$Trt)
temp_plot <- ggplot(posteriors_wp, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wp <- p$data[[1]][, c(1, 2)]
plot_data_wp$barriers <- 0
plot_data_wp$barriers[plot_data_wp$x < 1] <- 1
plot_data_wp$strength <- "Weak"
plot_data_wp$belief <- "Pessimistic"
plot_data_wp$label <- paste0("RR = ",
                             round(median(posteriors_wp$exptrt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_wp$exptrt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_wp$exptrt,
                                            0.975), 2), ")")

plot_data <- rbind(plot_data_sn, plot_data_mn, plot_data_wn,
                   plot_data_so, plot_data_mo, plot_data_wo,
                   plot_data_sp, plot_data_mp, plot_data_wp)

f_labels <-
  data.frame(belief = c(" Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic"),
             strength = c(" Strong", "Moderate", "Weak",
                          "Moderate", "Weak", " Strong",
                          "Weak", " Strong", "Moderate"),
             barriers = rep(0, 9),
             label = c(paste0("Median RR = ",
                              round(median(posteriors_so$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_so$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_so$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_mn$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_mn$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_mn$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_wp$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_wp$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_wp$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_mo$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_mo$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_mo$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_wn$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_wn$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_wn$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_sp$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_sp$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_sp$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_wo$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_wo$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_wo$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_sn$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_sn$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_sn$exptrt,
                                             0.975), 2), ")"),
                       paste0("Median RR = ",
                              round(median(posteriors_mp$exptrt), 2),
                              "\n", "95% CrI: (",
                              round(quantile(posteriors_mp$exptrt,
                                             0.025), 2), ", ",
                              round(quantile(posteriors_mp$exptrt,
                                             0.975), 2), ")")))

fig3 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  facet_grid(belief ~ strength, scales = "free") +
  labs(x = "Relative Risk", y = "") +
  scale_x_continuous(trans = "log", labels = seq(0.5, 3, by = 0.5),
                     breaks = seq(0.5, 3, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 3),
                  ylim = c(0, 2.4)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  #geom_segment(inherit.aes = FALSE,
  #data = subset(plot_data, x > 1/1.05 & x < 1.05)[
  #seq(1, nrow(subset(plot_data, x > 1/1.05 & x < 1.05)),
  # 5), ],
  #aes(x = x, y = 0, xend = x, yend = y)) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_bw() +
  scale_fill_brewer(guide = "none") +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1) +
  geom_text(data = f_labels, size = 2,
            aes(x = 2.1, y = 2, label = label))


#######################################################################
# Bayesian re-analysis of secondary outcome 1

# Flat prior
thapca_flat2 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
                    data = dat_secondary1,
                    prior = flat_prior,
                    family = "bernoulli", seed = 123)
posteriors <- insight::get_parameters(thapca_flat2, iterations = 10^5)

# Summarizes the flat prior analysis
summary(thapca_flat2)

# Some summaries used for plot legend
any_harm <- sum(posteriors$b_TreatRand > 0) / nrow(posteriors)
severe_harm <- sum(posteriors$b_TreatRand > log(1.25)) / nrow(posteriors)
rope <- sum(posteriors$b_TreatRand < log(1.1) &
              posteriors$b_TreatRand > log(1/1.1)) /
  nrow(posteriors)

# Neutral moderate prior (SD = 0.355)
m4 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
          data = dat_secondary1, prior = m1priors,
          family = "bernoulli", seed = 123)
posteriors_m4 <- insight::get_parameters(m4)[, 2]
any_harm_m4 <- sum(posteriors_m4 > 0) / length(posteriors_m4)
severe_harm_m4 <- sum(posteriors_m4 > log(1.25)) / length(posteriors_m4)

# Optimistic moderate prior (SD = 0.400)
m5 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
          data = dat_secondary1, prior = m2priors,
          family = "bernoulli", seed = 123)
posteriors_m5 <- insight::get_parameters(m5)[, 2]
any_harm_m5 <- sum(posteriors_m5 > 0) / length(posteriors_m5)
severe_harm_m5 <- sum(posteriors_m5 > log(1.25)) / length(posteriors_m5)

# Pessimistic weak prior (SD = 0.790)
m6 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
          data = dat_secondary1, prior = m3priors,
          family = "bernoulli", seed = 123)
posteriors_m6 <- insight::get_parameters(m6)[, 2]
any_harm_m6 <- sum(posteriors_m6 > 0) / length(posteriors_m6)
severe_harm_m6 <- sum(posteriors_m6 > log(1.25)) / length(posteriors_m6)
