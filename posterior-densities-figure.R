# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(flextable)
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
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
dat_primary$Trt <- abs(2 - dat_primary$TreatRand)

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


#######################################################################
# First do analysis on risk difference scale

# Optimistic priors
# Strong optimistic prior
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = 0.1)
so_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0.1, 0),
                                  scale = c(so_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_so <- insight::get_parameters(so_mod)

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = 0.1)
mo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0.1, 0),
                                  scale = c(mo_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mo <- insight::get_parameters(mo_mod)

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = 0.1)
wo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0.1, 0),
                                  scale = c(wo_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_wo <- insight::get_parameters(wo_mod)

# Neutral priors
# Strong neutral prior
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = -0.05,
                    priormean = 0)
sn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(sn_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_sn <- insight::get_parameters(sn_mod)

# Moderate neutral prior
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = -0.1,
                    priormean = 0)
mn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(mn_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mn <- insight::get_parameters(mn_mod)

# Weak neutral prior
wn_sd <- 5
wn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(wn_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_wn <- insight::get_parameters(wn_mod)

# Pessimistic priors
# Strong pessimistic prior
sp_sd <- so_sd
sp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(-0.1, 0),
                                  scale = c(sp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_sp <- insight::get_parameters(sp_mod)

# Moderate pessimistic prior
mp_sd <- mo_sd
mp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(mp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mp <- insight::get_parameters(mp_mod)

# Weak pessimistic prior
wp_sd <- wo_sd
wp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(wp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_wp <- insight::get_parameters(wp_mod)


# Make figure (change below to RD)
# Make data set with densities for each of the 9 priors
posteriors_sn$exptrt <- exp(posteriors_sn$Trt)
temp_plot <- ggplot(posteriors_sn, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn <- p$data[[1]][, c(1, 2)]
plot_data_sn$barriers <- 0
plot_data_sn$barriers[plot_data_sn$x < 1] <- 1
plot_data_sn$strength <- "Strong"
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
plot_data_wn$strength <- " Weak"
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
plot_data_so$strength <- "Strong"
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
plot_data_wo$strength <- " Weak"
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
plot_data_sp$strength <- "Strong"
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
plot_data_wp$strength <- " Weak"
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
# Now do analysis on relative risk scale