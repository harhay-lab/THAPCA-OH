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
o_mean <- 0.05
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = o_mean)
so_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(o_mean, 0),
                                  scale = c(so_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_so <- insight::get_parameters(so_mod)

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = o_mean)
mo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(o_mean, 0),
                                  scale = c(mo_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mo <- insight::get_parameters(mo_mod)

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = o_mean)
wo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(o_mean, 0),
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
p_mean <- -0.05
sp_sd <- so_sd
sp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(p_mean, 0),
                                  scale = c(sp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_sp <- insight::get_parameters(sp_mod)

# Moderate pessimistic prior
mp_sd <- mo_sd
mp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(p_mean, 0),
                                  scale = c(mp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mp <- insight::get_parameters(mp_mod)

# Weak pessimistic prior
wp_sd <- wo_sd
wp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(p_mean, 0),
                                  scale = c(wp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_wp <- insight::get_parameters(wp_mod)


# Prepare to make figure
# Make data set with densities for each of the 9 priors
temp_plot <- ggplot(posteriors_sn, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn <- p$data[[1]][, c(1, 2)]
plot_data_sn$Distribution <- "Posterior"

plot_data_sn2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_sn2$y <- dnorm(plot_data_sn2$x, 0, sn_sd)
plot_data_sn2$Distribution <- "Prior"

plot_data_sn <- rbind(plot_data_sn, plot_data_sn2)
plot_data_sn$barriers <- 0
plot_data_sn$barriers[plot_data_sn$x < 0] <- 1
plot_data_sn$strength <- "Strong"
plot_data_sn$belief <- "Neutral"
plot_data_sn$label <- paste0("RD = ",
                             round(median(posteriors_sn$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_sn$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_sn$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_mn, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mn <- p$data[[1]][, c(1, 2)]
plot_data_mn$Distribution <- "Posterior"

plot_data_mn2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_mn2$y <- dnorm(plot_data_mn2$x, 0, mn_sd)
plot_data_mn2$Distribution <- "Prior"

plot_data_mn <- rbind(plot_data_mn, plot_data_mn2)
plot_data_mn$barriers <- 0
plot_data_mn$barriers[plot_data_mn$x < 0] <- 1
plot_data_mn$strength <- "Moderate"
plot_data_mn$belief <- "Neutral"
plot_data_mn$label <- paste0("RD = ",
                             round(median(posteriors_mn$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_mn$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_mn$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_wn, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wn <- p$data[[1]][, c(1, 2)]
plot_data_wn$Distribution <- "Posterior"

plot_data_wn2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_wn2$y <- dnorm(plot_data_wn2$x, 0, wn_sd)
plot_data_wn2$Distribution <- "Prior"

plot_data_wn <- rbind(plot_data_wn, plot_data_wn2)
plot_data_wn$barriers <- 0
plot_data_wn$barriers[plot_data_wn$x < 0] <- 1
plot_data_wn$strength <- " Weak"
plot_data_wn$belief <- "Neutral"
plot_data_wn$label <- paste0("RD = ",
                             round(median(posteriors_wn$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_wn$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_wn$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_so, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_so <- p$data[[1]][, c(1, 2)]
plot_data_so$Distribution <- "Posterior"

plot_data_so2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_so2$y <- dnorm(plot_data_so2$x, o_mean, so_sd)
plot_data_so2$Distribution <- "Prior"

plot_data_so <- rbind(plot_data_so, plot_data_so2)
plot_data_so$barriers <- 0
plot_data_so$barriers[plot_data_so$x < 0] <- 1
plot_data_so$strength <- "Strong"
plot_data_so$belief <- " Optimistic"
plot_data_so$label <- paste0("RD = ",
                             round(median(posteriors_so$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_so$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_so$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_mo, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mo <- p$data[[1]][, c(1, 2)]
plot_data_mo$Distribution <- "Posterior"

plot_data_mo2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_mo2$y <- dnorm(plot_data_mo2$x, o_mean, mo_sd)
plot_data_mo2$Distribution <- "Prior"

plot_data_mo <- rbind(plot_data_mo, plot_data_mo2)
plot_data_mo$barriers <- 0
plot_data_mo$barriers[plot_data_mo$x < 0] <- 1
plot_data_mo$strength <- "Moderate"
plot_data_mo$belief <- " Optimistic"
plot_data_mo$label <- paste0("RD = ",
                             round(median(posteriors_mo$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_mo$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_mo$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_wo, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wo <- p$data[[1]][, c(1, 2)]
plot_data_wo$Distribution <- "Posterior"

plot_data_wo2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_wo2$y <- dnorm(plot_data_wo2$x, o_mean, wo_sd)
plot_data_wo2$Distribution <- "Prior"

plot_data_wo <- rbind(plot_data_wo, plot_data_wo2)
plot_data_wo$barriers <- 0
plot_data_wo$barriers[plot_data_wo$x < 0] <- 1
plot_data_wo$strength <- " Weak"
plot_data_wo$belief <- " Optimistic"
plot_data_wo$label <- paste0("RD = ",
                             round(median(posteriors_wo$Trt), 2),
                             "\n", "95% CrI: (",
                             round(quantile(posteriors_wo$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_wo$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_sp, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sp <- p$data[[1]][, c(1, 2)]
plot_data_sp$Distribution <- "Posterior"

plot_data_sp2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_sp2$y <- dnorm(plot_data_sp2$x, p_mean, sp_sd)
plot_data_sp2$Distribution <- "Prior"

plot_data_sp <- rbind(plot_data_sp, plot_data_sp2)
plot_data_sp$barriers <- 0
plot_data_sp$barriers[plot_data_sp$x < 0] <- 1
plot_data_sp$strength <- "Strong"
plot_data_sp$belief <- "Pessimistic"
plot_data_sp$label <- paste0("RD = ",
                             round(median(posteriors_sp$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_sp$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_sp$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_mp, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mp <- p$data[[1]][, c(1, 2)]
plot_data_mp$Distribution <- "Posterior"

plot_data_mp2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_mp2$y <- dnorm(plot_data_mp2$x, p_mean, mp_sd)
plot_data_mp2$Distribution <- "Prior"

plot_data_mp <- rbind(plot_data_mp, plot_data_mp2)
plot_data_mp$barriers <- 0
plot_data_mp$barriers[plot_data_mp$x < 0] <- 1
plot_data_mp$strength <- "Moderate"
plot_data_mp$belief <- "Pessimistic"
plot_data_mp$label <- paste0("RD = ",
                             round(median(posteriors_mp$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_mp$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_mp$Trt,
                                            0.975), 2), ")")

temp_plot <- ggplot(posteriors_wp, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wp <- p$data[[1]][, c(1, 2)]
plot_data_wp$Distribution <- "Posterior"

plot_data_wp2 <- data.frame(x = seq(-0.25, 0.25, 0.0005))
plot_data_wp2$y <- dnorm(plot_data_wp2$x, p_mean, wp_sd)
plot_data_wp2$Distribution <- "Prior"

plot_data_wp <- rbind(plot_data_wp, plot_data_wp2)
plot_data_wp$barriers <- 0
plot_data_wp$barriers[plot_data_wp$x < 0] <- 1
plot_data_wp$strength <- " Weak"
plot_data_wp$belief <- "Pessimistic"
plot_data_wp$label <- paste0("RD = ",
                             round(median(posteriors_wp$Trt), 2),
                             "\n", "95% CI: (",
                             round(quantile(posteriors_wp$Trt,
                                            0.025), 2), ", ",
                             round(quantile(posteriors_wp$Trt,
                                            0.975), 2), ")")

plot_data <- rbind(plot_data_sn, plot_data_mn, plot_data_wn,
                   plot_data_so, plot_data_mo, plot_data_wo,
                   plot_data_sp, plot_data_mp, plot_data_wp)

f_labels <-
  data.frame(belief = c(" Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic"),
             strength = c("Strong", "Moderate", " Weak",
                          "Moderate", " Weak", "Strong",
                          " Weak", "Strong", "Moderate"),
             Distribution = rep("Posterior", 9),
             barriers = rep(0, 9),
             label = c(paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_so$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_so$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_so$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_mn$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_mn$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_mn$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_wp$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_wp$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_wp$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_mo$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_mo$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_mo$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_wn$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_wn$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_wn$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_sp$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_sp$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_sp$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_wo$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_wo$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_wo$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_sn$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_sn$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_sn$Trt,
                                                       0.975)), ")"),
                       paste0("Median RD = ",
                              sprintf("%.2f", median(posteriors_mp$Trt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(posteriors_mp$Trt,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(posteriors_mp$Trt,
                                                       0.975)), ")")))

# Make figure for RD scale
plot_data$barriers[plot_data$Distribution == "Prior"] <- NA
fig1 <- ggplot(plot_data, aes(x = x, y = y, group = barriers,
                              lty = Distribution)) +
  facet_grid(belief ~ strength, scales = "free") +
  labs(x = "Risk Difference", y = "Density") +
  scale_x_continuous(labels = seq(-0.2, 0.2, by = 0.1),
                     breaks = seq(-0.2, 0.2, by = 0.1)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.25, 0.25),
                  ylim = c(0, 19)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_bw() +
  scale_fill_brewer(guide = "none") +
  geom_vline(xintercept = 0, color = "black",
             linetype = 1) +
  geom_line() +
  geom_text(data = f_labels, size = 2,
            aes(x = 0.15, y = 16, label = label))


#######################################################################
# Now do analysis on relative risk scale

# Optimistic priors
# Strong optimistic prior
o_mean <- log(1.25)
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = o_mean)
so_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(o_mean, 0),
                                  scale = c(so_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_so <- insight::get_parameters(so_mod)

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = o_mean)
mo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(o_mean, 0),
                                  scale = c(mo_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mo <- insight::get_parameters(mo_mod)

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = o_mean)
wo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(o_mean, 0),
                                  scale = c(wo_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wo <- insight::get_parameters(wo_mod)

# Neutral priors
# Strong neutral prior
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(1/1.5))
sn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(sn_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_sn <- insight::get_parameters(sn_mod)

# Moderate neutral prior
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.5))
mn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(mn_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mn <- insight::get_parameters(mn_mod)

# Weak neutral prior (SD = 5)
wn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(5, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wn <- insight::get_parameters(wn_mod)

# Pessimistic priors
# Strong pessimistic prior
p_mean <- log(1/1.25)
sp_sd <- so_sd
sp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(p_mean, 0),
                                  scale = c(sp_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_sp <- insight::get_parameters(sp_mod)

# Moderate pessimistic prior
mp_sd <- mo_sd
mp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(p_mean, 0),
                                  scale = c(mp_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mp <- insight::get_parameters(mp_mod)

# Weak pessimistic prior
wp_sd <- wo_sd
wp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(p_mean, 0),
                                  scale = c(wp_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wp <- insight::get_parameters(wp_mod)


# Prepare to make figure
# Make data set with densities for each of the 9 priors
posteriors_sn$exptrt <- exp(posteriors_sn$Trt)
temp_plot <- ggplot(posteriors_sn, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn <- p$data[[1]][, c(1, 2)]
plot_data_sn$Distribution <- "Posterior"

plot_data_sn2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_sn2$y <- dnorm(log(plot_data_sn2$x), 0, sn_sd)
plot_data_sn2$Distribution <- "Prior"

plot_data_sn <- rbind(plot_data_sn, plot_data_sn2)
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
plot_data_mn$Distribution <- "Posterior"

plot_data_mn2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_mn2$y <- dnorm(log(plot_data_mn2$x), 0, mn_sd)
plot_data_mn2$Distribution <- "Prior"

plot_data_mn <- rbind(plot_data_mn, plot_data_mn2)
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
plot_data_wn$Distribution <- "Posterior"

plot_data_wn2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_wn2$y <- dnorm(log(plot_data_wn2$x), 0, wn_sd)
plot_data_wn2$Distribution <- "Prior"

plot_data_wn <- rbind(plot_data_wn, plot_data_wn2)
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
plot_data_so$Distribution <- "Posterior"

plot_data_so2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_so2$y <- dnorm(log(plot_data_so2$x), o_mean, so_sd)
plot_data_so2$Distribution <- "Prior"

plot_data_so <- rbind(plot_data_so, plot_data_so2)
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
plot_data_mo$Distribution <- "Posterior"

plot_data_mo2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_mo2$y <- dnorm(log(plot_data_mo2$x), o_mean, mo_sd)
plot_data_mo2$Distribution <- "Prior"

plot_data_mo <- rbind(plot_data_mo, plot_data_mo2)
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
plot_data_wo$Distribution <- "Posterior"

plot_data_wo2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_wo2$y <- dnorm(log(plot_data_wo2$x), o_mean, wo_sd)
plot_data_wo2$Distribution <- "Prior"

plot_data_wo <- rbind(plot_data_wo, plot_data_wo2)
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
plot_data_sp$Distribution <- "Posterior"

plot_data_sp2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_sp2$y <- dnorm(log(plot_data_sp2$x), p_mean, sp_sd)
plot_data_sp2$Distribution <- "Prior"

plot_data_sp <- rbind(plot_data_sp, plot_data_sp2)
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
plot_data_mp$Distribution <- "Posterior"

plot_data_mp2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_mp2$y <- dnorm(log(plot_data_mp2$x), p_mean, mp_sd)
plot_data_mp2$Distribution <- "Prior"

plot_data_mp <- rbind(plot_data_mp, plot_data_mp2)
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
plot_data_wp$Distribution <- "Posterior"

plot_data_wp2 <- data.frame(x = seq(0.33, 3, 0.001))
plot_data_wp2$y <- dnorm(log(plot_data_wp2$x), p_mean, wp_sd)
plot_data_wp2$Distribution <- "Prior"

plot_data_wp <- rbind(plot_data_wp, plot_data_wp2)
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
                                            0.0975), 2), ")")

plot_data <- rbind(plot_data_sn, plot_data_mn, plot_data_wn,
                   plot_data_so, plot_data_mo, plot_data_wo,
                   plot_data_sp, plot_data_mp, plot_data_wp)

f_labels <-
  data.frame(belief = c(" Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic"),
             strength = c("Strong", "Moderate", " Weak",
                          "Moderate", " Weak", "Strong",
                          " Weak", "Strong", "Moderate"),
             Distribution = rep("Posterior", 9),
             barriers = rep(0, 9),
             label = c(paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_so$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_so$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_so$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_mn$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_mn$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_mn$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_wp$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_wp$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_wp$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_mo$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_mo$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_mo$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_wn$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_wn$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_wn$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_sp$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_sp$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_sp$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_wo$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_wo$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_wo$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_sn$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_sn$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_sn$exptrt,
                                                       0.975)), ")"),
                       paste0("Median RR = ",
                              sprintf("%.2f",
                                      median(posteriors_mp$exptrt)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f",
                                      quantile(posteriors_mp$exptrt,
                                                       0.025)), ", ",
                              sprintf("%.2f",
                                      quantile(posteriors_mp$exptrt,
                                                       0.975)), ")")))

# Make figure for RR scale
plot_data$barriers[plot_data$Distribution == "Prior"] <- NA
fig2 <- ggplot(plot_data, aes(x = x, y = y, group = barriers,
                              lty = Distribution)) +
  facet_grid(belief ~ strength, scales = "free") +
  labs(x = "Relative Risk", y = "Density") +
  scale_x_continuous(trans = "log", labels = seq(0.5, 3, by = 0.5),
                     breaks = seq(0.5, 3, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.33, 3.4),
                  ylim = c(0, 3.75)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_bw() +
  scale_fill_brewer(guide = "none") +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1) +
  geom_line() +
  geom_text(data = f_labels, size = 2,
            aes(x = 2.25, y = 3, label = label))
