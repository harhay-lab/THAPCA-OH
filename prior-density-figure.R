# Make figure displaying each of the prior distributions
# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(flextable)
library(insight)
library(bayestestR)
library(bayesplot)
library(cowplot)
library(brms)
library(broom)
library(multipanelfigure)
library(parameters)
library(sanon)
#library(rstanarm)

# Load trial data
dat <- read.csv("/Users/blette/Downloads/outcomes.csv")
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
dat_primary$Trt <- abs(2 - dat_primary$TreatRand)
dat_primary$AgeFactor <- as.factor(dat_primary$AgeGroup)

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

mod_primary_or <- glm(PrimaryEndpoint ~ Trt + AgeFactor,
                      data = dat_primary,
                      family = binomial(link = "logit"))

# First do each of the 9 standardized priors on the log(OR) scale
# Optimistic priors
# Strong optimistic prior
o_mean <- log(1.25)
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = o_mean)
ess_so <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
          (so_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, o_mean, so_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_so <- p$data[[1]][, c(1, 2)]
plot_data_so$Strength <- "Strong"
plot_data_so$Belief <- "Optimistic"

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = o_mean)
ess_mo <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
          (mo_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, o_mean, mo_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mo <- p$data[[1]][, c(1, 2)]
plot_data_mo$Strength <- "Moderate"
plot_data_mo$Belief <- "Optimistic"

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = o_mean)
ess_wo <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
          (wo_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, o_mean, wo_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wo <- p$data[[1]][, c(1, 2)]
plot_data_wo$Strength <- "Weak"
plot_data_wo$Belief <- "Optimistic"

# Neutral priors
# Strong neutral prior
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(1/1.5),
                    priormean = 0)
ess_sn <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
  (sn_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, 0, sn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn <- p$data[[1]][, c(1, 2)]
plot_data_sn$Strength <- "Strong"
plot_data_sn$Belief <- "Neutral"

# Moderate neutral prior
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.5),
                    priormean = 0)
ess_mn <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
  (mn_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, 0, mn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mn <- p$data[[1]][, c(1, 2)]
plot_data_mn$Strength <- "Moderate"
plot_data_mn$Belief <- "Neutral"

# Weak neutral prior
wn_sd <- 3
ess_wn <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
  (wn_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, 0, wn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wn <- p$data[[1]][, c(1, 2)]
plot_data_wn$Strength <- "Weak"
plot_data_wn$Belief <- "Neutral"

# Pessimistic priors
# Strong pessimistic prior
p_mean <- log(1/1.25)
sp_sd <- so_sd
ess_sp <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
  (sp_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, p_mean, sp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sp <- p$data[[1]][, c(1, 2)]
plot_data_sp$Strength <- "Strong"
plot_data_sp$Belief <- "Pessimistic"

# Moderate pessimistic prior
mp_sd <- mo_sd
ess_mp <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
  (mp_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, p_mean, mp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mp <- p$data[[1]][, c(1, 2)]
plot_data_mp$Strength <- "Moderate"
plot_data_mp$Belief <- "Pessimistic"

# Weak pessimistic prior
wp_sd <- wo_sd
ess_wp <- dim(dat_primary)[1] * summary(mod_primary_or)$coefficients[2, 2]^2 /
  (wp_sd^2)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, p_mean, wp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wp <- p$data[[1]][, c(1, 2)]
plot_data_wp$Strength <- "Weak"
plot_data_wp$Belief <- "Pessimistic"


# Combine each of the priors into one data set
plot_data <- rbind(plot_data_so, plot_data_mo, plot_data_wo,
                   plot_data_sn, plot_data_mn, plot_data_wn,
                   plot_data_sp, plot_data_mp, plot_data_wp)
plot_data$Strength <- factor(plot_data$Strength,
                             levels = c("Strong", "Moderate", "Weak"))

# Make full priors figure
f_labels <-
  data.frame(Belief = c("Optimistic", "Neutral", "Pessimistic"),
             Strength = "Strong",
             label = c(paste0("Strong ESS = ", sprintf("%.0f", ess_so), "\n",
                              "Moderate ESS = ", sprintf("%.0f", ess_mo), "\n",
                              "Weak ESS = ", sprintf("%.0f", ess_wo)),
                       paste0("Strong ESS = ", sprintf("%.0f", ess_sn), "\n",
                              "Moderate ESS = ", sprintf("%.0f", ess_mn), "\n",
                              "Weak ESS ~ 0"),
                       paste0("Strong ESS = ", sprintf("%.0f", ess_sp), "\n",
                              "Moderate ESS = ", sprintf("%.0f", ess_mp), "\n",
                              "Weak ESS = ", sprintf("%.0f", ess_wp))))
                        
p1 <- ggplot(plot_data, aes(x = x, y = y, lty = Strength)) +
  facet_wrap(vars(Belief), nrow = 3) +
  xlim(c(-1.5, 1.5)) +
  geom_line() +
  xlab("Log(OR)") +
  ylab("Density") +
  geom_vline(xintercept = 0, color = "black",
             linetype = 1) +
  geom_text(data = f_labels, size = 2,
            aes(x = 1, y = 2, label = label)) +
  theme_bw()


# Make additional panel for evidence-based priors
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

mod_primary_rr <- glm(PrimaryEndpoint ~ Trt + AgeFactor,
                      data = dat_primary,
                      family = binomial(link = "log"))

# Prior using Grandfelt 2021 meta-analysis and 50% weighting
grandfelt_mean <- log(1.2)
grandfelt_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.85),
                           priormean = grandfelt_mean)
dw <- 0.5
ess_gf <- dim(dat_primary)[1] * summary(mod_primary_rr)$coefficients[2, 2]^2 /
  sqrt(grandfelt_sd^2 / dw)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, grandfelt_mean,
                                sqrt(grandfelt_sd^2 / dw))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_gf <- p$data[[1]][, c(1, 2)]
plot_data_gf$Study <- "Grandfelt"

# Prior using TTM and TTM2 and 50% weighting
ttm_mean <- log(1)
ttm_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.9),
                     priormean = ttm_mean)
dw <- 0.5
ess_ttm <- dim(dat_primary)[1] * summary(mod_primary_rr)$coefficients[2, 2]^2 /
  sqrt(ttm_sd^2 / dw)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, ttm_mean,
                                sqrt(ttm_sd^2 / dw))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_ttm <- p$data[[1]][, c(1, 2)]
plot_data_ttm$Study <- "TTM"

# Prior using Hyperion with 50% weighting
hyp_mean <- log(1.8)
hyp_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.7),
                     priormean = hyp_mean)
dw <- 0.5
ess_hyp <- dim(dat_primary)[1] * summary(mod_primary_rr)$coefficients[2, 2]^2 /
  sqrt(hyp_sd^2 / dw)

temp_plot <-
  ggplot(data.frame(pri = rnorm(1e6, hyp_mean,
                                sqrt(hyp_sd^2 / dw))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_hyp <- p$data[[1]][, c(1, 2)]
plot_data_hyp$Study <- "Hyperion"

# Combine each of the priors into one data set
plot_data_ebp <- rbind(plot_data_gf, plot_data_ttm, plot_data_hyp)

# Make additional panel for evidence-based priors
p2 <- ggplot(plot_data_ebp, aes(x = x, y = y, lty = Study)) +
  xlim(c(-1, 2)) +
  geom_line() +
  xlab("Log(RB)") +
  ylab("Density") +
  geom_vline(xintercept = 0, color = "black",
             linetype = 1) +
  theme_bw()


# Combine panels into one figure
# Then Bayesian analyses on RD scale (first version of fig)
figure1 <- multi_panel_figure(columns = 2, rows = 1, width = 200,
                              height = 150)
figure1 <- fill_panel(figure1, p1)
figure1 <- fill_panel(figure1, p2)

# Output pdf of figure, dims may need to be changed
pdf("prior-density-figure.pdf", width = 8.5, height = 6.5)
figure1
dev.off()
