# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(insight)
library(bayestestR)
library(cowplot)
library(brms)
library(broom)
library(ggpubr)
library(parameters)
library(multipanelfigure)
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

mod_primary <- glm(PrimaryEndpoint ~ TreatRand + factor(AgeGroup),
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

# First figure panel
posteriors$exptrt <- exp(posteriors$Trt)
temp_plot <- ggplot(posteriors, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 1] <- 1
plot_data$barriers[plot_data$x < 1/1.25] <- 2

fig1 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
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
  annotate("text", x = 3.2, y = 0.9, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2),
                          "\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 3.05, xmax = 3.15, ymin = 0.943, ymax = 0.973,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 3.05, xmax = 3.15, ymin = 0.903, ymax = 0.933,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 2.95, xmax = 3.05, ymin = 0.903, ymax = 0.933,
           fill = "#3182BD") +
  annotate("rect", xmin = 3.05, xmax = 3.15, ymin = 0.863, ymax = 0.893,
           fill = "#3182BD") +
  annotate("segment", x = 3.06, xend = 3.06, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.08, xend = 3.08, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.10, xend = 3.10, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.12, xend = 3.12, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.14, xend = 3.14, y = 0.853, yend = 0.823) +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1)

#######################################################################
# Primary outcome on RD scale
thapca_flat2 <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                         data = dat_primary,
                         prior_intercept = normal(0, 100),
                         prior = normal(location = c(0, 0),
                                        scale = c(100, 100)),
                         family = gaussian(link = "identity"),
                         seed = 1234)

posteriors <- insight::get_parameters(thapca_flat2, iterations = 10^5)

# Some summaries used for plot legend
any_harm <- sum(posteriors$Trt < 0) / nrow(posteriors)
severe_harm <- sum(posteriors$Trt < -0.05) / nrow(posteriors)
rope <- sum(posteriors$Trt < 0.01 & posteriors$Trt > -0.01) /
  nrow(posteriors)

# Second figure panel
temp_plot <- ggplot(posteriors, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 0] <- 1
plot_data$barriers[plot_data$x < -0.05] <- 2

fig2 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Risk Difference", y = "") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(-0.1, 0.2, by = 0.05),
                     breaks = seq(-0.1, 0.2, by = 0.05)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.13, 0.3),
                  ylim = c(0, 9)) +
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
                 subset(plot_data, x > -0.01 & x < 0.01)[
                   c(1, 9, 24, 32), ],
               aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 0.165, y = 7.5, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) < 0.01",
                          #round(severe_harm, 2),
                          "\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 0.152, xmax = 0.160, ymin = 7.85, ymax = 8.15,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 0.152, xmax = 0.160, ymin = 7.5, ymax = 7.8,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 0.144, xmax = 0.152, ymin = 7.5, ymax = 7.8,
           fill = "#3182BD") +
  annotate("rect", xmin = 0.152, xmax = 0.160, ymin = 7.15, ymax = 7.45,
           fill = "#3182BD") +
  annotate("segment", x = 0.152, xend = 0.152, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.154, xend = 0.154, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.156, xend = 0.156, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.158, xend = 0.158, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.160, xend = 0.160, y = 6.8, yend = 7.07) +
  geom_vline(xintercept = 0, color = "black",
             linetype = 1)

#######################################################################
# Bayesian re-analysis of secondary outcome

# Flat prior
dat_secondary1$Trt <- abs(2 - dat_secondary1$TreatRand)
thapca_flat3 <- stan_glm(SurviveM12 ~ Trt + AgeGroup,
                         data = dat_secondary1,
                         prior_intercept = normal(0, 100),
                         prior = normal(location = c(0, 0),
                                        scale = c(100, 100)),
                         family = binomial(link = "log"), seed = 1234)

posteriors <- insight::get_parameters(thapca_flat3, iterations = 10^5)

# Some summaries used for plot legend
any_harm <- sum(exp(posteriors$Trt) < 1) / nrow(posteriors)
severe_harm <- sum(exp(posteriors$Trt) < 1/1.25) / nrow(posteriors)
rope <- sum(exp(posteriors$Trt) < 1.05 &
            exp(posteriors$Trt) > 1/1.05) / nrow(posteriors)

# Third figure panel
posteriors$exptrt <- exp(posteriors$Trt)
temp_plot <- ggplot(posteriors, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 1] <- 1
plot_data$barriers[plot_data$x < 1/1.25] <- 2

fig3 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Risk", y = "") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(0.5, 2.5, by = 0.5),
                     breaks = seq(0.5, 2.5, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 2.5),
                  ylim = c(0, 2)) +
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
                   c(1, 7, 21, 28), ],
               aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 1.8, y = 1.8, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) < 0.01",
                          #round(severe_harm, 2),
                          "\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 1.73, xmax = 1.78, ymin = 1.88, ymax = 1.95,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 1.73, xmax = 1.78, ymin = 1.8, ymax = 1.87,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 1.68, xmax = 1.73, ymin = 1.8, ymax = 1.87,
           fill = "#3182BD") +
  annotate("rect", xmin = 1.73, xmax = 1.78, ymin = 1.72, ymax = 1.79,
           fill = "#3182BD") +
  annotate("segment", x = 1.735, xend = 1.735, y = 1.64, yend = 1.71) +
  annotate("segment", x = 1.745, xend = 1.745, y = 1.64, yend = 1.71) +
  annotate("segment", x = 1.755, xend = 1.755, y = 1.64, yend = 1.71) +
  annotate("segment", x = 1.765, xend = 1.765, y = 1.64, yend = 1.71) +
  annotate("segment", x = 1.775, xend = 1.775, y = 1.64, yend = 1.71) +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1)

#######################################################################
# Secondary outcome on RD scale
thapca_flat4 <- stan_glm(SurviveM12 ~ Trt + AgeGroup,
                         data = dat_secondary1,
                         prior_intercept = normal(0, 100),
                         prior = normal(location = c(0, 0),
                                        scale = c(100, 100)),
                         family = gaussian(link = "identity"),
                         seed = 1234)

posteriors <- insight::get_parameters(thapca_flat4, iterations = 10^5)

# Some summaries used for plot legend
any_harm <- sum(posteriors$Trt < 0) / nrow(posteriors)
severe_harm <- sum(posteriors$Trt < -0.05) / nrow(posteriors)
rope <- sum(posteriors$Trt < 0.01 & posteriors$Trt > -0.01) /
  nrow(posteriors)

# Fourth figure panel
temp_plot <- ggplot(posteriors, aes(x = Trt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 0] <- 1
plot_data$barriers[plot_data$x < -0.05] <- 2

fig4 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Risk Difference", y = "") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(-0.1, 0.2, by = 0.05),
                     breaks = seq(-0.1, 0.2, by = 0.05)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.13, 0.3),
                  ylim = c(0, 9)) +
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
                 subset(plot_data, x > -0.01 & x < 0.01)[
                   c(1, 7, 19, 25), ],
               aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 0.165, y = 7.5, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2),
                          "\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 0.152, xmax = 0.160, ymin = 7.85, ymax = 8.15,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 0.152, xmax = 0.160, ymin = 7.5, ymax = 7.8,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 0.144, xmax = 0.152, ymin = 7.5, ymax = 7.8,
           fill = "#3182BD") +
  annotate("rect", xmin = 0.152, xmax = 0.160, ymin = 7.15, ymax = 7.45,
           fill = "#3182BD") +
  annotate("segment", x = 0.152, xend = 0.152, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.154, xend = 0.154, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.156, xend = 0.156, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.158, xend = 0.158, y = 6.8, yend = 7.07) +
  annotate("segment", x = 0.160, xend = 0.160, y = 6.8, yend = 7.07) +
  geom_vline(xintercept = 0, color = "black",
             linetype = 1)

# Combine figures
# Both figures exported as landscape PDFs with 5 in x 10 in dimensions
ggarrange(fig1, fig3, ncol = 2)
ggarrange(fig2, fig4, ncol = 2)

figure <- multi_panel_figure(columns = 2, rows = 1, width = 240,
                             height = 125)
figure <- fill_panel(figure, fig1)
figure <- fill_panel(figure, fig3)
figure

figure2 <- multi_panel_figure(columns = 2, rows = 1, width = 250,
                              height = 125)
figure2 <- fill_panel(figure2, fig2)
figure2 <- fill_panel(figure2, fig4)
figure2
