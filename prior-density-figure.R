# Make figure displaying each of the prior distributions
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

# First do each of the 9 standardized priors on the log(OR) scale
# Optimistic priors
# Strong optimistic prior
temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_sn), 0,
                                      sn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn2 <- p$data[[1]][, c(1, 2)]
plot_data_sn2$Distribution <- "Prior"