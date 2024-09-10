library(tidyverse)
library(copula)
library(gtable)
library(grid)
library(gridExtra)

# Loading the estimates of Lambda_phi from the simulation
load("./Simulation_Results/Simulation_MO_1000.Rda")
results_Lp <- results %>%
  filter(type %in% c("abs(x)^p"))

results_exp <- results %>%
  filter(type %in% c("exp(c * x) - 1"))

results_exp_abs <- results %>%
  filter(type %in% "exp(abs(c * x)) - 1")

# Loading the true values of Lambda_phi for all simulation scenarios
load("./Simulation_Results/True_Values.Rda")
limits_Lp <- limits %>%
  filter(type %in% c("abs(x)^p"))

limits_exp <- limits %>%
  filter(type %in% c("exp(c * x) - 1"))

limits_exp_abs <- limits %>%
  filter(type %in% "exp(abs(c * x)) - 1")

# Defining and transforming variables for nice Figures and merging true values of Lambda_phi into the frames containing the estimates
results_Lp <- results_Lp %>%
  mutate(pars = paste0("(", alpha, ", ", beta, ")")) %>%
  mutate(p = factor(p_c, levels = c(1, 2, 3), labels = c("p = 1", "p = 2", "p = 3"))) %>%
  mutate(n = factor(n, levels = c(10, 50, 100, 500, 1000, 5000, 10000), labels = c("n = 10", "n = 50", "n = 100", "n = 500", "n = 1.000", "n = 5.000", "n = 10.000")))

limits_Lp <- limits_Lp %>%
  mutate(p = factor(p_c, levels = c(1, 2, 3), labels = c("p = 1", "p = 2", "p = 3")))

results_Lp <- results_Lp %>%
  left_join(limits_Lp)

results_exp <- results_exp %>%
  mutate(pars = paste0("(", alpha, ", ", beta, ")")) %>%
  mutate(c = factor(p_c, levels = c("1/5", "1", "5"), labels = c("c = 1/5", "c = 1", "c = 5"))) %>%
  mutate(n = factor(n, levels = c(10, 50, 100, 500, 1000, 5000, 10000), labels = c("n = 10", "n = 50", "n = 100", "n = 500", "n = 1.000", "n = 5.000", "n = 10.000")))

limits_exp <- limits_exp %>%
  mutate(c = factor(p_c, levels = c("1/5", "1", "5"), labels = c("c = 1/5", "c = 1", "c = 5")))

results_exp <- results_exp %>%
  left_join(limits_exp) %>%
  rename(c_ = c)

results_exp_abs <- results_exp_abs %>%
  mutate(pars = paste0("(", alpha, ", ", beta, ")")) %>%
  mutate(c = factor(p_c, levels = c("1/5", "1", "5"), labels = c("c = 1/5", "c = 1", "c = 5"))) %>%
  mutate(n = factor(n, levels = c(10, 50, 100, 500, 1000, 5000, 10000), labels = c("n = 10", "n = 50", "n = 100", "n = 500", "n = 1.000", "n = 5.000", "n = 10.000")))

limits_exp_abs <- limits_exp_abs %>%
  mutate(c = factor(p_c, levels = c("1/5", "1", "5"), labels = c("c = 1/5", "c = 1", "c = 5")))

results_exp_abs <- results_exp_abs %>%
  left_join(limits_exp_abs) %>%
  rename(c_ = c)

my_pars <- unique(results_Lp %>% select(alpha, beta))

# Looping over all simulation settings and creating and saving the figures
for(i in 1:nrow(my_pars)){
  temp_copula <- moCopula(param = c(my_pars$alpha[i], my_pars$beta[i]))
  set.seed(1)
  temp_sample <- rCopula(n = 5000, copula = temp_copula) %>%
    as.data.frame() %>%
    rename(x = V1, y = V2)
  
  plot_example <- ggplot(temp_sample, aes(x = x, y = y)) +
    geom_point(size = 0.8) +
    ylab(NULL) +
    xlab(NULL) +
    theme_bw(base_size = 10) +
    facet_wrap(~1, labeller = label_bquote(alpha ~ " = " ~ .(my_pars$alpha[i]) ~ ", " ~ beta ~ " = " ~ .(my_pars$beta[i])))
  
  plot_Lp <- results_Lp %>%
    filter(alpha == my_pars$alpha[i], beta == my_pars$beta[i]) %>%
    ggplot(aes(x = n, y = lambda, fill = p)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw(base_size = 10) +
    geom_hline(aes(yintercept = limit, col = p), lty = 2) +
    ylim(0, 1) +
    ylab(NULL) +
    scale_fill_discrete(name = NULL) +
    facet_wrap(~1, labeller = label_bquote(varphi(x)~" = |x|"^p)) +
    scale_color_discrete(name = NULL) +
    theme(axis.title.x = element_blank())
  
  if(my_pars$alpha[i] == 1 & my_pars$beta[i] == 1){
  plot_Lp <- results_Lp %>%
    filter(alpha == my_pars$alpha[i], beta == my_pars$beta[i]) %>%
    ggplot(aes(x = n, y = lambda, fill = p, col = p)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw(base_size = 10) +
    geom_hline(aes(yintercept = limit, col = p), lty = 2) +
    ylim(0, 1) +
    ylab(NULL) +
    scale_fill_discrete(name = NULL) +
    facet_wrap(~1, labeller = label_bquote(varphi(x)~" = |x|"^p)) +
    scale_color_discrete(name = NULL) +
    theme(axis.title.x = element_blank())
  }
  
  plot_exp <- results_exp %>%
    filter(alpha == my_pars$alpha[i], beta == my_pars$beta[i]) %>%
    ggplot(aes(x = n, y = lambda, fill = c_)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw(base_size = 10) +
    geom_hline(aes(yintercept = limit, col = c_), lty = 2) +
    ylim(0, 1) +
    ylab(NULL) +
    scale_fill_discrete(name = NULL) +
    facet_wrap(~1, labeller = label_bquote(varphi(x) ~ " = " ~ e^{cx} - 1)) +
    scale_color_discrete(name = NULL) +
    theme(axis.title.x = element_blank()) +
    guides(color = "none", fill = "none")
  
  if(my_pars$alpha[i] == 1 & my_pars$beta[i] == 1){
  plot_exp <- results_exp %>%
    filter(alpha == my_pars$alpha[i], beta == my_pars$beta[i]) %>%
    ggplot(aes(x = n, y = lambda, fill = c_, col = c_)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw(base_size = 10) +
    geom_hline(aes(yintercept = limit, col = c_), lty = 2) +
    ylim(0, 1) +
    ylab(NULL) +
    scale_fill_discrete(name = NULL) +
    facet_wrap(~1, labeller = label_bquote(varphi(x) ~ " = " ~ e^{cx} - 1)) +
    scale_color_discrete(name = NULL) +
    theme(axis.title.x = element_blank()) +
    guides(color = "none", fill = "none")
  }
  
  plot_exp_abs <- results_exp_abs %>%
    filter(alpha == my_pars$alpha[i], beta == my_pars$beta[i]) %>%
    ggplot(aes(x = n, y = lambda, fill = c_)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw(base_size = 10) +
    geom_hline(aes(yintercept = limit, col = c_), lty = 2) +
    ylim(0, 1) +
    ylab(NULL) +
    scale_fill_discrete(name = NULL) +
    facet_wrap(~1, labeller = label_bquote(varphi(x) ~ " = " ~ e^{~"|"~cx~"|"} - 1)) +
    scale_color_discrete(name = NULL) +
    theme(axis.title.x = element_blank())
  
  if(my_pars$alpha[i] == 1 & my_pars$beta[i] == 1){
  plot_exp_abs <- results_exp_abs %>%
    filter(alpha == my_pars$alpha[i], beta == my_pars$beta[i]) %>%
    ggplot(aes(x = n, y = lambda, fill = c_, col = c_)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_bw(base_size = 10) +
    geom_hline(aes(yintercept = limit, col = c_), lty = 2) +
    ylim(0, 1) +
    ylab(NULL) +
    scale_fill_discrete(name = NULL) +
    facet_wrap(~1, labeller = label_bquote(varphi(x) ~ " = " ~ e^{~"|"~cx~"|"} - 1)) +
    scale_color_discrete(name = NULL) +
    theme(axis.title.x = element_blank())
  }
  
  
  temp_plot <- arrangeGrob(grobs = list(plot_example, plot_Lp, plot_exp, plot_exp_abs), widths = rep(1/9, times = 9), layout_matrix = rbind(c(1, 1, 1, 2, 2, 2, 2, 2, 2), c(3, 3, 3, 3, 4, 4, 4, 4, 4)))
  ggsave(plot = temp_plot, file = paste0("./Figures/MO_", my_pars$alpha[i], "_", my_pars$beta[i] ,".pdf"))
}