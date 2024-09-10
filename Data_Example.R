library(tidyverse)
source("./Core_Functions.R")

# Load Data
my_data <- read.csv("./Data_Example.csv")
names(my_data) <- c("ID", "age", "ethnicity", "diabetes", "diastole", "systole", "fat", "gestational_age", "pregnancies",
                    "fasting_glucose", "bmi", "gestational_age_at_birth", "delivery", "child_weight", "gestational_dm")

# Only use complete observations
my_data <- my_data %>%
  na.omit()

# Create Plot of Checkerboard Copula and conditional ecdfs with BMI as explanatory variable
ecbc <- ECBC(X = my_data$bmi, Y = my_data$child_weight, resolution = 9)

plot_data <- my_data %>%
  select(bmi, child_weight) %>%
  mutate(across(c(everything()), .fns = function(x) rank(x)/length(x)))

bmi_data_plot <- plot_density(ecbc) +
  theme_bw(base_size = 10) + 
  theme(axis.title.x = element_blank()) +
  facet_wrap(~1, labeller = label_bquote("ECBC for Birth Weight and BMI of the Mother")) +
  ylab(NULL) +
  ggtitle(NULL) +
  geom_point(data = plot_data, aes(x = bmi, y = child_weight), size = 0.5)

kernel_mat_bmi <- t(apply(ecbc, 1, FUN = density_to_kernel))
kernel_mat_bmi <- cbind(rep(0, times = 9), kernel_mat_bmi)
kernel_mat_bmi <- cbind(kernel_mat_bmi, rep(1, times = 9))

kernel_frame_bmi <- as.data.frame(kernel_mat_bmi) %>%
  mutate(x = 1:9) %>%
  gather(key = "key", value = "value", -x) %>%
  mutate(key = as.numeric(str_remove(key, "V"))) %>%
  mutate(key = case_when(key == 1 ~ 0,
                         key == 11 ~ 1,
                         TRUE ~ 1/18 + (key - 2)/9))

bmi_ecdf_plot <- ggplot(kernel_frame_bmi, aes(x = key, y = value, col = factor(x))) +
  geom_line() +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank()) +
  facet_wrap(~1, labeller = label_bquote("ECDFs for Birth Weight by BMI of the Mother")) +
  ylab(NULL) +
  scale_color_discrete(name = "BMI")

bmi_data_plot <- ggplotGrob(bmi_data_plot)
bmi_ecdf_plot <- ggplotGrob(bmi_ecdf_plot)

grid.arrange(bmi_data_plot, bmi_ecdf_plot, nrow = 1)
my_plot <- arrangeGrob(bmi_data_plot, bmi_ecdf_plot, nrow = 1)
ggsave(plot = my_plot, file = "./Figures/Data_example.pdf", height = 80, width = 200, units = "mm")

# Create Table with Age, Systole, Diastole, Fat, BMI and Fasting Glucose as explanatory variables and
# abs(x), exp(abs(x)) - 1 and exp(x) - 1 as choices for phi
continuous_vars <- c("age", "systole", "diastole", "fat", "bmi", "fasting_glucose")
phi_vec <- c("function(x) abs(x)", "function(x) exp(abs(x)) - 1", "function(x) exp(x) - 1")

argument_grid <- expand.grid(predictor = continuous_vars, phi_char = phi_vec) %>%
  mutate(phi_char = as.character(phi_char)) %>%
  mutate(predictor = as.character(predictor))

calculate_lambda <- function(predictor, phi_char){
  ecbc <- ECBC(X = unlist(my_data[, predictor]), Y = my_data$child_weight)
  phi <- eval(parse(text = phi_char))
  temp_lambda <- Lambda_phi(ecbc, phi = phi)
  return(data.frame(predictor = predictor, phi = phi_char, lambda = temp_lambda))
}

results <- pmap_dfr(.l = argument_grid, .f = calculate_lambda)

results <- results %>%
  spread(key = phi, value = lambda) %>%
  arrange(desc(`function(x) abs(x)`))

write.csv2(results, "./Results_Example.csv", row.names = F)