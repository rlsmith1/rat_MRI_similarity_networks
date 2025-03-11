
#######################################################################################

## Demonstrate what KL divergence is in terms of relationship between two distributions

#######################################################################################

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

## LIBRARIES
library(tidyverse)
library(LaplacesDemon)


## FUNCTIONS

# Compute KL divergence
kl_divergence <- function(mu1, sigma1, mu2, sigma2) {
    term1 <- log(sigma2 / sigma1)
    term2 <- (sigma1^2 + (mu1 - mu2)^2) / (2 * sigma2^2)
    return(term1 + term2 - 0.5)
}

# Generate data samples based on assigned parameters and target MIND similarity
f_generate_data <- function(mu, sigma, mind, seed = 123) {
    
    # Calculate target KL divergence based on target MIND
    target_kl <- (1 - mind)/mind
    
    # Optimization function to find mu2, sigma2
    optimize_kl <- function(params) {
        mu2 <- params[1]
        sigma2 <- abs(params[2])  # Ensure sigma2 is positive
        kl <- kl_divergence(mu1, sigma1, mu2, sigma2)
        return((kl - target_kl)^2)  # Minimize squared difference
    }
    
    # Find parameters that achieve the target KL divergence
    result <- optim(c(0, 1.5), optimize_kl, method = "L-BFGS-B")
    
    # Extract optimized parameters
    mu2_opt <- result$par[1]
    sigma2_opt <- abs(result$par[2])
    
    # Generate data samples
    set.seed(seed)
    N <- 5000
    data1 <- rnorm(N, mu1, sigma1)
    data2 <- rnorm(N, mu2_opt, sigma2_opt)
    
    # Combine into tibble
    df <- tibble(
        mind = paste0("MIND ~= ", mind),
        `P(x)` = data1,
        `Q(x)` = data2
    )
    return(df)
    
}

## SET PARAMETERS

# Base normal distribution (N1)
mu1 <- 5
sigma1 <- 2

## GENERATE DISTRIBUTIONS ACROSS MIND VALUES
df_kl_results <- map_dfr(
    .x = c(0.01, 0.5, 0.99),
    .f = ~ f_generate_data(mu = mu1, 
                           sigma = sigma1, 
                           mind = .x,
                           seed = 124
    )
)


## PLOT
df_kl_results %>% 
    pivot_longer(2:ncol(.), names_to = "distribution", values_to = "value") %>% 
    ggplot(aes(x = value)) +
    geom_density((aes(fill = distribution)), alpha = 0.5) +
    facet_wrap(vars(mind), scales = "free") +
    guides(fill = guide_legend(title = NULL)) +
    labs(title = "Kullback-Leibler Divergence") +
    theme_void() +
    theme(legend.position = c(0.945, 0.7),
          strip.text = element_text(size = 11, face = "bold"),
          plot.title = element_text(hjust = 0.5, margin = margin(b = 10))
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/example_KL_divergence", .x), width = 5, height = 1.8)
)

## PLOT P only
df_kl_results %>% 
    pivot_longer(2:ncol(.), names_to = "distribution", values_to = "value") %>% 
    mutate(distribution = factor(distribution, levels = c("Q(x)", "P(x)"))) %>% arrange(distribution) %>% 
    ggplot(aes(x = value)) +
    geom_density((aes(fill = distribution, color = distribution)), alpha = 0.5) +
    facet_wrap(vars(mind), scales = "free") +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = c("P(x)" = "#F8766D", "Q(x)" = "transparent")) +
    scale_color_manual(values = c("P(x)" = "black", "Q(x)" = "transparent")) +
    theme_void() +
    theme(legend.position = "none",
          strip.text = element_blank(),
          plot.title = element_text(hjust = 0.5, margin = margin(b = 10))
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/example_KL_divergence_Ponly", .x), width = 5, height = 1.8)
)



## GENERATE DUMMY NETWORK




