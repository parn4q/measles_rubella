##############################################################################
# This script contains the mapping information to use for analysis in Julia. #
# We will collect the country, latitude and longitutde of many countries to  #
# join to the other dataset.  We will conduct EDA here and use Julia for     #
# modeling.                                                                  #
##############################################################################


###
### Load in libraries and data
###

setwd('E:/TidyTuesday') # for windows
setwd('/Volumes/Lexar/TidyTuesday/Measles') # for mac

library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)
library(brms)
library(cmdstanr)

sf_use_s2(FALSE)


world = ne_countries(returnclass = "sf") |> select(sov_a3, label_x, label_y, geometry)

case_year = read.csv("cases_year.csv")

case_month = read.csv("cases_month.csv")

df = world |> left_join(case_year, join_by(sov_a3 == iso3))

df$year = as.factor(df$year)
df$region = as.factor(df$region)
df$country = as.factor(df$country)

df$log_total_population = log(df$total_population)



###
### EDA
###

na_counts_dplyr <- df |>
  summarise_all(~ sum(is.na(.)))

df = df |> na.omit()


# Plotting total population.  Nothing really comes out of it
df |> select(label_x, label_y, geometry, year, total_population) |>
  filter(year %in% c(2020:2025)) |>
  group_by_all() |>
  summarise(mean_population = mean(total_population)) |>
  ggplot() +
  geom_sf(aes(fill = mean_population), color = "gray30", size = 0.1) +
  scale_fill_viridis_c(option = "plasma", na.value = "lightgray", trans = "log10",
                       labels = scales::rescale()) +
  labs(title = "World Population by Country", fill = "Population") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(vars(year))


df |> select(year, region, country) |>
  filter(year == 2025 & region == 'AFRO') 

##########################################################################
# Countries that are missing may be truly missing from the case datasets #
# For example, Egypt appears in some years but not all                   #
##########################################################################


# Suspected Measles Rubella Cases by Country. 
df |> select(label_x, label_y, geometry, year, total_population, measles_total) |>
  filter(year %in% c(2020:2025)) |>
  group_by_all() |>
  summarise(mean_suspect = mean(measles_total)) |>
  ggplot() +
  geom_sf(aes(fill = mean_suspect), 
          color = "gray30", 
          size = 0.1) +
  scale_fill_viridis_c(option = "plasma", 
                       na.value = "lightgray") +
  labs(title = "Total Measles Cases by Country", 
       fill = "Measles Cases") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  facet_wrap(vars(year))




# Suspected Measles Rubella Cases by Country. 
df |> select(label_x, label_y, geometry, year, total_population, rubella_total) |>
  #filter(year %in% c(2018:2025)) |>
  group_by_all() |>
  summarise(mean_suspect = mean(rubella_total)) |>
  ggplot() +
  geom_sf(aes(fill = mean_suspect), color = "gray30", size = 0.1) +
  scale_fill_viridis_c(option = "plasma", na.value = "lightgray") +
  labs(title = "Total Rubella Cases by Country", fill = "Rubella Cases") +
  theme_minimal() +
  facet_wrap(vars(year))

df |> select(country, year, measles_total, rubella_total) |>
  group_by_all() |>
  summarise(mean_measles = mean(measles_total),
            mean_rubella = mean(rubella_total)) |>
  arrange(desc(measles_total))

df |> select(country, year, measles_total, rubella_total) |>
  group_by_all() |>
  summarise(mean_measles = mean(measles_total),
            mean_rubella = mean(rubella_total)) |>
  arrange(desc(rubella_total))

# correlation matrix

GGally::ggcorr(case_year |> select(where(is.numeric)),
               label = T, 
               size = 3,
               hjust = 0.95)


##############################################################################
# for both measles and Rubella total cases, the majority of the cases
# are Africa and Eurasia with some being in South America (primarily brazil) 
# Madagascar has the largest amount of measles in 2019, and 
# Poland had the largest about of rubella cases in 2013.
##############################################################################

##############################################################################################################
# some hypothesis:
# Could we predict the total number of measle AND rubella cases simulatneously (bivariate model?) for 2025
# or do it separately?
# Should we remove some countries that are not heavy on the either cases? 
# 
##############################################################################################################


# As population increases, does the number of measles and rubella increase?


cor(df$total_population, df$measles_total) # 46%
cor(df$total_population, df$rubella_total) # 50%

ggplot(data = df, aes(x = measles_total)) + 
  geom_histogram() +
  labs(title = "Frequency of the Total Number of Measles",
       x = NULL, 
       y = NULL
       ) +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(data = df, aes(total_population)) + 
  geom_histogram() + 
  scale_x_log10() +
  labs(title = "Frequency of the Log Total Population",
       x = NULL,
       y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = df, aes(x = measles_epi_linked)) + 
  geom_histogram() +
  labs(title = "Frequency of the Measles Epi Linked",
       x = NULL, 
       y = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5))



###
### Bayesian Model
###



prior = c(
          prior(normal(0, 1), class = 'b', coef = "log_total_population"),
          prior(normal(0, 10), class = "Intercept"),
          prior(normal(0,10), class = "Intercept", dpar = "zi"),
          prior(normal(0,10), class = "b", dpar = "zi")
)


# Poisson Model
## 1: Zero inflated but random effects on year and region
model = brm(bf(measles_total ~ log_total_population + (1| year/region/country),
               zi ~ 1 + log_total_population),
            data = df,
            family = zero_inflated_poisson(link = "log"),
            prior = prior,
            chains = 2,
            iter = 5000,
            warmup = 3000,
            cores = 2,
            backend = "cmdstanr",
            save_model = "zi_model.stan",
            control = list(adapt_delta = 0.99,
                           max_treedepth = 13)
)


summary(model)


plot(model)
pairs(model)

pp_check(model)

loo(model, model2)

pred = predict(model)

pred = as.data.frame(pred) |> mutate(measles_total = df$measles_total)

#write.csv(pred, "./poisson_pred.csv")

sqrt(mean((pred$measles_total - pred$Estimate)^2))

ggplot(pred) + geom_point(mapping = aes(Estimate, measles_total))


########################################################################################################################
# rmse for measles_total ~ log_total_population + (1|year/region)) (zero inflated on the intercept) is 6231.608
# prior = c(prior(normal(0,10), class = 'b'), prior(normal(0, 10), class = "Intercept"),prior(inv_gamma(5,1), class = 'sd', group = "year"))


#rmse for measles_total ~ log_total_population + (1|year) + (1|region/country) (zero inflated on the intercept) is ~6032.774. takes 80 minutes
# prior = c(prior(normal(0,10), class = 'b'),
#prior(normal(0, 10), class = "Intercept"),
#prior(inv_gamma(5,1), class = 'sd', group = "year"),
##prior(normal(0,10), class = "b", dpar = "zi"),
#prior(normal(0,10), class = "Intercept", dpar = "zi"))




# rmse for measles_total ~ log_total_population + region + (1|year) + (1|country), zi ~ 1) is ~6819.345.  

# rmse for measles_total ~ log_total_population + (1| year/region/country), zi ~ 1 + log_total_population) is ~330.8287



############################################################################################################################################

prior = c(
  prior(normal(0, 1), class = 'b', coef = "log_total_population"),
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0,10), class = "Intercept", dpar = "zi"),
  prior(normal(0,10), class = "b", dpar = "zi"),
  prior(gamma(2, 0.1), class = "shape"),
  prior(gamma(2, 0.1), class = "sd", group = "year")
)



# Negative Binomial Model
model = brm(bf(measles_total ~ log_total_population + (1| year/region/country),
               zi ~ 1 + log_total_population),
            data = df,
            family = zero_inflated_negbinomial(link = "log"),
            prior = prior,
            chains = 2,
            iter = 6000,
            warmup = 4000,
            init = 0,
            cores = 2,
            backend = "cmdstanr",
            save_model = "zi_model.stan",
            control = list(adapt_delta = 0.99,
                           max_treedepth = 14)
)


summary(model)
plot(model)

pairs(model)

pred = predict(model)

pred = as.data.frame(pred) |> mutate(measles_total = df$measles_total)

sqrt(mean((pred$measles_total - pred$Estimate)^2))




###
### bivariate Poisson
###

# Save the Stan file
writeLines(readLines("bivariate_poisson.stan"), "bivariate_poisson.stan")

# Compile model
mod <- cmdstan_model("bivariate_poisson.stan")

stan_data <- list(
  N = nrow(df),
  measles = df$measles_total,
  rubella = df$rubella_total,
  log_total_population = df$log_total_population
)


# Fit model
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 2,
  parallel_chains = 2,
  iter_sampling = 1000,
  iter_warmup = 1000
)

# Summarize
fit$summary()


fit$cmdstan_summary()

fit$diagnostic_summary()

# Extract draws
draws <- fit$draws()

# Convert to posterior::draws_array
draws_array <- as_draws_array(draws)

# Trace plot for key parameters
bayesplot::mcmc_trace(draws_array)


# Convert to data frame for sampling
draw_df <- as_draws_df(draws)

# Choose 500 posterior samples
n_samps <- 500
draw_subset <- draw_df %>% slice_sample(n = n_samps)

# Prepare prediction matrix
X <- df$log_total_population
n_obs <- length(X)
pred_y1 <- matrix(NA, nrow = n_samps, ncol = n_obs)
pred_y2 <- matrix(NA, nrow = n_samps, ncol = n_obs)

for (s in 1:n_samps) {
  for (i in 1:n_obs) {
    lambda1 <- exp(draw_subset$alpha1[s] + draw_subset$beta1[s] * X[i])
    lambda2 <- exp(draw_subset$alpha2[s] + draw_subset$beta2[s] * X[i])
    lambda3 <- exp(draw_subset$alpha3[s] + draw_subset$beta3[s] * X[i])
    
    # Simulate latent components
    a <- rpois(1, lambda1)
    b <- rpois(1, lambda2)
    c <- rpois(1, lambda3)
    
    pred_y1[s, i] <- a + c
    pred_y2[s, i] <- b + c
  }
}

# Mean prediction across posterior draws
pred_measles_mean <- apply(pred_y1, 2, mean)
pred_rubella_mean <- apply(pred_y2, 2, mean)

# 95% intervals
pred_measles_ci <- apply(pred_y1, 2, quantile, probs = c(0.025, 0.975))
pred_rubella_ci <- apply(pred_y2, 2, quantile, probs = c(0.025, 0.975))

# histograms
hist(pred_measles_mean, col = "lightblue", main = "Posterior Predicted Measles")
abline(v = mean(df$measles_total), col = "red")

# scatterplot
plot(df$measles_total, pred_measles_mean,
     xlab = "Observed Measles", ylab = "Predicted Measles")
abline(0, 1, col = "red")


draw_df <- as_draws_df(fit$draws())
X_mean <- mean(df$log_total_population)

corr_post <- with(draw_df, {
  lam1 <- exp(alpha1 + beta1 * X_mean)
  lam2 <- exp(alpha2 + beta2 * X_mean)
  lam3 <- exp(alpha3 + beta3 * X_mean)
  lam3 / sqrt((lam1 + lam3) * (lam2 + lam3))
})

hist(corr_post, main = "Posterior Correlation: Measles vs Rubella")



log_lik <- numeric(length = nrow(df))
for (n in 1:nrow(df)) {
  lam1 <- exp(draw_df$alpha1 + draw_df$beta1 * df$log_total_population[n])
  lam2 <- exp(draw_df$alpha2 + draw_df$beta2 * df$log_total_population[n])
  lam3 <- exp(draw_df$alpha3 + draw_df$beta3 * df$log_total_population[n])
  
  y1 <- df$measles_total[n]
  y2 <- df$rubella_total[n]
  log_p <- sapply(1:nrow(draw_df), function(s) {
    kmax <- min(y1, y2)
    terms <- sapply(0:kmax, function(k) {
      dpois(y1 - k, lam1[s], log = TRUE) +
        dpois(y2 - k, lam2[s], log = TRUE) +
        dpois(k, lam3[s], log = TRUE)
    })
    log_sum_exp(terms) - (lam1[s] + lam2[s] + lam3[s])
  })
  log_lik[n] <- mean(log_p)
}
mean(log_lik)


















