######################################################
####                                              ####
#### Empirical Bayesian spatial models using mgcv ####
####                                              ####
######################################################

## This tutorial gives a brief introduction to empirical Bayesian model-fitting using mgcv. 
# We will apply spatio-temporal models to a binary dengue outbreak indicator created
# using data from 2010 - 2020 from Rio de Janeiro state. Model estimation will be 
# carried out using INLA (a fully Bayesian approach) and mgcv (an EB allternative).

# We will compare the predicted probabilities and random effect estimates produced
# using these methods. The models' predictive abilities will also be compared using
# the Brier score and a receiver operator curve (ROC) analysis.


## Why empirical Bayes?
# Unlike fully Bayesian methods (e.g. INLA/MCMC), EB estimates priors from the data 
# rather than specifying them in advance. This can provide more flexibility, particularly
# in spatial models where the underlying structure may not be known.


#### Load packages and data from R script ####
source("")


#### Explore data ####
## Combine case data with the shapefile to plot the data
df_shp <- full_join(df, shp, 
                    by = "municip_code_ibge") %>% 
                      st_as_sf()


# Plot time series of dengue incidence rate (cases per 100,000 residents) per municipality 
# There is a wide range of values, with some explosive outbreaks, particularly 2013, 2015 and 2016.
ggplot(data = df) +
  geom_line(aes(x = year, y = DIR, colour = municip_code_ibge, 
                group = municip_code_ibge)) +
  scale_colour_viridis_c(name = "Municipality") +
  scale_x_continuous(name = "Year", breaks = seq(2010, 2020, by = 2)) +
  theme_bw()


# Let's look at DIR on a map to identify potential spatial patterns in the data
# Seem to be spatial clusters of high DIR in the northeast of the state and lower rates
# in the centre of the state. This spatial pattern should be incorporated into our model.
ggplot(data = df_shp) +
  geom_sf(aes(fill = DIR), lwd = .05) +
  # Produce a map per year
  facet_wrap(~ year) +
  # Transform onto a log scale to make patterns easier to identify
  scale_fill_viridis_c(name = "DIR", trans = "log1p",
                       breaks = c(0, 100, 300, 1000, 5000),
                       direction = 1) +
  # Remove coordinate grid
  theme_void()


## Let's look at the socioeconomic variables to build hypotheses about potential 
# driver of these outbreaks. First, let's create a df with just these variables:
df_socio <- df %>% 
  # Select just a single year as these are stationary across the time period
  filter(year == 2010) %>% 
  dplyr::select(municip_code_ibge, urban10:level18_num) %>% 
  # Convert level of influence to a factor with labels for each level
  mutate(regic18 = factor(level18_num, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre"))) %>% 
  # Join to shapefile to plot
  full_join(., shp, by = "municip_code_ibge") %>% 
  st_as_sf()


# % urbanisation
urban_map <- ggplot(data = df_socio) +
  geom_sf(aes(fill = urban10)) +
  scale_fill_distiller(name = "% urban", palette = "BuPu", direction = 1) +
  theme_void()


# % piped water
water_map <- ggplot(data = df_socio) +
  geom_sf(aes(fill = water_network)) +
  scale_fill_distiller(name = "% piped water", palette = "BuPu", direction = 1) +
  theme_void()


# % refuse collected
refuse_map <- ggplot(data = df_socio) +
  geom_sf(aes(fill = total_collected)) +
  scale_fill_distiller(name = "% refuse collect", palette = "BuPu", direction = 1) +
  theme_void()


# Levels of influence
connect_cols <- c("#355070", "#6d597a", "#b56576", "#e56b6f", "#eaac8b")
                  
regic_map <- ggplot(data = df_socio) +
  geom_sf(aes(fill = regic18)) +
  scale_fill_manual(name = "Level of \ninfluence", values = connect_cols) +
  theme_void()


## In Rio these socioeconomic variables are highly correlated. To avoid multicollinearity
# we will just include urbanisation as a measure of access to basic services and influence
socio_maps <- plot_grid(urban_map, water_map, refuse_map, regic_map, nrow = 2)
  


#### Fit an empirical Bayesian model (using mgcv) ####
# Include urbanisation and year as linear covariates, and a 2-d 
# thin plate spline applied to coordinates of the centroid of municipalities to 
# produce a spatial smooth surface
gam_model <- gam(outbreak ~ urban10 + fyear +
                   # spatial plane (using thin-plate spline)
                   s(lon, lat, bs = "tp"),
                 family = "binomial", data = st_drop_geometry(df_shp), 
                 # fit using restricted maximum likelihood
                 method = "REML")

summary(gam_model)


# Generic plot.gam function does not look good, use mgcViz package
plot(gam_model)

# Create mgcViz object
gam_viz <- getViz(gam_model)

# Plot the smoothed functions, apply exp transformation to show on the aOR scale
plot(sm(gam_viz, 1), tran = exp)


####  Fit a fully Bayesian model (in INLA) ####
## The most common spatial modelling approach in INLA assumes a Gaussian Markov random field structure
# (This means that regions are only expected to be correlated with neighbours)

# First, we must create a neighbourhood object 
# (any municipalities sharing a border are considered neighbours)
nb.rj <- poly2nb(shp, row.names = shp$municip_index)
# Attach index to nb object to allow matching between this and the df
names(nb.rj) <- attr(nb.rj, "region.id")


# Plot the neighbourhood object (lines indicate municipalities are connected)
nb_lines <- as(nb2lines(nb.rj, coords = st_centroid(shp)$geom), 'sf')

ggplot() +
  geom_sf(data = shp) +
  geom_sf(data = nb_lines)+
  theme_void()


# Convert into an INLA nb object
nb_inla <- nb2INLA("rj.graph", nb.rj)


# As INLA is a fully Bayesian approach, we can specify priors before model fitting
# We will use penalised complexity (PC) priors as these allow us to set meaningful, 
# interpretable priors

# Prior for phi (proportion of variance explained by the spatial structure)
phi.prior = list(prior = "pc", 
                 param = c(0.5, 2/3)) # P(phi < 0.5) = 2/3 

# Prior for random effect precision
prec.prior = list(prior = "pc.prec",
                  param = c(1, 0.01)) # P(st. dev > 1) = 0.01


# We are assuming a BYM2 spatial structure in our data, with the structured element 
# defined by our neighbourhood matrix. 
inla_model <- inla(outbreak ~ urban10 + fyear +
                     # spatial random effect
                     f(municip_index, model = "bym2", 
                       graph = "rj.graph",
                       # prior specification
                       hyper = list(prec = prec.prior,
                                    phi = phi.prior)),
                   data = st_drop_geometry(df_shp), family = "binomial", 
                   # as we will be sampling from the posterior, include the following
                   # to ensure these results are returned
                   control.compute = list(config = T))


summary(inla_model)







#### Compare estimated probabilities ####
# We can obtain predicted probabilities from the INLA model object, however for the
# GAM, we have to explicitly sample from the posterior distributions to obtain the
# predicted response.

## Step 1: Estimate the regression coefficients (betas):
# Function to obtain MVN random variates
rmvn <- function(n, mu, sig) { 
  L <- mroot(sig); m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}


# Specify number of simulations
n.sims <- 1000

# Simulate from beta posterior distributions
# Mean = coefficient estimates
# Variance = covariance matrix of the predictors (model$Vp)
betas_gam <- rmvn(n.sims, coef(gam_model), gam_model$Vp) 


## Step 2: Extract linear predictors from the model
model_mat <- predict(object = gam_model, 
                     # Feed through data used to fit the model
                     newdata = gam_model$model, 
                     # Return a matrix of linear predictors
                     type = "lpmatrix")

## Step 3: Multiple the simulated betas with linear predictors to produce an estimate 
# per municipality per year (per simulation)
# Apply probit transformation to convert this onto the probability scale:
# g(prob) = log(prob / (1 - prob)) = beta0 + beta1 + f(year) + f(space)
# prob = exp[g(prob)] / 1 + exp[g(prob)]
mean_prob_gam <- exp(model_mat %*% t(betas_gam)) /
  (1 + exp(model_mat %*% t(betas_gam)))


# Step 4: Take the mean of the 1000 simulations to return the predicted probability
pred_prob_gam <- apply(mean_prob_gam, 1, mean)


#### Compare predicted probabilities ####
# Combine predictions from INLA and mgcv into the same df
prob_pred <- data.table(municip_code_ibge = df$municip_code_ibge,
                        year = df$year,
                        # Keep observed outbreak variable
                        outbreak_obs = df$outbreak,
                        # Return predicted probability from the INLA model
                        pred_inla = inla_model$summary.fitted.values$mean,
                        # Predicted probability of an outbreak from the GAM model
                        pred_gam = pred_prob_gam)


# Plot predicted probabilities for INLA and GAM against each other
# Both models produce similar estimates, randomly scattered around the y = x line
ggplot(data = prob_pred) +
  geom_point(aes(x = pred_gam, y = pred_inla)) +
  # Add reference line (if models predictions were equal)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  theme_bw()




#### Compare spatial random effects ####
# As with predicted probabilities, we can obtain spatial random effects from the INLA
# model object, however mgcv requires simulation from the posterior distribution. 

## INLA
# Extract mean of posterior for spatial random effects (the first 92 rows relate to
# the spatially structured random effect, whereas the next 92 are a combination of 
# the structured and unstructured with the phi parameter)
inla_spat <- inla_model$summary.random$municip_index[(nrow(shp)+1):(nrow(shp)*2), 1:2] %>% 
  # Transform onto odds scale for easier interpretation
  transmute(municip_index = ID - 92, 
            spat_re = exp(mean)) 


## mgcv
# Follow the same process as before but just selecting beta and linear predictors 
# asociated with the spatial smooth
betas_gam_spat <- betas_gam[ ,substr(names(coef(gam_model)), 1, 5) == "s(lon"]
model_mat_spat <- model_mat[ ,substr(names(coef(gam_model)), 1, 5) == "s(lon"]


# Estimate spatial random effects + transform from log scale
spat_est_gam <- exp(model_mat_spat %*% t(betas_gam_spat))

# As all years are equal, just keep one of them (2010)
spat_est_gam <- spat_est_gam[gam_model$model$fyear == 2010, ] 

# Take average of simulations
mean_spat_gam <- apply(spat_est_gam, 1, mean)


# Join random effects to dataset, include  INLA estimates
spat_re <- data.table(municip_code_ibge = df[df$year == 2010,]$municip_code_ibge,
                      spat_re_gam = mean_spat_gam,
                      spat_re_inla = inla_spat$spat_re) %>%
  left_join(., shp, by = "municip_code_ibge") %>%
  st_as_sf()


## Plot random effects on a map to compare
# Return max random effect to set legend to same limits
max_re <- max(c(spat_re$spat_re_gam, spat_re$spat_re_inla))

# INLA spatial effects
inla_spat_map <- ggplot(data = spat_re) +
  geom_sf(aes(fill = spat_re_inla), lwd = .05) +
  scale_fill_viridis_c(name = "Spatial random \neffect") +
  expand_limits(fill = c(0, max_re)) +
  theme_void()

# GAM spatial effects
gam_spat_map <- ggplot(data = spat_re) +
  geom_sf(aes(fill = spat_re_gam), lwd = .05) +
  scale_fill_viridis_c(name = "Spatial random \neffect") +
  expand_limits(fill = c(0, max_re)) +
  theme_void()


# Obtain legend (to plot maps on the same grid)
spat_re_leg <- get_legend(inla_spat_map)

# Combine map onto the same grid
spat_re_maps <- plot_grid(inla_spat_map + theme(legend.position = "none"),
                          gam_spat_map +  theme(legend.position = "none"),
                          labels = c("INLA", "GAM"))

# Add the legend
spat_re_maps <- plot_grid(spat_re_maps, spat_re_leg, rel_widths = c(4, .3))

# There appears to be a similar pattern in the random effects, however the GAM
# is smoother 


## Compare the spatial random effect estimates using a scatterplot
# Generally similar results, however INLA estimates tend to be higher for larger
# random effects - potentially a result of the smoothness of the GAM model
ggplot(data = spat_re) +
  geom_point(aes(x = spat_re_gam, y = spat_re_inla)) +
  # Add a reference line showing equality
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "GAM", y = "INLA") +
  theme_bw()



#### Compare Brier score ####
## Calculate the Brier score (mean squared difference) for each model.
# The smaller the score, the better the model is doing of predicting the outcome.
brier_scores <- prob_pred %>% 
  mutate(sqdiff_inla = (pred_inla - outbreak_obs)^2,
         sqdiff_gam = (pred_gam - outbreak_obs)^2) %>% 
  group_by(municip_code_ibge) %>% 
  summarise(brier_inla = sum(sqdiff_inla),
            brier_gam = sum(sqdiff_gam)) %>% 
  ungroup() %>% 
  mutate(brier_diff = brier_inla - brier_gam) %>% 
  full_join(., shp, by = "municip_code_ibge") %>% 
  st_as_sf()


# Plot these Brier scores on a map
max_brier <- max(c(brier_scores$brier_inla, brier_scores$brier_gam))

brier_inla <- ggplot(data = brier_scores) +
  geom_sf(aes(fill = brier_inla), lwd = .05) +
  scale_fill_viridis_c(name = "Brier score") +
  expand_limits(fill = c(0, max_brier)) +
  theme_void()
  
brier_gam <- ggplot(data = brier_scores) +
  geom_sf(aes(fill = brier_gam), lwd = .05) +
  scale_fill_viridis_c(name = "Brier score") +
  expand_limits(fill = c(0, max_brier)) +
  theme_void() 

#  Plot maps on the same grid for easier comparison
brier_comp <- plot_grid(brier_inla + theme(legend.position = "none"), 
                        brier_gam + theme(legend.position = "bottom"),
                        labels = c("INLA", "GAM"))


# Plot the difference in Brier scores per municipality
# In general, INLA performs slightly better
brier_diff <- ggplot(data = brier_scores) +
  geom_sf(aes(fill = brier_diff), lwd = .05) +
  scale_fill_viridis_c(name = "Brier score") +
  theme_void() 

#### Compare ROC curves ####
roc_obj_inla <- roc(prob_pred$outbreak_obs, 
                    prob_pred$pred_inla,
                    auc = T, ci = T, plot = T)

roc_obj_gam <- roc(prob_pred$outbreak_obs, 
                   prob_pred$pred_gam,
                   auc = T, ci = T, plot = T)


roc_obj <- data.table(roc_inla_sens = roc_obj_inla$sensitivities,
                      roc_inla_spec = 1 - roc_obj_inla$specificities,
                      roc_gam_sens = roc_obj_gam$sensitivities,
                      roc_gam_spec = 1 - roc_obj_gam$specificities)


roc_curve <- ggplot(data = roc_obj) +
  # Add a line per threshold
  geom_line(aes(x = roc_inla_spec,  y = roc_inla_sens)) +
  geom_line(aes(x = roc_gam_spec,  y = roc_gam_sens), 
            linetype = "dashed", col = "red") +
  # Add reference line (y = x line represents chance)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "True negative rate", y = "True positive rate") +
  theme_light()

