#### Mixed model approach to SSF - from Muff et al 2019 ####

# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# install.packages("BMA")
# install.packages("ecospat")
# update.packages()

library(amt)
library(lubridate)
library(raster)
library(survival)
library(tidyverse)
library(car)
library(plyr)
library(dplyr)
library(maptools)
library(jtools)
library(BMA)
library(RColorBrewer)
#library(sen2r)

library(survival)
library(TwoStepCLogit)
library(INLA)
library(glmmTMB)
library(ecospat)
library(tictoc)
library(beepr)

options(timeout = 120)

# From 'RSF and SSF analysis of fisher' - Muff et al 2019 -----------------

T05 <- read_csv("data/CSV input data - dd_speed_6/T45505_dd_speed_6.csv") %>% filter(...1 != 1448)
T06 <- read_csv("data/CSV input data - dd_speed_6/T45506_dd_speed_6.csv") %>% mutate(id = 45506)
T07 <- read_csv("data/CSV input data - dd_speed_6/T45507_dd_speed_6.csv") %>% mutate(id = 45507)
T08 <- read_csv("data/CSV input data - dd_speed_6/T45508_dd_speed_6.csv") %>% mutate(id = 45508)
T09 <- read_csv("data/CSV input data - dd_speed_6/T45509_dd_speed_6.csv") %>% mutate(id = 45509)
T10 <- read_csv("data/CSV input data - dd_speed_6/T45510_dd_speed_6.csv") %>% mutate(id = 45510)
T11 <- read_csv("data/CSV input data - dd_speed_6/T45511_dd_speed_6.csv") %>% mutate(id = 45511)
T12 <- read_csv("data/CSV input data - dd_speed_6/T45512_dd_speed_6.csv") %>% mutate(id = 45512) %>% filter(...1 != 186)
T13 <- read_csv("data/CSV input data - dd_speed_6/T45513_dd_speed_6.csv") %>% mutate(id = 45513)
T14 <- read_csv("data/CSV input data - dd_speed_6/T45514_dd_speed_6.csv") %>% mutate(id = 45514)

T_all <- rbind(T06, T07, T08, T09, T10, T11, T12, T13, T14)

dat_all <- T_all %>% nest(-id)
dat_all$sex <- c("m", "m", "f", "f", "f", "f", "m", "m", "f")

T05t <- T05 %>% make_track(lon, lat, DateTime, crs = 4326) %>% transform_coords(crs_to = 2193)
plot(T05t$x_, T05t$y_)
summarize_sampling_rate(T05t)
T05t <- T05t %>% track_resample(rate = hours(2), tolerance = minutes(10)) # %>% 
  #filter_min_n_burst(min_n = 3) %>% steps_by_burst()

dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, lon, lat, DateTime, crs = 4326) %>% transform_coords(crs_to = 2193)
  }))

dat1 <- dat_all %>% mutate(dat_clean = map(trk, ~ {
  .x %>% track_resample(rate = hours(3), tolerance = minutes(10))
}))

# for boyce indices -------------------------------------------------------

# combining dfs into a list and mapping over a function to create projected tracks

all_dfs <- list(T05, T06, T07, T08, T09, T10, T11, T12, T13, T14)

track_function <- function(df) {
df %>% make_track(lon, lat, DateTime, crs = 4326) %>% transform_coords(crs_to = 2193)}

all_tracks <- map(all_dfs, track_function)


# Rasters -----------------------------------------------------------------

# using the INLA model requires using dummy variables (each habitat category is its own variable, with 1s and 0s, rather than a single variable for all habitat categories with e.g. 1-10)

DCChab <- raster("mapping/DCC BCNadj and LCDB.tif")
DCChab[DCChab>31] <- 0
plot(DCChab, col = terrain.colors(n = 8, rev = T))
names(DCChab) <- "DCChab"

shrubland <- DCChab %in% c(15, 16, 18)
names(shrubland) <- "shrubland"
plot(shrubland)
(cellStats(shrubland, stat = sum) * 25 * 25) / 10e6

# podocarp <- DCChab %in% c(22)
# plot(podocarp)
# names(podocarp) <- "podocarp"

# broadleaf <- DCChab %in% c(2, 3, 24)
# plot(broadleaf)
# names(broadleaf) <- "broadleaf"

# native forest combines podocarp and broadleaf
native_forest <- DCChab %in% c(22, 2, 3, 24)
names(native_forest) <- "native_forest"
plot(native_forest)
(cellStats(native_forest, stat = sum) * 25 * 25) / 10e6

exotic_conifers <- DCChab %in% c(7)
names(exotic_conifers) <- "exotic_conifers"
plot(exotic_conifers)
cellStats(exotic_conifers, stat = sum)

exotic_hardwoods <- DCChab %in% c(5, 6, 20)
names(exotic_hardwoods) <- "exotic_hardwoods"
plot(exotic_hardwoods)

agriculture <- DCChab %in% c(14)
names(agriculture) <- "agriculture"
plot(agriculture)

suburban <- DCChab %in% c(1)
names(suburban) <- "suburban"
plot(suburban)

other <- DCChab %in% c(4, 8, 9, 10, 11, 12, 13, 17, 19, 23, 25, 26, 27,28, 29, 30)
names(other) <- "other"
plot(other)

# grassland <- DCChab %in% c(14, 4, 8, 9, 10, 11, 12, 13, 17, 19, 23, 25, 26, 27,28, 29, 30)
# plot(grassland)
# names(grassland) <- "grassland"

water <- DCChab %in% c(0, 21)
names(water) <- "water"
plot(water)
cellStats(water, stat = sum)

# non_water <- water == 0
# plot(non_water)
# names(non_water) <- "non_water"

all_rasters <- raster::stack(DCChab, shrubland, native_forest, exotic_conifers, exotic_hardwoods, agriculture, suburban, other, water)

plot(all_rasters)


# creating random steps ---------------------------------------------------


ssf05 <- T05t %>% steps_by_burst() %>% 
  random_steps(n = 30,
               rand_ta = random_numbers(make_unif_distr(), n = 1e+05)) %>% 
  extract_covariates(all_rasters) %>% 
  mutate(
    y = as.numeric(case_),
    cos_ta_ = cos(ta_),
    log_sl_ = log(sl_))

T05_ssf <- tibble("id" = 45505, ssf05)


dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% steps_by_burst() %>% 
                      random_steps(n = 30,
                                   rand_ta = random_numbers(make_unif_distr(), n = 1e+05)) %>%
                      extract_covariates(all_rasters))) %>% 
  dplyr::select(id, stps) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    cos_ta_ = cos(ta_),
    log_sl_ = log(sl_))


all_ssf <- rbind(T05_ssf, dat_ssf)
dat_ssf <- all_ssf %>% mutate(id = as.numeric(factor(id)),
                              step_id = paste(id, step_id_, sep = "-"))


# remove random points that fell into water and remove NAs from DCChab

dat_ssf <- dat_ssf %>% filter(water != 1 & !is.na(DCChab))
dat_ssf_tbl <- as_tibble(dat_ssf)
dat_ssf_tbl %>% dplyr::group_by(id, step_id_) %>% dplyr::summarise(n=n())
dat_ssf <- dat_ssf_tbl %>% dplyr::group_by(id, step_id_) %>% dplyr::slice_head(n = 11) # only keep 10 random points per used point


# plot(water)
plot(DCChab, col = terrain.colors(n = 8, rev = T))
points(dat_ssf$x2_, dat_ssf$y2_)
points(dat_ssf$x1_, dat_ssf$y1_, col = "red")


# plotting - not essential to run models ----------------------------------

dat_ssf$idf <- as.factor(dat_ssf$id)

dat_ssf %>% filter(y == 0) %>% ggplot(aes(x = sl_, fill = idf)) +
  geom_density(alpha = 0.1, show.legend = F) +
  scale_x_continuous(limits = c(0,3000), name = "Step Length (m)") +
  scale_y_continuous(name = "Density") +
  theme_classic()

# ggsave("StepLengthDsitribution-available.png", width=150, height=90, units="mm", dpi = 300)

dat_ssf %>% filter(y == 0) %>% ggplot(aes(x = ta_, fill = idf)) +
  geom_density(alpha = 0.1, show.legend = F) +
  scale_x_continuous(name = "Turning Angle (radians)") +
  scale_y_continuous(name = "Density") +
  theme_classic()

# ggsave("TurningAngleDsitribution-used.png", width=150, height=90, units="mm", dpi = 300)

used <- dat_ssf %>% filter(y == 1)
available <- dat_ssf %>% filter(y == 0)
# write_csv(used, "usedSSF_points_30thMay21.csv")
# write_csv(available, "availableSSF_points_30thMay21.csv")

sum(is.na(dat_ssf$ta_))
# removes nas from cos_ta
dat_ssf <- dat_ssf %>% filter(!is.na(cos_ta_))

# required to run models

dat_ssf$id1 <- dat_ssf$id
dat_ssf$id2 <- dat_ssf$id
dat_ssf$id3 <- dat_ssf$id
dat_ssf$id4 <- dat_ssf$id
dat_ssf$id5 <- dat_ssf$id
dat_ssf$id6 <- dat_ssf$id
dat_ssf$id7 <- dat_ssf$id
dat_ssf$id8 <- dat_ssf$id
dat_ssf$id9 <- dat_ssf$id
dat_ssf$id10 <- dat_ssf$id


# reading/writing csv files ready for INLA --------------------------------

write_csv(dat_ssf, paste0("outputs/dat_ssf_INLA_ready_", Sys.Date(), ".csv"))
dat_ssf <- read_csv("outputs/dat_ssf_INLA_ready_2024-01-11.csv")

dat_ssf$sl_scaled <- scale(dat_ssf$sl_) # for some reason this disables the ability to write_csv, so do this after exporting/importing

used <- dat_ssf %>% filter(y == 1)
available <- dat_ssf %>% filter(y == 0)


# INLA mixed model --------------------------------------------------------

formula.random <- y ~ -1 + 
  
  native_forest + 
  exotic_conifers + 
  exotic_hardwoods + 
  agriculture + 
  suburban + 
  other + 
  # shrubland +
  sl_ + 
  log_sl_ + 
  cos_ta_ +
  
  f(step_id, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(id1, native_forest, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id2, exotic_conifers, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id3, exotic_hardwoods, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id4, agriculture, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id5, suburban, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id6, other, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  # f(id7, shrubland, values = 1:10, model = "iid",
  #   hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id8, sl_, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id9, log_sl_, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01)))) +
  f(id10, cos_ta_, values = 1:10, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(1,0.01))))
 
mean.beta <- 0
prec.beta <- 1e-4

tic()
m_all2 <- inla(formula.random, family = "Poisson", data = dat_ssf,
              control.fixed = list(
                mean = mean.beta,
                prec = list(default = prec.beta)
                )
              )
toc()
beep(sound = 2)

summary(m_all2)
m_all$summary.fixed
m_all$summary.fixed$mean # beta estimates (mean of posterior)
exp(m_all$summary.fixed$mean) # exponent of beta estimates (mean of posterior)
m_all$summary.fixed$mode # beta estimates (mode of posterior)
exp(m_all$summary.fixed$mode) # exponent of beta estimates (mode of posterior)
max(exp(m_all$summary.fixed$mean))
exp.beta <- exp(m_all$summary.fixed$mean)
m_all$summary.hyperpar

plot(m_all$marginals.fixed$native_forest, type = "l")
plot(m_all$marginals.fixed$exotic_conifers, type = "l")
plot(m_all$marginals.fixed$suburban, type = "l")

# calculate the posterior modes
native_forest_mode <- m_all$marginals.fixed$native_forest[which.max(m_all$marginals.fixed$native_forest[,2]),1]
exotic_conifers_mode <- m_all$marginals.fixed$exotic_conifers[which.max(m_all$marginals.fixed$exotic_conifers[,2]),1]
exotic_hardwoods_mode <- m_all$marginals.fixed$exotic_hardwoods[which.max(m_all$marginals.fixed$exotic_hardwoods[,2]),1]
agriculture_mode <- m_all$marginals.fixed$agriculture[which.max(m_all$marginals.fixed$agriculture[,2]),1]
suburban_mode <- m_all$marginals.fixed$suburban[which.max(m_all$marginals.fixed$suburban[,2]),1]
other_mode <- m_all$marginals.fixed$other[which.max(m_all$marginals.fixed$other[,2]),1]
sl_mode <- m_all$marginals.fixed$sl_[which.max(m_all$marginals.fixed$sl_[,2]),1]
log_sl_mode <- m_all$marginals.fixed$log_sl_[which.max(m_all$marginals.fixed$log_sl_[,2]),1]
cos_ta_mode <- m_all$marginals.fixed$cos_ta_[which.max(m_all$marginals.fixed$cos_ta_[,2]),1]

# create data frame of modes
covariate_modes <- c(native_forest_mode, exotic_conifers_mode, exotic_hardwoods_mode, agriculture_mode, suburban_mode, other_mode, sl_mode, log_sl_mode, cos_ta_mode)
covariate_modes <- data.frame("coef" = m_all$names.fixed, "posterior_mode" = covariate_modes)

# inla_model_sl_logsl_costa <- m_all
# m_all <- inla_model_sl_logsl_costa

saveRDS(m_all, file = paste0("outputs/INLA_model_move_params_", Sys.Date(), ".rds")) # to save a single object
m_all <- readRDS("outputs/INLA_model_move_params_2024-01-11.rds")

# save(m_ndvi, m_ndvip1, m_expndvi, m_quad_ndvi, file = "INLA_models_with_ndvi.RData") # to save multiple objects simultaneously
# load("INLA_models_with_ndvi.RData")

# I think just playing around here
# unscaled raster values
c(1,exp(m_all$summary.fixed$mode[c(1:6)])) / (1 + c(1,exp(m_all$summary.fixed$mode[c(1:6)]))) 
# scaled
(c(1,exp(m_all$summary.fixed$mode[c(1:6)])) / (1 + c(1,exp(m_all$summary.fixed$mode[c(1:6)])))) / max(c(1,exp(m_all$summary.fixed$mode[c(1:6)])) / (1 + c(1,exp(m_all$summary.fixed$mode[c(1:6)]))))

c(1,exp(m_all$summary.fixed$mode[c(1:6)])) / max(exp(m_all$summary.fixed$mode[c(1:6)]))


# creating prediction maps ------------------------------------------------

#' Source R functions for calculating posterior means 
#' and medians of the precisions
source("inla_emarginal.R")
source("inla_mmarginal.R")
inla_emarginal(m_all)
inla_mmarginal(m_all)


coef_rasters <- raster::stack(shrubland, native_forest, exotic_conifers, exotic_hardwoods, agriculture, suburban, other)


exp.beta <- c(1,exp(m_all$summary.fixed$mode[c(1:6)])) # use the mode of posterior

hk <- raster::raster(coef_rasters)
hk <- raster::setValues(hk, 0)

for (i in 1:7) {
  hk <- hk + coef_rasters[[i]] * exp.beta[[i]]
}

raster::plot(hk, col = brewer.pal(9, "Reds"))
hist(hk, breaks = 100)


# with slightly different scaling
for (i in 1:7) {
  hk <- hk + ((coef_rasters[[i]] * exp.beta[[i]])/(1 + coef_rasters[[i]] * exp.beta[[i]]))
}

raster::plot(hk, col = brewer.pal(9, "Reds"))
hist(hk, breaks = 100)



# scaling by linear stretch (from 0 to 1) -----------------------------------------------

hk_scaled <- hk/raster::maxValue(hk)
hk_scaled <- (hk-raster::minValue(hk))/(raster::maxValue(hk)-raster::minValue(hk)) # linear stretch
hist(hk_scaled)
raster::plot(hk_scaled, col = brewer.pal(9, "Reds"))

hist(getValues(hk_scaled), breaks = 100)


# writing rasters to file -------------------------------------------------

writeRaster(hk, "inla projection base_agri_other - 31stAugust21.tif")
writeRaster(hk_scaled, "inla projection base_agri_other - scaled to one - 2ndSept21.tif")

hk <- raster("inla projection base_agri_other - 31stAugust21.tif")
hk_scaled <- raster("inla projection base_agri_other - scaled to one - 2ndSept21.tif")




# Model validation --------------------------------------------------------

# jacknifing

ssf_list <- vector(mode = "list", length = 10)

# create a list of data frames that contain all individuals besides 1, sequentially for each id
for (i in seq_along(c(1:10))) {
  ssf_list[[i]] <- dat_ssf %>% filter(id != i)
}

names(ssf_list) <- c("T05", "T06", "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14")

# create function to run an SSF model on each object of the list
inla_function <- function(ssf_data) {
  inla(formula.random, family = "Poisson", data = ssf_data,
       control.fixed = list(
         mean = mean.beta,
         prec = list(default = prec.beta))) }

# run models (same as above, different data), this may take a while to run
inla_models <- map(ssf_list, inla_function)

summary(inla_models[[4]]) # look at a single model

saveRDS(inla_models, file = "INLA_models_jacknife.rds") # to save a single object
inla_models <- readRDS("INLA_models_jacknife.rds") # this model was used in the thesis



# jacknife rasters

hk_rasters <- vector(mode = "list", length = 10) # create empty list with 10 slots

for (j in seq_along(c(1:10))) {
  
exp.beta <- c(1,exp(inla_models[[j]]$summary.fixed$mode[c(1:6)]))
hk <- raster::raster(coef_rasters)
hk <- raster::setValues(hk, 0)

for (i in 1:7) {
  hk <- hk + ((coef_rasters[[i]] * exp.beta[[i]])/(1 + coef_rasters[[i]] * exp.beta[[i]]))
}
raster::plot(hk, col = brewer.pal(9, "Reds"))
hk_rasters[[j]] <- hk
}

saveRDS(hk_rasters, file = "all_hk_raster_files_minus_1_tag.rds") # to save a single object
hk_rasters <- readRDS("all_hk_raster_files_minus_1_tag.rds") # this model was used in the thesis



# assessing model predictions ---------------------------------------------

ecospat.boyce()
extract <- raster::extract

# T05_xy <- T05t[c(1,2)]

for(i in seq_along(all_tracks)){
  all_tracks[[i]]$t_ <- NULL
}

raster::extract(scaled_hk_rasters[[7]], all_tracks_xy[[7]]) # look at a raster

# T05_xy_df <- as.data.frame(T05_xy) # ecospat function requires a df

df_function <- function(df){
  as.data.frame(df)
}

all_tracks_xy <- map(all_tracks, df_function)

saveRDS(all_tracks_xy, "CBI_ready_all_tracks_xy.rds")
all_tracks_xy <- readRDS("CBI_ready_all_tracks_xy.rds")



# running continuous boyce index models -----------------------------------

names(all_tracks_xy) <- c("T05", "T06", "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14")
list2env(all_tracks_xy, .GlobalEnv)

boyce_function <- function(df) {
  ecospat.boyce(hk, df, window.w = 1/5, res = 100) }

# more code in the section below to run more Boyce indices

boyce_out_base2_rsf_10_100 <- map(all_tracks_xy, boyce_function)

ids <- c("T05", "T06", "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14")

boyce_indices_jacknife <- boyce_out_base2_rsf_10_100
#boyce_indices_jacknife <- boyce_indices_jacknife5_100

for(i in seq_along(boyce_indices_jacknife)){
  boyce_indices_jacknife[[i]]$id <- ids[[i]]
}

boyce_indices_jacknife_df <- bind_rows(boyce_indices_jacknife)

boyce_indices_jacknife_df %>% summarise(Spearman.avg = mean(Spearman.cor))

boyce_indices_jacknife_df %>% ggplot(aes(x = HS, y = F.ratio, colour = id)) +
  geom_line() +
  theme_classic()

summarise <- dplyr::summarise


boyce_summary <- boyce_indices_jacknife_df %>% group_by(HS) %>% summarise(SRC = mean(Spearman.cor), mean = mean(F.ratio), median = median(F.ratio), mode = mode(F.ratio), sd = sd(F.ratio), n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

boyce_summary %>% summarise(SRC = mean(SRC))

boyce_summary$res <- 100
boyce_summary_100 <- boyce_summary

boyce_summary_100 %>% ggplot(aes(x = HS, y = mean)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.25) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.15) +
  #geom_line(aes(x = HS, y = lower.ci), linetype = "dashed", alpha = 0.5) +
  #geom_line(aes(x = HS, y = upper.ci), linetype = "dashed", alpha = 0.5) +
  theme_classic()

setwd("~/Desktop/OneDrive - University of Otago/MSc - Scott Forrest/Figures, Plots and Images/SSF")
ggsave("jackknife boyce index - window5res100.png", width = 8, height = 5)


hist(getValues(hk_scaled), breaks = 100, main = "Habitat Suitability Values", xlab = "Habitat Suitability")


# if you want to try multiple window sizes and resolutions

# from running multiple Boyce indices with different windows and resolutions (just change model specification above and write as new object), not essential
boyce_summary_all_w5 <- rbind(boyce_summary_10, boyce_summary_20, boyce_summary_50, boyce_summary_100, boyce_summary_200)
boyce_summary_all_w5 %>% summarise(meanSRC = mean(SRC))
boyce_summary_all_w5$Res <- as.factor(boyce_summary_all_w5$res)

boyce_summary_all_w5 %>% ggplot(aes(x = HS, y = mean, colour = Res)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.25) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.1) +
  scale_x_continuous(n.breaks = 20) +
  scale_y_continuous("F.ratio") +
  theme_bw()

setwd("~/Desktop/OneDrive - University of Otago/MSc - Scott Forrest/Figures, Plots and Images/SSF")
ggsave("jackknife boyce index - window5 - 10-20-50-100-200.png", width = 8, height = 5)

write_csv(boyce_indices_jacknife_df, "cont_boyce_index_base_jacknife_1stSept.csv")


# jacknifing --------------------------------------------------------------

# boyce_function <- function(df) {
#   ecospat.boyce(hk, df, window.w = 1/10, res = 100) }

boyce_indices_jacknife5_200 <- vector(mode = "list", length = 10)

for(k in seq_along(c(1:10))) {
  boyce_indices_jacknife5_200[[k]] <- ecospat.boyce(scaled_hk_rasters[[k]], all_tracks_xy[[k]], window.w = 1/5, res = 200)
}

### return to code above to analyse jackknifed data

# scaling rasters to make the HS values line up and enable calculating the mean and CIs for plotting --------

max_rsf_values <- vector(mode = "numeric", length = 10)

for(k in seq(1:10)) {
  max_rsf_values[[k]] <- maxValue(hk_rasters[[k]])
}

mean(max_rsf_values)

scaled_hk_rasters <- vector(mode = "list", length = 10)

scaled_hk <- hk_rasters[[1]]/maxValue(hk_rasters[[1]]) # trialling for loop with [[1]]
plot(scaled_hk, col = brewer.pal(9, "Reds"))

for (l in seq(1:10)) {
  scaled_hk_rasters[[l]] <- hk_rasters[[l]]/maxValue(hk_rasters[[l]])
}

for (m in seq(1:10)) {
plot(scaled_hk_rasters[[m]], col = brewer.pal(9, "Reds"))
}
