# Multiple fixed-effect (individual-level) models --------------------------------------------

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
library(broom)
library(ggpubr)
library(patchwork)

# import data files and remove any obvious errors. The majority of erroneous locations were removed with the Shimada et al 2012 speed filter

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

# the track of T05 was generated separately as it had a sampling interval of 2 hours rather than 3 hours
T_all <- rbind(T06, T07, T08, T09, T10, T11, T12, T13, T14)

dat_all <- T_all %>% nest(data = -id)
dat_all$sex <- c("m", "m", "f", "f", "f", "f", "m", "m", "f")
dat_all$age <- c(10, 5, 1, 3, 2, 2, 3, 10, 8) # only 9 individuals - 05 is 1 year old


# the track of T05 was generated separately as it had a sampling interval of 2 hours rather than 3 hours

T05t <- T05 %>% make_track(lon, lat, DateTime, crs = 4326) %>% transform_coords(crs_to = 2193)
plot(T05t$x_, T05t$y_)
summarize_sampling_rate(T05t)

# regularise the track to 2 hour intervals
T05t <- T05t %>% track_resample(rate = hours(2), tolerance = minutes(10)) # %>% 
#filter_min_n_burst(min_n = 3) %>% steps_by_burst()


# the remainder of the tags

dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, lon, lat, DateTime, crs = 4326) %>% 
      transform_coords(crs_to = 2193)
  }))

dat1 <- dat_all %>% mutate(dat_clean = map(trk, ~ {
  .x %>% track_resample(rate = hours(3), tolerance = minutes(10))
}))



# Habitat rasters -----------------------------------------------------------------

# import raster
DCChab <- raster("mapping/DCC BCNadj and LCDB.tif")
DCChab[DCChab>31] <- 0
plot(DCChab, col = terrain.colors(n = 31, rev = T))

#podocarp and broadleaf as native forest - map used
rcl <- cbind(c(15,16,18, 22, 2, 3, 24, 7, 5, 6, 20, 14, 1, 31, 
               4, 8, 9, 10, 11, 13, 17, 19, 23, 25, 26, 27, 28, 29, 30, 12, 0, 21),
             c(1, 1,1, 2, 2, 2, 2, 3, 4,4,4, 5, 6,6, 
               7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8))

DCC <- reclassify(DCChab, rcl, right = NA)
plot(-DCC)
names(DCC) <- "habs"

# generating random steps -------------------------------------------------

# uniform dist for turning angle

ssf05 <- T05t %>% steps_by_burst() %>% 
  random_steps(n = 30, # 30 random steps, will be reduced to 10 random steps
                 rand_ta = random_numbers(make_unif_distr(), n = 1e+05)) %>% 
  extract_covariates(DCC) %>% 
  mutate(
    y = as.numeric(case_),
    cos_ta_ = cos(ta_),
    log_sl_ = log(sl_))

T05_ssf <- tibble("id" = 45505, ssf05)
T05_ssf$idf <- as.factor(T05_ssf$id)


dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% steps_by_burst() %>% 
                      random_steps(n = 30, 
                                   rand_ta = random_numbers(make_unif_distr(), n = 1e+05)) %>%
                      extract_covariates(DCC))) %>% 
  select(id, stps) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    cos_ta_ = cos(ta_),
    log_sl_ = log(sl_))


dat_ssf$idf <- as.factor(dat_ssf$id)

# adding all tags together
dat_ssf <- rbind(T05_ssf, dat_ssf)

plot(DCC == 8)

dat_ssf <- dat_ssf %>% filter(habs != 8) # if a random point lands in the water, remove it (non-habitat), this will make the used/random ratio unequal
dat_ssf_tbl <- as_tibble(dat_ssf)
dat_ssf_tbl %>% dplyr::group_by(id, step_id_) %>% dplyr::summarise(n=n())
dat_ssf <- dat_ssf_tbl %>% dplyr::group_by(id, step_id_) %>% dplyr::slice_head(n = 11) # only keep 10 random points per used point to equalise used/random ratio for all used points

dat_ssf %>% filter(case_ == F) %>% 
  ggplot(aes(x = ta_, fill = idf)) +
  geom_density(show.legend = F, alpha = 0.1) +
  theme_classic()

# csv file of locations of all tags with random points, reading for CLR (conditional logistic regression - i.e. SSF)
# write_csv(dat_ssf, "dat_ssf_CLR_ready_30thJune21.csv") 

# plotting turning angles, case_ is used or random
dat_ssf %>% ggplot(aes(x = ta_, fill = case_)) +
  geom_density(show.legend = F, alpha = 0.25) +
  scale_x_continuous("Turning Angle (radians)") +
  scale_y_continuous("Density") +
  theme_classic()

# ggsave("turn_angle_geom_density25thJune21.png", width=200, height=120, units="mm", dpi = 300)

# same for step length
dat_ssf %>% ggplot(aes(x = sl_, fill = case_)) +
  geom_density(show.legend = F, alpha = 0.25) +
  scale_x_continuous("Step Length (m)", limits = c(0,5000)) +
  scale_y_continuous("Density") +
  theme_classic()

# ggsave("step_length_geom_density25thJune21.png", width=200, height=120, units="mm", dpi = 300)


# plot used and actual steps in x/y coords
plot(-DCC)
points(dat_ssf$x2_, dat_ssf$y2_)
points(dat_ssf$x1_, dat_ssf$y1_, col = "red")


dat_ssf$habsF <- factor(dat_ssf$habs, labels = c("Kanuka", "Native Forest", "Exotic Conifers", "Exotic Hardwoods", "Agriculture", "Suburban", "Other"))



# plotting ----------------------------------------------------------------

# used vs available
gg_tag <- function (gg_tag) { gg_tag %>% 
    dplyr::group_by(case_, habsF) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr:: mutate(prop = n / sum(n), 
                   label = paste0(round(prop * 100, 1), "%")) %>% 
    ggplot(aes(habsF, prop, fill = case_, group=case_,label = label)) + 
    geom_col(position = position_dodge2()) +
    geom_text(size = 4, vjust = -0.25, position = position_dodge(width = 1)) +
    labs(x = "Habitat Type", y = "Proportion", fill = "case_")+
    scale_fill_brewer(palette = "Paired", name="", 
                      breaks=c("FALSE", "TRUE"), labels=c("Available", "Used")) +
    theme_classic() +
    theme(axis.text.x = element_text(margin = margin(t = 10, b = 15), angle = 0), legend.position = c(0.9, 0.9)) }

gg_tag(dat_ssf)
# ggsave("use-available-all_individuals25thJune21.png", width=200, height=120, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45505) %>% gg_tag()
# ggsave("use-available-45505.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45506) %>% gg_tag()
# ggsave("use-available-45506.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45507) %>% gg_tag()
# ggsave("use-available-45507.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45508) %>% gg_tag()
# ggsave("use-available-45508.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45509) %>% gg_tag()
# ggsave("use-available-45509.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45510) %>% gg_tag()
# ggsave("use-available-45510.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45511) %>% gg_tag()
# ggsave("use-available-45511.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45512) %>% gg_tag()
# ggsave("use-available-45512.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45513) %>% gg_tag()
# ggsave("use-available-45513.png", width=200, height=150, units="mm", dpi = 300)

dat_ssf %>% filter(id == 45514) %>% gg_tag()
# ggsave("use-available-45514.png", width=200, height=150, units="mm", dpi = 300)


write_csv(dat_ssf, file = paste0("outputs/dat_ssf_CLR_ready_", Sys.Date(), ".csv"))



# running SSF models -------------------------------------------------------------------

dat_ssf <- read_csv("outputs/dat_ssf_CLR_ready_2024-01-10.csv")

# dat_ssf$sl_scaled <- scale(dat_ssf$sl_)
dat_ssf$sl_log <- log(dat_ssf$sl_)
dat_ssf_nested <- dat_ssf %>% group_by(id) %>% nest() # nest models to run all at the same time

unique(dat_ssf$id)

# run models
dat_ssf_nested <- dat_ssf_nested %>%
  mutate(ssf = lapply(data, function(x) {
    x %>% fit_issf(case_ ~ habsF + sl_ + log_sl_ + cos_ta_ + strata(step_id_))
  } ))

dat_ssf_nested$ssf[2] # take a look at a single model
for(i in 1:10) print(dat_ssf_nested$ssf[i]) # take a look at all models

ids <- seq(1,10)

exp(dat_ssf_nested$ssf[[2]]$model$coefficients) # return coefficient (exponentiated), just to take a look

dat_ssf_nested$Sex <- c("m", "m", "m", "f", "f", "f", "f", "m", "m", "f")
dat_ssf_nested$Age <- c(1, 10, 5, 1, 3, 2, 2, 3, 10, 8)
dat_ssf_nested$Origin <- c("o", "o", "c", "o", "o", "o", "o", "c", "c","o")


# to create resource selection plot of all individuals together, relative to (ref category - kānuka), with confidence intervals

# pulls out all coefficients and adds to data frame with other kākā info
d2 <- dat_ssf_nested %>% mutate(coef = map(ssf, ~ broom::tidy(.x$model))) %>% 
  dplyr::select(id, Sex, coef, Age, Origin) %>% unnest(cols = c(coef)) %>% mutate(id = factor(id)) 

# change the ordering of the coefficients for plotting - with descending mean
d2$order <- rep(c(2, 3, 4, 5, 1, 6, 7, 8, 9), 10)
d2 <- d2 %>% drop_na()

d3 <- d2 %>% filter(grepl("habs", term)) %>%  dplyr::group_by(term) %>% 
  dplyr::summarise(
    mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
    ymax = mean + 1.96 * sd(estimate))

# change the ordering of the coefficients for plotting - with descending mean
d3$order <- c(5, 3, 4, 2, 6, 1)
d3 <- d3[order(d3$order),]

d4 <- d2 %>% filter(grepl("habs", term)) %>% filter(estimate > -5)
d4$id <- as.factor(d4$id)
# change the ordering
# d4$order <- rep(c(5, 3, 4, 2, 6, 1), 10)
d4 <- d4[order(d4$order),]

d5 <- d4 %>% filter(grepl("habs", term)) %>% dplyr::group_by(term) %>% 
  dplyr::summarise(
    mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
    ymax = mean + 1.96 * sd(estimate))

# change the ordering
d5$order <- c(5, 3, 4, 2, 6, 1)
d5 <- d5[order(d5$order),]
# d5$x <- 1:nrow(d5)
d5$term <- factor(d5$term)

coefficients <- d4
coefficients$x <- 1:nrow(coefficients)

# write_csv(coefficients, "SSF_coefficients_d4_25thJune21.csv")
# write_csv(d5, "SSF_coefficients_d5_25thJune21.csv")
# write_csv(d7, "SSF_coefficients_d7_25thJune21.csv")

# coefficients <- read_csv("SSF_coefficients_d4_25thJune21.csv")
# d5 <- read_csv("SSF_coefficients_d5_25thJune21.csv")
# d7 <- read_csv("SSF_coefficients_d7_25thJune21.csv")


ggplot(data = d4, 
       aes(x = reorder(term, order), y = estimate, group = id, col = Age, pch = Sex)) +
  scale_colour_viridis_c("Age", breaks = seq(0,10,2)) +
  geom_rect(data = d5, 
            aes(xmin = order - 0.4, xmax = order + 0.4, ymin = ymin, ymax = ymax), 
            inherit.aes = F,
            fill = "grey90") +
  geom_segment(data = d5, aes(x = order - 0.4, xend = order + 0.4, y = mean, yend = mean), 
               inherit.aes = F,
               size = 0.25) +
  geom_pointrange(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                  position = position_dodge(width = 0.7), size = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Land Cover Type", y = expression(Coefficient~values~(beta))) +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(t = 10, b = 10), angle = 0)) +
  scale_x_discrete(labels = 
                     c("Suburban", "Native Forest", "Ex. Conifers", "Ex. Hardwoods", "Agriculture", "Other"))
  
  # ggsave("Individual_coefficients_15thJune21.png", width=200, height=120, units="mm", dpi = 300)
ggsave(paste0("outputs/plots/Individual_coefficients_uniform_ta_", Sys.Date(), ".png"), width=180, height=120, units="mm", dpi = 300)


# coefficients %>% filter(term == "habsFNative Forest") %>% ggplot(aes(x = age, y = estimate, colour = origin)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = T) +
#   theme_classic()
# 
# ggsave("individual_SSF_esimates_lm_15thJune21.png", width=200, height=120, units="mm", dpi = 300)
# 
# ggarrange(nf_lm, ec_lm)
# ggsave("native_forest_exotic_conifers_lm_25thJune21.png", width=200, height=90, units="mm", dpi = 300)



# running linear models on coefficients of habitat categories -------------

# run linear model on coefficients of a single habitat category with age
native_forest_lm <- lm(estimate ~ Age, data = coefficients, subset = (term =="habsFNative Forest"))
summary(native_forest_lm)
# plot
nf_lm <- effect_plot(native_forest_lm, pred = Age, plot.points = T, interval = T) +
  labs(x = "Age", y = expression(paste("Estimate"))) +
  ggtitle("Native Forest") +
  geom_hline(yintercept =  0, linetype = "dashed") +
  theme_classic()

# run linear model on coefficients of a single habitat category with age
exotic_conifers_lm <- lm(estimate ~ Age, data = coefficients, subset = (term =="habsFExotic Conifers"))
summary(exotic_conifers_lm)
ec_lm <- effect_plot(exotic_conifers_lm, pred = Age, plot.points = T, interval = T) +
  labs(x = "Age", y = expression(paste("Estimate"))) +
  ggtitle("Exotic Conifers") +
  geom_hline(yintercept =  0, linetype = "dashed") +
  theme_classic()

# run linear model on coefficients of a single habitat category with age
exotic_hardwoods_lm <- lm(estimate ~ Age, data = coefficients, subset = (term =="habsFExotic Hardwoods"))
summary(exotic_hardwoods_lm)
eh_lm <- effect_plot(exotic_hardwoods_lm, pred = Age, plot.points = T, interval = T) +
  labs(x = "Age", y = expression(paste("Estimate"))) +
  ggtitle("Exotic Hardwoods") +
  geom_hline(yintercept =  0, linetype = "dashed") +
  theme_classic()

# run linear model on coefficients of a single habitat category with age
agriculture <- lm(estimate ~ Age, data = coefficients, subset = (term =="habsFAgriculture"))
summary(agriculture)
effect_plot(agriculture, pred = Age, plot.points = T, interval = T) +
  theme_classic()

# run linear model on coefficients of a single habitat category with age
suburban <- lm(estimate ~ Age, data = coefficients, subset = (term =="habsFSuburban"))
summary(suburban)
effect_plot(suburban, pred = Age, plot.points = T, interval = T) +
  theme_classic()

# run linear model on coefficients of a single habitat category with age
other <- lm(estimate ~ Age, data = coefficients, subset = (term =="habsFOther"))
summary(other)
effect_plot(other, pred = Age, plot.points = T, interval = T) +
  theme_classic()


nf_lm + ec_lm
ggsave(paste0("outputs/plots/linear_models_forest_conifers", Sys.Date(), ".png"), width=180, height=100, units="mm", dpi = 300)
