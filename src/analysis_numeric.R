rm(list = ls())
library(rstan)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(GGally)
library(plyr)
library(dplyr)
library(xtable)
library(tidyr)
library(randomForest)
library(Matching)
library(bayesplot)
library(scales)
source("src/old/functions_old.R")
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)
options(warn=-1)
# data prep ---------------------------

kenya <- read.csv("raw_data/ctovc_final.csv")


# X matrix ---------------------------

#kenya$head_age <- ifelse(kenya$head_age >=55, 1, 0) #1 if older
#based on matching or not, what hh codes are we using?
usehh <-unique(kenya$hhcode) #if using whole data set
#usehh <- unique(matched_rf$hhcode)  # if using matching scheme matched_rf
kenya <- subset(kenya, hhcode %in% usehh)

#what folder should the rstan RDS files be safed
subfolder <- "Numeric"#"Matched"#

full_X <- kenya%>%
  mutate(land = ihs_trans(land),
         initial_exp = ihs_trans(initial_exp))%>%#adjust skewness to avoid too much influence
  mutate_at(c("hhsize", "non_active", "head_age","head_educ","depend_ratio", "sex_ratio",
              "age_5under", "age_0617", "age_1834","land", "age_3549", "age_5064", "ovc",# "initial_exp",
              "adult_educ", "rooms"),
            gel_stand)%>% #use gelman standardization on continuous/numeric X columns
  mutate(treated_year = treated*year,#interactions
         treated_age = treated*head_age,
         year_age = year*head_age,
         treated_year_age = treated*year*head_age)%>% 
  dplyr::select(-c(hhcode, location, 
                   expenditure, food, micro, education1, education2, health, diversity, 
                   stunting, wasting, underweight, schooling,enrolled, illness1, illness2,
                   cal_pd,cal_micro ,cal_cereal,cal_rt,cal_pl,cal_anim,cal_fv,cal_oil,cal_sugar,cal_misc,
                   sample_nutrition, sample_education, initial_exp,
                   large,             small,            poultry,           animals, elderly ))  %>% #remove from X matrix
  dplyr::select(c(treated,head_age,year, treated_year, treated_age, year_age, treated_year_age), everything()) %>% #reorder
  as.matrix()

full_X_cntr <- kenya%>%
  mutate(land = ihs_trans(land),
         initial_exp = ihs_trans(initial_exp))%>%
  mutate_at(c("hhsize", "non_active", "head_age", "head_educ","depend_ratio", "sex_ratio",
              "age_5under", "age_0617", "age_1834","land", "age_3549", "age_5064", "ovc",# "initial_exp",
              "adult_educ", "rooms"),
            gel_stand)%>% #use gelman standardization on continuous/numeric X columns
  mutate(treated_year = 0,#treated*year,#interactions
         treated_age = treated*head_age,
         year_age = year*head_age,
         treated_year_age = 0)%>%#treated*year*head_sex)%>% 
  dplyr::select(-c(hhcode, location, 
                   expenditure, food, micro, education1, education2, health, diversity, 
                   stunting, wasting, underweight, schooling,enrolled, illness1, illness2,
                   cal_pd,cal_micro ,cal_cereal,cal_rt,cal_pl,cal_anim,cal_fv,cal_oil,cal_sugar,cal_misc,
                   sample_nutrition, sample_education, initial_exp,
                   large,             small,            poultry,           animals, elderly ))  %>% #remove from X matrix
  dplyr::select(c(treated,head_age,year, treated_year, treated_age, year_age, treated_year_age), everything()) %>% #reorder
  as.matrix()

#RUN STAN MODELS

#tto save space, only keep specified parameters; "lp__" is automatically included
#par_keep <- c("f_pred","f_pred_cntr", "c", "c_cntr","u", "u_cntr","beta",
#              "sigma_a", "sigma_v", "sigma_l", 
#              "gamma1", "gamma2", "phi", "rho", "mu_a", "alpha", "tau")
par_keep <- c("f_pred","f_pred_cntr", 
              #"c", "c_cntr",
              #"u", "u_cntr",
              "beta",
              "sigma_a", "sigma_v", "sigma_l", 
              "gamma1", "gamma2",
              "phi", "rho", "mu_a", "alpha", "tau")

#education X matrix for nutrition
nutrition_X <- full_X[kenya$sample_nutrition == 1,]
nutrition_X_cntr <- full_X_cntr[kenya$sample_nutrition == 1,]
nutrition_data <- kenya[kenya$sample_nutrition == 1,] 

#RUN MCMC - OR, read the RDS file in next section


# diversity ---------------------------
#standardize compared to control
#talk about impact in terms of standard deviation from control
#nov 14: diversity needs c0 = 1 in order to achieve good sampling performance
#comparing c0=1 to c0=10, the treatment effects appear similar
#feb 24 - lowest neff is 101
diversity_sampled <- do_sampling(y=y_stand(kenya$diversity, kenya),
                                 X=full_X, 
                                 X_cntr = full_X_cntr,
                                 hh_id=kenya$hhcode, loc_id=kenya$location,
                                 file = "src/selection_model2.stan", kappa = 1)
saveRDS(diversity_sampled, file = paste(subfolder,"diversity_sampled.rds", sep = "/"))


# nutrition - stunting ---------------------------
#2/20/25 stunting had 2 parameters w neff < 100 at katppa = 10
#rerunning with kappa = 1
stunting_sampled <- do_sampling(y=nutrition_data$stunting,
                                X=nutrition_X, 
                                X_cntr = nutrition_X_cntr,
                                hh_id=nutrition_data$hhcode, loc_id=nutrition_data$location, kappa = 1,
                                file = "src/selection_model2.stan" )
saveRDS(stunting_sampled, file =paste(subfolder, "stunting_sampled.rds", sep = "/"))

# nutrition - wasting ---------------------------
wasting_sampled <- do_sampling(y=nutrition_data$wasting,
                               X=nutrition_X, 
                               X_cntr = nutrition_X_cntr,
                               hh_id=nutrition_data$hhcode, loc_id=nutrition_data$location, kappa = 1,
                               file = "src/selection_model2.stan" )
saveRDS(wasting_sampled, file =paste(subfolder, "wasting_sampled.rds", sep = "/"))

# nutrition - underweight ---------------------------
underweight_sampled <- do_sampling(y=nutrition_data$underweight,
                                   X=nutrition_X, 
                                   X_cntr = nutrition_X_cntr,
                                   hh_id=nutrition_data$hhcode, loc_id=nutrition_data$location, kappa = 1,
                                   file = "src/selection_model2.stan" )
saveRDS(underweight_sampled, file =paste(subfolder, "underweight_sampled.rds", sep = "/"))


# nutrition - calories per day  ---------------------------
lower_bound <- quantile(kenya$cal_pd, 0.01, na.rm = TRUE)
upper_bound <- quantile(kenya$cal_pd, 0.99, na.rm = TRUE)

# Windsorize the variable
hist(kenya$cal_pd)
kenya$cal_pd_windsor <- pmin(pmax(kenya$cal_pd, lower_bound), upper_bound)
hist(kenya$cal_pd_windsor)
hist(log(kenya$cal_pd_windsor))

cal_pd_sampled <- do_sampling(y=log(kenya$cal_pd_windsor),
                              X=full_X, 
                              X_cntr = full_X_cntr,
                              hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 1,
                              file = "src/selection_model2.stan" )
saveRDS(cal_pd_sampled, file =paste(subfolder, "cal_pd_sampled.rds", sep = "/"))


# nutrition - calories from micro ---------------------------
lower_bound <- quantile(kenya$cal_micro, 0.01, na.rm = TRUE)# it's zero...
upper_bound <- quantile(kenya$cal_micro, 0.99, na.rm = TRUE)

# Windsorize the variable
kenya$cal_micro_windsor <- pmin(pmax(kenya$cal_micro, lower_bound), upper_bound)

hist(kenya$cal_micro_windsor)
hist(ihs_trans(kenya$cal_micro_windsor))

cal_micro_sampled <- do_sampling(y=ihs_trans(kenya$cal_micro_windsor),
                                 X=full_X, 
                                 X_cntr = full_X_cntr,
                                 hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 1,
                                 file = "src/selection_model2.stan" )
saveRDS(cal_micro_sampled, file =paste(subfolder, "cal_micro_sampled.rds", sep = "/"))

# nutrition - FRUITS & VEGGIES  ---------------------------
lower_bound <- quantile(kenya$cal_fv, 0.01, na.rm = TRUE)
upper_bound <- quantile(kenya$cal_fv, 0.99, na.rm = TRUE)

# Windsorize the variable
kenya$cal_fv_windsor <- pmin(pmax(kenya$cal_fv, lower_bound), upper_bound)
hist(kenya$cal_fv_windsor)
hist(ihs_trans(kenya$cal_fv_windsor))

cal_fv_sampled <- do_sampling(y=ihs_trans(kenya$cal_fv_windsor),
                              X=full_X, 
                              X_cntr = full_X_cntr,
                              hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 1,
                              file = "src/selection_model2.stan" )
saveRDS(cal_fv_sampled, file =paste(subfolder, "cal_fv_sampled.rds", sep = "/"))

# nutrition - ANIMAL PROTEIN  ---------------------------
lower_bound <- quantile(kenya$cal_anim, 0.01, na.rm = TRUE)
upper_bound <- quantile(kenya$cal_anim, 0.99, na.rm = TRUE)

# Windsorize the variable
kenya$cal_anim_windsor <- pmin(pmax(kenya$cal_anim, lower_bound), upper_bound)
hist(kenya$cal_anim_windsor)
hist(ihs_trans(kenya$cal_anim_windsor))

cal_anim_sampled <- do_sampling(y=ihs_trans(kenya$cal_anim_windsor),
                                X=full_X, 
                                X_cntr = full_X_cntr,
                                hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 1,
                                file = "src/selection_model2.stan" )
saveRDS(cal_anim_sampled, file =paste(subfolder, "cal_anim_sampled.rds", sep = "/"))


# read RDS files ---------------------
diversity_sampled <- readRDS(paste(subfolder,'diversity_sampled.rds', sep = "/"))
stunting_sampled <- readRDS(paste(subfolder,'stunting_sampled.rds', sep = "/"))
wasting_sampled <- readRDS(paste(subfolder,'wasting_sampled.rds', sep = "/"))
underweight_sampled <- readRDS(paste(subfolder,'underweight_sampled.rds', sep = "/"))
cal_fv_sampled <- readRDS(paste(subfolder,'cal_fv_sampled.rds', sep = "/"))
cal_anim_sampled <- readRDS(paste(subfolder,'cal_anim_sampled.rds', sep = "/"))
cal_micro_sampled <- readRDS(paste(subfolder,'cal_micro_sampled.rds', sep = "/"))
cal_pd_sampled <- readRDS(paste(subfolder,'cal_pd_sampled.rds', sep = "/"))

#-- DIAGNOSTICS ------------
#trace plots
stan_trace(diversity_sampled, par = c("phi","gamma1","gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))#need to redo
stan_trace(wasting_sampled, par = c("phi","gamma1","gamma2[4]", "gamma2[7]","sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__")) #one rogue chain - redo 11/29
stan_trace(underweight_sampled, par = c("phi","gamma1","gamma2[4]", "gamma2[7]","beta[4]", "beta[7]", "sigma_a", "sigma_v",  "beta[24]", "beta[7]","rho", "lp__", "f_pred[1]", "f_pred[10]"))
stan_trace(cal_fv_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))


#Summarise R-hat values, effective sample sizes
#sort by rhat, n_eff, see trouble spots
#matched : original
summary(diversity_sampled)$summary %>% View  #original:ok (low =93) 05/09/25 
summary(stunting_sampled)$summary %>%View   #original: good 05/09/2025     
summary(wasting_sampled)$summary %>%View    #original: good 05/09/2025
summary(underweight_sampled)$summary %>%View#original: good  05/09/2025
summary(cal_pd_sampled)$summary %>% View # original: good 05/12/2025
summary(cal_micro_sampled)$summary %>% View # original: good 05/12/2025
summary(cal_fv_sampled)$summary %>% View # original: NOT CONVERGED
summary(cal_anim_sampled)$summary %>% View # original: good 05/12/25

#treatment effects------

#capability set plots AGGREGATED
plot_cs_agg_numeric(sampled=underweight_sampled, y = nutrition_data$underweight, response = "underweight",data = nutrition_data)
ggsave(paste(subfolder,"underweight_CS_agg.pdf", sep = "/" ))



#LRF TO DO:
#finish numeric age runs
#update order of trt_gender and trt pots

###############################################

trt <- rbind(   
  data.frame(derivatives_numeric(full_X, diversity_sampled),    y = "dietary diversity"),
  data.frame(derivatives_numeric(nutrition_X, stunting_sampled),    y = "stunting"),
  data.frame(derivatives_numeric(nutrition_X, wasting_sampled),    y = "wasting"),
  data.frame(derivatives_numeric(nutrition_X, underweight_sampled),    y = "underweight"),
  data.frame(derivatives_numeric(full_X, cal_micro_sampled),    y = "micronutrients"),
  data.frame(derivatives_numeric(full_X, cal_pd_sampled),    y = "caloric intake"),
  data.frame(derivatives_numeric(full_X, cal_fv_sampled),    y = "fruits and vegetables"),
  data.frame(derivatives_numeric(full_X, cal_anim_sampled),    y = "animal products")
) %>% 
  mutate(y = factor(y, levels = c("caloric intake", "dietary diversity", "animal products",
                                  "fruits and vegetables", "micronutrients", "stunting", "wasting", "underweight")))

trt %>%
  group_by(y,type, age)%>% #add y (response category)
  dplyr::summarise(post_mean = mean(effect),
                   lower = quantile(effect, .025),
                   upper = quantile(effect, .975)) %>% 
  ggplot() +
  geom_line(aes(x = age, y = post_mean)) +
  geom_ribbon(aes(x = age, ymin = lower, ymax = upper ), alpha = I(.3))+
  facet_grid(type~y, scales = "free_y")  +
  labs(y = "Treatment Effect", x = "Head of Household Age") +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = 2)

ggsave("trt_numeric_age.pdf")
