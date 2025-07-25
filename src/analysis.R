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
source("src/functions.R")
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)
options(warn=-1)
# data prep ---------------------------

kenya <- read.csv("raw_data/ctovc_final_cg.csv")


#------------  Random forest propensity score matching -----------
#------------ covariate imbalance -----------

# SEE HOW MEAN DIFFERENCES, DISTRIBUTIONS CHANGE 

std_diffs_orig <- kenya_sub[,-1] %>%
  melt(id.var = "treated") %>%
  group_by(variable, treated)%>%
  summarise(meanvariable = mean(value),
            varvariable = var(value),
            n = length(value)) %>% 
  ungroup()%>%
  gather(variable2, value, -(treated:variable)) %>%
  unite(temp, variable2, treated) %>%
  spread(temp, value)%>%
  mutate(std_diff = (meanvariable_1-meanvariable_0)/sqrt(varvariable_0/2 + varvariable_1/2), #two variances
         variablef = reorder(variable, abs(std_diff)),
         method = "Original")

std_diffs_rf <- matched_rf[,-1] %>%
  melt(id.var = "treated") %>%
  group_by(variable, treated)%>%
  summarise(meanvariable = mean(value),
            varvariable = var(value),
            n = length(value)) %>% 
  ungroup()%>%
  gather(variable2, value, -(treated:variable)) %>%
  unite(temp, variable2, treated) %>%
  spread(temp, value)%>%
  mutate(std_diff = (meanvariable_1-meanvariable_0)/sqrt(varvariable_0/2 + varvariable_1/2), #two variances
         variablef = reorder(variable, abs(std_diff)),
         method = "Matched")

rbind(std_diffs_orig,
      std_diffs_rf)%>%
  subset(variable != "initial_exp")%>%
  ggplot() + 
  geom_point(aes(x = variablef, y = abs(std_diff),  fill = method), 
             colour = "black", shape=21, stroke=1.1, size = 4) + 
  scale_fill_manual("",values = c("black", "white")) +
  coord_flip() + labs(y = "Absolute Normalized Difference in Means", x = "X variable") + theme_bw()
ggsave( "Prop_matching.pdf", width = 9, height =7, units = "in" )



#TABLE 2  
std_diffs%>% 
  mutate_at(c('meanvariable_0', 'meanvariable_1','std_diff', 'lrsd', 'TCP'), round, 2) %>%
  mutate_at(c('varvariable_0', 'varvariable_1'), sqrt) %>%
  mutate(variablef = factor(variable, levels = c("caregiver"  ,"elderly"  ,"head_sex" ,"head_married" ,  
                                                 "head_educ"  , "head_sick" ,
                                                 "ovc"  ,"hhsize" , "age_5under" ,"age_0617",  "age_1834" ,  "age_3549",  "age_5064",
                                                 "depend_ratio"  , "sex_ratio" ,  "adult_educ" ,       "non_active"  , "initial_exp"  ,
                                                 "agriculture"   , "salaried"   ,"casual", "self" ,    "transfers" , "rooms" ,           
                                                 "land" ,  "livestock" , "water" , "bike" , "radio" ,  "phone" ,  
                                                 "mosquito" , "road" , "market" ,  "communication" , "electricity" ),
                            labels = 
                              c("Age of caregiver (years)",  
                                "Head over 64 years (yes=1) ",  
                                "Female head (yes=1) ",  
                                "Married head (yes=1)  ",  
                                "Education of head (years) ",  
                                "Head sick in last month (yes=1) ",  
                                "No.\ of OVC in household ",  
                                "No.\ household members ", 
                                "No.\ members under 6 years",  
                                "No.\ members 6-17 years",  
                                "No.\ members 18-34 years ",  
                                "No.\ members 35-49 years ",  
                                "No.\ members 50-64 years",  
                                "Proportion of members under 18 or over 64",  
                                "Proportion of female members",  
                                "No.\ of adults with at least 8 years of education ",  
                                "No.\ of adults not active in labor force",
                                "Monthly expenditures (KSh) ",  
                                "Main income source: Agriculture (yes=1) ",  
                                "Main income source: Salaried work (yes=1)",  
                                "Main income source: Casual labor (yes=1)",  
                                "Main income source: Self-employment (yes=1)",  
                                "Main income source: Transfers (yes=1)",  
                                "No.\ rooms in dwelling",  
                                "Land owned (acres)",  
                                "Livestock owner (yes=1)",  
                                "Water from unprotected source (yes=1)",  
                                "Owns bicycle (yes=1) ",  
                                "Owns radio (yes=1) ",  
                                "Owns phone (yes=1)",  
                                "Owns mosquito net (yes=1)",  
                                "Village has access to road (yes=1)",  
                                "Distance to market over 5 kilometers (yes=1) ",
                                "Majority of village has access to phone (yes=1) ",
                                "Majority of village has electricity (yes=1)")))%>%
  arrange(variablef)%>%
  dplyr::select(c(14,8,7,11,6,10,12,13,14,16  )) %>%
  xtable() %>% print(include.rownames=FALSE)


#------------  X matrix ---------------------------

median_age <- median(filter(kenya, sample_nutrition == 1)$caregiver)

kenya$caregiver <- ifelse(kenya$caregiver >=median_age, 1, 0) #1 if older
#based on matching or not, what hh codes are we using?
#usehh <-unique(kenya$hhcode) #if using whole data set
usehh <- unique(matched_rf$hhcode)  # if using matching scheme matched_rf
kenya <- subset(kenya, hhcode %in% usehh)

#what folder should the rstan RDS files be safed
subfolder <- "Matched"#"Original"#

full_X <- kenya%>%
  mutate(land = ihs_trans(land),
         initial_exp = ihs_trans(initial_exp))%>%#adjust skewness to avoid too much influence
  mutate_at(c("hhsize", "non_active", "head_educ","depend_ratio", "sex_ratio",
              "age_5under", "age_0617", "age_1834","land", "age_3549", "age_5064", "ovc",# "initial_exp",
              "adult_educ", "rooms"),
            gel_stand)%>% #use gelman standardization on continuous/numeric X columns
  mutate(treated_year = treated*year,#interactions
         treated_age = treated*caregiver,
         year_age = year*caregiver,
         treated_year_age = treated*year*caregiver)%>% 
  dplyr::select(-c(hhcode, location, head_age,
                    expenditure, food, micro, education1, education2, health, diversity, 
                    stunting, wasting, underweight, schooling,enrolled, illness1, illness2,
                   cal_pd,cal_micro ,cal_cereal,cal_rt,cal_pl,cal_anim,cal_fv,cal_oil,cal_sugar,cal_misc,
                    sample_nutrition, sample_education, initial_exp,
                    large,             small,            poultry,           animals, elderly ))  %>% #remove from X matrix
  dplyr::select(c(treated,caregiver,year, treated_year, treated_age, year_age, treated_year_age), everything()) %>% #reorder
  as.matrix()
#head_age
full_X_cntr <- kenya%>%
  mutate(land = ihs_trans(land),
         initial_exp = ihs_trans(initial_exp))%>%
  mutate_at(c("hhsize", "non_active", "head_educ","depend_ratio", "sex_ratio",
              "age_5under", "age_0617", "age_1834","land", "age_3549", "age_5064", "ovc",# "initial_exp",
              "adult_educ", "rooms"),
            gel_stand)%>% #use gelman standardization on continuous/numeric X columns
  mutate(treated_year = 0,#treated*year,#interactions
         treated_age = treated*caregiver,
         year_age = year*caregiver,
         treated_year_age = 0)%>%#treated*year*head_sex)%>% 
  dplyr::select(-c(hhcode, location, head_age,
                   expenditure, food, micro, education1, education2, health, diversity, 
                   stunting, wasting, underweight, schooling,enrolled, illness1, illness2,
                   cal_pd,cal_micro ,cal_cereal,cal_rt,cal_pl,cal_anim,cal_fv,cal_oil,cal_sugar,cal_misc,
                   sample_nutrition, sample_education, initial_exp,
                   large,             small,            poultry,           animals, elderly ))  %>% #remove from X matrix
  dplyr::select(c(treated,caregiver,year, treated_year, treated_age, year_age, treated_year_age), everything()) %>% #reorder
  as.matrix()

#RUN STAN MODELS

#to save space, only keep specified parameters; "lp__" is automatically included
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
lower_bound <- quantile(kenya$cal_micro, 0.01, na.rm = TRUE) # it's zero...
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
summary(diversity_sampled)$summary %>%as.data.frame() %>% arrange(n_eff) %>%   head(6)  #original:ok (low =93) 05/07/25 matched:good 6/11
summary(stunting_sampled)$summary %>%as.data.frame() %>% arrange(n_eff) %>%head(6)   #original: good 05/05/2025   matched: good 6/11/25  
summary(wasting_sampled)$summary %>%as.data.frame() %>% arrange(n_eff) %>%head(6)    #original: good 05/05/2025 matched: ok 6/11 lowest 63
summary(underweight_sampled)$summary %>%as.data.frame() %>% arrange(n_eff) %>%head(6)#original: good  05/05/2025 matched: good 6/11
summary(cal_pd_sampled)$summary %>% as.data.frame() %>% arrange(n_eff) %>%head(6) # original: good 05/07/2025 matched: good 6/11/25
summary(cal_micro_sampled)$summary %>%as.data.frame() %>% arrange(n_eff) %>% head(6) # original: good 05/07/2025 matched: good 6/11
summary(cal_fv_sampled)$summary %>% as.data.frame() %>% arrange(n_eff) %>%head(6) # original: good 05/07/25 matched: good 6/21
summary(cal_anim_sampled)$summary %>%as.data.frame() %>% arrange(n_eff) %>% head(6) # original: good 05/07/25 matched: good


# Treatment effects ---------------------------
#capability set plots
#grey errorbar = treated
#black line - control
plot_cs(sampled=diversity_sampled,       y=y_stand(kenya$diversity, kenya),      response =  "dietary diversity",       data = kenya, backtrans = "diversity")
ggsave(paste(subfolder,"diversity_CS.pdf", sep = "/" ))

plot_cs(sampled=underweight_sampled, y = nutrition_data$underweight, response = "underweight",data = nutrition_data)
ggsave(paste(subfolder,"underweight_CS.pdf", sep = "/" ))

plot_cs(sampled=wasting_sampled,       y=nutrition_data$wasting,      response =  "wasting",       data = nutrition_data)
ggsave(paste(subfolder,"wasting_CS.pdf", sep = "/" ))

plot_cs(sampled=stunting_sampled,       y=nutrition_data$stunting,      response =  "stunting",       data = nutrition_data)
ggsave(paste(subfolder,"stunting_CS.pdf", sep = "/" ))

plot_cs(sampled=cal_fv_sampled,       y=ihs_trans(kenya$cal_fv_windsor),      response =  "fruits and vegetables",       data = kenya, backtrans="IHS")
ggsave(paste(subfolder,"cal_fv_CS.pdf", sep = "/" ))

plot_cs(sampled=cal_anim_sampled,       y=ihs_trans(kenya$cal_anim_windsor),      response =  "animal products",       data = kenya, backtrans="IHS")
ggsave(paste(subfolder,"cal_animv_CS.pdf", sep = "/" ))

plot_cs(sampled=cal_micro_sampled,       y=ihs_trans(kenya$cal_micro_windsor),      response =  "micronutrients",       data = kenya, backtrans="IHS")
ggsave(paste(subfolder,"cal_micro_CS.pdf", sep = "/" ))

plot_cs(sampled=cal_pd_sampled,       y=log(kenya$cal_pd_windsor),      response =  "caloric intake",       data = kenya, backtrans="log")
ggsave(paste(subfolder,"cal_pd_CS.pdf", sep = "/" ))

#untransformed
plot_cs(sampled=cal_anim_sampled,       y=ihs_trans(kenya$cal_anim_windsor),      response =  "animal products",       data = kenya, backtrans=FALSE)
ggsave(paste(subfolder,"cal_animv_CS_UNTRANSFORMED.pdf", sep = "/" ))
#send sentence explaining droopiness - household random effects

plot_cs(sampled=cal_micro_sampled,       y=ihs_trans(kenya$cal_micro_windsor),      response =  "micronutrients",       data = kenya, backtrans=FALSE) 
ggsave(paste(subfolder,"cal_micro_CS_UNTRANSFORMED.pdf", sep = "/" ))

plot_cs(sampled=diversity_sampled,       y=y_stand(kenya$diversity, kenya),      response =  "dietary diversity",       data = kenya, backtrans = FALSE)
ggsave(paste(subfolder,"diversity_CS_UNTRANSFORMED.pdf", sep = "/" ))



#capability set plots AGGREGATED
plot_cs_agg(sampled=underweight_sampled, y = nutrition_data$underweight, response = "underweight",data = nutrition_data)
ggsave(paste(subfolder,"underweight_CS_agg.pdf", sep = "/" ))

plot_cs_agg(sampled=diversity_sampled,       y=y_stand(kenya$diversity, kenya),      response =  "dietary diversity",       data = kenya, backtrans = "diversity")
ggsave(paste(subfolder,"diversity_CS_agg.pdf", sep = "/" ))

plot_cs_agg(sampled=wasting_sampled,       y=nutrition_data$wasting,      response =  "wasting",       data = nutrition_data)
ggsave(paste(subfolder,"wasting_CS_agg.pdf", sep = "/" ))

plot_cs_agg(sampled=stunting_sampled,       y=nutrition_data$stunting,      response =  "stunting",       data = nutrition_data)
ggsave(paste(subfolder,"stunting_CS_agg.pdf", sep = "/" ))

plot_cs_agg(sampled=cal_fv_sampled,       y=ihs_trans(kenya$cal_fv_windsor),      response =  "fruits and vegetables",       data = kenya, backtrans="IHS")
ggsave(paste(subfolder,"cal_fv_CS_agg.pdf", sep = "/" ))

plot_cs_agg(sampled=cal_anim_sampled,       y=ihs_trans(kenya$cal_anim_windsor),      response =  "animal products",       data = kenya, backtrans="IHS")
ggsave(paste(subfolder,"cal_animv_CS_agg.pdf", sep = "/" ))

plot_cs_agg(sampled=cal_micro_sampled,       y=ihs_trans(kenya$cal_micro_windsor),      response =  "micronutrients",       data = kenya, backtrans="IHS")
ggsave(paste(subfolder,"cal_micro_CS_agg.pdf", sep = "/" ))

plot_cs_agg(sampled=cal_pd_sampled,       y=log(kenya$cal_pd_windsor),      response =  "caloric intake",       data = kenya, backtrans="log") 
ggsave(paste(subfolder,"cal_pd_CS_agg.pdf", sep = "/" ))


plot_cs_agg(sampled=cal_anim_sampled,       y=ihs_trans(kenya$cal_anim_windsor),      response =  "animal products",       data = kenya, backtrans=FALSE)
ggsave(paste(subfolder,"cal_animv_CS_agg_UNTRANSFORMED.pdf", sep = "/" ))
#send sentence explaining droopiness - household random effects

plot_cs_agg(sampled=cal_micro_sampled,       y=ihs_trans(kenya$cal_micro_windsor),      response =  "micronutrients",       data = kenya, backtrans=FALSE) 
ggsave(paste(subfolder,"cal_micro_CS_agg_UNTRANSFORMED.pdf", sep = "/" ))


plot_cs_agg(sampled=diversity_sampled,       y=y_stand(kenya$diversity, kenya),      response =  "dietary diversity",       data = kenya, backtrans = FALSE)
ggsave(paste(subfolder,"diversity_CS_agg_UNTRANSFORMED.pdf", sep = "/" ))

###############################################

trt <- rbind(   
   data.frame(derivatives(full_X, diversity_sampled),    y = "dietary diversity"),
   data.frame(derivatives(nutrition_X, stunting_sampled),    y = "stunting"),
   data.frame(derivatives(nutrition_X, wasting_sampled),    y = "wasting"),
   data.frame(derivatives(nutrition_X, underweight_sampled),    y = "underweight"),
  data.frame(derivatives(full_X, cal_micro_sampled),    y = "micronutrients"),
  data.frame(derivatives(full_X, cal_pd_sampled),    y = "caloric intake"),
  data.frame(derivatives(full_X, cal_fv_sampled),    y = "fruits and vegetables"),
  data.frame(derivatives(full_X, cal_anim_sampled),    y = "animal products")
   ) %>% 
  mutate(y = factor(y, levels = c("caloric intake", "dietary diversity", "animal products",
                                  "fruits and vegetables", "micronutrients", "stunting", "wasting", "underweight")[8:1]),
         type = factor(type, levels = c("Functioning", "Capabilities", "Choice")[c(3,2,1)]))
  


#summary stats - transformed scale
ss <- trt %>%
  group_by(y, by, type)%>%
  summarise(post_mean = mean(effect),
            post_prob_greater_0 = mean(effect > 0),
            cred_lower_025 = quantile(effect, .025),
            cred_lower_975 = quantile(effect, .975))
View(ss)

ss%>%
write.csv(paste(subfolder,"trt_effect_summary.csv", sep = "/"), row.names = FALSE)


#Overall plot

filter(trt, by == "Overall") %>%
  group_by(type, y, by)%>%
  dplyr::summarise(post_mean = mean(effect),
            lower = quantile(effect, .025),
            upper = quantile(effect, .975)) %>%
  ggplot()+
  geom_pointrange(aes(x = y,y = post_mean, ymin = lower, ymax = upper, colour = type), position = position_dodge(width=c(0.6,0.4))) +
  #facet_grid(.~type, scales = "free") + 
  theme_bw() + 
  geom_hline(aes(yintercept = 0)) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45)) +
  scale_colour_grey("")+
  labs(x = "", y = "Treatment Effect")+guides(color = guide_legend(reverse = TRUE))
ggsave(paste(subfolder, "/trt_effects_agg.pdf", sep = "/"), width = 8, height =7, units = "in" )



#By age

filter(trt, by != "Overall") %>%
  group_by(type, y, by)%>%
  dplyr::summarise(post_mean = mean(effect),
            lower = quantile(effect, .025),
            upper = quantile(effect, .975)) %>%
  ggplot()+
  geom_pointrange(aes(x = y,y = post_mean, ymin = lower, ymax = upper, colour = type), position = position_dodge(width=c(0.6,0.4))) +
  facet_grid(.~by) + 
  theme_bw() + 
  geom_hline(aes(yintercept = 0)) + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45)) +
  scale_colour_grey("")+guides(color = guide_legend(reverse = TRUE)) +
  labs(x = "", y = "Treatment Effect") 
ggsave(paste(subfolder, "trt_effects.pdf", sep = "/"), width = 9, height =7, units = "in" )


# make treatment table

df_formatted <- ss %>%
  mutate(value = sprintf("%.3f (%.2f)", post_mean, post_prob_greater_0))


df_wide <- df_formatted %>%
  dplyr::select(y, type, by, value) %>%
  pivot_wider(names_from = by, values_from = value)


df_wide <- df_wide %>%
  arrange(y, type) %>% 
  dplyr::select(y, type, Overall, Older, Younger)
library(knitr)
# Step 4: Print LaTeX table using kable
kable(df_wide, format = "latex", booktabs = TRUE, escape = FALSE)

