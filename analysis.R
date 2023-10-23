rm(list = ls())
library(rstan) #MCMC sampling
library(reshape2) #wide-to-long melting
library(tidyverse)
library(randomForest) #estimate propensity scores
library(Matching) #matching treatment, control groups
library(loo) #leave one out crossvalidation
source("functions.R")
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)
options(warn=-1)
RNGkind(sample.kind = "Rounding")
# Read in Data ---------------------------

kenya <- read.csv("ctovc_final2023.csv")

# Propensity Score Matching ---------------------------

#data frame of individual-level pre-treatment covariates 
kenya_sub <- kenya %>%  
  subset(year == 0, select = c("hhcode","treated","head_age","hhsize", "non_active", "head_educ","depend_ratio", "sex_ratio",
                               "ovc", "initial_exp",
                               "adult_educ", "rooms", "head_sex", "head_married", "head_sick",'elderly',
                               'agriculture', 'salaried', 'casual', 'self',
                               "age_5under", "age_0617", "age_1834", "age_3549", "age_5064",
                               'transfers', 'land', 'livestock', 'water', 'bike', 'radio', 'phone', 'mosquito', 'road',
                               'market', 'communication', 'electricity'))%>% `rownames<-`(seq_len(sum(kenya$year == 0)))

#random forest to estimate propensity scores
set.seed(34893460)
rf <- randomForest(as.factor(treated) ~ ., data = kenya_sub[,-1],
                   ntree = 1000,
                   proximity = TRUE,
                   seed=48394639)
pscores_rf <- predict(rf, type = "prob")[,2]
#match treatment, control
set.seed(9235728)
matches_rf <- Match(Y = NULL, Tr= 1-kenya_sub$treated, X = pscores_rf,replace = FALSE, caliper = 1)
#extract matched treated, control individuals
matched_rf <- rbind(kenya_sub[matches_rf$index.treated,],
                    kenya_sub[matches_rf$index.control,])

std_diffs_rf <- matched_rf[,-1]  %>%
  melt(id.var = "treated") %>%
  group_by(variable, treated)%>%
  summarise(meanvariable = mean(value),
            varvariable = var(value),
            n = length(value)) %>% 
  ungroup()%>%
  gather(variable2, value, -(treated:variable)) %>%
  unite(temp, variable2, treated) %>%
  spread(temp, value)%>%
  mutate(std_diff = (meanvariable_1-meanvariable_0)/sqrt(varvariable_0/2 + varvariable_1/2),
         variablef = reorder(variable, std_diff),
         method = "Matched")


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


rbind(std_diffs_orig,
      std_diffs_rf)%>%
  subset(variable != "initial_exp")%>%
  ggplot() + 
  geom_point(aes(x = variablef, y = abs(std_diff),  fill = method), 
             colour = "black", shape=21, stroke=1.1, size = 4) + 
  scale_fill_manual("",values = c("black", "white")) +
  coord_flip() + labs(y = "Absolute Normalized Difference in Means", x = "X variable") + theme_bw()
ggsave( "Prop_matching.pdf", width = 9, height =7, units = "in" )


# X matrix ---------------------------

#Original or Matched analysis
#note:results will be saved in a pre-existing file with this name
subfolder <- "Original"#"Matched"#

#If matching, subset kenya data
if (subfolder == "Original"){
  usehh <-unique(kenya$hhcode) #if using whole data set
}else if (subfolder == "Matched"){
  usehh <- unique(matched_rf$hhcode) #if using matching scheme matched_rf
}

kenya <- subset(kenya, hhcode %in% usehh) 


#create x matrix
full_X <- kenya%>%
  mutate(land = ihs_trans(land))%>%
  mutate_at(c("head_age","hhsize", "non_active", "head_educ","depend_ratio", "sex_ratio",
              "age_5under", "age_0617", "age_1834","land", "age_3549", "age_5064", "ovc",
              "adult_educ", "rooms"),
            gel_stand)%>% #use gelman standardization on continuous/numeric X columns
  mutate(treated_year = treated*year,#interactions
         treated_sex = treated*head_sex,
         year_sex = year*head_sex,
         treated_year_sex = treated*year*head_sex)%>% 
  dplyr::select(-c(hhcode, location, 
                   expenditure, food, micro, education1, education2, health, diversity, 
                   stunting, wasting, underweight, schooling,enrolled, illness1, illness2,
                   sample_nutrition, sample_education, initial_exp,
                   large,             small,            poultry,           animals ))  %>% #remove from X matrix
  dplyr::select(c(treated,head_sex,year, treated_year, treated_sex, year_sex, treated_year_sex), everything()) %>% #reorder
  as.matrix()

#counterfactual X matrix
full_X_cntr <- full_X
full_X_cntr[,"treated_year"] <- 0
full_X_cntr[,"treated_year_sex"] <- 0

#education X matrix for schooling
education_X <- full_X[kenya$sample_education == 1,]
education_X_cntr <- full_X_cntr[kenya$sample_education == 1,]
education_data <- kenya[kenya$sample_education == 1,] 

#education X matrix for nutrition
nutrition_X <- full_X[kenya$sample_nutrition == 1,]
nutrition_X_cntr <- full_X_cntr[kenya$sample_nutrition == 1,]
nutrition_data <- kenya[kenya$sample_nutrition == 1,] 


#---STAN MODEL PREP------------
#to save space, only keep specified parameters; "lp__" is automatically included
par_keep <- c("f_pred","f_pred_cntr", "c", "c_cntr","u", "u_cntr","beta",
              "sigma_a", "sigma_v", "sigma_l", 
              "gamma1", "gamma2", "phi", "rho", "mu_a", "alpha", "tau")


ggplot(data = kenya) +
  geom_histogram(aes(x = illness1))

#RUN MCMC - OR, read the RDS file in next section
# illness 1 ---------------------------
#very left skewed
illness1_sampled <- do_sampling(y=(kenya$illness1),
                            X=full_X, 
                            X_cntr = full_X_cntr,
                            hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 10,
                            file = "BSFA_model.stan")
saveRDS(illness1_sampled, file = paste(subfolder,"illness1_sampled.rds", sep = "/"))

# nutrition expenditure - diversity ---------------------------
#standardize compared to control
diversity_sampled <- do_sampling(y=y_stand(kenya$diversity, kenya),
                                 X=full_X, 
                                 X_cntr = full_X_cntr,
                                 hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 10,
                                 file = "BSFA_model.stan")
saveRDS(diversity_sampled, file = paste(subfolder,"diversity_sampled.rds", sep = "/"))



#ONLY DEFINED FOR A SUBSAMPLE
# nutrition - stunting ---------------------------
stunting_sampled <- do_sampling(y=(kenya$stunting),
                             X=full_X, 
                             X_cntr = full_X_cntr,
                             hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 10,
                             file = "BSFA_model.stan" )
saveRDS(stunting_sampled, file =paste(subfolder, "stunting_sampled.rds", sep = "/"))

# nutrition - wasting ---------------------------
wasting_sampled <- do_sampling(y=(kenya$wasting),
                                X=full_X, 
                                X_cntr = full_X_cntr,
                                hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 10,
                                file = "BSFA_model.stan" )
saveRDS(wasting_sampled, file =paste(subfolder, "wasting_sampled.rds", sep = "/"))

# nutrition - underweight ---------------------------
underweight_sampled <- do_sampling(y=(kenya$underweight),
                                X=full_X, 
                                X_cntr = full_X_cntr,
                                hh_id=kenya$hhcode, loc_id=kenya$location, kappa = 10,
                                file = "BSFA_model.stan" )
saveRDS(underweight_sampled, file =paste(subfolder, "underweight_sampled.rds", sep = "/"))

# education - schooling ---------------------------

schooling_sampled <- do_sampling(y=education_data$schooling,
                                 education_X, 
                                 X_cntr = education_X_cntr,
                                 hh_id=education_data$hhcode, 
                                 loc_id=education_data$location,kappa = 10,
                                 file = "BSFA_model.stan")

saveRDS(schooling_sampled, file = paste(subfolder,"schooling_sampled.rds", sep = "/"))


# OR, read existing RDS files ---------------------
diversity_sampled <- readRDS(paste(subfolder,'diversity_sampled.rds', sep = "/"))
food_sampled <- readRDS(paste(subfolder,'food_sampled.rds', sep = "/"))
micro_sampled <- readRDS(paste(subfolder,'micro_sampled.rds', sep = "/"))
schooling_sampled <- readRDS(paste(subfolder,'schooling_sampled.rds', sep = "/"))
health_sampled <- readRDS(paste(subfolder,'health_sampled.rds', sep = "/"))

#-- DIAGNOSTICS ------------
#trace plots
stan_trace(diversity_sampled, par = c("phi","gamma1","gamma2[5]", "gamma2[1]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))
stan_trace(education1_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))
stan_trace(micro_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))
stan_trace(schooling_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__")) 
stan_trace(food_sampled, par = c("phi","gamma1", "gamma2", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))
stan_trace(health_sampled, par = c("phi","gamma1","gamma2[4]", "gamma2[7]","sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__")) 

#Summarise R-hat values, effective sample sizes
#sort by rhat, n_eff, see trouble spots
summary(diversity_sampled)$summary %>%View  
summary(food_sampled)$summary %>%View       
summary(micro_sampled)$summary %>%View    
summary(health_sampled)$summary %>%View     
summary(schooling_sampled)$summary %>%View  



# Treatment effects ---------------------------
####Capability Set Plots
    #agg over gender
plot_cs_agg(sampled=food_sampled,     y=ihs_trans(kenya$food),          response =  "Food Expenditure",        data = kenya) 
ggsave(paste(subfolder,"food_comb_CS.pdf", sep = "/" ))
plot_cs_agg(sampled=micro_sampled,    y = ihs_trans(kenya$micro),       response = "Micronutrient Expenditure",data =kenya)
ggsave(paste(subfolder,"micro_comb_CS.pdf", sep = "/" ))
plot_cs_agg(sampled=diversity_sampled,y=y_stand(kenya$diversity, kenya),response =  "Dietary Diversity",       data = kenya)
ggsave(paste(subfolder,"diversity_comb_CS.pdf", sep = "/" ))
plot_cs_agg(sampled=schooling_sampled,y=education_data$schooling,       response =  "Schooling",               data = education_data)
ggsave(paste(subfolder,"schooling_comb_CS.pdf", sep = "/" ))
plot_cs_agg(sampled=health_sampled,   y=ihs_trans(kenya$health),        response =  "Health Expenditure",      data = kenya)
ggsave(paste(subfolder,"health_comb_CS.pdf", sep = "/" ))

    #agg over gender, backtransformed 
plot_cs_agg(sampled=health_sampled,     y=ihs_trans(kenya$health),    response =  "Health Expenditure",      data = kenya, backtrans = TRUE)
ggsave(paste(subfolder,"health_comb_origscale_CS.pdf", sep = "/" ))
plot_cs_agg(sampled=food_sampled,       y=ihs_trans(kenya$food),      response =  "Food Expenditure",        data = kenya, backtrans = TRUE)
ggsave(paste(subfolder,"food_comb_origscale_CS.pdf", sep = "/" ))
plot_cs_agg(sampled=micro_sampled,      y = ihs_trans(kenya$micro),   response = "Micronutrient Expenditure",data =kenya,  backtrans = TRUE)
ggsave(paste(subfolder,"micro_comb_origscale_CS.pdf", sep = "/" ))


    #by gender
plot_cs(sampled=food_sampled,     y=ihs_trans(kenya$food),          response =  "Food Expenditure",        data = kenya)
ggsave(paste(subfolder,"food_CS.pdf", sep = "/" ))
plot_cs(sampled=micro_sampled,    y=ihs_trans(kenya$micro),         response = "Micronutrient Expenditure",data =kenya)
ggsave(paste(subfolder,"micro_CS.pdf", sep = "/" ))
plot_cs(sampled=diversity_sampled,y=y_stand(kenya$diversity, kenya),response =  "Dietary Diversity",       data = kenya)
ggsave(paste(subfolder,"diversity_CS.pdf", sep = "/" ))
plot_cs(sampled=schooling_sampled,y=education_data$schooling,       response =  "Schooling",               data = education_data)
ggsave(paste(subfolder,"schooling_CS.pdf", sep = "/" ))
plot_cs(sampled=health_sampled,   y=ihs_trans(kenya$health),        response =  "Health Expenditure",      data = kenya)
ggsave(paste(subfolder,"health_CS.pdf", sep = "/" ))

    #by gender, backtransformed
plot_cs(sampled=health_sampled, y=ihs_trans(kenya$health),  response =  "Health Expenditure",      data = kenya, backtrans = TRUE)
ggsave(paste(subfolder,"health_origscale_CS.pdf", sep = "/" ))
plot_cs(sampled=food_sampled,   y=ihs_trans(kenya$food),    response =  "Food Expenditure",        data = kenya, backtrans = TRUE)
ggsave(paste(subfolder,"food_origscale_CS.pdf", sep = "/" ))
plot_cs(sampled=micro_sampled,  y = ihs_trans(kenya$micro), response = "Micronutrient Expenditure",data =kenya,  backtrans = TRUE)
ggsave(paste(subfolder,"micro_origscale_CS.pdf", sep = "/" ))

####Treatment effect (Eq (18)) calculations, plots

trt <- rbind(data.frame(derivatives(full_X, food_sampled),    y = "food"),
   data.frame(derivatives(full_X, micro_sampled),    y = "micronutrients"),
   data.frame(derivatives(full_X, diversity_sampled),    y = "diversity"),
   data.frame(derivatives(full_X, health_sampled),    y = "health"),
   data.frame(derivatives(education_X, schooling_sampled),    y = "schooling"))%>%
  mutate(y = factor(y, levels = c("schooling", "health", "education", 
                                  "diversity", "micronutrients", "food", "expenditure")))


#summary stats - original (ihs) scale
trt %>%
  group_by(y, by, type)%>%
  summarise(post_mean = mean(effect),
            post_prob_greater_0 = mean(effect > 0),
            cred_lower_025 = quantile(effect, .025),
            cred_lower_975 = quantile(effect, .975))%>%
write.csv(paste(subfolder,"trt_effect_summary.csv", sep = "/"), row.names = FALSE)

#on original scale - trt effect on potentials
trt %>%
subset(y %in% c("health", "education", "micronutrients", "food", "expenditure") &
         type %in% c("Potential"))%>%
  mutate(exp_effect = exp(effect)) %>%
  group_by(y, by, type) %>%
  summarise(post_mean = mean(exp_effect),
            post_prob_greater_1= mean(exp_effect > 1),
            cred_lower_025 = quantile(exp_effect, .025),
            cred_lower_975 = quantile(exp_effect, .975))%>%
  write.csv(paste(subfolder,"trt_effect_orig_scale_summary.csv", sep = "/"), row.names = FALSE)
  

#Overall plot
subset(trt, by == "Overall") %>%
  group_by(type, y, by)%>%
  summarise(post_mean = mean(effect),
            lower = quantile(effect, .025),
            upper = quantile(effect, .975)) %>%
  ggplot()+
  geom_pointrange(aes(x = y,y = post_mean, ymin = lower, ymax = upper, colour = type), position = position_dodge(width=c(0.6,0.4))) +
  theme_bw() + 
  geom_hline(aes(yintercept = 0)) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45)) +
  scale_colour_grey("")+
  labs(x = "", y = "Treatment Effect")+guides(color = guide_legend(reverse = TRUE))
ggsave(paste(subfolder, "trt_effects.pdf", sep = "/"), width = 8, height =7, units = "in" )


#By gender
subset(trt, by != "Overall") %>%
  mutate(by2 = factor(by, levels = c("Male", "Female")))%>%
  group_by(type, y, by2)%>%
  summarise(post_mean = mean(effect),
            lower = quantile(effect, .025),
            upper = quantile(effect, .975)) %>%
  ggplot()+
  geom_pointrange(aes(x = y,y = post_mean, ymin = lower, ymax = upper, colour = type), position = position_dodge(width=c(0.6,0.4))) +
  facet_grid(.~by2) + 
  theme_bw() + 
  geom_hline(aes(yintercept = 0)) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45)) +
  scale_colour_grey("")+guides(color = guide_legend(reverse = TRUE)) +
  labs(x = "", y = "Treatment Effect")
ggsave(paste(subfolder, "trt_effects_gender.pdf", sep = "/"), width = 9, height =7, units = "in" )

## LOO estimates of ELPD--------------

food<- loo(extract_log_lik(food_sampled))
micro<- loo(extract_log_lik(micro_sampled))
diversity<- loo(extract_log_lik(diversity_sampled))
schooling<- loo(extract_log_lik(schooling_sampled))
health<- loo(extract_log_lik(health_sampled))

#use compare() or loo_compare() to compare ELPD estimates

