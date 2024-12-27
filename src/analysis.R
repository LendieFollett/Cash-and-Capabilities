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
source("src/old/functions_old.R")
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)
options(warn=-1)
# data prep ---------------------------

kenya <- read.csv("raw_data/ctovc_final2023.csv")

# exploratory ---------------------------

kenya %>%
  mutate(treated_year = interaction(treated, year),
         expenditure = ihs_trans(expenditure),
         micro = ihs_trans(micro),
         education1 = ihs_trans(education1)) %>%
  ggduo(columnsY = c("expenditure", "micro", "education1", "diversity", "schooling", "illness1", "enrolled"),
        columnsX = c("treated_year","hhsize"),
        mapping = aes(colour = treated))

#Pre propensity/distance matching

kenya_sub <- kenya %>%  
  subset(year == 0, select = c("hhcode","treated","head_age","hhsize", "non_active", "head_educ","depend_ratio", "sex_ratio",
                               "ovc", "initial_exp",
                               "adult_educ", "rooms", "head_sex", "head_married", "head_sick",'elderly',
                               'agriculture', 'salaried', 'casual', 'self',
                               "age_5under", "age_0617", "age_1834", "age_3549", "age_5064",
                               'transfers', 'land', 'livestock', 'water', 'bike', 'radio', 'phone', 'mosquito', 'road',
                               'market', 'communication', 'electricity'))%>% `rownames<-`(seq_len(sum(kenya$year == 0)))

# RANDOM FOREST
#trying random forest because: (1) interactions are built in (2) no need for transformations / linearity (3) deals with NAs if we do this step pre imputation
rf <- randomForest(as.factor(treated) ~ ., data = kenya_sub[,-1],
                   ntree = 1000,
                   proximity = TRUE,
                   seed=48394639)
pscores_rf <- predict(rf, type = "prob")[,2]

# n <- names(data.frame( pscores_rf, rf$proximity))
# f <- as.formula(paste("treated ~", paste(n[!n %in% c("treated")], collapse = " + ")))
# matches_rf1 <- matchit(f, method = "nearest", data=data.frame(treated = kenya_sub$treated, pscores_rf, rf$proximity))
# dta_m <- match.data(matches_rf1)
#RF 
set.seed(9235728)
matches_rf <- Match(Y = NULL, Tr= 1-kenya_sub$treated, X = pscores_rf,replace = FALSE, caliper = 1)

matched_rf <- rbind(kenya_sub[matches_rf$index.treated,],
                    kenya_sub[matches_rf$index.control,])

#matched_rf <- kenya_sub[pscores_rf > quantile(pscores_rf, .2) & pscores_rf < quantile(pscores_rf, .8),] #didn't work as well

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

#------------ RANDOM FOREST 2-----------
#RANDOM FOREST 2 - use both propensity scores AND proximity matrix to match
#proximity matrix is based on the frequency with which observations fall in the same nodes in the trees
set.seed(9235728)
matches_rf2 <- Match(Y = NULL, Tr= 1-kenya_sub$treated, X = cbind(pscores_rf,  rf$proximity),
                     caliper = c(.5, rep(Inf, ncol(rf$proximity))),
                     replace = FALSE)

matched_rf2 <- rbind(kenya_sub[matches_rf2$index.treated,],
                     kenya_sub[matches_rf2$index.control,])
std_diffs_rf2 <- matched_rf2[,-1]  %>%
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
         method = "rf2")

#------------ RANDOM FOREST 3-----------
# RANDOM FOREST 3 - still use both propensity scores AND proximity matrix to match
# but algorithm aligns more closely with Zhao paper

prox_mat <- rf$proximity - diag(ncol(rf$proximity)) #diagonals are 0 (will be looking for max w/in each row) - O.W., bigger means 'closer'
prop_calip <- (pscores_rf %>% dist(diag = TRUE, upper = TRUE)%>%as.matrix())< .25*sd(pscores_rf) %>% as.numeric()#FALSE (0) if pscores too far apart
prop_prox <-prox_mat*prop_calip

prox_treat <- prop_prox[kenya_sub$treated == 1 & pscores_rf < .8 & pscores_rf > .2, kenya_sub$treated == 0] #rows for treated ONLY, columns for control ONLY
prox_control <- prop_prox[kenya_sub$treated == 0, kenya_sub$treated == 1] #rows for control ONLY, columns for treated ONLY
prox_treat %>% dim; prox_control %>% dim
treated.samples <- apply(prox_control, 1, function(x){which.max(x)[1]}) #for each of the 701 controls, which is the closest treated? i.e., treated indices we want to keep
control.samples <- apply(prox_treat, 1, function(x){which.max(x)[1]}) #for each of the 701 controls, which is the closest treated? i.e., treated indices we want to keep

for (idx in which(duplicated(control.samples))){
  if (prox_treat[idx,order(prox_treat[idx,], decreasing=TRUE)[2]] == 0){
    print("no match made")
    control.samples[idx] <- NA
  } else{
    control.samples[idx] <- order(prox_treat[idx,], decreasing=TRUE)[2]
  }
}

for (idx in which(duplicated(control.samples))){
  
  if (prox_treat[idx,order(prox_treat[idx,], decreasing=TRUE)[3]] == 0){
    print("no match made")
    control.samples[idx] <- NA
  } else{
    control.samples[idx] <- order(prox_treat[idx,], decreasing=TRUE)[3]
  }
}


matched_rf3 <- rbind(subset(kenya_sub, treated == 0)[unique(control.samples[!is.na(control.samples)]),],
                     subset(kenya_sub, treated == 1 & pscores_rf < .8 & pscores_rf > .2)[which(!is.na(control.samples)),])

std_diffs_rf3 <- matched_rf3[,-1]  %>%
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
         method = "rf3")




#------------ LOGISTIC REGRESSION-----------
lr <- glm(treated ~ ., data = kenya_sub[,-1],family = "binomial")

pscores_lr <- predict(lr, type = "response")
set.seed(9235728)
matches_lr <- Match(Y = NULL, Tr= 1-kenya_sub$treated, X = pscores_lr,replace = FALSE, caliper = 1)

matched_lr <- rbind(kenya_sub[matches_lr$index.treated,],
                    kenya_sub[matches_lr$index.control,])


std_diffs_lr <- matched_lr[,-1]  %>%
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
         method = "lr")


#------------ compare-----------

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
  #std_diffs_rf2,
  #std_diffs_rf3,
  #std_diffs_lr)%>%
  ggplot() + 
  #geom_line(aes(x = variablef, y = abs(std_diff), colour = method, group = method)) + 
  geom_point(aes(x = variablef, y = abs(std_diff),  fill = method), colour = "black", shape=21, stroke=1.5) + 
  scale_fill_manual(values = c("black", "white")) +
  coord_flip() + labs(y = "Absolute Standardized Difference in Means", x = "X variable") + theme_bw()



#TABLE 2  
merge(std_diffs_orig[,c(1,2,3,6,7,8)],std_diffs_rf[,c(1,8)],by = "variable")%>% 
  mutate_at(c('meanvariable_0', 'meanvariable_1', "std_diff.x", "std_diff.y"), round, 2) %>%
  mutate_at(c('varvariable_0', 'varvariable_1'), sqrt) %>%
  mutate(variablef = factor(variable, levels = c("head_age"  ,"elderly"  ,"head_sex" ,"head_married" ,  
                                                 "head_educ"  , "head_sick" ,
                                                 "ovc"  ,"hhsize" , "age_5under" ,"age_0617",  "age_1834" ,  "age_3549",  "age_5064",
                                                 "depend_ratio"  , "sex_ratio" ,  "adult_educ" ,       "non_active"  , "initial_exp"  ,
                                                 "agriculture"   , "salaried"   ,"casual", "self" ,    "transfers" , "rooms" ,           
                                                 "land" ,  "livestock" , "water" , "bike" , "radio" ,  "phone" ,  
                                                 "mosquito" , "road" , "market" ,  "communication" , "electricity" ),
                            labels = 
                              c("Age of head (years)",  
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
  dplyr::select(c(8, 2,3,4,5,6 )) %>%
  xtable() %>% print(include.rownames=FALSE)


#see whole distributions: continuous-ish
rbind(mutate(matched_lr, method = "lr"),
      mutate(matched_rf, method = "rf"),
      mutate(matched_rf2, method = "rf2"),
      mutate(matched_rf3, method = "rf3"),
      mutate(kenya_sub, method = "original"))%>%
  subset(select = c("depend_ratio", "sex_ratio", "head_age", "initial_exp","method", "hhcode", "treated", "ovc", "head_educ","adult_educ", "non_active", "hhsize", 
                    "age_0617", "age_5under", "age_1834", "age_5064","age_3549", "rooms", "land"))%>%
  melt( id.vars = c("method", "hhcode", "treated")) %>%
  mutate(treated = as.factor(treated))%>%
  ggplot() +
  geom_histogram(aes(x = value, fill = treated), position = "identity", alpha = I(.6)) +
  facet_grid(method~variable, scales = "free")

#see whole distributions: discrete
rbind(mutate(matched_lr, method = "lr"),
      mutate(matched_rf, method = "rf"),
      mutate(matched_rf2, method = "rf2"),
      mutate(kenya_sub, method = "original"))%>%
  subset(select = -c(depend_ratio, sex_ratio, head_age, initial_exp, ovc, head_educ, hhsize,age_0617, non_active,age_5under, age_1834, age_5064,age_3549,adult_educ, rooms, land))%>%
  melt( id.vars = c("method", "hhcode", "treated")) %>%
  mutate(treated = as.factor(treated),
         value = as.character(value))%>%
  ggplot() +
  geom_bar(aes(x = value, fill = treated), position = "dodge") +
  facet_grid(method~variable, scales = "free")


# X matrix ---------------------------

#based on matching or not, what hh codes are we using?
usehh <-unique(kenya$hhcode) #if using whole data set
#usehh <- unique(matched_rf$hhcode)  # if using matching scheme matched_rf
kenya <- subset(kenya, hhcode %in% usehh)

#what folder should the rstan RDS files be safed
subfolder <- "Original"#"Matched"#

full_X <- kenya%>%
  mutate(land = ihs_trans(land),
         initial_exp = ihs_trans(initial_exp))%>%#adjust skewness to avoid too much influence
  mutate_at(c("head_age","hhsize", "non_active", "head_educ","depend_ratio", "sex_ratio",
              "age_5under", "age_0617", "age_1834","land", "age_3549", "age_5064", "ovc",# "initial_exp",
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

full_X_cntr <- kenya%>%
  mutate(land = ihs_trans(land),
         initial_exp = ihs_trans(initial_exp))%>%
  mutate_at(c("head_age","hhsize", "non_active", "head_educ","depend_ratio", "sex_ratio",
              "age_5under", "age_0617", "age_1834","land", "age_3549", "age_5064", "ovc",# "initial_exp",
              "adult_educ", "rooms"),
            gel_stand)%>% #use gelman standardization on continuous/numeric X columns
  mutate(treated_year = 0,#treated*year,#interactions
         treated_sex = treated*head_sex,
         year_sex = year*head_sex,
         treated_year_sex = 0)%>%#treated*year*head_sex)%>% 
  dplyr::select(-c(hhcode, location, 
                   expenditure, food, micro, education1, education2, health, diversity, 
                   stunting, wasting, underweight, schooling,enrolled, illness1, illness2,
                   sample_nutrition, sample_education, initial_exp,
                   large,             small,            poultry,           animals ))  %>% #remove from X matrix
  dplyr::select(c(treated,head_sex,year, treated_year, treated_sex, year_sex, treated_year_sex), everything()) %>% #reorder
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
par_keep_lm <- c("f_pred", "beta", "sigma_a", "sigma_v", "sigma_l", "phi", "rho", "mu_a")

education_X <- full_X[kenya$sample_education == 1,]
education_X_cntr <- full_X_cntr[kenya$sample_education == 1,]
education_data <- kenya[kenya$sample_education == 1,] #use throughout education

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
diversity_sampled <- do_sampling(y=y_stand(kenya$diversity, kenya),
                                 X=full_X, 
                                 X_cntr = full_X_cntr,
                                 hh_id=kenya$hhcode, loc_id=kenya$location,
                                 file = "src/selection_model2.stan", kappa = 1)
saveRDS(diversity_sampled, file = paste(subfolder,"diversity_sampled.rds", sep = "/"))


# illness ---------------------------

illness_sampled <- do_sampling(y=ihs_trans(kenya$illness1),
                              full_X, 
                              X_cntr = full_X_cntr,
                              hh_id=kenya$hhcode, loc_id=kenya$location,
                              file = "src/selection_model2.stan")
saveRDS(illness_sampled, file = paste(subfolder,"illness_sampled.rds", sep = "/"))


# education - schooling ---------------------------

# no transformation
schooling_sampled <- do_sampling(y=education_data$schooling,
                                 education_X, 
                                 X_cntr = education_X_cntr,
                                 hh_id=education_data$hhcode, 
                                 loc_id=education_data$location, kappa = 10,
                                 file = "src/selection_model2.stan")

saveRDS(schooling_sampled, file = paste(subfolder,"schooling_sampled.rds", sep = "/"))


# nutrition - stunting ---------------------------
stunting_sampled <- do_sampling(y=nutrition_data$stunting,
                                X=nutrition_X, 
                                X_cntr = nutrition_X_cntr,
                                hh_id=nutrition_data$hhcode, loc_id=nutrition_data$location, kappa = 10,
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
                               hh_id=nutrition_data$hhcode, loc_id=nutrition_data$location, kappa = 10,
                               file = "src/selection_model2.stan" )
saveRDS(underweight_sampled, file =paste(subfolder, "underweight_sampled.rds", sep = "/"))


# read RDS files ---------------------
diversity_sampled <- readRDS(paste(subfolder,'diversity_sampled.rds', sep = "/"))
diversity_sampled_lm <- readRDS(paste(subfolder,'diversity_sampled_lm.rds', sep = "/"))
#diversity_sampled_lm2 <- readRDS(paste(subfolder,'diversity_sampled_lm2.rds', sep = "/"))

expenditure_sampled <- readRDS(paste(subfolder,'expenditure_sampled.rds', sep = "/"))
expenditure_sampled_lm <- readRDS(paste(subfolder,'expenditure_sampled_lm.rds', sep = "/"))
#expenditure_sampled_lm2 <- readRDS(paste(subfolder,'expenditure_sampled_lm2.rds', sep = "/"))

education1_sampled <- readRDS(paste(subfolder,'education1_sampled.rds', sep = "/"))
education1_sampled_lm <- readRDS(paste(subfolder,'education1_sampled_lm.rds', sep = "/"))
#education1_sampled_lm2 <- readRDS(paste(subfolder,'education1_sampled_lm2.rds', sep = "/"))

food_sampled <- readRDS(paste(subfolder,'food_sampled.rds', sep = "/"))
food_sampled_lm <- readRDS(paste(subfolder,'food_sampled_lm.rds', sep = "/"))
#food_sampled_lm2 <- readRDS(paste(subfolder,'food_sampled_lm2.rds', sep = "/"))

micro_sampled <- readRDS(paste(subfolder,'micro_sampled.rds', sep = "/"))
micro_sampled_lm <- readRDS(paste(subfolder,'micro_sampled_lm.rds', sep = "/"))
#micro_sampled_lm2 <- readRDS(paste(subfolder,'micro_sampled_lm2.rds', sep = "/"))

schooling_sampled <- readRDS(paste(subfolder,'schooling_sampled.rds', sep = "/"))
schooling_sampled_lm <- readRDS(paste(subfolder,'schooling_sampled_lm.rds', sep = "/"))
#schooling_sampled_lm2 <- readRDS(paste(subfolder,'schooling_sampled_lm2.rds', sep = "/"))

health_sampled <- readRDS(paste(subfolder,'health_sampled.rds', sep = "/"))
health_sampled_lm <- readRDS(paste(subfolder,'health_sampled_lm.rds', sep = "/"))
#health_sampled_lm2 <- readRDS(paste(subfolder,'health_sampled_lm2.rds', sep = "/"))


#-- DIAGNOSTICS ------------
#trace plots
stan_trace(diversity_sampled, par = c("phi","gamma1","gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))#need to redo
stan_trace(education1_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))
stan_trace(micro_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))
stan_trace(schooling_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__")) #need to redo 11/29
stan_trace(food_sampled, par = c("phi","gamma1", "gamma2", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))
stan_trace(wasting_sampled, par = c("phi","gamma1","gamma2[4]", "gamma2[7]","sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__")) #one rogue chain - redo 11/29
stan_trace(underweight_sampled, par = c("phi","gamma1","gamma2[4]", "gamma2[7]","beta[4]", "beta[7]", "sigma_a", "sigma_v",  "beta[24]", "beta[7]","rho", "lp__", "f_pred[1]", "f_pred[10]"))
stan_trace(stunting_sampled, par = c("phi","gamma1", "gamma2[4]", "gamma2[7]", "sigma_a", "sigma_v",  "beta[4]", "beta[7]","rho", "lp__"))


#Summarise R-hat values, effective sample sizes
#sort by rhat, n_eff, see trouble spots
#matched : original
summary(diversity_sampled)$summary %>%View  #good - c0 = 1 : ok
summary(food_sampled)$summary %>%View       #ok            : good
summary(education1_sampled)$summary %>%View #good          : good
summary(wasting_sampled)$summary %>%View      #good          : good
summary(health_sampled)$summary %>%View     #good          : good
summary(schooling_sampled)$summary %>%View  #ok c0 = 1     : ok
summary(expenditure_sampled)$summary %>%View#good          : ok (trace good)


# Treatment effects ---------------------------
#capability set plots
plot_cs(sampled=schooling_sampled, y=ihs_trans(kenya$schooling),response =  "schooling",       data = education_data,X = education_X)
ggsave(paste(subfolder,"schooling_CS.pdf", sep = "/" ))

plot_cs(sampled=underweight_sampled,y = (y_stand(nutrition_data$underweight,nutrition_data)), response = "underweight",data = nutrition_data, backtrans=FALSE)
ggsave(paste(subfolder,"underweight_CS.pdf", sep = "/" ))

plot_cs(sampled=diversity_sampled,       y=y_stand(kenya$diversity, kenya),      response =  "Dietary Diversity",       data = kenya)
ggsave(paste(subfolder,"diversity_CS.pdf", sep = "/" ))

plot_cs(sampled=schooling_sampled,       y=education_data$schooling,      response =  "schooling",       data = education_data)
ggsave(paste(subfolder,"schooling_CS.pdf", sep = "/" ))

plot_cs(sampled=wasting_sampled,       y=ihs_trans(nutrition_data$wasting),      response =  "wasting",       data = nutrition_data)
ggsave(paste(subfolder,"health_CS.pdf", sep = "/" ))

plot_cs(sampled=stunting_sampled,       y=nutrition_data$stunting,      response =  "stunting",       data = nutrition_data)


###############################################

trt <- rbind(   data.frame(derivatives(full_X, expenditure_sampled), y = "expenditure"),
   data.frame(derivatives(full_X, food_sampled),    y = "food"),
   data.frame(derivatives(full_X, micro_sampled),    y = "micronutrients"),
   data.frame(derivatives(full_X, diversity_sampled),    y = "diversity"),
   data.frame(derivatives(full_X, education1_sampled),    y = "education"),
   data.frame(derivatives(full_X, health_sampled),    y = "health"),
   data.frame(derivatives(education_X, schooling_sampled),    y = "schooling"))%>%
  mutate(y = factor(y, levels = c("schooling", "health", "education", 
                                  "diversity", "micronutrients", "food", "expenditure")))


#summary stats - transformed scale
trt %>%
  group_by(y, by, type)%>%
  summarise(post_mean = mean(effect),
            post_prob_greater_0 = mean(effect > 0),
            cred_lower_025 = quantile(effect, .025),
            cred_lower_975 = quantile(effect, .975))%>%
write.csv("trt_effect_summary.csv", row.names = FALSE)

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
  write.csv("trt_effect_orig_scale_summary.csv", row.names = FALSE)
  

#Overall plot

subset(trt, by == "Overall") %>%
  group_by(type, y, by)%>%
  summarise(post_mean = mean(effect),
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
ggsave(paste(subfolder, "trt_effects.pdf", sep = "/"), width = 8, height =7, units = "in" )

#gets better when remove inefficiency
ggplot(data = subset(trt,  by == "Overall"))+ 
  geom_violin(aes(x = y, y = effect, fill = type),draw_quantiles = .5) +
  theme_bw() + 
  geom_hline(aes(yintercept = 0)) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_grey() +
  labs(x = "", y = "Treatment Effect")





#By gender

subset(trt, by != "Overall") %>%
  group_by(type, y, by)%>%
  summarise(post_mean = mean(effect),
            lower = quantile(effect, .025),
            upper = quantile(effect, .975)) %>%
  ggplot()+
  geom_pointrange(aes(x = y,y = post_mean, ymin = lower, ymax = upper, colour = type), position = position_dodge(width=c(0.6,0.4))) +
  facet_grid(.~by, scales = "free") + 
  theme_bw() + 
  geom_hline(aes(yintercept = 0)) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45)) +
  scale_colour_grey("")+guides(color = guide_legend(reverse = TRUE)) +
  labs(x = "", y = "Treatment Effect")
ggsave(paste(subfolder, "trt_effects_gender.pdf", sep = "/"), width = 9, height =7, units = "in" )



ggplot(data = subset(trt, type == "Inefficiency"))+ 
  geom_violin(aes(x = y, y = effect, fill = type),draw_quantiles = .5) +
  facet_grid(.~by, scales = "free") + 
  theme_bw() + 
  geom_hline(aes(yintercept = 0)) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_grey()+
  labs(x = "", y = "Treatment Effect")

