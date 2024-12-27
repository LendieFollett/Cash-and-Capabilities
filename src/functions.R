###########################################################
#do_sampling()
#y = column of responses (including NA for missings, pre-transformed if needed)
#X = matrix of predictor variables, pre-standardized
#hh_id = household id vector (should be same dimension as y)
#loc_id = location id vector ("                           ")
###########################################################

do_sampling <- function(y, X, X_cntr, hh_id, loc_id,file){
  #model that takes missingness into consideration
  compiled_selection_model <- stan_model(file )
  
  f <- y
  index_miss <- which(is.na(f))
  N_miss <- length(index_miss)
  #Missingness
  R <- as.numeric(is.na(f))
  
  m1 <- lm(f~ X) #some of f elements are NA, X is complete
  imputed_f <- f
  #impute missing f values based on observed relationship with X
  imputed_f[index_miss] <- predict(m1, newdata = data.frame(X))[index_miss]
  
  #model linking P(missing) to y and X
  ## (needed imputed f to do this)
  ML <- glm(R ~ imputed_f , family = binomial)
  
  phi_prior_mean <- coef(ML)[1:2]
  phi_prior_var <- vcov(ML)
  
  #replace y values with missing indicator
  #these values will be replaced in Stan with parameter values to be sampled
  #(so they will get a capability set too)
  f[index_miss] <- -99

  data.list <- list(N = length(f), 
                    N2 = nrow(X_cntr),
                    F = length(unique(hh_id)), 
                    L = length(unique(loc_id)),
                    K = ncol(X), 
                    hh = as.numeric(as.factor(hh_id)),
                    loc = as.numeric(as.factor(loc_id)),
                    f = f, 
                    X = X,
                    X_cntr=X_cntr,
                    R = R,
                    N_miss = N_miss,
                    index_miss = index_miss,
                    phi_prior_mean = phi_prior_mean,
                    phi_prior_var = phi_prior_var,
                    c0 = 10)
  
  initc <- list(gamma1 = -5, 
                gamma2 = rep(0, 7),#length(coef(m1)[-1])), 
                mu_a = coef(m1)[1], 
                beta = coef(m1)[-1], 
                sigma_a = .5*mean(m1$residuals^2)^(.5), sigma_a2 = .5*mean(m1$residuals^2),
                sigma_v = .5*mean(m1$residuals^2)^(.5), sigma_v2 = .5*mean(m1$residuals^2),
                sigma_l = .5*mean(m1$residuals^2)^(.5), sigma_l2 =.5*mean(m1$residuals^2),
                u = .5*(max(imputed_f)-imputed_f) + .01, #approx. 50 percent of distance from  max
                rho = .9, 
                f_imp = rep(mean(imputed_f), N_miss),
                phi = phi_prior_mean)
  
  inits <- list(initc,initc,initc,initc)  
  
  #Feed ML results from logistic regression into phi priors in sampling()
   sm_sampled <- sampling(compiled_selection_model, 
                          data = data.list, 
                          chains = 1,
                          iter =1000,
                          warmup = 1000,
                          thin =2,
                          init = inits,
                          control = list(adapt_delta = 0.8)) 
  return(sm_sampled)
}

#Gelman standardization
gel_stand <- function(x){(x - mean(x))/(2*sd(x))}

#Inverse hyperbolic sine transformation
ihs_trans <- function(x){log(x + sqrt(x^2 + 1))}

#In terms of standard deviation of un-treated
y_stand <- function(x, data){
  x2 <-x[data$treated == 0]
  (x - mean(x2,na.rm = TRUE))/(sd(x2,na.rm = TRUE))
}

#sampled is stan sample object
#y is response variable - transfrmed
#response = character string of y variable (just for labels)
plot_cs <- function(sampled, y, response, data, X){
# Plot estimated capability sets using selection model
f_pred_frame <- rstan::extract(sampled, "f_pred")%>%
  data.frame()
trt_effects <- (rstan::extract(sampled, "beta[1]")[[1]]%*%t(rep(1, nrow(data))) +
                  rstan::extract(sampled, "beta[4]")[[1]]%*%t(rep(1, nrow(data))) +
                  rstan::extract(sampled, "beta[5]")[[1]]%*%t(X[,"head_sex"]) +
                  rstan::extract(sampled, "beta[7]")[[1]]%*%t(X[,"head_sex"]))%>%
  apply(1,function(x){(x*t(1-X[,"treated"]))})%>%t()#each column contains sampled treatment effect for that person*year


notrt_effects <- (rstan::extract(sampled, "beta[1]")[[1]]%*%t(rep(1, nrow(data))) +
                    rstan::extract(sampled, "beta[4]")[[1]]%*%t(rep(1, nrow(data))) +
                    rstan::extract(sampled, "beta[5]")[[1]]%*%t(X[,"head_sex"])+
                    rstan::extract(sampled, "beta[7]")[[1]]%*%t(X[,"head_sex"]))%>%
  apply(1, function(x){x*X[,"treated"]})%>%t()#each column contains sampled treatment effect for that person*year

#[kenya$year == 1,] because only want to look at treatment
#effects after treatment had time to take effect
my_index <- (data$year == 1)
p <- data.frame(f =y[my_index],
           trt = data$treated[my_index],
           missing = is.na(y[my_index]),
           sex = factor(data$head_sex[my_index], levels = c(0,1), labels = c("Male", "Female")),
           lb_trt =   apply((f_pred_frame+trt_effects)[,my_index], 2, "quantile", c(0.025)),
           ub_trt =   apply((f_pred_frame+trt_effects)[,my_index], 2, "quantile", c(0.975)),
           lb_notrt = apply((f_pred_frame-notrt_effects)[,my_index], 2, "quantile", c(0.025)),
           ub_notrt = apply((f_pred_frame-notrt_effects)[,my_index], 2, "quantile", c(0.975))) %>%
  arrange(ub_trt) %>%
  mutate(id = 1:sum(my_index)) %>%
  ggplot() +
  geom_errorbar(aes(x=id, ymin=lb_trt, ymax=ub_trt), colour = "grey50") +
  geom_line(aes(x=id, y=lb_notrt  )) +
  geom_line(aes(x=id, y=ub_notrt )) +
  labs(x = "Household ID", y = response)+
  geom_point(aes(x = id, y =f),size = .1)+
  facet_wrap(.~sex)+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_grey()
print(p)
}

#sampled is stan sample object
#y is response variable - transfrmed
#response = character string of y variable (just for labels)
plot_cs <- function(sampled, y, response, data, X){
  
  f_pred_trt <-  rstan::extract(sampled, "f_pred")%>%
    data.frame()
  f_pred_cntr <-  rstan::extract(sampled, "f_pred_cntr")%>%
    data.frame()
  
  d <-  data %>%
    subset(year == 1 & treated == 1) 
  
  data.frame(f =y[data$year == 1& data$treated == 1],
             missing = is.na(y[data$year == 1& data$treated == 1]),
             sex = factor(d$head_sex, levels = c(0,1), labels = c("Male", "Female")),
             lb_trt =   apply((f_pred_trt)[,data$year == 1& data$treated == 1], 2, "quantile", c(0.10)),
             ub_trt =   apply((f_pred_trt)[,data$year == 1& data$treated == 1], 2, "quantile", c(0.95)),
             lb_notrt = apply((f_pred_cntr)[,data$year == 1& data$treated == 1], 2, "quantile", c(0.10)),
             ub_notrt = apply((f_pred_cntr)[,data$year == 1& data$treated == 1], 2, "quantile", c(0.95))) %>%
    arrange(sex, ub_notrt)%>%
    mutate(id = c(1:sum(d$head_sex == 0), 1:sum(d$head_sex == 1))) %>%
    #subset(missing == FALSE)%>%
    ggplot() +
    geom_errorbar(aes(x=id, ymin=lb_trt, ymax=ub_trt), colour = "grey50") +
    geom_line(aes(x=id, y=lb_notrt  )) +
    geom_line(aes(x=id, y=ub_notrt )) +
    labs(x = "Household ID", y = response)+
    geom_point(aes(x = id, y =f),size = .5,alpha = I(.3), colour  ="black")+
    facet_wrap(.~sex, scales = "free_x")+
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_grey() #+scale_y_continuous(limits = c(-4, 4))
  
}

#b1"treated"+b2"head_sex"+b3"year"+b4"treated_year"+b5"treated_sex"+b6"year_sex" +b7"treated_year_sex"
#female trt effect: B4 + B7
#male trt effect:   B4
get_trt_effects <- function(sampled, ca = TRUE){
  x <- rstan::extract(sampled, par = c("beta"))[[1]] %>%
    as.data.frame()%>%
    rename_all(funs(colnames(full_X)))%>%
    mutate(female = treated_year + treated_year_sex,
           male = treated_year) %>%
    dplyr::select(c(female, male))%>%
    melt()
  if(ca){
  y <- rstan::extract(sampled, par = c("gamma2"))[[1]] %>%
    as.data.frame()%>%
    rename_all(funs(colnames(full_X)))%>%
    mutate(female = treated_year + treated_year_sex,
           male = treated_year) %>%
    dplyr::select(c(female, male))%>%
    melt()
  d <- rbind(
    data.frame(x, type = "capability"),
    data.frame(y, type = "inefficiency")
  )
  }else{
   d <- data.frame(x, type = "functioning") 
  }
  return(d)
  #%>% #uncomment to use errorbar plot
  #ddply(.(variable), summarise,
  #    mean = mean(value),
  #    ub = quantile(value, .975),
  #    lb = quantile(value, .0275)) 
}

get_trt_effects2 <- function(sampled, ca = TRUE){
  x <- rstan::extract(sampled, par = c("beta"))[[1]] %>%
    as.data.frame()%>%
    rename_all(funs(colnames(full_X[,-24])))%>%
    mutate(female = treated_year + treated_year_sex,
           male = treated_year) %>%
    dplyr::select(c(female, male))%>%
    melt()
  if(ca){
    y <- rstan::extract(sampled, par = c("gamma2"))[[1]] %>%
      as.data.frame()%>%
      rename_all(funs(colnames(full_X[,-24])))%>%
      mutate(female = treated_year + treated_year_sex,
             male = treated_year) %>%
      dplyr::select(c(female, male))%>%
      melt()
    d <- rbind(
      data.frame(x, type = "capability"),
      data.frame(y, type = "inefficiency")
    )
  }else{
    d <- data.frame(x, type = "functioning") 
  }
  return(d)
  #%>% #uncomment to use errorbar plot
  #ddply(.(variable), summarise,
  #    mean = mean(value),
  #    ub = quantile(value, .975),
  #    lb = quantile(value, .0275)) 
}


cleanObject <- function(object,pars){
  pars <- c(pars,'lp__')
  nn <- paste0('^',pars,'(\\[|$)',collapse="|")
  ids <-  grep(nn,  object@sim$fnames_oi)
  ids.2 <- which(names(object@par_dims) %in% pars)
  for(i in 1:4){
    a <- attributes(object@sim$samples[[i]])
    x <- object@sim$samples[[i]][ids]
    for(j in c('names','inits','mean_pars'))
      a[[j]] <- a[[j]][ids]
    attributes(x) <- a
    object@sim$samples[[i]] <- x
  }
  object@par_dims <- object@par_dims[ids.2]
  object@sim$dims_oi <-   object@sim$dims_oi[ids.2]  
  object@sim$pars_oi<- object@sim$pars_oi[ids.2]
  object@sim$fnames_oi <-  object@sim$fnames_oi[ids]
  object@sim$n_flatnames <- length(object@sim$fnames_oi)
  object
}





derivatives <- function(full_X, samples){ 
  trt_mean <- 1
  full_X_mean <- c(1,apply(full_X, 2, mean))
  full_X_mean[c("treated", "head_sex", "year", 
                "treated_year", "treated_sex", "year_sex", 
                "treated_year_sex")] <- c(trt_mean, full_X_mean["head_sex"], 1,
                                          trt_mean*1, trt_mean*full_X_mean["head_sex"],full_X_mean["head_sex"]*1,
                                          full_X_mean["head_sex"]*1*trt_mean)
  
  female_X_mean <- c(1,apply(full_X, 2, mean))
  female_X_mean[c("treated", "head_sex", "year", 
                  "treated_year", "treated_sex", "year_sex", 
                  "treated_year_sex")] <- c(trt_mean, 1, 1,
                                            trt_mean*1, trt_mean*1,1*1,
                                            1*1*trt_mean)
  
  male_X_mean <- c(1,apply(full_X, 2, mean))
  male_X_mean[c("treated", "head_sex", "year", 
                "treated_year", "treated_sex", "year_sex", 
                "treated_year_sex")] <- c(trt_mean, 0, 1,
                                          trt_mean*1, trt_mean*0,0*1,
                                          0*1*trt_mean)
  
  
  gamma_samps <- do.call("cbind", rstan::extract(samples, par = c("gamma1","gamma2")))
  rho_samps <- rstan::extract(samples, par = c("rho"))[[1]]
  
  
  #INEFFICIENCY effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%full_X_mean[1:8])}) #exp(X*gamma) = lambda(1)
  temp2 <- gamma_samps[,5] + gamma_samps[,8]*full_X_mean["head_sex"] #d/db X*gamma
  dudb <- rho_samps*temp1*temp2
  dudb_trim <- dudb[dudb < quantile(dudb, .975) & dudb > quantile(dudb, .025)]
  #female effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%female_X_mean[1:8])}) #exp(X*gamma)
  temp2 <- gamma_samps[,5]  + gamma_samps[,8]*1 #d/db X*gamma
  dudb_female <- temp1*temp2*rho_samps
  dudb_female_trim <- dudb_female[dudb_female < quantile(dudb_female, .975) & dudb_female > quantile(dudb_female, .025)]
  #male effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%male_X_mean[1:8])}) #exp(X*gamma)
  temp2 <- gamma_samps[,5] + gamma_samps[,8]*0 #d/db X*gamma
  dudb_male <- temp1*temp2*rho_samps
  dudb_male_trim <- dudb_male[dudb_male < quantile(dudb_male, .975) & dudb_male > quantile(dudb_male, .025)]
  rm(gamma_samps)
  rm(rho_samps)
  
  #CAPABILITY EFFECT
  
  beta_samps <- do.call("cbind", rstan::extract(samples, par = c("beta")))
  #overall effect
  dcdb <- beta_samps[,4] + beta_samps[,7]*full_X_mean["head_sex"]
  
  #female effect
  dcdb_female <- beta_samps[,4] + beta_samps[,7]*1
  
  #male effect
  dcdb_male <- beta_samps[,4] + beta_samps[,7]*0
  rm(beta_samps)
  
  #FUNCTIONING effect 
  dydb <- dcdb - dudb
  dydb_female <- dcdb_female - dudb_female
  dydb_male <- dcdb_male - dudb_male

  
  all <- rbind(data.frame(effect = dydb, type = "Functionings", by = "Overall"), 
               data.frame(effect = dydb_female, type = "Functionings", by = "Female"), 
               data.frame(effect = dydb_male, type = "Functionings", by = "Male"), 
               
               data.frame(effect = dudb, type = "Inefficiencies", by = "Overall"), 
               data.frame(effect = dudb_female, type = "Inefficiencies", by = "Female"), 
               data.frame(effect = dudb_male, type = "Inefficiencies", by = "Male"), 
               
               data.frame(effect = dcdb, type = "Capabilities", by = "Overall"), 
               data.frame(effect = dcdb_female, type = "Capabilities", by = "Female"), 
               data.frame(effect = dcdb_male, type = "Capabilities", by = "Male"))
  return(all)
}


derivatives <- function(full_X, samples){ 
  full_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  full_X_mean[c("treated", "head_sex", "year", 
                "treated_year", "treated_sex", "year_sex", 
                "treated_year_sex")] <- c(1, full_X_mean["head_sex"], 1,
                                          1*1, 1*full_X_mean["head_sex"],full_X_mean["head_sex"]*1,
                                          full_X_mean["head_sex"]*1*1)
  
  female_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  female_X_mean[c("treated", "head_sex", "year", 
                  "treated_year", "treated_sex", "year_sex", 
                  "treated_year_sex")] <- c(1, 1, 1,
                                            1*1, 1*1,1*1,
                                            1*1*1)
  
  male_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  male_X_mean[c("treated", "head_sex", "year", 
                "treated_year", "treated_sex", "year_sex", 
                "treated_year_sex")] <- c(1, 0, 1,
                                           1*1, 1*0,0*1,
                                           0*1*1)
  
  
  gamma_samps <- do.call("cbind", rstan::extract(samples, par = c("gamma1","gamma2")))
  rho_samps <- rstan::extract(samples, par = c("rho"))[[1]]
  
  
  #INEFFICIENCY effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%full_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%full_X_mean[c(5,8)]})) # 1/exp(g2 + g3*.6244)
  dudb <- temp1*(temp2-1)*rho_samps
  dudb_trim <- dudb[dudb < quantile(dudb, .975) & dudb > quantile(dudb, .025)]
  #female effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%female_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%female_X_mean[c(5,8)]})) # 1/exp(g2 + g3*1)
  dudb_female <- temp1*(temp2-1)*rho_samps
  dudb_female_trim <- dudb_female[dudb_female < quantile(dudb_female, .975) & dudb_female > quantile(dudb_female, .025)]
  #male effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%male_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%male_X_mean[c(5,8)]})) # 1/exp(g2 + g3*0)
  dudb_male <- temp1*(temp2-1)*rho_samps
  dudb_male_trim <- dudb_male[dudb_male < quantile(dudb_male, .975) & dudb_male > quantile(dudb_male, .025)]
  rm(gamma_samps)
  rm(rho_samps)
  
  #CAPABILITY EFFECT
  
  beta_samps <- do.call("cbind", rstan::extract(samples, par = c("beta")))
  #overall effect
  dcdb <- beta_samps[,4] + beta_samps[,7]*full_X_mean["head_sex"]
  #female effect
  dcdb_female <- beta_samps[,4] + beta_samps[,7]*1
  #male effect
  dcdb_male <- beta_samps[,4] + beta_samps[,7]*0
  rm(beta_samps)
  #FUNCTIONING effect 
  dydb <-        dcdb        + dudb
  dydb_female <- dcdb_female + dudb_female
  dydb_male <-   dcdb_male   + dudb_male
  
  
  all <- rbind(data.frame(effect = dydb, type = "Functioning", by = "Overall"), 
               data.frame(effect = dydb_female, type = "Functioning", by = "Female"), 
               data.frame(effect = dydb_male, type = "Functioning", by = "Male"), 
               
               data.frame(effect = dudb, type = "Deviation", by = "Overall"), 
               data.frame(effect = dudb_female, type = "Deviation", by = "Female"), 
               data.frame(effect = dudb_male, type = "Deviation", by = "Male"), 
               
               data.frame(effect = dcdb, type = "Potential", by = "Overall"), 
               data.frame(effect = dcdb_female, type = "Potential", by = "Female"), 
               data.frame(effect = dcdb_male, type = "Potential", by = "Male"))
  return(all)
}


orig_scale_u_effect <- function(full_X, samples){
  gamma_samps <- do.call("cbind", rstan::extract(samples, par = c("gamma1","gamma2")))
  rho_samps <- rstan::extract(samples, par = c("rho"))[[1]]

  full_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  full_X_mean[c("treated", "head_sex", "year", 
                "treated_year", "treated_sex", "year_sex", 
                "treated_year_sex")] <- c(1, full_X_mean["head_sex"], 1,
                                          1*1, 1*full_X_mean["head_sex"],full_X_mean["head_sex"]*1,
                                          full_X_mean["head_sex"]*1*1)
  
  female_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  female_X_mean[c("treated", "head_sex", "year", 
                  "treated_year", "treated_sex", "year_sex", 
                  "treated_year_sex")] <- c(1, 1, 1,
                                            1*1, 1*1,1*1,
                                            1*1*1)
  
  male_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  male_X_mean[c("treated", "head_sex", "year", 
                "treated_year", "treated_sex", "year_sex", 
                "treated_year_sex")] <- c(1, 0, 1,
                                          1*1, 1*0,0*1,
                                          0*1*1)
  

  lambda1 <- exp(apply(gamma_samps, 1, function(x){x%*%full_X_mean}))
  lambda1_female <- exp(apply(gamma_samps, 1, function(x){x%*%female_X_mean}))  
  lambda1_male <- exp(apply(gamma_samps, 1, function(x){x%*%male_X_mean}))
  
  full_X_mean["treated_year"] <-  full_X_mean["treated_year_sex"] <- 0
  female_X_mean["treated_year"] <-  female_X_mean["treated_year_sex"] <- 0  
  male_X_mean["treated_year"] <-  male_X_mean["treated_year_sex"] <- 0  
  
  lambda0 <- exp(apply(gamma_samps, 1, function(x){x%*%full_X_mean}))
  lambda0_female <- exp(apply(gamma_samps, 1, function(x){x%*%female_X_mean}))  
  lambda0_male <- exp(apply(gamma_samps, 1, function(x){x%*%male_X_mean}))
  n <- length(lambda1)
  z <- rbinom(n, size = 1,prob = rho_samps)
  u_effect <- exp(z*(rexp(n, 1/lambda0) - rexp(n, 1/lambda1)))
  u_effect_male <- exp(z*(rexp(n, 1/lambda0_male) - rexp(n, 1/lambda1_male)))
  u_effect_female <- exp(z*(rexp(n, 1/lambda0_female) - rexp(n, 1/lambda1_female)))
  
  full_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  full_X_mean[c("treated", "head_sex", "year", 
                "treated_year", "treated_sex", "year_sex", 
                "treated_year_sex")] <- c(1, full_X_mean["head_sex"], 1,
                                          1*1, 1*full_X_mean["head_sex"],full_X_mean["head_sex"]*1,
                                          full_X_mean["head_sex"]*1*1)
  temp <- exp(rho_samps*lambda1*(exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%full_X_mean[c(5,8)]})) -1))
  
  
  }



