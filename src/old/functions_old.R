
#Gelman standardization
gel_stand <- function(x){(x - mean(x))/(2*sd(x))}

#Inverse hyperbolic sine transformation
ihs_trans <- function(x){log(x + sqrt(x^2 + 1))}

#In terms of standard deviation of un-treated
y_stand <- function(x, data){
  x2 <-x[data$treated == 0]
  (x - mean(x2,na.rm = TRUE))/(sd(x2,na.rm = TRUE))
}

##########do_sampling()-------------
#y = column of responses (including NA for missings, pre-transformed if needed)
#X = matrix of predictor variables, pre-standardized
#hh_id = household id vector (should be same dimension as y)
#loc_id = location id vector ("                           ")

do_sampling <- function(y, X, X_cntr, hh_id, loc_id,kappa,file){
  #compile BSFA model incorporating missingness mechanism 
  compiled_selection_model <- stan_model(file)
  
  f <- y
  index_miss <- which(is.na(f))
  N_miss <- length(index_miss)
  #Missingness - empirical Bayes
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
  f[index_miss] <-imputed_f[index_miss] #-99
  todo = which(X[,"treated"] == 1 & X[,"year"] == 1)
  data.list <- list(N = length(f), 
                    N_subset = length(todo),
                    subset_indices = todo,
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
                    kappa = kappa)
  #set reasonable initial values
  initc <- list(gamma1 = -1, 
                gamma2 = rep(0, 7),
                mu_a = coef(m1)[1], 
                beta = coef(m1)[-1], 
                sigma_a = .5*mean(m1$residuals^2)^(.5), sigma_a2 = .5*mean(m1$residuals^2),
                sigma_v = .5*mean(m1$residuals^2)^(.5), sigma_v2 = .5*mean(m1$residuals^2),
                sigma_l = .5*mean(m1$residuals^2)^(.5), sigma_l2 =.5*mean(m1$residuals^2),
                u =.5*(max(imputed_f)-imputed_f) + .01, #approx. 50 percent of distance from  max
                rho = .6, 
                f_imp = rep(mean(imputed_f), N_miss),
                phi = phi_prior_mean)
  
  inits <- list(initc,initc,initc,initc)  
  
  #Sample with rstan
  sm_sampled <- sampling(compiled_selection_model, 
                         data = data.list, 
                         chains = 4,
                         iter =15000,
                         warmup = 5000,
                         pars = par_keep,
                         include=TRUE,
                         init = inits,
                         control = list(adapt_delta = 0.95)) 
  return(sm_sampled)
}  #Sample with rstan
  #sm_sampled <- sampling(compiled_selection_model, 
  ##                       data = data.list, 
  #                       chains = 2,
  #                       iter =1000,
  #                       warmup = 500,
  #                       pars = par_keep,
  #                       include=TRUE,
  #                       init = inits,
  #                       control = list(adapt_delta = 0.95)) 

##########plot_cs()-----------
#sampled = stanfit object
#y = response variable - transformed if applicable
#response = character string of y variable (just for labels)
#data = kenya data frame
#backtrans = logical TRUE if on original (untransformed) scale

plot_cs <- function(sampled, y, response, data, backtrans = FALSE){
  #draws from posterior predictive (CS)
  f_pred_trt <-  rstan::extract(sampled, "f_pred")%>%
    data.frame() 
  #draws from counterfactual posterior predictive (CS)
  f_pred_cntr <-  rstan::extract(sampled, "f_pred_cntr")%>%
    data.frame()
  if (backtrans == TRUE) {
    f_pred_trt <- exp(f_pred_trt)
    f_pred_cntr<- exp(f_pred_cntr)
    y <- exp(y)
  }
  
  
  todo = which(data[,"treated"] == 1 & data[,"year"] == 1)
  d <-  data[todo,]
  
  data.frame(f =y[todo],
             missing = is.na(y[todo]),
             sex = factor(d$head_sex, levels = c(0,1), labels = c("Male", "Female")),
             lb_trt =   apply(f_pred_trt , 2, "quantile", c(0.05)),
             ub_trt =   apply(f_pred_trt , 2, "quantile", c(0.95)),
             lb_notrt = apply(f_pred_cntr, 2, "quantile", c(0.05)),
             ub_notrt = apply(f_pred_cntr, 2, "quantile", c(0.95))) %>%
    arrange(sex, ub_notrt)%>%
    mutate(id = c(1:sum(d$head_sex == 0), 1:sum(d$head_sex == 1))) %>%
    subset(missing == FALSE)%>%
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

#########plot_cs_agg()-------
#same arguments as plot_cs()

plot_cs_agg <- function(sampled, y, response, data,  backtrans = FALSE){
  
  f_pred_trt <-  rstan::extract(sampled, "f_pred")%>%
    data.frame() 
  f_pred_trt <- f_pred_trt[,data$year == 1 & data$treated == 1]
  f_pred_cntr <-  rstan::extract(sampled, "f_pred_cntr")%>%
    data.frame()
  f_pred_cntr <- f_pred_cntr[,data$year == 1 & data$treated == 1]
  
  if (backtrans == TRUE) {
    f_pred_trt <- exp(f_pred_trt)
    f_pred_cntr<- exp(f_pred_cntr)
    y <- exp(y)
  }
  
  
  
  d <-  data %>%
    subset(year == 1 & treated == 1) 
  
  data.frame(f =y[data$year == 1 & data$treated == 1],
             missing = is.na(y[data$year == 1 & data$treated == 1]),
             sex = factor(d$head_sex, levels = c(0,1), labels = c("Male", "Female")),
             lb_trt =   apply(f_pred_trt , 2, "quantile", c(0.05)),
             ub_trt =   apply(f_pred_trt , 2, "quantile", c(0.95)),
             lb_notrt = apply(f_pred_cntr, 2, "quantile", c(0.05)),
             ub_notrt = apply(f_pred_cntr, 2, "quantile", c(0.95))) %>%
    arrange( ub_notrt)%>%
    mutate(id = c(1:length(d$head_sex))) %>%
    subset(missing == FALSE)%>%
    ggplot() +
    geom_errorbar(aes(x=id, ymin=lb_trt, ymax=ub_trt), colour = "grey50") +
    geom_line(aes(x=id, y=lb_notrt  )) +
    geom_line(aes(x=id, y=ub_notrt )) +
    labs(x = "Household ID", y = response)+
    geom_point(aes(x = id, y =f),size = .5,alpha = I(.3), colour  ="black")+
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_grey() #+scale_y_continuous(limits = c(-4, 4))
  
}

#########derivatives() : posterior draws from quantities in Eq (18)-------
#full_X: X matrix
#samples = stanfit object

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
  
  
  #CHOICE effect
  #overall effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%full_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%full_X_mean[c(5,8)]})) # 1/exp(g2 + g3*.6244)
  dudb <- temp1*(temp2-1)*rho_samps

  #female effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%female_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%female_X_mean[c(5,8)]})) # 1/exp(g2 + g3*1)
  dudb_female <- temp1*(temp2-1)*rho_samps
  
  #male effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%male_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%male_X_mean[c(5,8)]})) # 1/exp(g2 + g3*0)
  dudb_male <- temp1*(temp2-1)*rho_samps

  rm(gamma_samps)
  rm(rho_samps)
  
  #CAPABILITIES effect
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
               
               data.frame(effect = dudb, type = "Choice", by = "Overall"), 
               data.frame(effect = dudb_female, type = "Choice", by = "Female"), 
               data.frame(effect = dudb_male, type = "Choice", by = "Male"), 
               
               data.frame(effect = dcdb, type = "Capabilities", by = "Overall"), 
               data.frame(effect = dcdb_female, type = "Capabilities", by = "Female"), 
               data.frame(effect = dcdb_male, type = "Capabilities", by = "Male"))
  return(all)
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



