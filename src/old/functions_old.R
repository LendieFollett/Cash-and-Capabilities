
#Gelman standardization
gel_stand <- function(x){(x - mean(x))/(2*sd(x))}

#Inverse hyperbolic sine transformation
ihs_trans <- function(x){log(x + sqrt(x^2 + 1))}
hs_trans <- function(x){(exp(x) - exp(-x))/2}

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
                         control = list(adapt_delta = 0.9)) 
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
  if (backtrans == "log") {
    f_pred_trt <- exp(f_pred_trt)
    f_pred_cntr<- exp(f_pred_cntr)
    y <- exp(y)
  } else if(backtrans == "IHS"){
    f_pred_trt <- hs_trans(f_pred_trt)
    f_pred_cntr<- hs_trans(f_pred_cntr)
    y <- hs_trans(y) 
  }
  
  todo = which(data[,"treated"] == 1 & data[,"year"] == 1)
  d <-  data[todo,]
  
  data.frame(f =y[todo],
             missing = is.na(y[todo]),
             age = factor(d$caregiver, levels = c(0,1), labels = c("Younger", "Older")),
             lb_trt =   apply(f_pred_trt , 2, "quantile", c(0.05)),
             ub_trt =   apply(f_pred_trt , 2, "quantile", c(0.95)),
             lb_notrt = apply(f_pred_cntr, 2, "quantile", c(0.05)),
             ub_notrt = apply(f_pred_cntr, 2, "quantile", c(0.95))) %>%
    arrange(age, ub_notrt)%>%
    mutate(id = c(1:sum(d$caregiver == 0), 1:sum(d$caregiver == 1))) %>%
    subset(missing == FALSE)%>%
    ggplot() +
    geom_errorbar(aes(x=id, ymin=lb_trt, ymax=ub_trt), colour = "grey50") +
    geom_line(aes(x=id, y=lb_notrt  )) +
    geom_line(aes(x=id, y=ub_notrt )) +
    labs(x = "Household ID", y = response)+
    geom_point(aes(x = id, y =f),size = .5,alpha = I(.3), colour  ="black")+
    facet_wrap(.~age, scales = "free_x")+
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_grey() + scale_y_continuous(labels = label_comma())
  
}


#########plot_cs_agg()-------
#same arguments as plot_cs()

plot_cs_agg <- function(sampled, y, response, data,  backtrans = FALSE){
  
  #draws from posterior predictive (CS)
  f_pred_trt <-  rstan::extract(sampled, "f_pred")%>%
    data.frame() 
  #draws from counterfactual posterior predictive (CS)
  f_pred_cntr <-  rstan::extract(sampled, "f_pred_cntr")%>%
    data.frame()
  if (backtrans == "log") {
    f_pred_trt <- exp(f_pred_trt)
    f_pred_cntr<- exp(f_pred_cntr)
    y <- exp(y)
  } else if(backtrans == "IHS"){
    f_pred_trt <- hs_trans(f_pred_trt)
    f_pred_cntr<- hs_trans(f_pred_cntr)
    y <- hs_trans(y) 
  }
  
  
  todo = which(data[,"treated"] == 1 & data[,"year"] == 1)
  d <-  data[todo,]
  
  data.frame(f =y[todo],
             missing = is.na(y[todo]),
             age = factor(d$caregiver, levels = c(0,1), labels = c("Younger", "Older")),
             lb_trt =   apply(f_pred_trt , 2, "quantile", c(0.05)),
             ub_trt =   apply(f_pred_trt , 2, "quantile", c(0.95)),
             lb_notrt = apply(f_pred_cntr, 2, "quantile", c(0.05)),
             ub_notrt = apply(f_pred_cntr, 2, "quantile", c(0.95))) %>%
    arrange( ub_notrt)%>%
    mutate(id = c(1:length(d$caregiver))) %>%
    subset(missing == FALSE)%>%
    ggplot() +
    geom_errorbar(aes(x=id, ymin=lb_trt, ymax=ub_trt),colour = "grey50") +#colour = "grey50"
    geom_line(aes(x=id, y=lb_notrt  )) +
    geom_line(aes(x=id, y=ub_notrt )) +
    labs(x = "Household ID", y = response)+
    geom_point(aes(x = id, y =f),size = .5,alpha = I(.3), colour  ="black")+
    scale_colour_grey() + scale_y_continuous(labels = label_comma())
  

}

plot_cs_agg_numeric <- function(sampled, y, response, data,  backtrans = FALSE){
  
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
             head_age = d$head_age,
             lb_trt =   apply(f_pred_trt , 2, "quantile", c(0.05)),
             ub_trt =   apply(f_pred_trt , 2, "quantile", c(0.95)),
             lb_notrt = apply(f_pred_cntr, 2, "quantile", c(0.05)),
             ub_notrt = apply(f_pred_cntr, 2, "quantile", c(0.95))) %>%
    arrange( ub_notrt)%>% #changed from ub_notrt
    mutate(id = c(1:length(d$head_age))) %>%
    subset(missing == FALSE)%>%
    ggplot() +
    geom_errorbar(aes(x=id, ymin=lb_trt, ymax=ub_trt, colour = ub_notrt)) +
    geom_line(aes(x=id, y=lb_notrt  )) +
    geom_line(aes(x=id, y=ub_notrt )) +
    labs(x = "Household ID", y = response)+
    geom_point(aes(x = id, y =f),size = .5,alpha = I(.3), colour  ="black")+
    scale_colour_grey() + scale_y_continuous(labels = label_comma())
  
  
}

#########derivatives() : posterior draws from quantities in Eq (18)-------
#full_X: X matrix
#samples = stanfit object

derivatives <- function(full_X, samples){ 
  full_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  full_X_mean[c("treated", "caregiver", "year", 
                "treated_year", "treated_age", "year_age", 
                "treated_year_age")] <- c(1, full_X_mean["caregiver"], 1,
                                          1*1, 1*full_X_mean["caregiver"],full_X_mean["caregiver"]*1,
                                          full_X_mean["caregiver"]*1*1)
  
  older_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  older_X_mean[c("treated", "caregiver", "year", 
                 "treated_year", "treated_age", "year_age", 
                 "treated_year_age")] <- c(1, 1, 1,
                                            1*1, 1*1,1*1,
                                            1*1*1)
  
  younger_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  younger_X_mean[c("treated", "caregiver", "year", 
                   "treated_year", "treated_age", "year_age", 
                   "treated_year_age")] <- c(1, 0, 1,
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
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%older_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%older_X_mean[c(5,8)]})) # 1/exp(g2 + g3*1)
  dudb_older <- temp1*(temp2-1)*rho_samps
  
  #male effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%younger_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%younger_X_mean[c(5,8)]})) # 1/exp(g2 + g3*0)
  dudb_younger <- temp1*(temp2-1)*rho_samps

  rm(gamma_samps)
  rm(rho_samps)
  
  #CAPABILITIES effect
  beta_samps <- do.call("cbind", rstan::extract(samples, par = c("beta")))
  #overall effect
  dcdb <- beta_samps[,4] + beta_samps[,7]*full_X_mean["caregiver"]
  #older effect
  dcdb_older <- beta_samps[,4] + beta_samps[,7]*1
  #younger effect
  dcdb_younger <- beta_samps[,4] + beta_samps[,7]*0
  rm(beta_samps)
  
  #FUNCTIONING effect 
  dydb <-        dcdb        + dudb
  dydb_older <- dcdb_older + dudb_older
  dydb_younger <-   dcdb_younger   + dudb_younger
  
  
  all <- rbind(data.frame(effect = dydb, type = "Functioning", by = "Overall"), 
               data.frame(effect = dydb_older, type = "Functioning", by = "Older"), 
               data.frame(effect = dydb_younger, type = "Functioning", by = "Younger"), 
               
               data.frame(effect = dudb, type = "Choice", by = "Overall"), 
               data.frame(effect = dudb_older, type = "Choice", by = "Older"), 
               data.frame(effect = dudb_younger, type = "Choice", by = "Younger"), 
               
               data.frame(effect = dcdb, type = "Capabilities", by = "Overall"), 
               data.frame(effect = dcdb_older, type = "Capabilities", by = "Older"), 
               data.frame(effect = dcdb_younger, type = "Capabilities", by = "Younger"))
  return(all)
}



derivatives_numeric <- function(full_X, samples){ 
  full_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
  full_X_mean[c("treated", "head_age", "year", 
                "treated_year", "treated_age", "year_age", 
                "treated_year_age")] <- c(1, full_X_mean["head_age"], 1,
                                          1*1, 1*full_X_mean["head_age"],full_X_mean["head_age"]*1,
                                          full_X_mean["head_age"]*1*1)
 
  gamma_samps <- do.call("cbind", rstan::extract(samples, par = c("gamma1","gamma2")))
  rho_samps <- rstan::extract(samples, par = c("rho"))[[1]]
  beta_samps <- do.call("cbind", rstan::extract(samples, par = c("beta")))
  
  all <- list()
  j = 0
for (a in seq(from = -1, to = 1, by = .25)) {  
  j = j + 1
  print(a)
age_X_mean <- c(1,apply(full_X, 2, mean))[1:8]
age_X_mean[c("treated", "head_age", "year", 
                 "treated_year", "treated_age", "year_age", 
                 "treated_year_age")] <- c(1, a, 1,
                                           1*1, 1*a,a*1,
                                           a*1*1)
  
  #CHOICE effect
  temp1 <- apply(gamma_samps, 1, function(x){exp(x%*%age_X_mean)}) #exp(X*gamma) = lambda(1)
  temp2 <- exp(-apply(gamma_samps[,c(5,8)], 1, function(x){x%*%age_X_mean[c(5,8)]})) # 1/exp(g2 + g3*0)
  dudb_younger <- temp1*(temp2-1)*rho_samps
  
  #CAPABILITIES effect
  #younger effect
  dcdb_younger <- beta_samps[,4] + beta_samps[,7]*a
  
  #FUNCTIONING effect 
  dydb_younger <-   dcdb_younger   + dudb_younger
  
  
  all[[j]] <- rbind(
               data.frame(effect = dydb_younger, type = "Functioning", age = a), 

               data.frame(effect = dudb_younger, type = "Choice", age = a), 
               
               data.frame(effect = dcdb_younger, type = "Capabilities", age = a))
  
}
  
  all <-do.call(rbind, all) 
  
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



