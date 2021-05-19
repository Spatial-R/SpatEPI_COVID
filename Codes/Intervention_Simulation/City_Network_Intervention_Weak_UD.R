library(dplyr)
library(pomp)
library(reshape2)
library(spatPomp)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggsci)
library(stringi)
library(tidyr)
library(ggrepel)
library(readr)

load("Process_Data/Cities_Spatial.RData")
source("Codes/City_Network_Sim_Function.R")
time_range <- seq(7, 150, by = 1)

measles_rprocess <- Csnippet('
                             double beta,foi;
                             // double dw,sus_del;
                             double rate[8], trans[8];
                             double *S = &S1;
                             double *E = &E1;
                             double *I = &I1;
                             double *R = &R1;
                             double *J = &J1;
                             double *P = &P1;
                             double *C = &C1;
                             double powVec[U];
                             const double *pop = &pop1;
                             //int obstime = 0;
                             int u,v;
                             
                             int kl = floor(t);
                             
                             // transmission rate
                             // beta = R0*(gamma+mu);
                             
                             // pre-computing this saves substantial time
                             
                             for (u = 0 ; u < U ; u++) {
                             powVec[u] = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                             }
                             
                             for (u = 0 ; u < U ; u++) {
                             
                             if(kl < 15){
                             beta = R0IN*gamma/(1-omega+omega*kappa);
                             } else {
                             beta = R0IN*gamma/(1-omega+omega*kappa)*contreduce[u];
                             } 
                             
                             
                             // expected force of infection
                             foi = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                             
                             
                             for (v=0; v < U ; v++) {
                             if(v != u){
                             if(kl > 15){
                             foi += g * (v_by_g[v+U][u] * powVec[v] - (v_by_g[u+U][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                             } else {
                             foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                             }
                             }
                             }
                             if(foi < 0) foi = 0;
                             // white noise (extrademographic stochasticity)
                             
                             double qm;
                             
                             if(kl > 15){
                             kl = 15;
                             }                             
                             
                             if (u == (targetp-1)){
                             qm = 1- (exp(quawh*kl+log(qualow))/(exp(quawh*kl+log(qualow))+1));
                             } else {
                             qm = 1- (exp(qua*kl+log(qualow))/(exp(qua*kl+log(qualow))+1));
                             }
                             
                             rate[0] = beta*foi;
                             rate[1] = sigma*(1-qm)*omega;    // hOSP
                             rate[2] = sigma*(qm)*omega;       // HOUSE
                             rate[3] = sigma*(1-omega);
                             rate[4] = gamma;  // for the exposed to quarten 
                             rate[5] = gamma;		       // recovery for quarten   
                             rate[6] = gamma;		       // recovery for quarten                     
                             
                             // transitions between classes
                             
                             trans[0]  = (S[u]*rate[0]*dt);
                             trans[1]  = (E[u]*rate[1]*dt);
                             trans[2]  = (E[u]*rate[2]*dt);
                             trans[3]  = (E[u]*rate[3]*dt);
                             trans[4]  = (J[u]*rate[4]*dt);
                             trans[5]  = (P[u]*rate[5]*dt);
                             trans[6]  = (I[u]*rate[6]*dt);
                             
                             S[u] += -trans[0];
                             E[u] +=  trans[0] - trans[1] - trans[2] - trans[3];
                             J[u] +=  trans[1] - trans[4];
                             P[u] +=  trans[2] - trans[5];
                             I[u] +=  trans[3] - trans[6];
                             R[u] +=  trans[4] + trans[5] + trans[6];
                             C[u] += trans[5]+trans[6];           // death
                             
                             if(S[u] < 0) S[u] = 0;
                             if(E[u] < 0) E[u] = 0;
                             if(J[u] < 0) J[u] = 0;
                             if(I[u] < 0) I[u] = 0;
                             if(P[u] < 0) P[u] = 0;
                             if(R[u] < 0) R[u] = 0;
                             }
                             ')


case_all_1 <- data.frame(expand.grid(day = time_range,city = sort(unique(case_all$city))))
case_all_1$cases <- NA; case_all_1$city <- as.character(case_all_1$city)

pop_list <- lapply(sort(unique(case_all$city)),function(data){
  dat_tem <- filter(pop_all,city == data)
  dat_tem <- data.frame(day = time_range, city = data, pop = dat_tem[1,"pop"])
  return(dat_tem)
})
pop_all <- bind_rows(pop_list)


spatPomp(case_all_1,
         units = "city",
         times = "day",
         t0 = min(case_all$day)-1,
         unit_statenames = measles_unit_statenames,
         covar = pop_all,
         tcovar = "day",
         rprocess=euler(measles_rprocess, delta.t=dt),
         # skeleton=vectorfield(measles_skel),
         accumvars = c(paste0("C",1:U)),
         paramnames=measles_paramnames,
         covarnames=measles_covarnames,
         globals=measles_globals,
         rinit=measles_rinit,
         partrans=parameter_trans(log=c("alpha","g","R0IN","kappa"),
                                  logit = c("psi","qua","sigma","gamma","quawh",
                                            "qualow","rho","omega",
                                            paste0("S",1:U,"_0"),
                                            paste0("E",1:U,"_0"),
                                            paste0("J",1:U,"_0"),
                                            paste0("P",1:U,"_0"),
                                            paste0("I",1:U,"_0"),
                                            paste0("R",1:U,"_0"))),
         dmeasure=measles_dmeasure,
         unit_emeasure=measles_unit_emeasure,
         unit_mmeasure=measles_unit_mmeasure,
         unit_vmeasure=measles_unit_vmeasure,
         rmeasure=measles_rmeasure,
         unit_dmeasure=measles_unit_dmeasure) -> origin_model


s1_pos <- which(names(parameter) == "S1_0"); 
lgh_pos <- which(names(parameter) == "loglik")

estimate_p <- unlist(lapply(measles_unit_statenames[1:(length(measles_unit_statenames)-1)], function(data){
  test_p <- parameter[(s1_pos+which(substr(names(parameter)[s1_pos:lgh_pos],1,1) == data)-1)]
  names(test_p) <- paste(data,test_col,"_0",sep = "")
  return(test_p)
}))


fixed_p <- unlist(lapply(measles_unit_statenames[1:(length(measles_unit_statenames)-1)], function(data){
  fix_col <- c(1:ncol(flow_all))[-c(test_col)]
  if(!length(fix_col) == 0){
    if(data == "S"){
      test_p <- rep(1,length(fix_col))
    } else {
      test_p <- rep(1e-8,length(fix_col))
    }
    names(test_p) <- paste(data,fix_col,"_0",sep = "")
  } else {
    test_p <- NULL
  }
  return(test_p)
}))

parameter_all <- unlist(c(fixed_p,estimate_p))
parameter_all <- parameter_all[unlist(lapply(measles_IVPnames, function(data)which(names(parameter_all) == data)))]
parameter_all <- c(parameter_all,unlist(c(parameter[1:(s1_pos-1)])))


dat_origin <-  simulated_timeseries(model = origin_model,parameter = parameter_all,
                                    case_origin = case_all_1)

dat_origin_day <- data.frame(summarise(group_by(dat_origin,day),
                                       sim = floor(sum(median)),low = floor(sum(low)),high = floor(sum(high))))


###################################################################################################
################################ if wuhan have not block down   ###################################
###################################################################################################

measles_rprocess_fc <- Csnippet('
                                double beta,foi;
                                // double dw,sus_del;
                                double rate[8], trans[8];
                                double *S = &S1;
                                double *E = &E1;
                                double *I = &I1;
                                double *R = &R1;
                                double *J = &J1;
                                double *P = &P1;
                                double *C = &C1;
                                double powVec[U];
                                const double *pop = &pop1;
                                //int obstime = 0;
                                int u,v;
                                
                                int kl = floor(t);
                                
                                // transmission rate
                                // beta = R0*(gamma+mu);
                                
                                // pre-computing this saves substantial time
                                
                                for (u = 0 ; u < U ; u++) {
                                powVec[u] = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                                }
                                
                                for (u = 0 ; u < U ; u++) {
                                
                                if(kl < 15){
                                beta = R0IN*gamma/(1-omega+omega*kappa);
                                } else {
                                beta = R0IN*gamma/(1-omega+omega*kappa)*contreduce[u];
                                } 
                                
                                // expected force of infection
                                foi = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                                
                                
                                for (v=0; v < U ; v++) {
                                if(v != u){
                                if(kl > 15){
                                foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                                } else {
                                foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                                }
                                }
                                }
                                if(foi < 0) foi = 0;
                                // white noise (extrademographic stochasticity)
                                
                                double qm;
                                
                                if(kl > 15){
                                kl = 15;
                                }      
                                
                                if (u == (targetp-1)){
                                qm = 1- (exp(quawh*kl+log(qualow))/(exp(quawh*kl+log(qualow))+1));
                                } else {
                                qm = 1- (exp(qua*kl+log(qualow))/(exp(qua*kl+log(qualow))+1));
                                }
                                
                                rate[0] = beta*foi;
                                rate[1] = sigma*(1-qm)*omega;    // hOSP
                                rate[2] = sigma*(qm)*omega;       // HOUSE
                                rate[3] = sigma*(1-omega);
                                rate[4] = gamma;  // for the exposed to quarten 
                                rate[5] = gamma;		       // recovery for quarten   
                                rate[6] = gamma;		       // recovery for quarten                     
                                
                                // transitions between classes
                                
                                trans[0]  = (S[u]*rate[0]*dt);
                                trans[1]  = (E[u]*rate[1]*dt);
                                trans[2]  = (E[u]*rate[2]*dt);
                                trans[3]  = (E[u]*rate[3]*dt);
                                trans[4]  = (J[u]*rate[4]*dt);
                                trans[5]  = (P[u]*rate[5]*dt);
                                trans[6]  = (I[u]*rate[6]*dt);
                                
                                S[u] += -trans[0];
                                E[u] +=  trans[0] - trans[1] - trans[2] - trans[3];
                                J[u] +=  trans[1] - trans[4];
                                P[u] +=  trans[2] - trans[5];
                                I[u] +=  trans[3] - trans[6];
                                R[u] +=  trans[4] + trans[5] + trans[6];
                                C[u] += trans[5]+trans[6];           // death
                                
                                if(S[u] < 0) S[u] = 0;
                                if(E[u] < 0) E[u] = 0;
                                if(J[u] < 0) J[u] = 0;
                                if(I[u] < 0) I[u] = 0;
                                if(P[u] < 0) P[u] = 0;
                                if(R[u] < 0) R[u] = 0;
                                }
                                ')


spatPomp(case_all_1,
         units = "city",
         times = "day",
         t0 = min(case_all$day)-1,
         unit_statenames = measles_unit_statenames,
         covar = pop_all,
         tcovar = "day",
         rprocess=euler(measles_rprocess_fc, delta.t=dt),
         # skeleton=vectorfield(measles_skel),
         accumvars = c(paste0("C",1:U)),
         paramnames=measles_paramnames,
         covarnames=measles_covarnames,
         globals=measles_globals,
         rinit=measles_rinit,
         partrans=parameter_trans(log=c("alpha","g","R0IN"),
                                  logit = c("psi","qua","sigma","gamma","quawh","qualow","rho","omega","kappa",
                                            paste0("S",1:U,"_0"),
                                            paste0("E",1:U,"_0"),
                                            paste0("J",1:U,"_0"),
                                            paste0("P",1:U,"_0"),
                                            paste0("I",1:U,"_0"),
                                            paste0("R",1:U,"_0"))),
         dmeasure=measles_dmeasure,
         unit_emeasure=measles_unit_emeasure,
         unit_mmeasure=measles_unit_mmeasure,
         unit_vmeasure=measles_unit_vmeasure,
         rmeasure=measles_rmeasure,
         unit_dmeasure=measles_unit_dmeasure) -> origin_model_fc



Sys.setlocale("LC_TIME", "English") # Windows


dat_fc <-  simulated_timeseries(model = origin_model_fc,parameter = parameter_all,
                                case_origin = case_all_1)

dat_fc_day <- data.frame(summarise(group_by(dat_fc,day),
                                   sim = floor(sum(median)),low = floor(sum(low)),high = floor(sum(high))))


####################################  fengcheng + contact  ####################################

measles_rprocess_fcct <- Csnippet('
                                  double beta,foi;
                                  // double dw,sus_del;
                                  double rate[8], trans[8];
                                  double *S = &S1;
                                  double *E = &E1;
                                  double *I = &I1;
                                  double *R = &R1;
                                  double *J = &J1;
                                  double *P = &P1;
                                  double *C = &C1;
                                  double powVec[U];
                                  const double *pop = &pop1;
                                  //int obstime = 0;
                                  int u,v;
                                  
                                  int kl = floor(t);
                                  
                                  // transmission rate
                                  // beta = R0*(gamma+mu);
                                  
                                  // pre-computing this saves substantial time
                                  
                                  for (u = 0 ; u < U ; u++) {
                                  powVec[u] = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                                  }
                                  
                                  for (u = 0 ; u < U ; u++) {
                                  
                                  if(kl < 15){
                                  beta = R0IN*gamma/(1-omega+omega*kappa);
                                  } else {
                                  beta = R0IN*gamma/(1-omega+omega*kappa);
                                  } 
                                  
                                  // expected force of infection
                                  foi = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                                  
                                  
                                  for (v=0; v < U ; v++) {
                                  if(v != u){
                                  if(kl > 15){
                                  foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                                  } else {
                                  foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                                  }
                                  }
                                  }
                                  if(foi < 0) foi = 0;
                                  // white noise (extrademographic stochasticity)
                                  
                                  double qm;
                                  
                                  if(kl > 15){
                                  kl = 15;
                                  }

                                  if (u == (targetp-1)){
                                  qm = 1- (exp(quawh*kl+log(qualow))/(exp(quawh*kl+log(qualow))+1));
                                  } else {
                                  qm = 1- (exp(qua*kl+log(qualow))/(exp(qua*kl+log(qualow))+1));
                                  }
                                  
                                  rate[0] = beta*foi;
                                  rate[1] = sigma*(1-qm)*omega;    // hOSP
                                  rate[2] = sigma*(qm)*omega;       // HOUSE
                                  rate[3] = sigma*(1-omega);
                                  rate[4] = gamma;  // for the exposed to quarten 
                                  rate[5] = gamma;		       // recovery for quarten   
                                  rate[6] = gamma;		       // recovery for quarten                     
                                  
                                  // transitions between classes
                                  
                                  trans[0]  = (S[u]*rate[0]*dt);
                                  trans[1]  = (E[u]*rate[1]*dt);
                                  trans[2]  = (E[u]*rate[2]*dt);
                                  trans[3]  = (E[u]*rate[3]*dt);
                                  trans[4]  = (J[u]*rate[4]*dt);
                                  trans[5]  = (P[u]*rate[5]*dt);
                                  trans[6]  = (I[u]*rate[6]*dt);
                                  
                                  S[u] += -trans[0];
                                  E[u] +=  trans[0] - trans[1] - trans[2] - trans[3];
                                  J[u] +=  trans[1] - trans[4];
                                  P[u] +=  trans[2] - trans[5];
                                  I[u] +=  trans[3] - trans[6];
                                  R[u] +=  trans[4] + trans[5] + trans[6];
                                  C[u] += trans[6]+trans[5];           // true incidence
                                  
                                  if(S[u] < 0) S[u] = 0;
                                  if(E[u] < 0) E[u] = 0;
                                  if(J[u] < 0) J[u] = 0;
                                  if(I[u] < 0) I[u] = 0;
                                  if(P[u] < 0) P[u] = 0;
                                  if(R[u] < 0) R[u] = 0;
                                  }
                                  ')


spatPomp(case_all_1,
         units = "city",
         times = "day",
         t0 = min(case_all$day)-1,
         unit_statenames = measles_unit_statenames,
         covar = pop_all,
         tcovar = "day",
         rprocess=euler(measles_rprocess_fcct, delta.t=dt),
         # skeleton=vectorfield(measles_skel),
         accumvars = c(paste0("C",1:U)),
         paramnames=measles_paramnames,
         covarnames=measles_covarnames,
         globals=measles_globals,
         rinit=measles_rinit,
         partrans=parameter_trans(log=c("alpha","g","R0IN","kappa"),
                                  logit = c("psi","qua","sigma","gamma","quawh",
                                            "qualow","rho","omega",
                                            paste0("S",1:U,"_0"),
                                            paste0("E",1:U,"_0"),
                                            paste0("J",1:U,"_0"),
                                            paste0("P",1:U,"_0"),
                                            paste0("I",1:U,"_0"),
                                            paste0("R",1:U,"_0"))),
         dmeasure=measles_dmeasure,
         unit_emeasure=measles_unit_emeasure,
         unit_mmeasure=measles_unit_mmeasure,
         unit_vmeasure=measles_unit_vmeasure,
         rmeasure=measles_rmeasure,
         unit_dmeasure=measles_unit_dmeasure) -> origin_model_fcct


dat_fcct <-  simulated_timeseries(model = origin_model_fcct,parameter = parameter_all,
                                  case_origin = case_all_1)


dat_fcct_day <- data.frame(summarise(group_by(dat_fcct,day),
                                     sim = floor(sum(median)),low = floor(sum(low)),
                                     high = floor(sum(high))))


####################################  fengcheng + contact  ####################################



measles_rprocess_ct <- Csnippet('
                                double beta,foi;
                                // double dw,sus_del;
                                double rate[8], trans[8];
                                double *S = &S1;
                                double *E = &E1;
                                double *I = &I1;
                                double *R = &R1;
                                double *J = &J1;
                                double *P = &P1;
                                double *C = &C1;
                                double powVec[U];
                                const double *pop = &pop1;
                                //int obstime = 0;
                                int u,v;
                                
                                int kl = floor(t);
                                
                                // transmission rate
                                // beta = R0*(gamma+mu);
                                
                                // pre-computing this saves substantial time
                                
                                for (u = 0 ; u < U ; u++) {
                                powVec[u] = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                                }
                                
                                for (u = 0 ; u < U ; u++) {
                                
                                if(kl < 15){
                                beta = R0IN*gamma/(1-omega+omega*kappa);
                                } else {
                                beta = R0IN*gamma/(1-omega+omega*kappa);
                                } 
                                
                                // expected force of infection
                                foi = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                                
                                
                                for (v=0; v < U ; v++) {
                                if(v != u){
                                if(kl > 15){
                                foi += g * (v_by_g[v+U][u] * powVec[v] - (v_by_g[u+U][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                                } else {
                                foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                                }
                                }
                                }
                                if(foi < 0) foi = 0;
                                // white noise (extrademographic stochasticity)
                                
                                double qm;
                                
                                if(kl > 15){
                                kl = 15;
                                }      
                                
                                if (u == (targetp-1)){
                                qm = 1- (exp(quawh*kl+log(qualow))/(exp(quawh*kl+log(qualow))+1));
                                } else {
                                qm = 1- (exp(qua*kl+log(qualow))/(exp(qua*kl+log(qualow))+1));
                                }
                                
                                rate[0] = beta*foi;
                                rate[1] = sigma*(1-qm)*omega;    // hOSP
                                rate[2] = sigma*(qm)*omega;       // HOUSE
                                rate[3] = sigma*(1-omega);
                                rate[4] = gamma;  // for the exposed to quarten 
                                rate[5] = gamma;		       // recovery for quarten   
                                rate[6] = gamma;		       // recovery for quarten                     
                                
                                // transitions between classes
                                
                                trans[0]  = (S[u]*rate[0]*dt);
                                trans[1]  = (E[u]*rate[1]*dt);
                                trans[2]  = (E[u]*rate[2]*dt);
                                trans[3]  = (E[u]*rate[3]*dt);
                                trans[4]  = (J[u]*rate[4]*dt);
                                trans[5]  = (P[u]*rate[5]*dt);
                                trans[6]  = (I[u]*rate[6]*dt);
                                
                                S[u] += -trans[0];
                                E[u] +=  trans[0] - trans[1] - trans[2] - trans[3];
                                J[u] +=  trans[1] - trans[4];
                                P[u] +=  trans[2] - trans[5];
                                I[u] +=  trans[3] - trans[6];
                                R[u] +=  trans[4] + trans[5] + trans[6];
                                C[u] += trans[5]+trans[6];           // death
                                
                                if(S[u] < 0) S[u] = 0;
                                if(E[u] < 0) E[u] = 0;
                                if(J[u] < 0) J[u] = 0;
                                if(I[u] < 0) I[u] = 0;
                                if(P[u] < 0) P[u] = 0;
                                if(R[u] < 0) R[u] = 0;
                                }
                                ')


spatPomp(case_all_1,
         units = "city",
         times = "day",
         t0 = min(case_all_1$day)-1,
         unit_statenames = measles_unit_statenames,
         covar = pop_all,
         tcovar = "day",
         rprocess=euler(measles_rprocess_ct, delta.t=dt),
         # skeleton=vectorfield(measles_skel),
         accumvars = c(paste0("C",1:U)),
         paramnames=measles_paramnames,
         covarnames=measles_covarnames,
         globals=measles_globals,
         rinit=measles_rinit,
         partrans=parameter_trans(log=c("alpha","g","R0IN","kappa"),
                                  logit = c("psi","qua","sigma","gamma",
                                            "quawh","qualow","rho","omega",
                                            paste0("S",1:U,"_0"),
                                            paste0("E",1:U,"_0"),
                                            paste0("J",1:U,"_0"),
                                            paste0("P",1:U,"_0"),
                                            paste0("I",1:U,"_0"),
                                            paste0("R",1:U,"_0"))),
         dmeasure=measles_dmeasure,
         unit_emeasure=measles_unit_emeasure,
         unit_mmeasure=measles_unit_mmeasure,
         unit_vmeasure=measles_unit_vmeasure,
         rmeasure=measles_rmeasure,
         unit_dmeasure=measles_unit_dmeasure) -> origin_model_ct


dat_ct <-  simulated_timeseries(model = origin_model_ct,parameter = parameter_all,
                                case_origin = case_all_1)

dat_ct_day <- data.frame(summarise(group_by(dat_ct,day),
                                   sim = floor(sum(median)),low = floor(sum(low)),
                                   high = floor(sum(high))))


dat_ct_day$type <- "Travel restrication";
dat_fcct_day$type <- "No social distancing and \n no travel restrication";
dat_fc_day$type <- "Social distancing"
dat_origin_day$type <- "Social distancing and \n travel restrication"

dat_day <- rbind(dat_ct_day,dat_fcct_day,dat_fc_day,dat_origin_day)

write.csv(dat_day,file = "Process_Data/City_Network_Intervention_Weak_UD.csv",row.names = F)

save(origin_model,origin_model_ct,origin_model_fc,origin_model_fcct,
     file = "E:/nCOV/Process_Data/City_Network_Intervention_Weak_UD.RData")


