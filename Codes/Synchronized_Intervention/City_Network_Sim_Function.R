

simulated_timeseries <- function(model = origin_model_1,
                                 parameter = parameter_all,
                                 case_origin = case_all,target_name){
  
  sims <- simulate(pomp(model),param = unlist(parameter),nsim = 1,format = "data.frame")
  case_simu <- sims[,c(1,which(stri_detect_fixed(names(sims),"cases")))]
  case_simu_low <- aggregate(case_simu[,-1],by = list(day = case_simu$day),quantile,probs = 0.05)
  case_simu_median <- aggregate(case_simu[,-1],by = list(day = case_simu$day),quantile,probs = 0.5)
  case_simu_high <- aggregate(case_simu[,-1],by = list(day = case_simu$day),quantile,probs = 0.95)
  names(case_simu_low)[-1] <- names(case_simu_median)[-1] <- names(case_simu_high)[-1]  <- 
    target_name[as.numeric(gsub("cases","",names(case_simu_low)[-1]))]
  case_simu_low_1 <- melt(case_simu_low,id= "day");names(case_simu_low_1) <- c("day","city","low")
  case_simu_median_1 <- melt(case_simu_median,id= "day");names(case_simu_median_1) <- c("day","city","median")
  case_simu_high_1 <- melt(case_simu_high,id= "day");names(case_simu_high_1) <- c("day","city","high")
  case_plot <- merge(case_simu_high_1,case_simu_low_1,by = c("day","city"))
  case_plot <- merge(case_plot,case_simu_median_1,by = c("day","city"))
  case_plot <- merge(case_plot,case_origin,by = c("day","city"))
  case_plot$day <- case_plot$day + as.Date("2020-01-09")
  return(case_plot)
}


measles_dmeasure <- Csnippet("
                             const double *C = &C1;
                             const double *cases = &cases1;
                             double m,v;
                             double tol = 1e-300;
                             int u;
                             int kl = floor(t);
                             double rhom;
                             
                             lik = 0;
                             for (u = 0; u < U; u++) {
                             rhom =  rho;
                             m = rhom*C[u];
                             v = m*(1.0-rhom+psi*psi*m);
                             
                             if (cases[u] > 0.0) {
                             lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0)+tol);
                             } else {
                             lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
                             }
                             }
                             if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
                             ")

measles_rmeasure <- Csnippet("
                             const double *C = &C1;
                             double *cases = &cases1;
                             double m,v;
                             double tol = 1.0e-300;
                             int u;
                             int kl = floor(t);
                             double rhom;
                             
                             for (u = 0; u < U; u++) {
                             rhom =  rho;
                             m = rhom*C[u];
                             v = m*(1.0-rhom+psi*psi*m);
                             
                            // cases[u] = rnorm(m,sqrt(v)+tol);
                               cases[u] = C[u];
                             if (cases[u] > 0.0) {
                              cases[u] = cases[u];
                             } else {
                              cases[u] = 0.0;
                             }
                             }
                             ")

measles_unit_dmeasure <- Csnippet('
                                  // consider adding 1 to the variance for the case C = 0
                                  int kl = floor(t);
                                  double m,v;
                                  double rhom;
                                  
                                  rhom =  rho;
                                  m = rhom*C;
                                  v = m*(1.0-rhom+psi*psi*m);
                                  
                                  double tol = 1e-300;
                                  if (cases > 0.0) {
                                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                                  } else {
                                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                                  }
                                  if(give_log) lik = log(lik);
                                  ')

measles_unit_emeasure <- Csnippet("
                                  int kl = floor(t);
                                  double rhom;
                                  rhom =  rho;
                                  ey = rhom*C;
                                  ")

measles_unit_vmeasure <- Csnippet("
                                  //consider adding 1 to the variance for the case C = 0
                                  double m;
                                  int kl = floor(t);
                                  double rhom;
                                  rhom =  rho;
                                  m = rhom*C;
                                  vc = m*(1.0-rhom+psi*psi*m);
                                  ")

measles_unit_mmeasure <- Csnippet("
                                  double binomial_var;
                                  int kl = floor(t);
                                  double m;
                                  double rhom;
                                  rhom =  rho;
                                  m = rhom*C;
                                  binomial_var = rhom*(1-rhom)*C;
                                  
                                  if(vc > binomial_var) {
                                  M_psi = sqrt(vc - binomial_var)/m;
                                  }
                                  ")



measles_rinit <- Csnippet("
                          double *S = &S1;
                          double *E = &E1;
                          double *J = &J1;
                          double *I = &I1;
                          double *P = &P1;
                          double *R = &R1;
                          double *C = &C1;
                          const double *S_0 = &S1_0;
                          const double *E_0 = &E1_0;
                          const double *J_0 = &J1_0;
                          const double *I_0 = &I1_0;
                          const double *P_0 = &P1_0;
                          const double *R_0 = &R1_0;
                          
                          const double *pop = &pop1;
                          double m;
                          int u;
                          for (u = 0; u < U; u++) {
                          m = pop[u]/(S_0[u]+E_0[u]+J_0[u]+I_0[u]+R_0[u]+P_0[u]);
                          S[u] = nearbyint(m*S_0[u]);
                          E[u] = nearbyint(m*E_0[u]);
                          J[u] = nearbyint(m*J_0[u]);
                          I[u] = nearbyint(m*I_0[u]);
                          P[u] = nearbyint(m*P_0[u]);
                          R[u] = nearbyint(m*R_0[u]);
                          C[u] = 0;
                          }
                          ")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

