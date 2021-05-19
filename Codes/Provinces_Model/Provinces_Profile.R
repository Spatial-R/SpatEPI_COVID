library(dplyr)
#library(pomp,lib = "/home/spatialr/R/x86_64-pc-linux-gnu-library/3.5")
library(pomp)
library(reshape2)
#library(spatPomp,lib = "/home/spatialr/R/x86_64-pc-linux-gnu-library/3.5")
library(spatPomp)
library(ggplot2)
library(stringi)
library(tidyr)
library(readr)

U <- 18

measles_unit_statenames <- c('S','E','I',"R","J","P",'C')
measles_statenames <- paste0(rep(measles_unit_statenames,each=U),1:U)
measles_IVPnames <- paste0(measles_statenames[1:((length(measles_unit_statenames)-1)*U)],"_0")

measles_RPnames <- c("alpha","gamma","sigma","psi","g","R0IN",
                     "qua","quawh","qualow","rho","kappa","omega","quatar")
measles_paramnames <- c(measles_RPnames,measles_IVPnames)
measles_covarnames <- paste0(rep(c("pop"),each=U),1:U)


measles_rprocess <- Csnippet('
                             double beta,foi;
                             // double dw,sus_del;
                             double rate[7], trans[7];
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
                             foi += g * (v_by_g[v+U][u] * powVec[v] - flowall[2][u]* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                             } else {
                             foi += g * (v_by_g[v][u] * powVec[v] - flowall[1][u]* powVec[u]) /(S[u]+I[u]+R[u]+E[u]+P[u]);
                             }
                             }
                             }
                             if(foi < 0) foi = 0;
                             // white noise (extrademographic stochasticity)
                             
                             double qm;
                             
                             if (u == (12-1)){
                             qm = 1- (exp(quawh*kl+log(qualow))/(exp(quawh*kl+log(qualow))+1));
                             } else if (u == (U-1)) {
                             qm = 1- (exp(quatar*kl+log(qualow))/(exp(quatar*kl+log(qualow))+1));
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
                             
                             trans[0]  = nearbyint(S[u]*rate[0]*dt);
                             trans[1]  = nearbyint(E[u]*rate[1]*dt);
                             trans[2]  = nearbyint(E[u]*rate[2]*dt);
                             trans[3]  = nearbyint(E[u]*rate[3]*dt);
                             trans[4]  = nearbyint(J[u]*rate[4]*dt);
                             trans[5]  = nearbyint(P[u]*rate[5]*dt);
                             trans[6]  = nearbyint(I[u]*rate[6]*dt);
                             
                             S[u] += -trans[0];
                             E[u] +=  trans[0] - trans[1] - trans[2] - trans[3];
                             J[u] +=  trans[1] - trans[4];
                             P[u] +=  trans[2] - trans[5];
                             I[u] +=  trans[3] - trans[6];
                             R[u] +=  trans[4] + trans[5] + trans[6];
                             C[u] += trans[4];           // true incidence
                             
                             if(S[u] < 0) S[u] = 0;
                             if(E[u] < 0) E[u] = 0;
                             if(J[u] < 0) J[u] = 0;
                             if(I[u] < 0) I[u] = 0;
                             if(P[u] < 0) P[u] = 0;
                             if(R[u] < 0) R[u] = 0;
                             }
                             ')


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
                             
                             cases[u] = rnorm(m,sqrt(v)+tol);
                             if (cases[u] > 0.0) {
                             cases[u] = nearbyint(cases[u]);
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

for (pro_n in c("tianjin","shanxi","shan-xi","hainan","guizhou","yunnan")){


### "shanghai","chongqing","jiangxi","jiangsu","anhui","fujian"， "zhejiang","shanghai","shandong","jiangsu","anhui"


load(paste("Data/Cities_Spatial_",pro_n,".RData",sep = ""))

spatPomp(case_all_final,
         units = "city",
         times = "day",
         t0 = min(case_all_final$day)-1,
         unit_statenames = measles_unit_statenames,
         covar = pop_all_final,
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
                                            "qualow","rho","omega","quatar",
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

simlu <- read.csv(paste("Hope_UHO_PRO_IVP_L_",pro_n,".csv",sep=""),header=T) %>% arrange(desc(loglik));
#simlu <- simlu[simlu$loglik > (max(simlu$loglik) - 5),c(1:29)]
simlu <- simlu[1:10,c(1:(ncol(simlu) - 2))]
guess.param <- apply(simlu,2, function(data)quantile(data,probs = c(0.05,0.95)))
measles_params <- data.frame(guess.param)


library(doSNOW)
mcors <- 120
random.number <- 1000
Np <- 3000
Nmif <- 100
rw.sd_rp <- 0.02
rw.sd_ivp <- 0.02

cl <- makeCluster(mcors,type = "SOCK"); registerDoSNOW(cl)
set.seed(998468235L,kind = "L'Ecuyer")
mcor.set <- list(preschedule = FALSE,set.seed = TRUE)
parameter.low  <- unlist(measles_params[1,-13]);
parameter.high <- unlist(measles_params[2,-13]);
#min_tar <- min(simlu$quatar) - 0.3;
#max_tar <- max(simlu$quatar) + 0.3
min_tar <- 0.3; max_tar <- 0.6
if(min_tar < 0) min_tar <- 0
if(max_tar > 1) max_tar <- 1
startparameter <- profileDesign(quatar = seq(from=min_tar,to = max_tar,length = 20),
                                  lower = parameter.low, upper = parameter.high, nprof = 15)
  
  global.result <- foreach(p = iter(startparameter,"row"), .combine=rbind,
                           .packages = c("magrittr","spatPomp","pomp",
                                         "stats","methods","dplyr","tidyr"),.errorhandling = "remove",
                           .inorder = FALSE, .options.multicore = mcor.set) %dopar%  {
                             
                             
                             mif2(origin_model,params = unlist(p),
                                  Np = Np,
                                  Nmif = Nmif,
                                  #islands = 10,
                                  tol =  1e-300,
                                  cooling.type = "geometric",cooling.fraction.50 = 0.5,
                                  rw.sd = rw.sd(S18_0 = ivp(rw.sd_rp),
                                                E18_0 = ivp(rw.sd_rp),
                                                I18_0 = ivp(rw.sd_rp),
                                                J18_0 = ivp(rw.sd_rp),
                                                P18_0 = ivp(rw.sd_rp),
                                                R18_0 = ivp(rw.sd_rp),
                                                psi = rw.sd_rp
                                  ) 
                             ) %>%  mif2() -> mf
                             pf <- replicate(10, pfilter(mf, Np = Np))
                             ll <- sapply(pf,logLik)
                             ll <- logmeanexp(ll, se = TRUE)
                             data.frame(as.list(coef(mf)),
                                        loglik = ll[1],
                                        loglik.se = ll[2])
                           }
  write.csv(global.result,file=paste("Hope_Qua_NF_",pro_n,".csv",sep=""),row.names = F)
  
  stopCluster(cl)
  registerDoSEQ()
  
}