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

load("Spat_Pomp_big_covar_Update.RData")

measles_unit_statenames <- c('S','E','I',"R","J","P",'C')
measles_statenames <- paste0(rep(measles_unit_statenames,each=U),1:U)
measles_IVPnames <- paste0(measles_statenames[1:((length(measles_unit_statenames)-1)*U)],"_0")

measles_RPnames <- c("alpha","gamma","sigma","psi","g","R0IN",
                     "qua","quawh","qualow","rho","kappa","omega")
measles_paramnames <- c(measles_RPnames,measles_IVPnames)
measles_covarnames <- paste0(rep(c("pop"),each=U),1:U)


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
                             powVec[u] = pow((kappa*I[u]+P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                             }
                             
                             for (u = 0 ; u < U ; u++) {
                             
                             if(kl < 15){
                                beta = R0IN*gamma/(1-omega+omega*kappa);
                             } else {
                                beta = R0IN*gamma/(1-omega+omega*kappa)*contreduce[u];
                             } 
                             
                             // expected force of infection
                             foi = pow((kappa*I[u]+P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                             
                             
                             for (v=0; v < U ; v++) {
                             if(v != u){
                             if(kl > 15){
                             foi += g * (v_by_g[v+U][u] * powVec[v] - (v_by_g[u+U][v]+flowno[2][u])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]);
                             } else {
                             foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v]+flowno[1][u])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]);
                             }
                             }
                             }
                             if(foi < 0) foi = 0;
                             // white noise (extrademographic stochasticity)
                             
                             double qm;
                             
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

spatPomp(ncov_cases,
         units = "city",
         times = "day",
         t0 = min(ncov_cases$day)-1,
         unit_statenames = measles_unit_statenames,
         covar = dat_pop_1,
         tcovar = "day",
         rprocess=euler(measles_rprocess, delta.t=dt),
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
         unit_dmeasure=measles_unit_dmeasure) -> origin_model


simlu <- read.csv("Hope_UHQ_2.csv",header=T) %>% arrange(desc(loglik));
#simlu <- simlu[simlu$loglik > (max(simlu$loglik) - 5),c(1:29)]
simlu <- simlu[1:20,c(1:(ncol(simlu) -4))]
guess.param <- apply(simlu,2, function(data)quantile(data,probs = c(0.05,0.95)))
measles_params <- data.frame(guess.param)
measles_params$R0IN <- c(1.2,3.5)
measles_params$gamma <- c(1/10.5,1/3.5)
measles_params$sigma <- c(1/10.2,1/5.2)
measles_params$g <- c(9.7,9.7)

library(doSNOW)
mcors <- 80
random.number <- 1000
Np <- 3000
Nmif <- 100
rw.sd_rp <- 0.02
rw.sd_ivp <- 0.02


fix_col <- c(1:12)[!measles_paramnames[c(1:12)] %in% c("alpha","g")]
rw_text <- paste(c(paste(measles_paramnames[fix_col],"=",rw.sd_rp,sep = ""),
                   paste(measles_paramnames[-c(1:12)],"=ivp(",rw.sd_ivp,")",sep = "")),collapse = ",")

rw_text <- paste("rw.sd(",rw_text,")",sep="")

parameter.rw.sd <- eval(parse(text = rw_text))
### (length(names(parameter.rw.sd)) - 3): length(names(parameter.rw.sd))

cl <- makeCluster(mcors,type = "SOCK"); registerDoSNOW(cl)
set.seed(998468235L,kind = "L'Ecuyer")
mcor.set <- list(preschedule = FALSE,set.seed = TRUE)

global.result <- foreach(i = 1:random.number, .combine=rbind,
                         .packages = c("magrittr","spatPomp","pomp",
                                       "stats","methods","dplyr","tidyr"),.errorhandling = "remove",
                         .inorder = FALSE, .options.multicore = mcor.set) %dopar%  {
                           
                           start_params <- apply(measles_params,2,function(data)runif(1,data[1],data[2]))
                           
                           #origin_model@params <- start_params
                           
                           mif2(origin_model,params = start_params,
                                Np = Np,
                                Nmif = Nmif,
                                #islands = 10,
                                tol =  1e-300,
                                cooling.type = "geometric",cooling.fraction.50 = 0.5,
                                rw.sd = parameter.rw.sd 
                           ) %>%  mif2() -> mf
                           pf <- replicate(10, pfilter(mf, Np = Np))
                           ll <- sapply(pf,logLik)
                           ll <- logmeanexp(ll, se = TRUE)
                           data.frame(as.list(coef(mf)),
                                      loglik = ll[1],
                                      loglik.se = ll[2])
                         }

global.result  <- arrange(global.result,desc(loglik))
write.csv(global.result,file=paste("Hope_UPQ_NF_IVP_K",".csv",sep = ""),row.names = F)

stopCluster(cl)
registerDoSEQ()

global.result.1 <- filter(global.result,!is.nan(alpha))

parameter <- arrange(global.result.1,desc(loglik))[1,]

sims <- simulate(pomp(origin_model),param = unlist(parameter),nsim = 10,format = "data.frame")
case_simu <- sims[,c(1,which(stri_detect_fixed(names(sims),"cases")))]
case_simu_low <- aggregate(case_simu[,-1],by = list(day = case_simu$day),quantile,probs = 0.05)
case_simu_median <- aggregate(case_simu[,-1],by = list(day = case_simu$day),quantile,probs = 0.5)
case_simu_high <- aggregate(case_simu[,-1],by = list(day = case_simu$day),quantile,probs = 0.95)
names(case_simu_low)[-1] <- names(case_simu_median)[-1] <- names(case_simu_high)[-1]  <- unique(ncov_cases$city)
case_simu_low_1 <- melt(case_simu_low,id= "day");names(case_simu_low_1) <- c("day","city","low")
case_simu_median_1 <- melt(case_simu_median,id= "day");names(case_simu_median_1) <- c("day","city","median")
case_simu_high_1 <- melt(case_simu_high,id= "day");names(case_simu_high_1) <- c("day","city","high")
case_plot <- merge(case_simu_high_1,case_simu_low_1,by = c("day","city"))
case_plot <- merge(case_plot,case_simu_median_1,by = c("day","city"))
case_plot <- merge(ncov_cases,case_plot,by = c("day","city"))
case_plot$day <- case_plot$day + as.Date("2020-01-09")
Sys.setlocale("LC_TIME", "German") # Windows

ggplot(case_plot,aes(x = day,y = median,ymin = low,ymax = high))+
  geom_ribbon(alpha = 0.2,fill = "#EE6363") +
  geom_line(aes(y = cases),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = median),color = "red",size = 0.8,alpha =0.9) +
  #theme_bw(base_family = "serif",base_size = 11) +
  facet_wrap(~city,scales = "free_y")+
  scale_x_date(date_breaks = "5 days", date_labels = "%b-%d")+
  theme_classic(base_size = 18,base_family = "serif") + xlab("") + ylab("Cases")

ggsave(paste("Hope_UPQ_NF_IVP_K",".pdf",sep = ""),width = 22,height = 16,units = "cm",dpi = 300)


