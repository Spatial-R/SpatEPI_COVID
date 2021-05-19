library(dplyr)
library(pomp)
library(reshape2)
library(spatPomp)
library(ggplot2)
library(stringi)
library(tidyr)
library(doSNOW)


mcors <- 80
cl <- makeCluster(mcors,type = "SOCK"); registerDoSNOW(cl)
set.seed(998468235L,kind = "L'Ecuyer")
mcor.set <- list(preschedule = FALSE,set.seed = TRUE)

load("Cities_Spatial_All.RData")
global.result.1 <- read.csv("Hope_UPQ_NF_IVP_K1.csv",header=T) %>% arrange(desc(loglik));

#search_target <- data.frame(expand.grid(delay = c(1,4,7,10),city = c(0,10,20,30),sims = c(1:50)))
search_target <- data.frame(expand.grid(delay = seq(0,20,2),city_num = c(0,10,20,30,40,50),sims = c(1:50),
                             city_name = c("haikou","lasa","wulumuqi")))

measles_unit_statenames <- c('S','E','I',"R","J","P",'C')
measles_statenames <- paste0(rep(measles_unit_statenames,each=U),1:U)
measles_IVPnames <- paste0(measles_statenames[1:((length(measles_unit_statenames)-1)*U)],"_0")

measles_RPnames <- c("alpha","gamma","sigma","psi","g","R0IN",
                     "qua","quawh","qualow","rho","kappa","omega")
measles_paramnames <- c(measles_RPnames,measles_IVPnames)
measles_covarnames <- paste0(rep(c("pop"),each=U),1:U)

test_dat <- data.frame(col = test_col,city = names(flow_all)[test_col])

foreach(i = 1:(nrow(search_target)), #.combine=rbind,
        .packages = c("magrittr","spatPomp","dplyr","tidyr",
                      "reshape2","stringi","base"),.errorhandling = "stop",
        .inorder = FALSE, .options.multicore = mcor.set) %dopar%  {
          
          source("City_Network_Sim_Function.R")
          
          k <- search_target[i,"delay"];j <- search_target[i,"city_num"];
          city_tem <- as.character(search_target[i,"city_name"])
          pos_wuhan <- match(city_tem,names(flow_all))
          flow_dat_wuhan <- flow_all[pos_wuhan,1:ncol(flow_all)]
          flow_wuhan <- data.frame(flow = as.numeric(flow_dat_wuhan),city = names(flow_all))
          flow_wuhan <- arrange(flow_wuhan,desc(flow))
          flow_wuhan_17 <- as.character(flow_wuhan[1:16,2])
          test_col_1 <- match(c(flow_wuhan_17,city_tem),names(flow_all))
          
          pos_wuhan_0 <- match("wuhan",names(flow_all))
          flow_dat_wuhan_0 <- flow_all[pos_wuhan_0,test_col]
          flow_dat_wuhan_0 <- data.frame(t(flow_dat_wuhan_0))
          flow_dat_wuhan_0$id <- 1:17;
          flow_dat_wuhan_0 <- arrange(flow_dat_wuhan_0,desc(X268))
          flow_dat_wuhan_0$id_1 <- test_col_1
          flow_dat_wuhan_0 <- arrange(flow_dat_wuhan_0,id)
          test_col_2 <- flow_dat_wuhan_0$id_1
          
          flow_C_rows <- apply(flow_all,1,to_C_array)
          flow_C_array <- to_C_array(flow_C_rows)
          
          city_num <- sample((1+k):(7+k),size = ncol(flow_all), replace = T) + 15
          
          if(!(j == 0)){
            other <- match(flow_wuhan[1:j,"city"],names(flow_all))
            city_num[c(other)] <- 15 + k
          }
          city_num[c(pos_wuhan)] <- 15+k
          
          time_point_array <- to_C_array(city_num)
          #test_col <- match(flow_wuhan[1:17,"city"],names(flow_all))
          #test_col <- pos_wuhan
          v_by_g_C <- Csnippet(paste0("const double v_by_g[",nrow(flow_all),"][",U,"] = ",
                                      flow_C_array,"; \n ",
                                      "const int targetp = ",pos_wuhan,"; \n ",
                                      "const double contreduce[",U,"] = ",flow_CT_array,";",
                                      "const double inte_time[",U,"] = ",time_point_array,";"
                                      #  "const double contact[",nrow(contact_all),"][",U,"] = ",
                                      #  contact_C_array,";"
          ))
          
          measles_globals <- Csnippet(
            paste0(v_by_g_C)
          )
          
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
                             
                             if(kl < inte_time[u]){
                             beta = R0IN*gamma/(1-omega+omega*kappa);
                             } else {
                             beta = R0IN*gamma/(1-omega+omega*kappa)*contreduce[u];
                             } 
                             
                             // expected force of infection
                             foi = pow((I[u]+kappa*P[u])/(S[u]+I[u]+R[u]+E[u]+P[u]),alpha);
                             
                             
                             for (v=0; v < U ; v++) {
                             if(v != u){
                             if(kl > inte_time[u]){
                             foi += g * (v_by_g[v+U][u] * powVec[v] - (v_by_g[u+U][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]);
                             } else {
                             foi += g * (v_by_g[v][u] * powVec[v] - (v_by_g[u][v])* powVec[u]) /(S[u]+I[u]+R[u]+E[u]);
                             }
                             }
                             }
                             if(foi < 0) foi = 0;
                             // white noise (extrademographic stochasticity)
                             
                             double qm,jk;
                             
                             if( kl < 15) {
                              jk  = kl;                                  
                             } else if ( kl <= inte_time[u]) {
                              jk = 15;
                             } else {
                              jk = 15 + (kl - inte_time[u]);
                             }
                              
                               
                             if (u == (targetp-1)){
                             qm = 1- (exp(quawh*kl+log(qualow))/(exp(quawh*kl+log(qualow))+1));
                             } else {
                             qm = 1- (exp(qua*jk+log(qualow))/(exp(qua*jk+log(qualow))+1));
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
                             C[u] +=  trans[4] + trans[5] + trans[6];           // true incidence
                             
                             if(S[u] < 0) S[u] = 0;
                             if(E[u] < 0) E[u] = 0;
                             if(J[u] < 0) J[u] = 0 ;
                             if(I[u] < 0) I[u] = 0;
                             if(P[u] < 0) P[u] = 0;
                             if(R[u] < 0) R[u] = 0;
                             }
                             ')
          
          
          case_all_1 <- data.frame(expand.grid(day = 7:90,city = sort(unique(case_all$city))))
          case_all_1$cases <- NA; case_all_1$city <- as.character(case_all_1$city)
          
          pop_list <- lapply(sort(unique(case_all$city)),function(data){
            dat_tem <- filter(pop_all,city == data) 
            dat_tem <- data.frame(day = 7:90, city = data, pop = dat_tem[1,"pop"])
            return(dat_tem)
          })
          pop_all <- bind_rows(pop_list)
          
          
          spatPomp(case_all_1,
                   units = "city",
                   times = "day",
                   t0 = min(case_all$day)-1,
                   unit_statenames = measles_unit_statenames,
                   covar = pop_all, 
                   rprocess=euler(measles_rprocess, delta.t=dt),
                   unit_accumvars = c("C"),
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
                   eunit_measure=measles_unit_emeasure,
                   munit_measure=measles_unit_mmeasure,
                   vunit_measure=measles_unit_vmeasure,
                   rmeasure=measles_rmeasure,
                   dunit_measure=measles_unit_dmeasure) -> origin_model
          
          parameter <- arrange(global.result.1,desc(loglik))[1,];
          s1_pos <- which(names(parameter) == "S1_0");
          lgh_pos <- which(names(parameter) == "loglik")
          
          estimate_p <- unlist(lapply(measles_unit_statenames[1:(length(measles_unit_statenames)-1)], function(data){
            test_p <- parameter[(s1_pos+which(substr(names(parameter)[s1_pos:lgh_pos],1,1) == data)-1)]
            names(test_p) <- paste(data,test_col_2,"_0",sep = "")
            return(test_p)
          }))
          
          fixed_p <- unlist(lapply(measles_unit_statenames[1:(length(measles_unit_statenames)-1)], function(data){
            fix_col <- c(1:ncol(flow_all))[-c(test_col_2)]
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
                                              case_origin = case_all_1,target_name = names(flow_all))
          #dat_origin <- filter(dat_origin,!city == city_tem)
          #dat_origin_day <- data.frame(summarise(group_by(dat_origin,day),
          #                                       sim = floor(sum(median)),low = floor(sum(low)),high = floor(sum(high))))
          #result_fin <- rbind(result_fin,dat_origin_day)
          write.csv(dat_origin,file = paste0("Result_5/Synchronized_NWH_",city_tem,"&",j,"&",k,"&",i,".csv"),row.names = F)
        }
stopCluster(cl)
registerDoSEQ()
