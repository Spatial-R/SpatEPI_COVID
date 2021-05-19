library(pomp)
library(ggsci)
library(dplyr)
library(stringi)
library(viridis)
library(cowplot)
library(ggplot2)

cols <- pal_npg("nrc",alpha = 0.6)(9); scales::show_col(cols)
time_range <- seq(7, 150, by = 1)
source("Codes/Single_Region_Dynamic_Model.R")
Sys.setlocale("LC_TIME", "English") # Windows

load("Process_Data/Spat_Pomp_big_covar_Update.RData")

res_tem <- read.csv("E:/nCOV/Server/Last/Hope_UPQ_NF_IVP_K1.csv") %>% arrange(desc(loglik))

param_col <- match(c("quawh","R0IN","sigma","gamma","psi","g",
                     "omega","kappa","rho","qualow"),names(res_tem))
initial_state_col <- which(stri_detect_fixed(names(res_tem),"12_0"))

pre_de <- res_tem[2,c(param_col,initial_state_col)]
names(pre_de)[11:16] <- c("S_0","E_0","IU_0","R_0","IQ_0","IH_0")
pre_de$socidst <- 1;


#plot(1/(exp(-(as.numeric(pre_de["quawh"])*c(7:34)+log(as.numeric(pre_de["qualow"]))))+1))

res_fin <- data.frame()

for (j in c(seq(0.01,0.99,0.01))){
  
  pre_de[match("socidst",names(pre_de))] <- 1-j
  print(paste0("& socidst=",j))
  
  sims_total_strong <- simulate(model_total_strong,params = pre_de,nsim = 100,format = "data.frame")
  case_total_strong <- sims_total_strong[,c(1,match(c("cases"),names(sims_total_strong)))]
  sims_confirm_strong <- simulate(model_confirm_strong,params = pre_de,nsim = 100,format = "data.frame")
  case_confirm_strong <- sims_confirm_strong[,c(1,match(c("cases"),names(sims_confirm_strong)))]
  sims_ud_strong <- simulate(model_ud_strong,params = pre_de,nsim = 100,format = "data.frame")
  case_ud_strong <- sims_ud_strong[,c(1,match(c("cases"),names(sims_ud_strong)))]
  case_total_strong_1 <- data.frame(summarise(group_by(case_total_strong,time),total_low = quantile(cases,0.05),
                                              total_mid = quantile(cases,0.5),total_high = quantile(cases,0.95)))
  case_confirm_strong_1 <- data.frame(summarise(group_by(case_confirm_strong,time),confirm_low = quantile(cases,0.05),
                                                confirm_mid = quantile(cases,0.5),confirm_high = quantile(cases,0.95)))
  case_ud_strong_1 <- data.frame(summarise(group_by(case_ud_strong,time),ud_low = quantile(cases,0.05),
                                           ud_mid = quantile(cases,0.5),ud_high = quantile(cases,0.95)))
  
  sims_total_weak <- simulate(model_total_weak,params = pre_de,nsim = 100,format = "data.frame")
  case_total_weak <- sims_total_weak[,c(1,match(c("cases"),names(sims_total_weak)))]
  sims_confirm_weak <- simulate(model_confirm_weak,params = pre_de,nsim = 100,format = "data.frame")
  case_confirm_weak <- sims_confirm_weak[,c(1,match(c("cases"),names(sims_confirm_weak)))]
  sims_ud_weak <- simulate(model_ud_weak,params = pre_de,nsim = 100,format = "data.frame")
  case_ud_weak <- sims_ud_weak[,c(1,match(c("cases"),names(sims_ud_weak)))]
  case_total_weak_1 <- data.frame(summarise(group_by(case_total_weak,time),total_low = quantile(cases,0.05),
                                            total_mid = quantile(cases,0.5),total_high = quantile(cases,0.95)))
  case_confirm_weak_1 <- data.frame(summarise(group_by(case_confirm_weak,time),confirm_low = quantile(cases,0.05),
                                              confirm_mid = quantile(cases,0.5),confirm_high = quantile(cases,0.95)))
  case_ud_weak_1 <- data.frame(summarise(group_by(case_ud_weak,time),ud_low = quantile(cases,0.05),
                                         ud_mid = quantile(cases,0.5),ud_high = quantile(cases,0.95)))
  
  simu_strong <- data.frame(cbind(case_total_strong_1,case_confirm_strong_1[,-1],case_ud_strong_1[,-1]))
  simu_weak <- data.frame(cbind(case_total_weak_1,case_confirm_weak_1[,-1],case_ud_weak_1[,-1]))
  
  simu_strong$type <- "strong"; simu_weak$type <- "weak"; 
  simu_exm <- rbind(simu_strong,simu_weak)
  simu_exm$socidst <- j
  res_fin <- rbind(res_fin,simu_exm)
}

res_fin <- filter(res_fin,time > 7)

wuhan_case <- filter(ncov_cases,city == "wuhan")
wuhan_case <- mutate(wuhan_case,time = as.Date("2020-01-09")+ day)
full_date <- data.frame(time = unique(res_fin$time))
wuhan_case <- merge(wuhan_case[,-c(1,2)],full_date,by.x = "time",all.y = T)

res_fin_1 <- mutate(res_fin,time = as.Date("2020-01-9") + time)
#res_fin_2 <- filter(res_fin_1, time > as.Date("2020-1-16") & time < as.Date("2020-03-30"))
res_fin_1 <- mutate(res_fin_1,type = factor(type,levels = c("weak","strong"),labels = c("Weak","Strong")))
res_fin_2 <- filter(res_fin_1, time > as.Date("2020-1-16") & time < as.Date("2020-03-30"))


peak_list <- lapply(c("Strong","Weak"),function(id){
  dat_tem <- filter(res_fin_2,type == id)
  soc_list <- lapply(unique(dat_tem$socidst),function(sodid){
    dat_tem1 <- filter(dat_tem,socidst == sodid)
    peak_pos <-  which.max(dat_tem1$total_mid)
    data.frame(peak_time = dat_tem1[peak_pos,"time"],
               peak_int = dat_tem1[peak_pos,"total_mid"],
               socidst = sodid)
  })
  soc_dat <- bind_rows(soc_list)
  soc_dat$type <- id
  return(soc_dat)
})
peak_dat <- bind_rows(peak_list)
peak_dat <- mutate(peak_dat,type = factor(type,levels = c("Weak","Strong")))

fig_peak <- ggplot(data=peak_dat,aes(x = socidst, y = peak_time,group = factor(type),
                                     color = factor(type))) + geom_point() +
  scale_y_date() +
  theme_classic(base_family = "serif",base_size = 12) +
  xlab("Strength of social distancing") +
  ylab("Peak time") +
  scale_color_manual(values = rev(cols[c(8,4)]),guide = "none")


fig_int <- ggplot(data=peak_dat,aes(x = socidst, y = log(peak_int),group = factor(type),
                                    color = factor(type))) + geom_point() + 
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Strength of social distancing") +
  ylab("Peak intensity \n (Nature logarithms scale)") +
  scale_color_manual(values = rev(cols[c(8,4)]),guide = "none") +
  theme(legend.position = c(0.85,0.8)) +
  scale_y_continuous(breaks = seq(11,12.5,0.3))


dat_sum <- data.frame(summarise(group_by(res_fin_2,socidst,type),total = sum(total_mid)))

fig_sum <- ggplot(data=dat_sum,aes(x = socidst, y = log(total),group = factor(type),
                                   color = factor(type))) + geom_point() + 
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Strength of social distancing") +
  ylab("Epidemic size \n (Nature logarithms scale)") +
  scale_y_continuous(breaks = seq(13,16,0.5))+
  scale_color_manual(values = rev(cols[c(8,4)]),guide = guide_legend(title = "Infection isolation")) +
  theme(legend.position = c(0.70,0.9),legend.background = element_blank())

fig_merge <- plot_grid(fig_peak,fig_int,fig_sum,labels = c("C","D","E"),ncol = 3)



fig_total_list <- lapply(c("Weak","Strong"),function(id){
  
  (fig_sac_1 <- ggplot() +
     geom_line(data = res_fin_2[res_fin_2$type == id,],
               aes(x = time, y = log(total_mid),color = socidst,group = factor(socidst)),size = 1.2) + 
     #geom_line(data = wuhan_case,aes(x = time, y = log(cases)),color = "red")+
     theme_classic(base_family = "serif",base_size = 12) +
     xlab("Time (days)") + ylab("Daily new COVID-19 cases \n (Nature logarithms scale)") +
     #scale_x_date(date_breaks = "20 days",date_labels = "%b-%d")+
     scale_x_date(breaks = c(as.Date("2020-02-01"),as.Date("2020-03-01"),as.Date("2020-04-01")),
                  labels = c("Feb","Mar","April"))+
     scale_color_viridis(option="B") + 
     scale_y_continuous(limits = c(0,16),breaks = seq(0,16,3))+
     geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
     guides(color = guide_colourbar(title = "Strength of social distancing",
                                    title.position = "right", 
                                    title.theme = element_text(angle = -90), 
                                    title.hjust = 0.5, 
                                    barwidth = 1, 
                                    barheight = 10))) #+
     #theme(legend.position = c(0.85,0.75))+
     #geom_segment(aes(x = as.Date("2020-1-31"),xend = as.Date("2020-01-25"),
      #                y = 14,yend = 14),arrow= arrow(length = unit(0.3,"cm")),linetype = 1)+
     #annotate("text",x=as.Date("2020-02-27"),y = 14.2,label=""))
  #### The time point where social distancing \n measures were implemented
  
  return(fig_sac_1)
})



fig_confirm_list <- lapply(c("Weak","Strong"),function(id){
  res_fin_3 <- mutate(res_fin_2,confirm_mid = ifelse(confirm_mid < 1, 0, log(confirm_mid)))
  (fig_sac_1 <- ggplot() +
      geom_line(data = res_fin_3[res_fin_3$type == id,],
                aes(x = time, y = confirm_mid,color = socidst,group = factor(socidst)),size = 1.2) + 
      #geom_line(data = wuhan_case,aes(x = time, y = log(cases)),color = "red")+
      theme_classic(base_family = "serif",base_size = 12) +
      xlab("Time (days)") + ylab("Daily new confirmed cases \n (Nature logarithms scale)") +
      scale_x_date(breaks = c(as.Date("2020-02-01"),as.Date("2020-03-01"),as.Date("2020-04-01")),
                   labels = c("Feb","Mar","April"))+
      scale_color_viridis(option="B") + 
      scale_y_continuous(limits = c(0,16),breaks = seq(0,16,3))+
      geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
      guides(color = guide_colourbar(title = "Strength of social distancing",
                                     title.position = "right", 
                                     title.theme = element_text(angle = -90), 
                                     title.hjust = 0.5, 
                                     barwidth = 1, 
                                     barheight = 10)) +
      theme(legend.position = c(0.85,0.75)))#+
      #geom_segment(aes(x = as.Date("2020-1-31"),xend = as.Date("2020-01-25"),
       #                y = 14,yend = 14),arrow= arrow(length = unit(0.3,"cm")),linetype = 1)+
      #annotate("text",x=as.Date("2020-02-18"),y = 14.2,label=""))
  
  return(fig_sac_1)
})




fig_ud_list <- lapply(c("Strong","Weak"),function(id){
  res_fin_3 <- mutate(res_fin_2,ud_mid = ifelse(ud_mid < 1, 0, log(ud_mid)))
  (fig_sac_1 <- ggplot() +
      geom_line(data = res_fin_3[res_fin_3$type == id,],
                aes(x = time, y = ud_mid,color = socidst,group = factor(socidst)),size = 1.2) + 
      #geom_line(data = wuhan_case,aes(x = time, y = log(cases)),color = "red")+
      theme_classic(base_family = "serif",base_size = 12) +
      xlab("Time (days)") + ylab("Daily new undocumented cases \n (Nature logarithms scale)") +
      scale_x_date(breaks = c(as.Date("2020-02-01"),as.Date("2020-03-01"),as.Date("2020-04-01")),
                   labels = c("Feb","Mar","April"))+
      scale_color_viridis(option="B") + 
      scale_y_continuous(limits = c(0,16),breaks = seq(0,16,3))+
      geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
      guides(color = guide_colourbar(title = "Strength of social distancing",
                                     title.position = "right", 
                                     title.theme = element_text(angle = -90), 
                                     title.hjust = 0.5, 
                                     barwidth = 1, 
                                     barheight = 10)) +
      theme(legend.position = c(0.85,0.75))) #+
      #geom_segment(aes(x = as.Date("2020-1-31"),xend = as.Date("2020-01-25"),
      #                 y = 14,yend = 14),arrow= arrow(length = unit(0.3,"cm")),linetype = 1)+
      #annotate("text",x=as.Date("2020-02-17"),y = 14.2,label=""))
  
  return(fig_sac_1)
})


fig_total <- plot_grid(fig_total_list[[1]]+guides(color = "none"),fig_total_list[[2]],
                       labels = c("A","B"),rel_widths = c(0.45,0.55))
fig_total1 <- plot_grid(fig_total,fig_merge,ncol = 1)

ggsave(filename = "Figures/Epi_total.pdf",fig_total1,width = 23,height = 18,units = "cm",dpi = 300)

fig_confirm <- plot_grid(fig_confirm_list[[2]]+guides(color = "none"),fig_confirm_list[[1]],labels = c("A","B"))
ggsave(filename = "Figures/Epi_confirm.pdf",fig_confirm,width = 26,height = 14,units = "cm",dpi = 300)

fig_ud <- plot_grid(fig_ud_list[[2]]+guides(color = "none"),fig_ud_list[[1]],labels = c("A","B"))
ggsave(filename = "Figures/Epi_undocumented.pdf",fig_ud,width = 26,height = 14,units = "cm",dpi = 300)




time_range_1 <- seq(7,250,1)

pomp(data = data.frame(cases = NA, time = time_range_1),
     times = "time",
     t0 = 1, 
     rprocess = euler(step.fun=seir_rprocess_total_end,delta.t=1/2),
     rinit = initlz,
     dmeasure = seir_dmeasure,
     rmeasure = seir_rmeasure_total,
     partrans=parameter_trans(log=c("R0IN"),
                              logit = c("psi","sigma","gamma","quawh","omega","kappa","qualow",
                                        "socidst","dettem",
                                        paste0("S","_0"),
                                        paste0("E","_0"),
                                        paste0("IU","_0"),
                                        paste0("IQ","_0"),
                                        paste0("IH","_0"),
                                        paste0("R","_0"))),
     accumvars = c("C","err"),
     statenames = names_seir,
     paramnames = c(para_seir,"dettem")) ->  model_total_end




res_end <- data.frame(); pre_de$dettem <- 0.2

for (i in c(seq(0.01,0.99,0.05))){
for (j in c(seq(0.01,0.99,0.05))){
  
  pre_de[match("socidst",names(pre_de))] <- 1-j
  pre_de$dettem <- i
  print(paste0("& socidst=",j,"; detect =",i))
  
  sims_total_strong <- simulate(model_total_end,params = pre_de,nsim = 10,format = "data.frame")
  case_total_strong <- sims_total_strong[,c(1,match(c("cases"),names(sims_total_strong)))]
  case_total_strong_1 <- data.frame(summarise(group_by(case_total_strong,time),total_low = quantile(cases,0.05),
                                              total_mid = quantile(cases,0.5),total_high = quantile(cases,0.95)))
  
  case_total_strong_2 <- case_total_strong_1[-(1:which.max(case_total_strong_1$total_mid)),]
  end_date <- case_total_strong_2[which(case_total_strong_2$total_mid < 50)[1],"time"]
  case_end <- data.frame(date = end_date,detect = i, socidst = j)
  res_end <- rbind(res_end,case_end)
}
}


res_end <- mutate(res_end,date1  = ifelse(date > 150,NA,date))

(fig_end <- ggplot(data = res_end,aes(y = detect,x = socidst)) +
    geom_tile(aes(fill = date1)) +
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = 5,
                         name = "Distance between end date of the epidemic", 
                         guide = guide_colorbar(title.position = "right", 
                                                title.theme = element_text(angle = -90), 
                                                title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 18)) +
    scale_y_continuous(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) + 
    scale_x_continuous(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) +
    ylab("Strength of contact tracing & testing") + xlab("Strength of social distancing"))





dat_total_strong <- data.frame(summarise(group_by(res_fin_2[res_fin_2$type == "Strong",],socidst),cases = sum(total_mid)))
dat_total_weak <- data.frame(summarise(group_by(res_fin_2[res_fin_2$type == "Weak",],socidst),cases = sum(total_mid)))
dat_confirm_strong <- data.frame(summarise(group_by(res_fin_2[res_fin_2$type == "Strong",],socidst),cases = sum(confirm_mid)))
dat_confirm_weak <- data.frame(summarise(group_by(res_fin_2[res_fin_2$type == "Weak",],socidst),cases = sum(confirm_mid)))
dat_ud_strong <- data.frame(summarise(group_by(res_fin_2[res_fin_2$type == "Strong",],socidst),cases = sum(ud_mid)))
dat_ud_weak <- data.frame(summarise(group_by(res_fin_2[res_fin_2$type == "Weak",],socidst),cases = sum(ud_mid)))

(dat_ud_weak[99,2]-dat_ud_strong[99,2])/dat_ud_strong[99,2]


total_strong_target <- filter(dat_total_strong,socidst > 0.825 & socidst < 0.84)$cases 
total_weak_target <- filter(dat_total_weak,socidst > 0.825 & socidst < 0.84)$cases 
(total_weak_target - total_strong_target)/total_weak_target

ud_strong_target <- filter(dat_ud_strong,socidst > 0.825 & socidst < 0.84)$cases 
ud_weak_target <- filter(dat_ud_weak,socidst > 0.825 & socidst < 0.84)$cases 
(ud_weak_target - ud_strong_target)/ud_weak_target

ud_strong_target <- filter(dat_ud_strong,socidst > 0.825 & socidst < 0.84)$cases 
ud_weak_target <- filter(dat_ud_weak,socidst > 0.825 & socidst < 0.84)$cases 
(ud_weak_target - ud_strong_target)/ud_weak_target


total_perc_list <- lapply(unique(dat_total_strong$socidst),function(id){
  strong_target <- filter(dat_total_strong,socidst == id)$cases 
  weak_target <- filter(dat_total_weak,socidst == id)$cases 
  perc_tem <- (strong_target)/weak_target
  data.frame(perc = perc_tem,socidt = id)
})
total_perc <- bind_rows(total_perc_list)


confirm_perc_list <- lapply(unique(dat_total_strong$socidst),function(id){
  strong_target <- filter(dat_confirm_strong,socidst == id)$cases 
  weak_target <- filter(dat_confirm_weak,socidst == id)$cases 
  perc_tem <- (strong_target)/weak_target
  data.frame(perc = perc_tem,socidt = id)
})
confirm_perc <- bind_rows(confirm_perc_list)


ud_perc_list <- lapply(unique(dat_total_strong$socidst),function(id){
  strong_target <- filter(dat_ud_strong,socidst == id)$cases 
  weak_target <- filter(dat_confirm_weak,socidst == id)$cases 
  perc_tem <- (strong_target)/weak_target
  data.frame(perc = perc_tem,socidt = id)
})
ud_perc <- bind_rows(ud_perc_list)



fig_total_perc <- ggplot(data = total_perc,aes(x = socidt, y = 1- perc)) + geom_point()+
  theme_classic(base_family = "serif",base_size = 12) +
  xlab("Strength of social distancing") +
  # ylab(expression(frac("Total caese in Strong CCT","Total caese in the weak CCT")))
  scale_y_continuous(limits = c(0,0.6))+
  ylab("Percentage of infection averted \n by the infection isolation")
ggsave(filename = "Figures/Perc_Total.pdf",fig_total_perc,width = 14,height = 12,units = "cm",dpi = 300)


fig_confirm_perc <- ggplot(data = confirm_perc,aes(x = socidt, y = perc)) + geom_point()+
  theme_classic(base_family = "serif",base_size = 12) +
  xlab("Strength of social distancing") +
  ylab(expression(frac("Confirmed caese in the strong infection isolation",
                       "Confirmed caese in the weak infection isolation")))
ggsave(filename = "Figures/Perc_Confirm.pdf",fig_confirm_perc,width = 14,height = 12,units = "cm",dpi = 300)


fig_ud_perc <- ggplot(data = ud_perc,aes(x = socidt, y = perc)) + geom_point()+
  theme_classic(base_family = "serif",base_size = 12) +
  xlab("Strength of social distancing") +
  ylab(expression(frac("Confirmed caese in the strong infection isolation",
                       "Confirmed caese in the weak infection isolation")))
ggsave(filename = "Figures/Perc_ud.pdf",fig_ud_perc,width = 14,height = 12,units = "cm",dpi = 300)




pre_de_imp <- pre_de; 
pre_de_imp$socidst <- 1; pre_de_imp$quawh <- 0

pre_de_imp[11:16] <- c(1,0,10/1e07,0,0,0)
names(pre_de_imp)[11:16] <- c("S_0","E_0","IU_0","R_0","IQ_0","IH_0")

res_fin_1 <- data.frame()

for (i in c(seq(0.1,0.9,0.1))){
  for (j in seq(0.1,1.0,0.1)){
    pre_de_imp[10] <- i/(1-i);
    pre_de_imp[17] <- j
    print(paste0("rho=",i,"& socidst=",j))
    sims <- simulate(model_total_weak,params = pre_de_imp,nsim = 10,format = "data.frame")
    case_simu <- sims[,c(1,which(stri_detect_fixed(names(sims),"cases")))]
    case_simu_exm <- aggregate(case_simu[,-1],by = list(time = case_simu$time),quantile,probs = 0.5)
    names(case_simu_exm) <- c("time","cases")
    case_simu_exm$rho <- i; case_simu_exm$socidst <- 1-j
    res_fin_1 <- rbind(res_fin_1,case_simu_exm)
  }
}
res_fin_1 <- filter(res_fin_1,time > 1 & time < 100)
res_fin_1$rho <- paste("CTT = ",res_fin_1$rho,sep = "")
res_fin_1$time <- as.Date("2020-01-16") + res_fin_1$time


(fig_sac <- ggplot(data = res_fin_1,aes(x = time, y = cases,color = socidst,group = factor(socidst))) +
    geom_line(size = 1.2) + facet_wrap(~rho)+#,scales = "free_y") +
    theme_bw(base_family = "serif",base_size = 12) +
    #scale_y_log10()+
    xlab("Time(days)") + ylab("Daily new COVID-19 cases") +
    scale_color_viridis(option="C") +
    guides(color = guide_colourbar(title = "Strength of social distancing",
                                   title.position = "right", 
                                   title.theme = element_text(angle = -90), 
                                   title.hjust = 0.5, 
                                   barwidth = 1, 
                                   barheight = 12)))


ggsave(filename = "Figures/Imported_con_trac.pdf",fig_sac,width = 16,height = 12,units = "cm",dpi = 300)




res_fin_2 <- data.frame()

for (i in c(seq(0.1,0.9,0.05))){
  for (j in seq(0.1,1.0,0.05)){
    pre_de_imp[10] <- i/(1-i);
    pre_de_imp[17] <- j
    print(paste0("rho=",i,"& socidst=",j))
    sims <- simulate(model_total_weak,params = pre_de_imp,nsim = 10,format = "data.frame")
    case_simu <- sims[,c(1,which(stri_detect_fixed(names(sims),"cases")))]
    case_simu_exm <- aggregate(case_simu[,-1],by = list(time = case_simu$time),quantile,probs = 0.5)
    names(case_simu_exm) <- c("time","cases")
    case_simu_exm$rho <- i; case_simu_exm$socidst <- 1-j
    res_fin_2 <- rbind(res_fin_2,case_simu_exm)
  }
}
res_fin_2 <- filter(res_fin_2,time > 1)

res_plot <- data.frame(summarise(group_by(res_fin_2,rho,socidst),cases = log(sum(cases))))

(fig_sc_sd <- ggplot(data = res_plot,aes(x = rho,y = socidst)) +
    geom_tile(aes(fill = cases)) +
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         #midpoint = 6,
                         name = "Total confirmed cases (Log-transformation)", 
                         guide = guide_colorbar(title.position = "right", 
                                                title.theme = element_text(angle = -90), 
                                                title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 18)) +
    scale_y_continuous(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) + 
    scale_x_continuous(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) +
    xlab("Strength of contact tracing & testing") + ylab("Strength of social distancing"))

ggsave(filename = "Figures/Imported_con_trac_heatmap.pdf",fig_sc_sd,width = 18,height = 16,units = "cm",dpi = 300)




#############################################################################################
####################################### For R0 ##############################################
#############################################################################################

Rt <- function(params){
  
  beta0 <- params["R0IN"]*params["gamma"]/(1-params["omega"] + params["omega"]*params["kappa"])
  rt_tem <- ((1 - params["omega"])*beta0*params["contact"]+ params["omega"]*params["kappa"]*
               (1 - params["zeta"])*beta0*params["contact"])/params["gamma"]
  return(rt_tem)
}


params <- c(omega = 0.93,R0IN = 1.71,kappa = 1.22,gamma = 1/5.08)

params <- data.frame(omega = c(0.84,0.96),R0IN = c(1.27,2.13),kappa = c(0.65,1.65),gamma = c(1/6.03,1/4.09))


res_fin <- data.frame()

for (i in seq(0,1,0.05)){
  for (j in seq(0,1,0.05)){
    tem <- c(contact = i, zeta = j)
    res_list <- lapply(1:100,function(rm){
      params_1 <- apply(params,2,function(id)runif(1,id[1],id[2]))
      params_tem <- c(tem,params_1)
      res_tem <- data.frame(rt = Rt(params = params_tem),ssd = 1 - i, ctt = j)
      return(res_tem)
    })
    res_rm <- bind_rows(res_list)
    res_rm <- data.frame(summarise(group_by(res_rm,ssd,ctt),rt = mean(rt)))
    res_fin <- rbind(res_fin,res_rm)
  }
}


rt_1_p_list <- lapply(unique(res_fin$ssd),function(id){
  #res_fin <- mutate(res_fin,rt = rt*0.7)
  dat_tem <- filter(res_fin, ssd == id & rt < 1.0)
  dat_tem_1 <- arrange(dat_tem,ctt)[1,]
  return(dat_tem_1)
})
rt_1_p <- bind_rows(rt_1_p_list)

(fig_1 <- ggplot(data = res_fin,aes(x = as.character(ssd),y = as.character(ctt))) +
    geom_tile(aes(fill = rt)) +
    scale_fill_gradient2(low = "blue", high = "red", na.value = "transparent", 
                         midpoint = 1,
                         name = expression(R[t]), 
                         guide = guide_colorbar(#title.position = "right", 
                           #title.theme = element_text(angle = -90), 
                           title.hjust = 0.1, 
                           barwidth = 1, 
                           barheight = 14)) +
    geom_point(data = rt_1_p,aes(x = as.character(ssd),y = as.character(ctt)),size =0.5)+
    theme_bw(base_size = 12,base_family = "serif")+
    theme(legend.margin = margin(0,0,0,0,unit = "cm"))+
    scale_y_discrete(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) + 
    scale_x_discrete(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) +
    ylab("Strength of infection isolation") + xlab("Strength of social distancing"))

ggsave(filename = "Figures/Rt.pdf",fig_1,width = 14,height = 12,units = "cm",dpi = 300)


fig_merge <- plot_grid(fig_peak,fig_int,fig_sum,labels = c("A","B","C"),ncol = 3)

fig_all <- plot_grid(fig_merge,fig_1,labels = c("","D"),ncol = 1)
ggsave(filename = "Figures/Rt_all.pdf",fig_all,width = 20,height = 20,units = "cm",dpi = 300)



fig_merge_final <- plot_grid(fig_total_list[[2]]+guides(color = "none"),fig_total_list[[1]],fig_total_perc,fig_1,
                             ncol = 2,labels = c("A","","B","C"))

ggsave(filename = "Figures/Control_Strategy.pdf",fig_merge_final,width = 24,height = 22,units = "cm",dpi = 300)
