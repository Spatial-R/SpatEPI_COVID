library(dplyr)
library(ggplot2)
library(cowplot)
library(ggsci)

re_run <- TRUE

global.result.1 <- read.csv("E:\\nCOV\\Server\\Last\\Hope_UPQ_NF_IVP_K1.csv",header=T) %>% arrange(desc(loglik));
parameter <- arrange(global.result.1,desc(loglik))[1,];
#parameter$qua <- parameter$quawh

if(isTRUE(re_run)){
  source("Codes/City_Network_Intervention_Strong_Total.R")
  print("City_Network_Intervention_Strong_Total")
  source("Codes/City_Network_Intervention_Strong_Confirm.R")
  print("City_Network_Intervention_Strong_Confirm")
  source("Codes/City_Network_Intervention_Strong_UD.R")
  print("City_Network_Intervention_Strong_UD")
  source("Codes/City_Network_Intervention_Weak_Total.R")
  print("City_Network_Intervention_Weak_Total")
  source("Codes/City_Network_Intervention_Weak_Confirm.R")
  print("City_Network_Intervention_Weak_Confirm")
  source("Codes/City_Network_Intervention_Weak_UD.R")
  print("City_Network_Intervention_Weak_UD")
}



cols <- pal_npg("nrc",alpha = 0.6)(9); scales::show_col(cols)
Sys.setlocale(category="LC_ALL",locale = "English")

dat_strong_total <- read.csv(file = "Process_Data/City_Network_Intervention_Strong_Total.csv",stringsAsFactors = F)
dat_weak_total <- read.csv(file = "Process_Data/City_Network_Intervention_Weak_Total.csv",stringsAsFactors = F)
dat_weak_total <- filter(dat_weak_total, day < as.Date("2020-04-01"))
dat_strong_total <- filter(dat_strong_total, day < as.Date("2020-04-01"))

dat_weak_total$ctt <- "Weak infection isolation";dat_strong_total$ctt <- "Strong infection isolation"

dat_total_plot <- rbind(dat_weak_total,dat_strong_total)
dat_total_plot <- mutate(dat_total_plot,day = as.Date(day),type = factor(type,
                                                                         levels = c("No social distancing and \n no travel restrication",
                                                                                    "Travel restrication","Social distancing",
                                                                                    "Social distancing and \n travel restrication")))
dat_total_plot <- filter(dat_total_plot,day < as.Date("2020-04-01"))
dat_total_plot[,c(2:4)] <- apply(dat_total_plot[,c(2:4)],2,log)


fig_strong_total <- ggplot(data = dat_total_plot[dat_total_plot$ctt == "Strong infection isolation",],
                   aes(x = day, y = sim, group = factor(type),color = factor(type))) + 
  geom_pointrange(aes(ymin = low,ymax = high),fatten =1,
                  position = position_dodge(width = 0.4)) +
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 12,base_family = "serif") +
  scale_x_date() +
  scale_y_continuous(limits = c(0,17),breaks = seq(0,167,3))+
  scale_color_manual(values = c(cols[c(3,9,8)],"black"),guide = guide_legend(title = "")) +
  xlab("") + ylab("Daily new COVID-19 cases (Nature logarithm scale)") +
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  theme(legend.position = c(0.25,0.86),
        legend.key.height = unit(1,"cm"))

fig_weak_total <- ggplot(data = dat_total_plot[dat_total_plot$ctt == "Weak infection isolation",],
                           aes(x = day, y = sim, group = factor(type),color = factor(type))) + 
  geom_pointrange(aes(ymin = low,ymax = high),fatten =1,
                  position = position_dodge(width = 0.4)) +
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  theme_classic(base_size = 12,base_family = "serif") +
  scale_x_date() +
  scale_y_continuous(limits = c(0,17),breaks = seq(0,17,3))+
  scale_color_manual(values = c(cols[c(3,9,8)],"black"),guide = guide_legend(title = "")) +
  xlab("") + ylab("Daily new COVID-19 cases (Nature logarithm scale)") +
  theme(legend.position = c(0.25,0.86),
        legend.key.height = unit(1,"cm"))

fig_total <- plot_grid(fig_weak_total,fig_strong_total,labels = c("A","B"))

ggsave(filename = "Figures/City_Network_Intervention_Total.pdf",
       fig_total,width = 24,height = 16,units = "cm",dpi = 300) 




fig_me1 <- ggplot(data = dat_total_plot,
                  aes(x = day, y = sim, group = factor(type),color = factor(type),shape = factor(ctt))) + 
  geom_point(size = 0.8,position = "jitter")+
  scale_shape_manual(values = c(2,1),guide = guide_legend(title = "Infection isolation"))+
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 10,base_family = "serif") +
  scale_x_date() +
  scale_y_continuous(limits = c(0,17),breaks = seq(0,16,3))+
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  theme(legend.spacing.y = unit(1.5,"cm"))+
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),
                     guide = guide_legend(title = "Simulation scaenarios",
                                                                      override.aes = list(size =1.2))) +
  xlab("") + ylab("Daily new COVID-19 cases \n (Nature logarithm scale)")



dat_total_plot_1 <- mutate(dat_total_plot,ctt = gsub(" infection isolation","",ctt))


fig_me1 <- ggplot(data = dat_total_plot_1,
                  aes(x = day, y = sim, group = factor(type),color = factor(type),shape = factor(ctt))) + 
  geom_point(size = 0.8,position = "jitter")+
  scale_shape_manual(values = c(2,1),guide = guide_legend(title = "Infection isolation",order = 1))+
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 10,base_family = "serif") +
  scale_x_date() +
  theme(legend.spacing.y = unit(0.5,"cm"))+
  scale_y_continuous(limits = c(0,17),breaks = seq(0,16,3))+
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),
                     guide = guide_legend(title = "Simulation scenarios",
                                          override.aes = list(size =1.2))) +
  xlab("") + ylab("Daily new COVID-19 cases \n (Nature logarithm scale)")


library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)

fig_legend1 <- ggplot(data = dat_total_plot_1,
                  aes(x = day, y = sim, group = factor(type),color = factor(type),shape = factor(ctt))) + 
  geom_point(size = 0.8,position = "jitter")+
  theme(legend.spacing.y = unit(0,"cm"),legend.key = element_rect(colour = "transparent",fill = "white"))+
  scale_shape_manual(values = c(2,1),guide = "none")+
  theme_classic(base_size = 10,base_family = "serif")+ 
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),
                     guide = guide_legend(title = "Simulation scenarios",
                                          override.aes = list(size =1.2))) 

fig_legend1 <-  gtable_filter(ggplot_gtable(ggplot_build(fig_legend1)), "guide-box") 
  

fig_legend2 <- ggplot(data = dat_total_plot_1,
                      aes(x = day, y = sim, group = factor(type),color = factor(type),shape = factor(ctt))) + 
  geom_point(size = 0.8,position = "jitter")+
  scale_shape_manual(values = c(2,1),guide = guide_legend(title = "Infection isolation",order = 1))+
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 10,base_family = "serif") +
  theme(legend.spacing.y = unit(0,"cm"),legend.key = element_rect(colour = "transparent",fill = "white"))+
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none")

fig_legend2 <-  gtable_filter(ggplot_gtable(ggplot_build(fig_legend2)), "guide-box") 

  
  
fig_me1 <- ggplot(data = dat_total_plot_1,
                  aes(x = day, y = sim, group = factor(type),color = factor(type),shape = factor(ctt))) + 
  geom_point(size = 0.8,position = "jitter")+
  scale_shape_manual(values = c(2,1),guide = guide_legend(title = "Infection isolation",order = 1))+
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 10,base_family = "serif") +
  scale_x_date() +
  theme(legend.spacing.y = unit(0,"cm"),legend.key = element_rect(colour = "transparent",fill = "white"))+
  scale_y_continuous(limits = c(0,17),breaks = seq(0,16,3))+
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),
                     guide = guide_legend(title = "Simulation scenarios",
                                          override.aes = list(size =1.2))) +
  xlab("") + ylab("Daily new COVID-19 cases \n (Nature logarithm scale)")





fig_final <- ggplot(data = dat_total_plot_1,
                  aes(x = day, y = sim, group = factor(type),color = factor(type),shape = factor(ctt))) + 
  geom_point(size = 0.8,position = "jitter")+
  scale_shape_manual(values = c(2,1),guide = "none")+
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 11,base_family = "serif") +
  scale_x_date() +
  theme(legend.spacing.y = unit(0,"cm"))+
  scale_y_continuous(limits = c(0,18),breaks = seq(0,18,3))+
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none") +
  xlab("") + ylab("Daily new COVID-19 cases \n (Nature logarithm scale)")


plotNew <- fig_final + 
  annotation_custom(grob = fig_legend1, xmin = as.Date("2020-02-02"), 
                    xmax = as.Date("2020-02-10"), ymin = 14, ymax = 18)+
  annotation_custom(grob = fig_legend2, xmin = as.Date("2020-02-20"), 
                    xmax = as.Date("2020-03-10"), ymin = 16, ymax = 18)


ggsave(filename = "Figures/City_Network_Intervention_Total_1.pdf",
       plotNew,width = 16,height = 14,units = "cm",dpi = 300) 



##############################################################################################
################################  Undocumented    ############################################
##############################################################################################


dat_strong_ud <- read.csv(file = "Process_Data/City_Network_Intervention_Strong_UD.csv",stringsAsFactors = F)
dat_weak_ud <- read.csv(file = "Process_Data/City_Network_Intervention_Weak_UD.csv",stringsAsFactors = F)

dat_weak_ud$ctt <- "Weakly CTT";dat_strong_ud$ctt <- "Strong  CTT"

dat_ud_plot <- rbind(dat_weak_ud,dat_strong_ud)
dat_ud_plot <- mutate(dat_ud_plot,day = as.Date(day),type = factor(type,
                                                                         levels = c("No social distancing and \n no travel restrication",
                                                                                    "Travel restrication","Social distancing",
                                                                                    "Social distancing and \n travel restrication")))
dat_ud_plot <- filter(dat_ud_plot,day < as.Date("2020-04-01"))
dat_ud_plot[,c(2:4)] <- apply(dat_ud_plot[,c(2:4)],2,log)


fig_strong_ud <- ggplot(data = dat_ud_plot[dat_ud_plot$ctt == "Strong  CTT",],
                           aes(x = day, y = sim, group = factor(type),color = factor(type))) + 
  #geom_pointrange(aes(ymin = low,ymax = high),fatten = 1,position = position_dodge(width = 0.4)) +
  geom_point(size = 0.8,position = "jitter") +
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 12,base_family = "serif") +
  scale_x_date() +
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  scale_y_continuous(limits = c(0,14),breaks = seq(0,14,3))+
  scale_color_manual(values = c(cols[c(3,9,8)],"black"),guide = guide_legend(title = "")) +
  #scale_color_manual(values = c(cols[c(3,9,8)],"black"),guide = "none") 
  xlab("") + ylab("") +
  theme(legend.position = c(0.75,0.85),
        legend.key.height = unit(0.5,"cm"))


fig_weak_ud <- ggplot(data = dat_ud_plot[dat_ud_plot$ctt == "Weakly CTT",],
                         aes(x = day, y = sim, group = factor(type),color = factor(type))) + 
  #geom_pointrange(aes(ymin = low,ymax = high),fatten = 1,position = position_dodge(width = 0.4)) +
  geom_point(size = 0.8,position = "jitter") +
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 12,base_family = "serif") +
  scale_x_date() +
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  scale_y_continuous(limits = c(0,14),breaks = seq(0,14,3))+
  scale_color_manual(values = c(cols[c(3,9,8)],"black"),guide = "none") +
  xlab("") + ylab("Daily new undocumented COVID-19 cases \n (Nature logarithm scale)") +
  theme(legend.position = c(0.7,0.86),
        legend.key.height = unit(0.5,"cm"))

fig_ud <- plot_grid(fig_weak_ud,fig_strong_ud,labels = c("A","B"))

ggsave(filename = "Figures/City_Network_Intervention_UD.pdf",
       fig_ud,width = 20,height = 12,units = "cm",dpi = 300) 



##############################################################################################
###################################  Confirmed    ############################################
##############################################################################################



dat_strong_confirm <- read.csv(file = "Process_Data/City_Network_Intervention_Strong_Confirm.csv",stringsAsFactors = F)
dat_weak_confirm <- read.csv(file = "Process_Data/City_Network_Intervention_Weak_Confirm.csv",stringsAsFactors = F)

dat_weak_confirm$ctt <- "Weakly CTT";dat_strong_confirm$ctt <- "Strong  CTT"

dat_confirm_plot <- rbind(dat_weak_confirm,dat_strong_confirm)
dat_confirm_plot <- mutate(dat_confirm_plot,day = as.Date(day),type = factor(type,
                                                                   levels = c("No social distancing and \n no travel restrication",
                                                                              "Travel restrication","Social distancing",
                                                                              "Social distancing and \n travel restrication")))
dat_confirm_plot <- filter(dat_confirm_plot,day < as.Date("2020-04-01"))
dat_confirm_plot[,c(2:4)] <- apply(dat_confirm_plot[,c(2:4)],2,log)

fig_strong_confirm <- ggplot(data = dat_confirm_plot[dat_confirm_plot$ctt == "Strong  CTT",],
                        aes(x = day, y = sim, group = factor(type),color = factor(type))) + 
  geom_pointrange(aes(ymin = low,ymax = high),fatten = 1,
                  position = position_dodge(width = 0.4)) +
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 12,base_family = "serif") +
  scale_x_date() +
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  scale_y_continuous(limits = c(0,17),breaks = seq(0,167,3))+
  scale_color_manual(values = c(cols[c(3,9,8)],"black"),guide = guide_legend(title = "")) +
  xlab("") + ylab("Number of COVID-19 cases (Nature logarithm scale)") +
  theme(legend.position = c(0.25,0.86),
        legend.key.height = unit(1,"cm"))

fig_weak_confirm <- ggplot(data = dat_confirm_plot[dat_confirm_plot$ctt == "Weakly CTT",],
                      aes(x = day, y = sim, group = factor(type),color = factor(type))) + 
  geom_pointrange(aes(ymin = low,ymax = high),fatten = 1,
                  position = position_dodge(width = 0.4)) +
  #geom_line(data = dat_cases,aes(x = date, y = cases),color = "black") +
  theme_classic(base_size = 12,base_family = "serif") +
  scale_x_date() +
  geom_vline(xintercept = as.Date("2020-01-24"),linetype = 2)+
  scale_y_continuous(limits = c(0,17),breaks = seq(0,167,3))+
  scale_color_manual(values = c(cols[c(3,9,8)],"black"),guide = guide_legend(title = "")) +
  xlab("") + ylab("Number of COVID-19 cases (Nature logarithm scale)") +
  theme(legend.position = c(0.25,0.86),
        legend.key.height = unit(1,"cm"))

fig_confirm<- plot_grid(fig_weak_confirm,fig_strong_confirm,labels = c("A","B"))

ggsave(filename = "Figures/City_Network_Intervention_Confirm.pdf",
       fig_confirm,width = 24,height = 16,units = "cm",dpi = 300) 


dat_day$low <- ifelse(dat_day$low < 1,1,dat_day$low)

dat_weak_total <- filter(dat_weak_total,day > as.Date("2020-01-23"))
dat_strong_total <- filter(dat_strong_total,day > as.Date("2020-01-23"))

dat_sd_weak <- apply(dat_weak_total[dat_weak_total$type == "Social distancing",2:4],2,sum)
dat_sdtr_weak <- apply(dat_weak_total[dat_weak_total$type == "Social distancing and \n travel restrication",2:4],2,sum)
dat_n_weak <- apply(dat_weak_total[dat_weak_total$type == "No social distancing and \n no travel restrication",2:4],2,sum)
dat_tr_weak <- apply(dat_weak_total[dat_weak_total$type == "Travel restrication",2:4],2,sum)

dat_sd_strong <- apply(dat_strong_total[dat_strong_total$type == "Social distancing",2:4],2,sum)
dat_sdtr_strong <- apply(dat_strong_total[dat_strong_total$type == "Social distancing and \n travel restrication",2:4],2,sum)
dat_n_strong <- apply(dat_strong_total[dat_strong_total$type == "No social distancing and \n no travel restrication",2:4],2,sum)
dat_tr_strong <- apply(dat_strong_total[dat_strong_total$type == "Travel restrication",2:4],2,sum)


(dat_n_weak[1] - dat_sd_weak[1])/dat_n_weak[1]           ### social distancing 
(dat_n_weak[1] - dat_sd_strong[1])/dat_n_weak[1]         ### social distancing and contact tracing
(dat_n_weak[1] - dat_n_strong[1])/dat_n_weak[1]          ### contact tracing 

(dat_n_weak[1] - dat_sd_weak[1])/dat_n_weak[1]           ### social distancing 
(dat_sd_strong[1] - dat_sdtr_strong[1])/dat_n_strong[1]    ### social distancing and contact tracing
(dat_tr_strong[1] - dat_sdtr_strong[1])/dat_n_strong[1]    ### social distancing and contact tracing

(dat_n_weak[1] - dat_n_strong[1])/dat_n_weak[1]          ### contact tracing 


(dat_n_weak[1] - dat_tr_weak[1])/dat_n_weak[1]           ### social distancing 
(dat_n_strong[1] - dat_tr_strong[1])/dat_n_strong[1]     ### social distancing and contact tracing
(dat_n_weak[1] - dat_n_strong[1])/dat_n_weak[1]          ### contact tracing 


(dat_tr_weak[1] - dat_sdtr_weak[1])/dat_n_weak[1]
(dat_sd_weak[1] - dat_sdtr_weak[1])/dat_n_weak[1]


(dat_n_strong[1] - dat_sd_strong[1])/dat_n_strong[1]
(dat_n_strong[1] - dat_sdtr_strong[1])/dat_n_strong[1]
(dat_n_strong[1] - dat_tr_strong[1])/dat_n_strong[1]

(dat_n_weak[1] - dat_n_strong[1])/dat_n_weak[1]              #### contact tracing
(dat_n_weak[1] - dat_sdtr_strong[1])/dat_n_weak[1]           #### contact tracing


(dat_tr_weak[1] - dat_sdtr_weak[1])/dat_n_weak[1]            #### social distancing 
(dat_sd_weak[1] - dat_sdtr_weak[1])/dat_n_weak[1]            #### 
(dat_tr_strong[1]  - dat_sdtr_strong[1])/dat_n_strong[1]     #### relarity
(dat_sd_strong[1]  - dat_sdtr_strong[1])/dat_n_strong[1]     #### travel



(dat_sdtr_weak[1]  -  dat_sdtr_strong[1])/dat_sdtr_weak[1]   #### infection isolation
(dat_sd_weak[1]    -  dat_sd_strong[1])/dat_sd_weak[1]       #### social distancing
(dat_tr_weak[1]    -  dat_tr_strong[1])/dat_tr_weak[1]       #### with travel restriction
(dat_n_weak[1]     -  dat_n_strong[1])/dat_n_weak[1]         #### no other interventions



which.max(dat_strong_total[dat_strong_total$type == "Social distancing",2])
which.max(dat_strong_total[dat_strong_total$type == "No social distancing and \n no travel restrication",2])

which.max(dat_weak_total[dat_weak_total$type == "Social distancing",2])
which.max(dat_weak_total[dat_weak_total$type == "No social distancing and \n no travel restrication",2])




#############################################################################################
####################################  region number #########################################
#############################################################################################

library(ggspatial)

china_map <- sf::st_read("Data/Maps/bou2_4p.shp")
sf::st_crs(china_map) <- "+proj=longlat +datum=WGS84 +no_defs"
china_map$NAME <- iconv(china_map$NAME,"gb18030","UTF-8")
dat_coord <- read.csv("Data/Demo/Citycode.csv",header = T,stringsAsFactors = F)


dat_origin_city <- data.frame(summarise(group_by(dat_origin,city),
                                        sim = floor(sum(median)),low = floor(sum(low)),high = floor(sum(high))))
dat_fc_city <- data.frame(summarise(group_by(dat_fc,city),
                                    sim = floor(sum(median)),low = floor(sum(low)),high = floor(sum(high))))
dat_fcct_city <- data.frame(summarise(group_by(dat_fcct,city),
                                      sim = floor(sum(median)),low = floor(sum(low)),high = floor(sum(high))))
dat_ct_city <- data.frame(summarise(group_by(dat_ct,city),
                                    sim = floor(sum(median)),low = floor(sum(low)),high = floor(sum(high))))


sum_dat_origin <- merge(dat_origin_city,dat_coord,by.x = "city",by.y = "name",all.x = T)
sum_dat_origin <- filter(sum_dat_origin,!city %in% test_name)
sum_dat_origin <- mutate(sum_dat_origin,sim = ifelse(sim < 1 ,NA,sim))


sum_dat_fc <- merge(dat_fc_city,dat_coord,by.x = "city",by.y = "name",all.x = T)
sum_dat_fc <- filter(sum_dat_fc,!city %in% test_name)
sum_dat_fc <- mutate(sum_dat_fc,sim = ifelse(sim < 1 ,NA,sim))


sum_dat_fcct <- merge(dat_fcct_city,dat_coord,by.x = "city",by.y = "name",all.x = T)
sum_dat_fcct <- filter(sum_dat_fcct,!city %in% test_name)
sum_dat_fcct <- mutate(sum_dat_fcct,sim = ifelse(sim < 1 ,NA,sim))


sum_dat_ct <- merge(dat_ct_city,dat_coord,by.x = "city",by.y = "name",all.x = T)
sum_dat_ct <- filter(sum_dat_ct,!city %in% test_name)
sum_dat_ct <- mutate(sum_dat_ct,sim = ifelse(sim < 1 ,NA,sim))



fig_origin <- ggplot(data = china_map) +
  geom_sf(fill = "gray88") +
  coord_sf()+
  geom_point(data = sum_dat_origin,aes(x = lat,y = lon,size = sim),color = cols[8])+
  scale_size_continuous(guide = guide_legend(title = "Cases"),range = c(0.1,6)) +
  #scale_color_manual(values = rev(cols[-2]),guide = F)+
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.35, "in"), pad_y = unit(0.3, "in"),
                         height = unit(0.3, "in"), width = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = c(0.9,0.2))


fig_fcct <- ggplot(data = china_map) +
  geom_sf(fill = "gray88") +
  coord_sf()+
  geom_point(data = sum_dat_fcct,aes(x = lat,y = lon,size = sim),color = cols[8])+
  scale_size_continuous(guide = guide_legend(title = "Cases"),range = c(0.1,6)) +
  #scale_color_manual(values = rev(cols[-2]),guide = F)+
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.35, "in"), pad_y = unit(0.3, "in"),
                         height = unit(0.3, "in"), width = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = c(0.9,0.2))




fig_fc <- ggplot(data = china_map) +
  geom_sf(fill = "gray88") +
  coord_sf()+
  geom_point(data = sum_dat_fc,aes(x = lat,y = lon,size = sim),color = cols[8])+
  scale_size_continuous(guide = guide_legend(title = "Cases"),range = c(0.1,6)) +
  #scale_color_manual(values = rev(cols[-2]),guide = F)+
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.35, "in"), pad_y = unit(0.3, "in"),
                         height = unit(0.3, "in"), width = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = c(0.9,0.2))


fig_ct<- ggplot(data = china_map) +
  geom_sf(fill = "gray88") +
  coord_sf()+
  geom_point(data = sum_dat_ct,aes(x = lat,y = lon,size = sim),color = cols[8])+
  scale_size_continuous(guide = guide_legend(title = "Cases"),range = c(0.1,6)) +
  #scale_color_manual(values = rev(cols[-2]),guide = F)+
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.35, "in"), pad_y = unit(0.3, "in"),
                         height = unit(0.3, "in"), width = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = c(0.9,0.2))

map_merge <- plot_grid(fig_ct,fig_fcct,fig_origin,fig_fc,labels = c("A","B","C","D"))

ggsave(filename = "Figures/City_Net_Intervention_Strong_Map_Total.pdf",map_merge,width = 22,height = 18,units = "cm",dpi = 300 )

