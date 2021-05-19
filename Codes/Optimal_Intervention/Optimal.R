library(dplyr)
library(ggplot2)

Rt <- function(params){
 rt_tem <- ((1 - params["omega"])*params["R0IN"]*params["contact"]+ params["omega"]*params["kappa"]*
    (1 - params["zeta"])*params["R0IN"]*params["contact"])/params["gamma"]
 return(rt_tem)
}

params <- c(omega = 0.93,R0IN = 1.7,kappa = 1.2,gamma = 1/5)

res_fin <- data.frame()

for (i in seq(0,1,0.01)){
  for (j in seq(0,1,0.01)){
    tem <- c(contact = i, zeta = j)
    params_tem <- c(tem,params)
    res_tem <- data.frame(rt = Rt(params = params_tem),ssd = 1 - i, ctt = j)
    res_fin <- rbind(res_fin,res_tem)
  }
}


rt_1_p_list <- lapply(unique(res_fin$ssd),function(id){
  #res_fin <- mutate(res_fin,rt = rt*0.7)
  dat_tem <- filter(res_fin, ssd == id & rt < 1.1)
  dat_tem_1 <- arrange(dat_tem,ctt)[1,]
  return(dat_tem_1)
})
rt_1_p <- bind_rows(rt_1_p_list)

fig_1 <- ggplot(data = res_fin,aes(x = as.character(ssd),y = as.character(ctt))) +
  geom_tile(aes(fill = rt)) +
  scale_fill_gradient2(low = "blue", high = "red", na.value = "transparent", 
                       midpoint = 1,
                       name = "Rt", 
                       guide = guide_colorbar(#title.position = "right", 
                                              #title.theme = element_text(angle = -90), 
                                              title.hjust = 0.1, 
                                              barwidth = 1, 
                                              barheight = 18)) +
  geom_point(data = rt_1_p,aes(x = as.character(ssd),y = as.character(ctt)),
             size =0.5)+
  scale_y_discrete(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) + 
  scale_x_discrete(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) +
  ylab("The strength of contact tracing & testing") + xlab("The strength of social distancing")




est_beta <- read.csv("Processed_Data/Estimated_beta0.csv")
est_beta <- filter(est_beta,!is.na(growrate))
iso_name <- as.character(unique(est_beta$iso3))
est_beta$growrate <- est_beta$growrate/max(est_beta$growrate)

res_fin <- mutate(res_fin,ssd = as.character(ssd), ctt = as.character(ctt))

week_fin_list <-  lapply(iso_name, function(iso_id){
    dat_tem <- filter(est_beta,iso3 == iso_id)
    print(iso_id)
    ssd_list <- lapply(unique(res_fin$ssd),function(id){
    week_list <- lapply(1:nrow(dat_tem), function(grow_id){
          res_fin_1 <- mutate(res_fin,rt = rt*dat_tem[grow_id,"growrate"])
          dat_tem_1 <- filter(res_fin_1, ssd == id & rt < 1)
          dat_tem_2 <- arrange(dat_tem_1,ctt)[1,]
          dat_tem_2$week <- grow_id
          return(dat_tem_2)
        })
       week_res <- bind_rows(week_list)
       return(week_res)
      })
    ssd_res <- bind_rows(ssd_list)
    ssd_res$iso3 <- iso_id
    return(ssd_res)
})

week_fin <- bind_rows(week_fin_list)

week_fin_ssd0 <- filter(week_fin,ssd == "0")
week_fin_ssd0$iso3 <- factor(week_fin_ssd0$iso3,levels = c(dat_coord_1$iso3),labels = c(1:length(dat_coord_1$iso3)))

equa_pos <-   which.min(abs(dat_coord_1$lat - 0))
north_tropic <-  which.min(abs(dat_coord_1$lat - 23.4))
south_tropic <-  which.min(abs(dat_coord_1$lat + 23.4))
north_template <-  which.min(abs(dat_coord_1$lat - 66.5))
south_template <-  which.min(abs(dat_coord_1$lat + 66.5))

template_region0 <- filter(week_fin_ssd0,iso3 %in% (dat_coord_1[dat_coord_1$lat > 23.4,"iso3"]))

template_region0 <- filter(week_fin_ssd0,iso3 %in% 
                             (dat_coord_1[dat_coord_1$lat > 23.4 & dat_coord_1$lat < 66.5,"iso3"]))
template_region0 <- filter(week_fin_ssd0,iso3 %in% 
                             (dat_coord_1[dat_coord_1$lat > 0 & dat_coord_1$lat < 23.4,"iso3"]))

ggplot(data = week_fin_ssd0,aes(x = factor(week), y = as.numeric(ctt))) +
  geom_boxplot()




ggplot(data = template_region0,aes(x = as.character(ssd),y = as.character(ctt))) +
  geom_tile(aes(fill = rt)) +
  scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                       # midpoint = 3.5,
                       name = "Rt", 
                       guide = guide_colorbar(#title.position = "right", 
                         #title.theme = element_text(angle = -90), 
                         title.hjust = 0.1, 
                         barwidth = 1, 
                         barheight = 18)) +
  geom_point(data = rt_1_p,aes(x = as.character(ssd),y = as.character(ctt)),
             size =0.5)+
  scale_y_discrete(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) + 
  scale_x_discrete(expand = c(0, 0),breaks = seq(0.1,0.9,0.1)) +
  ylab("The strength of contact tracing & testing") + xlab("The strength of social distancing")


(fig_tem_ah <- ggplot(data = week_fin_ssd0, aes(y = iso3, x = week)) + 
    geom_tile(aes(fill =as.numeric(ctt))) + 
    geom_hline(yintercept = c(yintercept= c(equa_pos,north_tropic,south_tropic)),linetype = 2)+
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = 0.35,
                         name = "", 
                         guide = guide_colorbar(title.position = "right", 
                                                title.theme = element_text(angle = -90), 
                                                title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 26)) + 
    scale_y_discrete(expand = c(0, 0),breaks = c(equa_pos,north_tropic,south_tropic),
                     labels =c("Equator","Tropic of Cancer","Tropic of Capricorn")) + 
    #geom_point(data = max_date,aes(x = week, y = iso3))+
    scale_x_continuous(breaks = seq(1,53,2),expand = c(0, 0)) +
    labs(x = "Time (Weeks)",y ="Counrties & regions") + 
    theme_bw(base_family = "serif",base_size = 12) + 
    theme(legend.position = "right", 
          legend.margin = margin(0.5, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -9), 
          #axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = .5), 
          plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
          plot.title = element_text(vjust = 5, hjust = .5)))


week_fin_ssd1 <- filter(week_fin,ssd == "0.1")
week_fin_ssd1$iso3 <- factor(week_fin_ssd1$iso3,levels = c(dat_coord_1$iso3),labels = c(1:length(dat_coord_1$iso3)))


(fig_tem_ah <- ggplot(data = week_fin_ssd1, aes(y = iso3, x = week)) + 
    geom_tile(aes(fill =as.numeric(ctt))) + 
    geom_hline(yintercept = c(yintercept= c(equa_pos,north_tropic,south_tropic)),linetype = 2)+
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = 0.5,
                         name = "Growth rate", 
                         guide = guide_colorbar(title.position = "right", 
                                                title.theme = element_text(angle = -90), 
                                                title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 26)) + 
    scale_y_discrete(expand = c(0, 0),breaks = c(equa_pos,north_tropic,south_tropic),
                     labels =c("Equator","Tropic of Cancer","Tropic of Capricorn")) + 
    #geom_point(data = max_date,aes(x = week, y = iso3))+
    scale_x_continuous(breaks = seq(1,53,2),expand = c(0, 0)) +
    labs(x = "Time (Weeks)",y ="Counrties & regions") + 
    theme_bw(base_family = "serif",base_size = 12) + 
    theme(legend.position = "right", 
          legend.margin = margin(0.5, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -9), 
          #axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = .5), 
          plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
          plot.title = element_text(vjust = 5, hjust = .5)))
