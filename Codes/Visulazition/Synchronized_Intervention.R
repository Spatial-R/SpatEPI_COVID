#remove(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringi)
library(ggsci)
library(pomp)
library(Hmisc)
library(patchwork)
load("Other.RData")
load("Merge.RData")
text_size <- 8

Sys.setlocale("LC_TIME","en_US.UTF-8")
#Sys.setenv(LANG = "en_US.UTF-8")
source("City_Network_Sim_Function.R")
dat_top_list_1 <- lapply(unique(plot_dat_top$delay),function(id){
  dat_tem <- filter(plot_dat_top,delay == id)
  dat_tem <- filter(plot_dat_top,as.Date(day) > (as.Date("2020-01-23")+as.numeric(as.character(id))))
  return(dat_tem)
})
dat_top_2 <- bind_rows(dat_top_list_1)
names(dat_top_2)[5] <- "sim"
plot_dat_top <- data.frame(summarise(group_by(dat_top_2,top,delay,region),sim1 = mean(sim),sd = sd(sim)))
plot_dat_top_1 <- data.frame(summarise(group_by(plot_dat_top,top,delay,region),sim1 = sum(sim1)))
plot_dat_top_1 <- tidyr::spread(plot_dat_top_1,top, sim1)
plot_dat_top_1 <- arrange(plot_dat_top_1,as.numeric(as.character(delay)))
1 - plot_dat_top_1[,4:(ncol(plot_dat_top_1))]/plot_dat_top_1[,3]

plot_dat_top_1$delay <- as.numeric(plot_dat_top_1$delay) + as.Date("2020-01-23")
plot_dat_top_2 <- plot_dat_top_1
plot_dat_top_2[,4:(ncol(plot_dat_top_2))] <- 1 - plot_dat_top_2[,4:(ncol(plot_dat_top_2))]/plot_dat_top_2[,3]
plot_dat_top_2 <- plot_dat_top_2[,-3]
plot_dat_top_2 <- reshape::melt(plot_dat_top_2,id = c("delay","region"))
plot_dat_top_2 <- mutate(plot_dat_top_2,
                         #delay = factor(delay,levels = seq(0,20,2)),
                         variable = factor(variable,levels = seq(10,50,10)) )
plot_dat_top_select <- filter(plot_dat_top_2, region %in% c("wuhan","beijing","chengdu"))
plot_dat_top_select$region <- capitalize(plot_dat_top_select$region)


plot_list <- lapply(c("Wuhan","Chengdu","Beijing"),function(id){
  plot_dat_top_tem <- filter(plot_dat_top_select,region == id)
  
  if(id == "Beijing"){
    x_tem <- xlab("")
  } else {
    x_tem <- xlab("Intervention Time")
  }
  
  ggplot(data = plot_dat_top_tem,aes(x = delay, y = variable)) + 
    geom_tile(aes(fill = value)) +
    x_tem + ylab("Number of synchronized regions") +
    
    ggtitle(id)+
    scale_x_date(breaks = seq(as.Date("2020-01-23"),as.Date("2020-02-12"),by = "4 days"),
                 labels = c("Jan-23","Jan-27","Jan-31","Feb-04","Feb-08","Feb-12"),
                 expand = c(0,0)) +
    scale_fill_gradient2(low="#2c7fb8", high="red",na.value = "white",midpoint = 0.25,
                         name = "Proportion of infection averted",
                         limits = c(0,0.5),
                         guide = guide_colorbar(title.position = "top", 
                                                title.theme = element_text(angle = 0,size = text_size), 
                                                title.hjust = 0.5, 
                                                barwidth = 10, 
                                                barheight = 0.4)) +
    theme(legend.position = "bottom",text = element_text(size = text_size),
          legend.box.margin = margin(-12.5,0,2,0)
          )
})
#Sys.setlocale("LC_MESSAGES", 'en_GB.UTF-8')
p_no_legend <- lapply(plot_list, function(x) x + theme(legend.position = "none"))
legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom"))
p_grid <- p_no_legend[[1]] + p_no_legend[[2]] + p_no_legend[[3]] + plot_layout(ncol = 1)
plot.out <- cowplot::plot_grid(p_grid, legend, ncol = 1, nrow = 2, 
                               rel_heights = c(40,1), rel_widths = c(6,0.9))
ggsave(filename = "Region_select.pdf",width = 30,height = 40,units = "cm",dpi = 300)


################################# select the less 

library(statnet)
load("Cities_Spatial_All.RData")
network_dat <- flow_all[1:U,1:U]
network_dat$city <-  names(flow_all)
flow_netwwork <- reshape::melt(network_dat,id = "city")
names(flow_netwwork) <- c("source","target","weight")
flow_netwwork <- mutate(flow_netwwork,source = as.character(source),
                        target = as.character(target))
nodes <- data.frame(id = names(flow_all))
network <- igraph::graph_from_data_frame(vertices = nodes,d = flow_netwwork,
                                         directed = T)
network_statnet <- intergraph::asNetwork(network)

(degree_weighted <- myDegree(network_statnet, gmode = "digraph", ignore.eval = FALSE, attrname = "weight",rescale = T))

degree_dat <- data.frame(degree = degree_weighted,city = nodes$id)

degree_dat_1 <- filter(degree_dat,city %in% unique(plot_dat_top_2$region))
degree_dat_1 <- arrange(degree_dat_1,degree)
degree_dat_1 <- mutate(degree_dat_1,city = factor(degree_dat_1$city,levels = degree_dat_1$city))

cut_off <- as.numeric(quantile(degree_dat_1$degree,probs = c(0,0.25,0.75,1.00)))
degree_dat_1 <- mutate(degree_dat_1,group = cut(degree,breaks = cut_off,include.lowest = T))
degree_dat_1$city <- as.character(degree_dat_1$city)
degree_dat_1$city <- capitalize(degree_dat_1$city)
degree_dat_1$city <- factor(degree_dat_1$city,levels = degree_dat_1$city)
save(degree_dat_1,file = "degree.RData")

fig_degree_city <- ggplot(data = degree_dat_1,aes(x = factor(city), y = degree)) +
  geom_point(aes(color = factor(group))) + coord_flip() + 
  geom_hline(yintercept = cut_off[-1],linetype = 2)+
  scale_color_manual(values = cols[c(3,8,9)],guide = NULL)+
  theme_classic(base_family = "serif",base_size = text_size) +
  xlab("Capital cities") + ylab("Interconnectivity") +
  annotate("text",x = c(25,25,25),y = c(0.0025,0.009,0.018),label = c("Small \n interconnectivity",
                                                                      "Middle \n interconnectivity","Large \n interconnncetivity"))
ggsave(fig_degree_city,filename = "Region_Degree.pdf",width = 20,height = 15,units = "cm",dpi = 300)

############################  test the connection with effectiveness  #################################

lm_list <- lapply(unique(plot_dat_top_2$delay),function(id){
  dat_tem <- filter(plot_dat_top_2,delay == id)
  reg_list <- lapply(unique(dat_tem$variable),function(id_tem){
    dat_tem1 <- filter(dat_tem,variable == id_tem)
    dat_tem2 <- merge(dat_tem1,degree_dat,by.x = "region",by.y = "city")
    lm_res <- summary(lm(value ~ degree,data = dat_tem2))
    data.frame(var = id_tem, delay = id,est = lm_res$coefficients[2,1],sd = lm_res$coefficients[2,2],
               r2 = lm_res$r.squared)
  })
  reg_dat <- bind_rows(reg_list)
  return(reg_dat)
})

lm_dat <- bind_rows(lm_list)

lm_dat <- mutate(lm_dat,low = est - 1.96*sd, high = est + 1.96*sd)

plot_dat_top_3 <- mutate(plot_dat_top_2,variable_1  = paste("n = ",variable,sep = ""))


fig_eff <- ggplot(data = plot_dat_top_3,
                  aes(x = reorder(format(delay,'%b-%d'),delay), y = value)) +
  geom_boxplot(aes(fill = factor(variable_1)),outlier.size = 0.5) + 
  theme_classic(base_family = "serif",base_size = text_size) +
  #scale_x_date()+
  scale_x_discrete(breaks = c("Jan-23","Jan-27","Jan-31","Feb-04","Feb-08","Feb-12"),
               labels = c("Jan-23","Jan-27","Jan-31","Feb-04","Feb-08","Feb-12"))+
  scale_fill_manual(values = cols[c(3:5,9,8)])+ 
  xlab("Intervention time") + ylab("Proportion of infection averted") +
  guides(fill = guide_legend(title = "",
                             keyheight = 0.8,
                             title.position = "right",title.theme = element_text(angle = -90,size = text_size)))+
  theme(legend.position = c(0.85,0.8)) 


lm_dat <- mutate(lm_dat,var = paste("n = ",var,sep = ""))

(fig_degree <- ggplot(data = lm_dat,aes(x = delay, y = est)) +
  geom_pointrange(aes(ymin = low, ymax = high,group = factor(var),color = factor(var)),
                  position = position_dodge(width = 1.8),size = 0.5,fatten = 1) +
  theme_classic(base_family = "serif",base_size = text_size) + 
  scale_x_date(breaks = seq(as.Date("2020-01-23"),as.Date("2020-02-12"),by = "4 days"),
                     labels = c("Jan-23","Jan-27","Jan-31","Feb-04","Feb-08","Feb-12")) +
  scale_color_manual(values = cols[c(3:5,9,8)],
                     guide = guide_legend(title = "",title.position = "right",
                                          keyheight = 0.8,
                                          title.theme = element_text(angle = -90,size = text_size))) +
  xlab("Intervention time") + ylab("Estimated coefficient") +
  geom_hline(yintercept = 0,linetype = 2) +
  theme(legend.position = c(0.85,0.2)))

#ggsave(filename = "Region_Degree.pdf",width = 20,height = 15,units = "cm",dpi = 300)


cols <- pal_npg("nrc",alpha = 0.6)(9); scales::show_col(cols)


fig_two <- plot_grid(fig_degree,fig_eff,labels = c("(b)","(c)"),label_fontfamily = "serif",label_fontface = "plain",ncol = 1)
fig_merge <- plot_grid(plot.out,fig_two,labels = c("(a)",""),label_fontfamily = "serif",label_fontface = "plain",ncol = 2)

ggsave(fig_merge,filename = "merge.tiff",width = 16,height = 16,units = "cm",dpi = 600)

