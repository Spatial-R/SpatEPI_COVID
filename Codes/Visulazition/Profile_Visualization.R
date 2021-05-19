library(dplyr)
library(ggplot2)
library(ggsci)
library(stringi)
library(cowplot)
library(reshape2)
cols <- pal_npg("nrc",alpha = 0.6)(9); scales::show_col(cols)


#########################################################################################
################# Paper: Monte Carlo profile confidence intervals #######################
################ https://arxiv.org/abs/1612.02710  ######################################
#########################################################################################

mcap <- function(lp,parameter,confidence = 0.95,lambda = 0.75,Ngrid = 1000){
  smooth_fit <- loess(lp ~ parameter,span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0,length(parameter))
  weight[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp,b=parameter,a=-parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  var_b <- m["b","b"]
  var_a <- m["a","a"]
  cov_ab <- m["a","b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence,df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp,parameter=parameter,confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(
         parameter=parameter_grid,
         smoothed=smoothed_loglik,
         quadratic=predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
       ),
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), se=sqrt(se_total_squared)
  )
}


first_up <- function(data){
  substr(data,1,1) <- toupper(substr(data,1,1))
  return(data)
}




pro_dir <- "E:/nCOV/Server/Last/Hubei_Profile"

file_dir <- list.files(path = pro_dir,pattern = "csv",full.names = T)
file_dir <- file_dir[stri_detect_fixed(file_dir,"Hope_UHQ_NF_IVP_K")]

q_res_list <- lapply(file_dir,function(id){
  
  dat.tem <- read.csv(id,stringsAsFactors = F)
  #dat.tem %>% 
  #  subset(is.finite(loglik) & nfail.max == 0,-c(nfail.max,nfail.min)) -> dat.tem
  region_id <- gsub(pro_dir,"",id);region_id <- gsub(").csv","",region_id);
  region_id <- gsub("/Hope_UHQ_NF_IVP_K","",region_id);
  region_id <- substr(region_id,2,nchar(region_id))
  names(dat.tem)[which(names(dat.tem) == region_id)] <- "target"
  
  dat.tem %>%
    plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.pred
  #dat.pred <- dat.tem
  dat.pred <- filter(dat.pred,loglik > -1700)
  para.tem <- mcap(lp =dat.pred$loglik,parameter = dat.pred$target)
  plot(dat.pred$loglik~dat.pred$target)
  
  res_tem <- data.frame(low = para.tem$ci[1],
                        high = para.tem$ci[2],
                        median = para.tem$mle,
                        param = region_id)
})

seq_res <- do.call("rbind",q_res_list)

#############################################################################################
#########################  Parameter estimation for Hubei province ##########################
#############################################################################################


pro_dir <- "E:/nCOV/Server/Last/Hubei_Profile"

file_dir <- list.files(path = pro_dir,pattern = "csv",full.names = T)
file_ms <- file_dir[stri_detect_fixed(file_dir,"Hope_UHQ_NF_IVP_K")]

ms_res_list <- lapply(file_ms,function(id){
  
  dat.tem <- read.csv(id,stringsAsFactors = F)
  dat.tem %>% 
    subset(is.finite(loglik)) -> dat.tem
  region_id <- gsub(pro_dir,"",id);region_id <- gsub(").csv","",region_id);
  region_id <- gsub("/Hope_UHQ_NF_IVP_K","",region_id);
  region_id <- stri_split_fixed(region_id,"(")[[1]]
  param <- region_id[2];
  names(dat.tem)[which(names(dat.tem) == param)] <- "target"
  #dat.tem$target <- as.character(dat.tem$target)
  dat.tem %>%
    plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.pred
  #dat.pred <- dat.tem
  dat.pred$parm <- param
  return(dat.pred)
})

ms_res <- bind_rows(ms_res_list)



Profile.Plot <- function(dataset,smooth.type = F,profile.figure.dir = NULL,
                         alpha = 0.95, p.value, seq.long,mcap.t =T,
                         parm.name = c("omega","kappa","qua","quawh","qualow","R0IN",
                                       "sigma","gamma")){
  
  
  cols <- pal_lancet("lanonc",alpha = 0.5)(8); #show_col(cols)
  
  type.name <- unique(dataset$type)
  
  if(is.null(profile.figure.dir)){
    profile.figure.dir <- getwd()
  } 
  
  
  plot.total <- list(); max_dat <- data.frame()
  
  for (i in (1:length(parm.name))){
    
    dat.all.1 <- filter(dataset, (parm == parm.name[i]) & loglik > -1900 )  
  
    # if(parm.name[i] == "omega"){
    #   dat.all.1 <- filter(dat.all.1,loglik > -1740 & target < 0.96)
    # } else if (parm.name[i] == "kappa"){
    #   dat.all.1 <- filter(dat.all.1,target < 0.8 & loglik > - 1730)
    # } else if (parm.name[i] == "qualow"){
    #   dat.all.1 <- filter(dat.all.1,target < 0.8)
    # } else {
    #   dat.all.1 <- dat.all.1
    # }
    # 
    if(isTRUE(mcap.t)){
      max.x  <-  mcap(dat.all.1$loglik,parameter = dat.all.1$target)$mle;
      dat.ci <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$ci
      dat.plot <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$fit
      names(dat.plot) <- c("x","y","z")
      x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
      dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
      lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
      
    } else{
      
      fit <- locfit.robust(loglik~target,data = dat.all.1,family = "qrgauss",alpha = alpha);
      dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
      dat.plot <- data.frame(y = predict(fit,dat.all.1$target),x = dat.all.1$target)
      
      if(isTRUE(smooth.type)){
        max.x <- dat.plot[which.max(dat.plot$y),"x"]
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
        x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
        x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
        n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
        n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
        x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
      }else{
        max.x  <-  dat.all.1[which.max(dat.all.1$loglik),"target"];
        lik.95 <-  max(dat.all.1$loglik) - 0.5*qchisq(p = p.value,df = 1)
        x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
        x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
        n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
        n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
        x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2];
      }}
    
    # 
    # if(parm.name[i] == "qualow"){
    #   title_text <- paste("MLE: ",fill.character(round(max.x,2),6)," (",
    #                       fill.character(round(x.value.1,2),6),", ",
    #                       fill.character(round(x.value.2,2),6),")",sep = "")
    # } else {
    #   title_text <-  paste("MLE: ",fill.character(round(max.x,2),4)," (",
    #                        fill.character(round(x.value.1,2),4),", ",
    #                        fill.character(round(x.value.2,2),4),")",sep = "")
    # }
    
    max_dat <- rbind(max_dat,data.frame(parm = parm.name[i],
                                        max = max.x,low = x.value.1,
                                        high = x.value.2))
    
    if (all(dat.all.1$parm == "quawh")) {
      #  xlab.label <- xlab(expression(paste("AH on TR: ", eta[AH],sep = "")))
      xlab.label <- expression(paste("Growth rate of CTT in Wuhan (",zeta[wh],")",sep = ""))
    } else if(all(dat.all.1$parm == "qua")){
      # xlab.label <- xlab(expression(paste("MT on TR: ",eta[MT],sep = "")))
      xlab.label <- expression(paste("Growth rate of CTT in Hubei* (",zeta[e],")",sep = ""))
    } else if(all(dat.all.1$parm == "qualow")){
      # xlab.label <- xlab(expression(paste(" MT on RR: ",psi[MT],sep = "")))
      xlab.label <- expression(paste("The initial strength of CTT (",zeta[0],")",sep = ""))
    } else if(all(dat.all.1$parm == "kappa")){
      # xlab.label <- xlab(expression(paste(" RH on RR: ", psi[RH],sep = "")))
      xlab.label <- expression(paste("Relative transmission rate of IH compared with IU (",
                                     kappa,")",sep = ""))
    } else if(all(dat.all.1$parm == "sigma")){
      xlab.label <- expression(paste("Latency period (",sigma,")",sep = ""))
    } else if(all(dat.all.1$parm == "gamma")){
      xlab.label <- expression(paste("Infectious period (",gamma,")",sep = ""))
    } else if(all(dat.all.1$parm == "R0IN")){
      xlab.label <- expression(paste("Basic reprodduction number (",R[0],")",sep = ""))
    } else {
      xlab.label <- expression(paste("Symptomatic rate of COVID-19 (",omega,")",sep = ""))
    }
    
    plot.total[[i]] <- ggplot(data=dat.plot,aes(x = x,y = y)) + geom_line(linetype = 1,alpha = 0.5) + 
      theme_classic(base_family = "serif",base_size = 10) + 
      geom_pointrange(data=dat.all.1,aes(x = target,y = loglik,ymin = LL,ymax = UL),colour = cols[3],
                      fill="yellow",size=0.2,linetype = 1) +
      geom_hline(yintercept=lik.95,colour="black",linetype = 2,alpha = 0.5) +
      geom_vline(xintercept = c(x.value.1,x.value.2,max.x),colour=cols[2],linetype = c(4,4,1)) +
      #labs(title = title_text) + 
      xlab(xlab.label) + ylab("Profile log likehood")
  }
  fig_mer <- plot_grid(plotlist = plot.total,
                       labels = LETTERS[1:length(parm.name)],ncol = 2)
  
  ggsave(paste(profile.figure.dir,"/Profile.pdf",sep = ""),
         dpi = 300,width = 18,height = 16,units = "cm")
  return(list(fig = fig_mer,result = max_dat))
}

fig_profile <- Profile.Plot(dataset = ms_res,p.value = 0.95,profile.figure.dir = "Figures")

est_qe <- 1/(exp(-(fig_profile$result[3,2]*c(7:34)+log(fig_profile$result[5,2])))+1)
est_wh <- 1/(exp(-(fig_profile$result[4,2]*c(7:34)+log(fig_profile$result[5,2])))+1)

1/(exp(-(fig_profile$result[4,2]*c(7:42)+log(fig_profile$result[5,2])))+1)
as.Date("2020-01-09") + 42

res_plot <- data.frame(qwh = est_wh,qeh = est_qe,date = as.Date("2020-01-16")+(0:length(7:33)))
res_plot <- melt(res_plot,id = c("date"))
res_plot$variable <- factor(res_plot$variable,levels = c("qwh","qeh"),
                            labels = c("Wuhan","Hubei excluding Wuhan"))

(fig_hubei <- ggplot(data=res_plot,aes(x = date, y = value,color = factor(variable)))+
  geom_point() + geom_line() + 
  theme_classic(base_size = 12,base_family = "serif") +
  scale_color_manual(values = cols[c(1,3)],guide = guide_legend(title = "")) +
  xlab("") + ylab("Strength of infection isolation") +
  theme(legend.position = c(0.2,0.9)))

ggsave("Figures/hubei_contact_tracing.pdf",fig_hubei,width = 14,height = 12,units = 'cm',dpi = 300)


#############################################################################################
#########################  Parameter estimation for specific province #######################
#############################################################################################


prov_id <- c("guangdong","henan","zhejiang","anhui","chongqing","jiangsu","hunan",
             "jiangxi","shandong","sichuan","heilongjiang","shanghai","beijing")

pro_dir <- "E:/nCOV/Server/Provinces/Profile_1"

file_dir <- list.files(path = pro_dir,pattern = "csv",full.names = T)

res_list <- lapply(file_dir,function(id){
  
  dat.tem <- read.csv(id,stringsAsFactors = F)
  dat.tem %>% 
    subset(is.finite(loglik)) -> dat.tem
  dat.tem %>%
    plyr::ddply(~quatar,subset,loglik == max(loglik)) -> dat.pred
  dat.pred <- filter(dat.pred,loglik > -2000)
  region_id <- gsub(pro_dir,"",id);region_id <- gsub(".csv","",region_id);
  region_id <- gsub("/Hope_Qua_NF_","",region_id);
  
  para.tem <- mcap(lp =dat.pred$loglik,parameter = dat.pred$quatar)
  # plot(dat.pred$loglik~dat.pred$quatar)
  res_tem <- data.frame(low = para.tem$ci[1],
                        high = para.tem$ci[2],
                        median = para.tem$mle,
                        region = region_id)
  return(res_tem)
})

res_fin_pro <- bind_rows(res_list)
res_fin_pro <- filter(res_fin_pro,region %in% prov_id)

prov_id[!prov_id %in% res_fin_pro$region]

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

res_fin_pro$region <- firstup(res_fin_pro$region)
res_fin_pro <- arrange(res_fin_pro,median)
res_fin_pro$region <- factor(res_fin_pro$region,levels = c(res_fin_pro$region))

(fig_pro_1 <- ggplot(data = res_fin_pro,aes(x = region, y = median)) +
  geom_pointrange(aes(ymin = low,ymax = high)) +
  theme_classic(base_family = "serif",base_size = 12) +
  xlab("Provinces / municipalities")+
 # ylab(expression(paste("Growth rate of infection isolation (",zeta[s],")",sep = ""))) +
  ylab(expression(paste("Growth coefficient of infection isolation","",sep = ""))) +
  scale_y_continuous(limits = c(0.3,0.9))+
  coord_flip())

ggsave(filename = "Figures/zeta.pdf",fig_pro_1,width = 14,height = 14,units = "cm",dpi = 300 )


q_res_list <- lapply(as.character(unique(res_fin_pro$region)),function(id){
  dat.tem <- filter(res_fin_pro,region == id)
  # res_tem_list <- lapply(1:10,function(data){
  #   1/(exp(-(dat.tem[data,"quatar"]*c(7:36)+log(dat.tem[data,"qualow"])))+1)
  # })
  # res_tem <- bind_cols(res_tem_list)
  # est_q <- apply(res_tem,1,mean)
  est_q <- 1/(exp(-(dat.tem$median*c(7:34)+log(0.0003)))+1)
  data.frame(q = est_q,region = id,date = as.Date("2020-01-16")+(0:length(7:33)))
})
q_res <- bind_rows(q_res_list)

# Sys.setlocale("LC_TIME", "English") # Windows


filter(q_res,date == as.Date("2020-02-12"))

q_res$region <- factor(q_res$region,levels = levels(res_fin_pro$region))

fig_pro_2 <- ggplot(data =  q_res,aes(x = date, y = region)) +
  geom_tile(aes(fill = q)) + 
  theme_bw(base_family = "serif",base_size = 12) +
  scale_fill_gradient2(low = "white", high = "indianred1", na.value = "transparent", 
                       guide = guide_colorbar(title.position = "right", 
                                              title = expression(paste("Strength of infection isolation","",sep = "")),
                                              title.theme = element_text(angle = -90), 
                                              title.hjust = 0.5, 
                                              barwidth = 1, 
                                              barheight = 22)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.margin = margin(0,0,0,-0.1,unit = "cm"))+
  scale_x_date(expand = c(0, 0))+xlab("") + ylab("Provinces / municipalities")

fig_pro <- cowplot::plot_grid(fig_pro_1,fig_pro_2,labels = c("A","B"))

ggsave(filename = "Figures/zeta_heatmap.pdf",fig_pro_2,width = 14,height = 14,units = "cm",dpi = 300 )
ggsave(filename = "Figures/zeta_pro.pdf",fig_pro,width = 26,height = 14,units = "cm",dpi = 300 )


##############################################################################################
######################################## Cities   ############################################
##############################################################################################


city_dir <- "E:/nCOV/Server/Cities/Results/Profile_1"

file_dir <- list.files(path = city_dir,pattern = "csv",full.names = T)

res_list <- lapply(file_dir,function(id){
  dat.tem <- read.csv(id,stringsAsFactors = F)
  dat.tem %>% 
    subset(is.finite(loglik)) -> dat.tem
  dat.tem %>%
    plyr::ddply(~quatar,subset,loglik == max(loglik)) -> dat.pred
  dat.pred <- filter(dat.pred, loglik > -2000)
  region_id <- gsub(city_dir,"",id);region_id <- gsub(".csv","",region_id);
  region_id <- gsub("/Hope_Qua_3_","",region_id);
  
  para.tem <- mcap(lp =dat.pred$loglik,parameter = dat.pred$quatar)
  
  res_tem <- data.frame(low = para.tem$ci[1],
                        high = para.tem$ci[2],
                        median = para.tem$mle,
                        region = region_id)
  return(res_tem)
})

res_fin_city <- bind_rows(res_list)
res_fin_city$region <- firstup(res_fin_city$region)
res_fin_city <- arrange(res_fin_city,median)
res_fin_city$region <- factor(res_fin_city$region,levels = c(res_fin_city$region))



fig_city_1 <- ggplot(data = res_fin_city,aes(x = region, y = median)) +
  geom_pointrange(aes(ymin = low,ymax = high)) +
  theme_classic(base_family = "serif",base_size = 12) +
  xlab("Cities") + ylab(expression(zeta[s])) +
  scale_y_continuous(limits = c(0.3,0.9))+
  ylab(expression(paste("Growth coefficient of infection isolation","",sep = ""))) +
  coord_flip()



city_res_list <- lapply(unique(res_fin_city$region),function(id){
  
  dat.tem <- filter(res_fin_city,region == id)
  # res_tem_list <- lapply(1:10,function(data){
  #   1/(exp(-(dat.tem[data,"quatar"]*c(7:36)+log(dat.tem[data,"qualow"])))+1)
  # })
  # res_tem <- bind_cols(res_tem_list)
  # est_q <- apply(res_tem,1,mean)
  est_q <- 1/(exp(-(dat.tem$median*c(7:34)+log(0.0003)))+1)
  data.frame(q = est_q,region = id,date = as.Date("2020-01-16")+(0:length(7:33)))
})
city_res <- bind_rows(city_res_list)

Sys.setlocale("LC_TIME", "English") # Windows

fig_city_2 <- ggplot(data =  city_res,aes(x = date, y = region)) +
  geom_tile(aes(fill = q)) + 
  theme_bw(base_family = "serif",base_size = 12) +
  scale_fill_gradient2(low = "white", high = "indianred1", na.value = "transparent", 
                       guide = guide_colorbar(title.position = "right", 
                                              title = expression(paste("Strength of infection isolation","",sep = "")),
                                              title.theme = element_text(angle = -90), 
                                              title.hjust = 0.5, 
                                              barwidth = 1, 
                                              barheight = 11)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.margin = margin(0,0,0,-0.1,unit = "cm"))+
  scale_x_date(expand = c(0, 0))+xlab("") + ylab("Cities")



fig_city <- cowplot::plot_grid(fig_city_1,fig_city_2,labels = c("A","B"))

ggsave(filename = "Figures/zeta_city_heatmap.pdf",fig_city_2,width = 14,height = 14,units = "cm",dpi = 300 )
ggsave(filename = "Figures/zeta_city.pdf",fig_city,width = 24,height = 14,units = "cm",dpi = 300 )


fig_city <- cowplot::plot_grid(fig_city_1,fig_city_2,labels = c("B",""))
fig_pro <- cowplot::plot_grid(fig_pro_1,fig_pro_2,labels = c("A",""))


fig_merge <- cowplot::plot_grid(fig_pro,fig_city,ncol = 1,rel_heights = c(0.65,0.35))
ggsave(filename = "Figures/zeta_merge.pdf",fig_merge,width = 24,height = 20,units = "cm",dpi = 300 )


res_fin <- rbind(res_fin_pro,res_fin_city)

write.csv(res_fin,file = "Eestimated_Q.csv",row.names = F)

