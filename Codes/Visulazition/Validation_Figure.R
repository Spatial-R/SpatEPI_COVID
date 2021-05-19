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




pro_dir <- "E:/nCOV/Server/Final_2/Model_Results"

file_dir <- list.files(path = pro_dir,pattern = "csv",full.names = T)
file_dir <- file_dir[stri_detect_fixed(file_dir,"Hope_UHQ_NF")]

q_res_list <- lapply(file_dir,function(id){
  
  dat.tem <- read.csv(id,stringsAsFactors = F)
  #dat.tem %>% 
  #  subset(is.finite(loglik) & nfail.max == 0,-c(nfail.max,nfail.min)) -> dat.tem
  region_id <- gsub(pro_dir,"",id);region_id <- gsub(").csv","",region_id);
  region_id <- gsub("/Hope_UHQ_NF","",region_id);
  region_id <- substr(region_id,2,nchar(region_id))
  names(dat.tem)[which(names(dat.tem) == region_id)] <- "target"
  
  dat.tem %>%
    plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.pred
  #dat.pred <- dat.tem
  dat.pred <- filter(dat.pred,loglik > -1100)
  para.tem <- mcap(lp =dat.pred$loglik,parameter = dat.pred$target)
  plot(dat.pred$loglik~dat.pred$target)
  
  res_tem <- data.frame(low = para.tem$ci[1],
                        high = para.tem$ci[2],
                        median = para.tem$mle,
                        param = region_id)
})

seq_res <- do.call("rbind",q_res_list)

seq_pre <- data.frame(param = seq_res$param,
                      origin = c(0.28,0.3,0.1,0.2,0.03,0.6,2.3,0.19)) 

seq_all <- merge(seq_res,seq_pre,by = "param")
seq_all$param <- factor(seq_all$param,levels = c("sigma","gamma","R0IN","omega","kappa",
                                                 "qua","quawh","qualow"))

fig_l <- ggplot(data = seq_all,aes(x = param)) + 
  geom_pointrange(aes( y = median, ymin = low, ymax = high)) +
  geom_point(aes(y = origin),color = "red") +
  theme_classic(base_size = 10,base_family = "serif") +
  xlab('') + ylab("Estimation of parameter values") +
  scale_x_discrete(labels = c("sigma" = expression(sigma),
                              "gamma" = expression(gamma),
                              "R0IN" = expression(R[0]^"'"),
                              "omega" = expression(omega),
                              "kappa" = expression(kappa),
                              "qua" = expression(zeta[e]),
                              "quawh" = expression(zeta[w]),
                              "qualow" = expression(zeta[0])
  ))




pro_dir <- "E:/nCOV/Server/Final_2/Model_Results_E"

file_dir <- list.files(path = pro_dir,pattern = "csv",full.names = T)
file_dir <- file_dir[stri_detect_fixed(file_dir,"Hope_UHQ_NF_EQ")]

q_res_list <- lapply(file_dir,function(id){
  
  dat.tem <- read.csv(id,stringsAsFactors = F)
  #dat.tem %>% 
  #  subset(is.finite(loglik) & nfail.max == 0,-c(nfail.max,nfail.min)) -> dat.tem
  region_id <- gsub(pro_dir,"",id);region_id <- gsub(").csv","",region_id);
  region_id <- gsub("/Hope_UHQ_NF_EQ","",region_id);
  region_id <- substr(region_id,2,nchar(region_id))
  names(dat.tem)[which(names(dat.tem) == region_id)] <- "target"
  
  dat.tem %>%
    plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.pred
  #dat.pred <- dat.tem
  dat.pred <- filter(dat.pred,loglik > -1100)
  para.tem <- mcap(lp =dat.pred$loglik,parameter = dat.pred$target)
  plot(dat.pred$loglik~dat.pred$target)
  
  res_tem <- data.frame(low = para.tem$ci[1],
                        high = para.tem$ci[2],
                        median = para.tem$mle,
                        param = region_id)
})

seq_res <- do.call("rbind",q_res_list)

seq_pre <- data.frame(param = seq_res$param,
                      origin = c(0.28,0.3,0.1,0.2,0.03,0.2,2.3,0.19)) 

seq_all <- merge(seq_res,seq_pre,by = "param")
seq_all$param <- factor(seq_all$param,levels = c("sigma","gamma","R0IN","omega","kappa",
                                                 "qua","quawh","qualow"))

fig_e <- ggplot(data = seq_all,aes(x = param)) + 
  geom_pointrange(aes( y = median, ymin = low, ymax = high)) +
  geom_point(aes(y = origin),color = "red") +
  theme_classic(base_size = 10,base_family = "serif") +
  xlab('') + ylab("Estimation of parameter values") +
  scale_x_discrete(labels = c("sigma" = expression(sigma),
                              "gamma" = expression(gamma),
                              "R0IN" = expression(R[0]^"'"),
                               "omega" = expression(omega),
                              "kappa" = expression(kappa),
                              "qua" = expression(zeta[e]),
                              "quawh" = expression(zeta[w]),
                              "qualow" = expression(zeta[0])
                           ))

fig_valid <- cowplot::plot_grid(fig_e,fig_l,labels = c("A","B"),ncol = 1)

ggsave2(fig_valid,filename ="Figures/validation.pdf",width = 13,height = 14,units = "cm")
