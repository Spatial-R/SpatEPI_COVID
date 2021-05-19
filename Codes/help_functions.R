geo_latlon <- function(data,add_var="address",api_key = "ukBIiSwlUoEj0PUMKS98PXWzPFX5c4YG"){
  dat.tem <- data
  url.base <- "http://api.map.baidu.com/geocoding/v3/?address="
  url.true <- paste("&output=xml&ak=",api_key,"&callback=showLocation",sep = "")
  dat.list <- lapply(1:nrow(dat.tem), function(id){
    url.new <- paste(url.base,dat.tem[id,add_var],url.true, sep = "")
    url.result <- try(getURL(url.new))
    longlat <- unlist(stri_match_all_regex(url.result,"[0-9]+[.]*[0-9]*[<>]"))[c(2:3)]
    longlat <- data.frame(t(longlat));
    longlat[,1:2] <- apply(longlat[,1:2],2,as.character)
    return(longlat)
  })
  lt <- bind_rows(dat.list)
  lt[,1:2] <- apply(lt[,1:2],2,function(data){
    lt.tem <- as.numeric(gsub("<","",data));
    return(lt.tem)
  })
  names(lt) <- c("lat","lon")
  dat.final <- cbind(dat.tem,lt); return(dat.final)
}

firstup <- function(data){
substr(data,1,1) <- toupper(substr(data,1,1))
data
}


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

