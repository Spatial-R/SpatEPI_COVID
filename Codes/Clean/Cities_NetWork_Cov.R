###########################################################################################################
##########################################   NetWork for Cities  ##########################################
###########################################################################################################

library(tidyr)
library(dplyr)
library(stringi)
library(reshape2)
library(pomp)
library(zoo)
library(openxlsx)
library(spatPomp)
library(readr)
Sys.setlocale(category="LC_ALL")
select_city <- F

date_target <- seq.Date(as.Date("2020-01-01"),as.Date("2020-04-01"),by = "day")

flow_dat <- read_csv("Data/Flow_20200214/20200101-20200208_out.csv");flow_dat$flow <- NULL
flow_dat_1 <- read.csv("Data/Flow_20200214/20200209_out.csv",encoding = "GB18030",stringsAsFactors = F)
flow_dat_2 <- read.csv("Data/Flow_20200214/20200210_out.csv",encoding = "GB18030",stringsAsFactors = F)
flow_dat_3 <- read_csv("Data/Flow_20200214/20200211_out.csv")
flow_dat_4 <- read_csv("Data/Flow_20200214/20200212_out.csv")

flow_dat_5 <- rbind(flow_dat,flow_dat_1,flow_dat_2,flow_dat_3,flow_dat_4)
flow_dat_5 <- data.frame(flow_dat_5)

flow_dat_5 <- mutate(flow_dat_5,flow = qianxi_index*qianxi_value/100,
                     date = as.Date(as.character(date),"%Y%m%d"))

######################################################################################################
#########################################  Merge the cities ##########################################
######################################################################################################

pop_dat <- read.csv("Data/Demo/Citycode.csv",header = T,stringsAsFactors = F)
GDP_dat <- read.csv("Data/Demo/GDP.csv",header = T,stringsAsFactors = F)

pop_dat <- merge(pop_dat,GDP_dat[,-c(1,2)],by = "city_code")

final_city_name <- intersect(unique(c(flow_dat_5$from_city,flow_dat_5$to_city)),pop_dat$city)
pop_dat_1 <- filter(pop_dat,city %in% final_city_name);which(duplicated(pop_dat_1$city))

flow_city_1 <- filter(flow_dat_5,from_city %in% final_city_name & to_city %in% final_city_name)

flow_city_1 <- mutate(flow_city_1,type = ifelse(date > as.Date("2020-01-24"),1,0))

flow_city_2 <- data.frame(summarise(group_by(flow_city_1,from_city,to_city,type),flow = mean(flow)))

flow_city_3 <- merge(flow_city_2,pop_dat_1[,c(2,4)],by.x = "from_city",by.y = "city",all.x = T)
flow_city_4 <- merge(flow_city_3,pop_dat_1[,c(2,4)],by.x = "to_city",by.y = "city",all.x = T)
flow_city_5 <- flow_city_4[,-c(1,2)];names(flow_city_5) <- c("type","flow","from_city","to_city")

flow_city_5_0 <- filter(flow_city_5,type == 0)[,-1];flow_city_5_1 <- filter(flow_city_5,type == 1)[,-1]

flow_city_5_0 <- spread(flow_city_5_0,to_city,flow);flow_city_5_1 <- spread(flow_city_5_1,to_city,flow)

identical(sort(flow_city_5_0$from_city),sort(names(flow_city_5_0)))
identical(sort(flow_city_5_1$from_city),sort(names(flow_city_5_1)))

combn_name_0 <- intersect(flow_city_5_0$from_city,names(flow_city_5_0))
combn_name_1 <- intersect(flow_city_5_1$from_city,names(flow_city_5_1))
combn_fin <- intersect(combn_name_0,combn_name_1)

col_keep_0 <- which(names(flow_city_5_0) %in% combn_fin)
row_keep_0 <- which(flow_city_5_0$from_city %in% combn_fin)
flow_city_5_0_1 <- flow_city_5_0[row_keep_0,c(1,col_keep_0)]


flow_manuc <- function(data,na_fill = "0",gravity_full = T,gdp_adj = T){
  
  dat_tem <- data

  if (na_fill == "min"){
    for ( i in c(1:nrow(dat_tem))){
      dat_tem[i,which(is.na(dat_tem[i,]))] <- min(dat_tem[i,-1],na.rm = T)
      dat_tem[i,(i+1)] <- NA
    }
  } else {
   #### we used the gravity model to complete the whole flow
    
    if(isTRUE(gravity_full)){
      ana_tem <- reshape2::melt(dat_tem,id = "from_city")
      ana_tem <- filter(ana_tem,!(from_city == variable))
      if(isTRUE(gdp_adj)){
        ana_tem_1 <- merge(ana_tem,pop_dat_1[,c(4,7:10)],by.x = "from_city",by.y = "name")
        names(ana_tem_1) <- c("origin","target","flow","pop_origin","lat_origin","lon_origin","gdp_origin")
        ana_tem_2 <- merge(ana_tem_1,pop_dat_1[,c(4,7:10)],by.x = "target",by.y = "name")
        names(ana_tem_2) <- c("target","origin","flow","pop_origin","lat_origin","lon_origin","gdp_origin",
                              "pop_target","lat_target","lon_target","gdp_target")
        
        min_dat <- min(ana_tem_2[ana_tem_2$flow > 0,"flow"]) 
        
        ana_tem_3 <- mutate(ana_tem_2,
                            distance = log(sqrt((lat_target - lat_origin)^2 + (lon_target - lon_origin)^2)),
                            pop_target = log(as.numeric(pop_target)),
                            pop_origin = log(as.numeric(pop_origin)),
                            gdp_target = log(as.numeric(gdp_target)),
                            gdp_origin = log(as.numeric(gdp_origin)))
        ana_tem_fu <- filter(ana_tem_3,!is.na(flow))
        ana_tem_na <- filter(ana_tem_3,is.na(flow))
        ana_tem_fu <- mutate(ana_tem_fu,flow  = ifelse(flow == 0, log(min_dat),log(flow)))
        gravity_model <- lm(flow~pop_target+distance+pop_origin+gdp_target+gdp_origin,data = ana_tem_fu)
        
      } else {
        ana_tem_1 <- merge(ana_tem,pop_dat_1[,c(4,7:9)],by.x = "from_city",by.y = "name")
        names(ana_tem_1) <- c("origin","target","flow","pop_origin","lat_origin","lon_origin")
        ana_tem_2 <- merge(ana_tem_1,pop_dat_1[,c(4,7:9)],by.x = "target",by.y = "name")
        names(ana_tem_2) <- c("target","origin","flow","pop_origin","lat_origin","lon_origin",
                              "pop_target","lat_target","lon_target")
        
        min_dat <- min(ana_tem_2[ana_tem_2$flow > 0,"flow"]) 
        
        ana_tem_3 <- mutate(ana_tem_2,
                            distance = log(sqrt((lat_target - lat_origin)^2 + (lon_target - lon_origin)^2)),
                            pop_target = log(as.numeric(pop_target)),
                            pop_origin = log(as.numeric(pop_origin)))
        ana_tem_fu <- filter(ana_tem_3,!is.na(flow))
        ana_tem_na <- filter(ana_tem_3,is.na(flow))
        ana_tem_fu <- mutate(ana_tem_fu,flow  = ifelse(flow == 0, log(min_dat),log(flow)))
        gravity_model <- lm(flow~pop_target+distance+pop_origin,data = ana_tem_fu)
      }

      ana_tem_na$flow <- predict(gravity_model,ana_tem_na)
      ana_tem_fin <- rbind(ana_tem_fu,ana_tem_na)[,c(1:3)]; ana_tem_fin$flow <- exp(ana_tem_fin$flow)
      ana_tem_fin <- spread(ana_tem_fin,target,flow)
      
      if(identical(ana_tem_fin$origin,names(ana_tem_fin)[-1])){
         ### check for the row and col name
         print("The flow data is complete in the gravity model")
      }
      
      if(identical(ana_tem_fin$origin,names(dat_tem)[-1])){
        ### check for the row and col name
        print("The flow data is completely matched with the input data matrix")
      }
      
      check_na <- unlist(lapply(1:nrow(ana_tem_fin),function(id){
        which(is.na(ana_tem_fin[id,-1]))
      }))
      
      if(identical(1:nrow(ana_tem_fin),check_na)){
        ### check for the row and col name
        print("The predicted flow data in the full gravity model is completely full")
      }
      
      for (i in c(1:nrow(dat_tem))){
       # print(as.character(dat_tem[i,1]))
        dat_tem_1 <- (dat_tem[i,-1]); na_pos <- which(is.na(dat_tem[i,-1]))
        pred_dat <- as.numeric(ana_tem_fin[i,(na_pos+1)]);
        min_value <- min(as.numeric(dat_tem_1),na.rm = T)
        pred_dat <- ifelse(pred_dat > min_value,min_value,pred_dat)
        dat_tem[i,(1+na_pos)] <- pred_dat
    }
    } else {
      
    for (i in c(1:nrow(dat_tem))){
      print(as.character(dat_tem[i,1]))
      dat_tem_1 <- data.frame(t(dat_tem[i,-1])); dat_tem_1$city <- row.names(dat_tem_1)
      names(dat_tem_1)[1] <- "X2"
      if(isTRUE(gdp_adj)){
        dat_tem_2 <- merge(dat_tem_1,pop_dat_1[,c(4,7:10)],by.x = "city",by.y = "name",all.x = T)
        dat_target <- filter(pop_dat_1,name == dat_tem[i,1])
        dat_tem_2 <- mutate(dat_tem_2,
                            distance = sqrt((dat_target$lat - lat)^2 + (dat_target$lon - lon)^2),
                            X2 = log(X2), distance = log(distance),pop = log(pop),gdp = log(GDP))
        dat_tem_ana <- filter(dat_tem_2,!is.na(X2))
        dat_tem_na <- filter(dat_tem_2,is.na(X2))
        gravity_model <- lm(X2~pop+distance+gdp,data = dat_tem_ana)
        dat_tem_na$X2 <- predict(gravity_model,newdata = dat_tem_na)
      } else {
        dat_tem_2 <- merge(dat_tem_1,pop_dat_1[,c(4,7:9)],by.x = "city",by.y = "name",all.x = T)
        dat_target <- filter(pop_dat_1,name == dat_tem[i,1])
        dat_tem_2 <- mutate(dat_tem_2,
                            distance = sqrt((dat_target$lat - lat)^2 + (dat_target$lon - lon)^2),
                            X2 = log(X2), distance = log(distance),pop = log(pop))
        dat_tem_ana <- filter(dat_tem_2,!is.na(X2))
        dat_tem_na <- filter(dat_tem_2,is.na(X2))
        gravity_model <- lm(X2~pop+distance,data = dat_tem_ana) 
        dat_tem_na$X2 <- predict(gravity_model,newdata = dat_tem_na)
      }
      dat_tem_na$X2 <- predict(gravity_model,newdata = dat_tem_na)
      dat_tem_na$X2 <- ifelse(dat_tem_na$X2 > min(dat_tem_ana$X2),min(dat_tem_ana$X2),dat_tem_na$X2)
      dat_final <- rbind(dat_tem_ana,dat_tem_na)
      dat_final <- arrange(dat_final,city)
      if(identical(dat_final$city,dat_tem_2$city)){
        dat_tem[i,-1] <- exp(dat_final$X2)
        dat_tem[i,(i+1)] <- NA
      }
    }
    }
  }
  
  fin_check <- unlist(lapply(1:nrow(dat_tem),function(id){
    which(is.na(dat_tem[id,-1]))
  }))
  
  if(identical(1:nrow(dat_tem),fin_check)){
    ### check for the row and col name
    print("Predication is done and the ouput file is correlated")
  }
  
  return(dat_tem)
}

flow_city_5_0_1 <- flow_manuc(data = flow_city_5_0_1,na_fill = "gravity",gravity_full = T)
flow_city_5_0_1[,-1] <- apply(flow_city_5_0_1[,-1], 2,function(data)ifelse(is.na(data),0,data))
identical(flow_city_5_0_1$from_city,names(flow_city_5_0_1)[-1]) ### check for the row and col name

col_keep_1 <- which(names(flow_city_5_1) %in% combn_fin)
row_keep_1 <- which(flow_city_5_1$from_city %in% combn_fin)
flow_city_5_1_1 <- flow_city_5_1[row_keep_1,c(1,col_keep_1)]
flow_city_5_1_1 <- flow_manuc(data = flow_city_5_1_1,na_fill = "gravity",gravity_full = T)
flow_city_5_1_1[,-1] <- apply(flow_city_5_1_1[,-1], 2,function(data)ifelse(is.na(data),0,data))

identical(flow_city_5_1_1$from_city,names(flow_city_5_1_1)[-1]) ### check for the row and col name
identical(flow_city_5_1_1$from_city,flow_city_5_0_1$from_city) ### check two dataset



flow_city_5_0_1[,-1] <- apply(flow_city_5_0_1[,-1],2,round,2)
flow_city_5_1_1[,-1] <- apply(flow_city_5_1_1[,-1],2,round,2)


###############################################################################################
####################################### Population ############################################
###############################################################################################

pop_dat_1 <- filter(pop_dat_1,name %in% flow_city_5_1_1$from_city)
dat_pop_list <- lapply(unique(date_target),function(data){
  pop_dat_1$date <- data
  pop_pos <- match(c("name","pop","date"),names(pop_dat_1))
  return(pop_dat_1[,c(pop_pos)])
})
pop_dat_2 <- bind_rows(dat_pop_list)


contact <- read.csv("Data/Contact/all_city.csv",header = T,stringsAsFactors = F)
contact_1 <- filter(contact,city_code %in% pop_dat_1$city_code)
contact_1 <- merge(contact_1,pop_dat_1[,c(1,4)],by = "city_code",all.x = T)

contact_2 <- mutate(contact_1,date = as.Date(as.character(date),"%Y%m%d"),
                     type = ifelse(date > as.Date("2020-01-24"),1,0))
contact_3 <- spread(contact_2[,c(2,5,6)],name,interflow)

full_date <- data.frame(date = seq.Date(range(contact_3$date)[1],range(contact_3$date)[2],by="day"))
contact_3_0 <- merge(contact_3,full_date,by = "date",all.x = T)

contact_3_0[,-1] <- apply(contact_3_0[,-1],2,function(data)na.approx(data))


contact_p_list <- lapply(2:ncol(contact_3_0),function(data){
  dat_tem <- contact_3_0[,data]/max(contact_3_0[,data])
})

contact_p <- data.frame(bind_cols(contact_p_list));names(contact_p) <- names(contact_3_0)[-1]
contact_4_0 <- contact_3_0;
contact_4_0[,-1] <- contact_p

# cot_tem <- data.frame()
# for (i in c(2:(ncol(contact_4_0) -1 ))){
#   print(i)
#   for (j in c(i:ncol(contact_4_0))){
#     cor_test <- cor.test(contact_4_0[,i],contact_4_0[,j])
#     cor_tem <- data.frame(cor = cor_test$estimate,i =i,j =j)
#     cot_tem <- rbind(cot_tem,cor_tem)
#   }
# }


identical(names(flow_city_5_1_1)[-1],names(contact_3_0)[-1]) ### check for the row and col name


# contact_2 <- data.frame(summarise(group_by(contact_2,name,type),flow = mean(interflow)))
# contact_2_0 <- filter(contact_2, type == 0)[,-2];contact_2_1 <- filter(contact_2, type == 1)[,-2]


pro_name_list <- list.files("E:/nCOV/Data/Cases/",full.names = T)


case_list <- lapply(pro_name_list,function(data){
  print(data)
  dat_tem <- tryCatch(read_csv(paste(data,"/Cases.csv",sep = "")))
  if(!any(stri_detect_fixed(dat_tem$X1,"月"))){
    dat_tem <- tryCatch(read.csv(paste(data,"/Cases.csv",sep = "")))
  }
  dat_tem <- data.frame(dat_tem);names(dat_tem)[1] <- "date"
  dat_tem$date <- gsub("月","-",dat_tem$date);dat_tem$date <- gsub("日","",dat_tem$date)
  dat_tem$date <- gsub("/","-",dat_tem$date)
  dat_tem$date <- as.Date(paste("2020-",dat_tem$date,sep = ""))
  dat_tem_1 <- arrange(dat_tem,date)
  if(ncol(dat_tem) == 2){
    zero_pos <- which.min(dat_tem_1[,2] > 0) -1
  } else {
  zero_pos <- as.numeric(apply(dat_tem_1[,-1],2,function(data)which.min(data > 0)) -1)
  }
  for (i in c(2:ncol(dat_tem))) {
    dat_tem[zero_pos[i-1],i] <- 0
  }
  
  case_list <- lapply(2:ncol(dat_tem_1),function(data){
    c(0,diff(dat_tem_1[,data]))
  })
  
  dat_cases_1 <- bind_cols(case_list); dat_tem_1[,-1] <- dat_cases_1
  if(ncol(dat_tem) == 2){
    dat_tem_1[,2] <- ifelse(is.na(dat_tem_1[,2]),0,dat_tem_1[,2])
  } else {
    dat_tem_1[,-1] <- apply(dat_tem_1[,-1],2,function(data) ifelse(is.na(data),0,data)) 
  }
  dat_tem_2 <- melt(dat_tem_1,id = "date")
  names(dat_tem_2) <- c("date","city","cases")
  return(dat_tem_2)
  })

case_fin <- bind_rows(case_list)

case_code <- data.frame(read_csv("Data/Demo/Normalized_citycode.csv"))

for (i in  (1:nrow(case_code))){ 
  case_code[i,3] <- ifelse(is.na(case_code[i,3]),case_code[i,2],case_code[i,3])
}

case_fin <- mutate(case_fin,city = ifelse(city == "永城市","商丘市",
                                          ifelse(city == "宿松","安庆",
                                                  ifelse(city == "公主岭市","长春市",
                                                  ifelse(city == "邓州市","南阳市",
                                                  ifelse(city == "滑县","安阳市",
                                                  ifelse(city == "长垣市","新乡市",
                                                  ifelse(city == "梅河口市","通化市",
                                                  ifelse(city == "宁东","银川市",
                                                  ifelse(city == "杨凌示范区","咸阳市",
                                                  ifelse(city == "赣江新区","南昌市",
                                                  ifelse(city == "韩城市","渭南市",
                                                  ifelse(city == "济源示范区","新乡市",city)))))))))))))

case_fin <- data.frame(summarise(group_by(case_fin,date,city),cases = sum(cases)))
  
case_fin_1 <- merge(case_fin,case_code[,-2],by.x = "city",by.y = "cases_city",all.x = T)
case_fin_2 <- merge(case_fin_1,pop_dat_1[,c(1,4)],by = "city_code",all.x = T)[,c(3,5,4)]
ko <- case_fin_2[which(is.na(case_fin_2$name)),]; sum(ko$cases)
write.csv(case_fin_2,file = "Case_all.csv",row.names = F)

unique(case_fin_2$name)

case_full <- data.frame(expand.grid(date = seq.Date(range(case_fin_2$date)[1],range(case_fin_2$date)[2],by ="day"),
                         name = sort(pop_dat_1$name)))
case_fin_3 <- merge(case_full,case_fin_2,by = c("date","name"),all.x = T)
case_fin_3$cases <- ifelse(is.na(case_fin_3$cases),0,case_fin_3$cases)
case_fin_3 <- arrange(case_fin_3,date,name)
case_fin_4 <- spread(case_fin_3,name,cases)

identical(names(flow_city_5_1_1)[-1],names(case_fin_4)[-1]) ### check for the row and col name


########################################  Model data  constructue ##############################

contact_all <- contact_3_0;
contact_all[,-1] <- apply(contact_all[,-1],2,round,2)
contact_all <- mutate(contact_all,day = date - as.Date("2020-01-15"))
contact_all <- filter(contact_all,day > 0)[,-c(1,ncol(contact_all))]

contact_all_p <- contact_all;
contact_all_p_list <- lapply(1:ncol(contact_all),function(data){
  round(contact_all[,data]/(max(contact_all[,data])),2)
})
contact_all_p <- data.frame(bind_cols(contact_all_p_list))
names(contact_all_p) <- names(contact_all)



pop_all <- pop_dat_2;
pop_all <- mutate(pop_all,day = as.numeric(date - as.Date("2020-01-09")))
pop_all <- filter(pop_all,day > 6)[,-c(3)];names(pop_all)[1] <- c("city")
pop_all <- arrange(pop_all[,c(3,1,2)],city,day)

flow_all <- rbind(flow_city_5_0_1[,-1],flow_city_5_1_1[,-1])

case_all <- case_fin_4
case_all <- mutate(case_all,day = as.numeric(date - as.Date("2020-01-09")))
case_all <- filter(case_all,day > 6)[,-c(1)]
case_all <- melt(case_all,id = "day"); names(case_all)[c(2:3)] <- c("city","cases")

case_all <- mutate(case_all,cases=ifelse(cases < 0,0,cases),city = as.character(city))

pop_dat_1 <- arrange(pop_dat_1,name)

if(identical(names(flow_all),pop_dat_1$name)){
  print("population dataset is right")
  pop_target <- pop_dat_1$pop
}


#######################################  Select cities ########################################


test_name <- sort(as.character(filter(pop_dat,Name == "hubei")$name))

hubei_col <- which(names(flow_city_5_0_1) %in% test_name)
flow_city_5_0_1_hubei <- filter(flow_city_5_0_1,from_city %in% test_name)[,-c(hubei_col)]
flow_city_5_1_1_hubei <- filter(flow_city_5_1_1,from_city %in% test_name)[,-c(hubei_col)]

pre_dat <- apply(flow_city_5_0_1_hubei[,-1],1,sum);
after_dat <- apply(flow_city_5_1_1_hubei[,-1],1,sum);
names(after_dat) <- names(pre_dat) <- flow_city_5_0_1_hubei$from_city
other_dat <- rbind(pre_dat,after_dat)
save(other_dat,file = "Flow_Other_Cities.RData")

if(isTRUE(select_city)){
  
  target_city <- filter(pop_dat_1,pop > 30591000)$name
  target_city_1 <- union(target_city,test_name)
  test_col <- unlist(lapply(test_name, function(data) which(target_city_1 == data)))
  
  target__col <- unlist(lapply(target_city_1, function(data) which(names(flow_city_5_0_1)[-1] == data)))
  target__col_1 <- unlist(lapply(target_city_1, function(data) which(names(contact_all) == data)))
  identical(target__col,target__col_1)
  
  contact_all <- contact_all[,target__col]
  contact_all_p <- contact_all_p[,target__col]

  case_all <- filter(case_all,city %in% target_city_1)
  flow_all <- flow_all[c(target__col,target__col+ncol(flow_all)),target__col]

  pop_all <- filter(pop_all,city %in% target_city_1)
} else {
  test_col <- unlist(lapply(test_name, function(data) which(names(flow_all) == data)))
}


if(length(test_col) == length(test_name)){
  print("All the dataset was in the flow")
}

#################################################################################################
###################################### Matrix ###################################################
#################################################################################################


to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
wuhan_pos <- which(names(flow_all) == "wuhan")
flow_dat_wuhan <- flow_all$wuhan[1:(nrow(flow_all)/2)]
flow_dat_p <- round(flow_dat_wuhan/(sum(flow_dat_wuhan)),4)
flow_C_wuhan_array <- to_C_array(flow_dat_p)



contact_C_rows <- apply(contact_all,1,to_C_array)
contact_C_array <- to_C_array(contact_C_rows)

contact_CP_rows <- apply(contact_all_p,1,to_C_array)
contact_CP_array <- to_C_array(contact_CP_rows)


flow_C_rows <- apply(flow_all,1,to_C_array)
flow_C_array <- to_C_array(flow_C_rows)

contact_all_per <- contact_all_p[-c(1:as.numeric(as.Date("2020-01-24") - as.Date("2020-01-08"))),]
contact_all_per <- round(apply(contact_all_per, 2, mean),3)
flow_CT_array <- to_C_array(contact_all_per)

U <- dim(flow_all)[2];dt <- 1/2

v_by_g_C <- Csnippet(paste0("const double v_by_g[",nrow(flow_all),"][",U,"] = ",
                            flow_C_array,"; \n ",
                            "const int targetp = ",wuhan_pos,"; \n ",
                           # "const double contactp[",nrow(contact_all_p),"][",U,"] = ",
                           # contact_CP_array,"; \n ",
                            "const double contreduce[",U,"] = ",
                           flow_CT_array,";"
                          #  "const double contact[",nrow(contact_all),"][",U,"] = ",
                          #  contact_C_array,";"
))

measles_globals <- Csnippet(
  paste0("const int U = ",U,"; \n ", v_by_g_C)
)

save(measles_globals,U,flow_all,contact_all,pop_all,case_all,dt,test_col,test_name,
     pop_target,contact_all_p,pop_dat_1,v_by_g_C,
     file = "Process_Data/Cities_Spatial.RData")
