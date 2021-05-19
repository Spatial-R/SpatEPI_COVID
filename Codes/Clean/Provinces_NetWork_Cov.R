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

setwd("Server")
select_city <- F

date_target <- seq.Date(as.Date("2020-01-01"),as.Date("2020-04-01"),by = "day")

flow_dat <- read_csv("../Data/Flow_20200214/20200101-20200208_out.csv");flow_dat$flow <- NULL
flow_dat_1 <- read.csv("../Data/Flow_20200214/20200209_out.csv",stringsAsFactors = F); 
flow_dat_2 <- read.csv("../Data/Flow_20200214/20200210_out.csv",stringsAsFactors = F)
flow_dat_3 <- read_csv("../Data/Flow_20200214/20200211_out.csv")
flow_dat_4 <- read_csv("../Data/Flow_20200214/20200212_out.csv")

flow_dat_5 <- rbind(flow_dat,flow_dat_1,flow_dat_2,flow_dat_3,flow_dat_4)
flow_dat_5 <- data.frame(flow_dat_5)

flow_dat_5 <- mutate(flow_dat_5,flow = qianxi_index*qianxi_value/100,   ###*9.735272727272728
                     date = as.Date(as.character(date),"%Y%m%d"))


######################################################################################################
#########################################  Merge the cities ##########################################
######################################################################################################

pop_dat <- read.csv("../Data/Demo/Citycode.csv",header = T,stringsAsFactors = F)

for (i in c(1:nrow(pop_dat))){
  pop_dat[i,"city1"] <- ifelse(pop_dat[i,6] == "hubei",pop_dat[i,4],pop_dat[i,6])
}

pop_dat <- filter(pop_dat,!city1 %in% c("taiwan","hongkong","macao"))
final_city_name <- intersect(unique(c(flow_dat_5$from_city,flow_dat_5$to_city)),pop_dat$city)
pop_dat_1 <- filter(pop_dat,city %in% final_city_name);which(duplicated(pop_dat_1$city))

flow_city_1 <- filter(flow_dat_5,from_city %in% final_city_name & to_city %in% final_city_name)

flow_city_1 <- mutate(flow_city_1,type = ifelse(date > as.Date("2020-01-24"),1,0))

flow_city_2 <- data.frame(summarise(group_by(flow_city_1,from_city,to_city,type),flow = mean(flow)))

flow_city_3 <- merge(flow_city_2,pop_dat_1[,c(2,10)],by.x = "from_city",by.y = "city",all.x = T)
flow_city_4 <- merge(flow_city_3,pop_dat_1[,c(2,10)],by.x = "to_city",by.y = "city",all.x = T)
flow_city_5 <- flow_city_4[,-c(1,2)];names(flow_city_5) <- c("type","flow","from_city","to_city")
flow_city_5 <- data.frame(summarise(group_by(flow_city_5,type,from_city,to_city),flow = sum(flow,na.rm = T)))
#flow_city_5 <- filter(flow_city_5,!(from_city %in% c("macao","hongkong") | to_city  %in% c ("macao","hongkong")))

flow_city_5_0 <- filter(flow_city_5,type == 0)[,-1];flow_city_5_1 <- filter(flow_city_5,type == 1)[,-1]

flow_city_5_0 <- spread(flow_city_5_0,to_city,flow);flow_city_5_1 <- spread(flow_city_5_1,to_city,flow)

identical(sort(flow_city_5_0$from_city),sort(names(flow_city_5_0)[-1]))
identical(sort(flow_city_5_1$from_city),sort(names(flow_city_5_1)[-1]))

combn_name_0 <- intersect(flow_city_5_0$from_city,names(flow_city_5_0))
combn_name_1 <- intersect(flow_city_5_1$from_city,names(flow_city_5_1))
combn_fin <- intersect(combn_name_0,combn_name_1)

col_keep_0 <- which(names(flow_city_5_0) %in% combn_fin)
row_keep_0 <- which(flow_city_5_0$from_city %in% combn_fin)
flow_city_5_0_1 <- flow_city_5_0[row_keep_0,c(1,col_keep_0)]
flow_city_5_0_1[,-1] <- apply(flow_city_5_0_1[,-1], 2,function(data)ifelse(is.na(data),0,data))
identical(flow_city_5_0_1$from_city,names(flow_city_5_0_1)[-1]) ### check for the row and col name

col_keep_1 <- which(names(flow_city_5_1) %in% combn_fin)
row_keep_1 <- which(flow_city_5_1$from_city %in% combn_fin)
flow_city_5_1_1 <- flow_city_5_1[row_keep_1,c(1,col_keep_1)]
flow_city_5_1_1[,-1] <- apply(flow_city_5_1_1[,-1], 2,function(data)ifelse(is.na(data),0,data))
identical(flow_city_5_1_1$from_city,names(flow_city_5_1_1)[-1]) ### check for the row and col name

identical(flow_city_5_1_1$from_city,flow_city_5_0_1$from_city) ### check two dataset


flow_city_5_0_1[,-1] <- apply(flow_city_5_0_1[,-1],2,floor)
flow_city_5_1_1[,-1] <- apply(flow_city_5_1_1[,-1],2,floor)

diag(flow_city_5_0_1[,-1]) <- 0;diag(flow_city_5_1_1[,-1]) <- 0;

#flow_city_5_0_1 <- data.frame(flow_city_5_0_1);flow_city_5_1_1 <- data.frame(flow_city_5_1_1)

###############################################################################################
####################################### Population ############################################
###############################################################################################

pop_dat_1 <- filter(pop_dat_1,city1 %in% flow_city_5_1_1$from_city)
pop_dat_2 <- data.frame(summarise(group_by(pop_dat_1,city1),pop = sum(pop)))

dat_pop_list <- lapply(unique(date_target),function(data){
  pop_dat_2$date <- data
  return(pop_dat_2)
})
pop_dat_3 <- bind_rows(dat_pop_list)


contact <- read.csv("../Data/Contact/all_city.csv",header = T,stringsAsFactors = F)
contact_1 <- filter(contact,city_code %in% pop_dat_1$city_code)
contact_1 <- merge(contact_1,pop_dat_1[,c(1,10)],by = "city_code",all.x = T)

contact_2 <- mutate(contact_1,date = as.Date(as.character(date),"%Y%m%d"),
                    type = ifelse(date > as.Date("2020-01-24"),1,0))
contact_2 <- filter(contact_2,date > as.Date("2020-1-08"))
contact_2 <- data.frame(summarise(group_by(contact_2,date,city1),interflow = mean(interflow)))
contact_3 <- spread(contact_2,city1,interflow)

full_date <- data.frame(date = seq.Date(range(contact_3$date)[1],range(contact_3$date)[2],by="day"))
contact_3_0 <- merge(contact_3,full_date,by = "date",all.x = T)

contact_3_0[,-1] <- apply(contact_3_0[,-1],2,function(data)zoo::na.approx(data))


contact_p_list <- lapply(2:ncol(contact_3_0),function(data){
  dat_tem <- contact_3_0[,data]/max(contact_3_0[,data])
})

contact_p <- data.frame(bind_cols(contact_p_list));names(contact_p) <- names(contact_3_0)[-1]
contact_4_0 <- contact_3_0;contact_4_0[,-1] <- contact_p

identical(names(flow_city_5_1_1)[-1],names(contact_3_0)[-1]) ### check for the row and col name

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

case_code <- data.frame(read_csv("../Data/Demo/Normalized_citycode.csv"))

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
case_fin_2 <- merge(case_fin_1,pop_dat_1[,c(1,10)],by = "city_code",all.x = T)[,c(3,5,4)]
case_fin_2 <- filter(case_fin_2,!is.na(city1))
#write.csv(case_fin_2,file = "Case_all.csv",row.names = F)

case_fin_2 <- data.frame(summarise(group_by(case_fin_2,date,city1),cases = sum(cases)))
unique(case_fin_2$city1)

case_full <- data.frame(expand.grid(date = seq.Date(range(case_fin_2$date)[1],range(case_fin_2$date)[2],by ="day"),
                                    city1 = unique(pop_dat_1$city1)))
case_fin_3 <- merge(case_full,case_fin_2,by = c("date","city1"),all.x = T)
case_fin_3$cases <- ifelse(is.na(case_fin_3$cases),0,case_fin_3$cases)
case_fin_3 <- data.frame(summarise(group_by(case_fin_3,date,city1),cases = sum(cases)))
case_fin_3 <- arrange(case_fin_3,date,city1)
case_fin_4 <- spread(case_fin_3,city1,cases)

identical(names(flow_city_5_1_1)[-1],sort(names(case_fin_4)[-1])) ### check for the row and col name


########################################  Model data  constructue ##############################

contact_all <- contact_3_0;
contact_all[,-1] <- apply(contact_all[,-1],2,round,2)
contact_all <- mutate(contact_all,day = date - as.Date("2020-01-08"))
contact_all <- filter(contact_all,day > 0)[,-c(1,ncol(contact_all))]

contact_all_p <- contact_all;
contact_all_p_list <- lapply(1:ncol(contact_all),function(data){
  round(contact_all[,data]/(max(contact_all[,data])),2)
})
contact_all_p <- data.frame(bind_cols(contact_all_p_list))
names(contact_all_p) <- names(contact_all)


pop_all <- pop_dat_3;
pop_all <- mutate(pop_all,day = as.numeric(date - as.Date("2020-01-09")))
pop_all <- filter(pop_all,day > 6)[,-c(3)];names(pop_all)[1] <- c("city")
pop_all <- arrange(pop_all[,c(3,1,2)],city,day)

flow_all <- rbind(flow_city_5_0_1[,-1],flow_city_5_1_1[,-1])

case_all <- case_fin_4
case_all <- mutate(case_all,day = as.numeric(date - as.Date("2020-01-09")))
case_all <- filter(case_all,day > 6)[,-c(1)]
case_all <- melt(case_all,id = "day"); names(case_all)[c(2:3)] <- c("city","cases")

case_all <- mutate(case_all,cases=ifelse(cases < 0,0,cases),city = as.character(city))

pop_dat_2 <- arrange(pop_dat_2,city1)

if(identical(names(flow_all),pop_dat_2$city1)){
  print("population dataset is right")
  pop_target <- pop_dat_1$pop
}


#################################################################################################
###################################### Matrix ###################################################
#################################################################################################


to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
wuhan_pos <- which(names(flow_all) == "wuhan") 

test_name  <-  sort(pop_dat_1[pop_dat_1$Name=="hubei","name"])
test_col <- unlist(lapply(test_name,function(data) which(names(flow_all) == data)))

flow_all <- apply(flow_all,2,function(data)ifelse(is.na(data),0,data))
flow_all <- data.frame(flow_all)

flow_C_rows <- apply(flow_all,1,to_C_array)
flow_C_array <- to_C_array(flow_C_rows)


contact_kl <- contact_all_p[-c(1:as.numeric(as.Date("2020-01-24") - as.Date("2020-01-08"))),]
contact_kl <- round(apply(contact_kl, 2, mean),3)

contact_kl_array <- to_C_array(contact_kl)

U <- dim(flow_all)[2];dt <- 1/2

v_by_g_C <- Csnippet(paste0("const double v_by_g[",nrow(flow_all),"][",U,"] = ",
                            flow_C_array,"; \n ",  
                            "const int targetp = ",wuhan_pos,"; \n ",
                            #"const double contactp[",nrow(contact_final),"][",U,"] = ",
                            #contact_CP_array,"; \n ",
                            # "const double casep[",length(case_p),"] = ",contact_cp_array, "; \n ",
                            "const double contreduce[",length(contact_kl),"] = ",contact_kl_array, "; "
                            #  "const double contact[",nrow(contact_all),"][",U,"] = ",
                            #  contact_C_array,";"
))

measles_globals <- Csnippet(
  paste0("const int U = ",U,"; \n ", v_by_g_C)
)

test_name 

save(measles_globals,U,flow_all,contact_all_p,case_all,dt,
     pop_dat_1,test_col,test_name,pop_all,
     file = paste("Data/Cities_Spatial_Prov",".RData",sep = ""))
