
MSPQ_tidy <- function(df, data_name){

if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
library(pacman)
p_load(plotly, robustbase, TestDataImputation, dplyr, ggplot2, imputeTS, installr, corrplot, stringr, magrittr, ez, lubridate, tidyverse)

  
  #-----delete Plot.1 column (repeated)-----
  r_full <- nrow(df)
  c_full <- ncol(df)
  delete <- which(names(df) == "Plot.1")
  if(!is.empty(delete)){df <- df[,-delete]}
  rm(delete)
  #-----Replacing null for NA-----
  cat("Replacing null for NA\n")
  
  for(i in 1:ncol(df)){
    if(any(df[,i]=="null", na.rm = T)){
      df[,i] <- as.numeric(gsub(x = df[,i], pattern = "null", replacement = NA))
    }
  }
  rm(i)
  #-----------------Calculating Phi index--------------------------------
  cat("Calculating Phi Index\n")
  
  df %<>% 
    mutate(phi_index = Phi2/(PhiNPQ+PhiNO))
  
  
  #-----Formating dates and creating time (morning/afternoon) column-----
  cat("Formating dates and creating time (morning/afternoon) column\n")
  
  df$time %<>%  as.character
  
  if(grepl(x = df$time[1], pattern = ":.*AM|PM&AM|PM")){
    
    x <- list()
    
    for(i in 1:length(df$time)){
      
      x[[i]] <- unlist(strsplit(df$time[i], " "))[-2]
      x[[i]] <- paste0(x[[i]][1]," ", x[[i]][2])
    }
    
    x <- unlist(x)
    df$time <- x
    rm(x,i)
    
    names(df)[names(df) == "time"] <- "date"
    
    df <- separate(df,date,into = c("date", "time"), sep = " ")
    
    df$date <- as.Date(df$date, format = "%m/%d/%Y")
    df$time <- as.factor(df$time)
  } else {
    stop("there is not AM/PM indicator in column time and/or hour is missing, check out first\n")
  }
  
  #-----Discarding rows with issues-----
  
  cat("Discarding rows with issues\n")
  
  issues_index <- which(!is.na(df$Issues))
  r_issue <- length(issues_index)
  if(r_issue==0){issues_index <- c()}
  
  #-----Subsetting factors and character columns out from df (NOT to be included into final df)-----
  
  cat("Subsetting factors and character columns out from df (NOT to be included into final df)\n") 
  factor_index <-c() 
  
  for(i in 1:ncol(df)){
    if(is.factor(df[,i])){
      factor_index[i] <- i
    }else if(is.character(df[,i])){
      factor_index[i] <- i
    }else if(length(which(is.na(df[,i]))) == nrow(df) || length(unique(df[,i])) == 1 ||length(unique(df[,i])) == 2){
      factor_index[i] <- i  
    }else if(length(unique(df[,i])) == 3 || length(which(is.na(df[,i]))) >= round(nrow(df)*0.5,digits = 0)){
      factor_index[i] <- i
    } else if (is.Date(df[,i])){
      factor_index[i] <- i
    }
    
  }
  rm(i)
  factor_index <- factor_index[!is.na(factor_index)]
  factor_index <- c(factor_index,which(names(df)=="ID"))
  #-----Applying control structures-----
  cat("Applying control structures\n")
  
  SoV_index <- which(names(df)=="date"): (which(names(df)=="absorbance_420")-1)
  factor_index <- factor_index[!factor_index %in% SoV_index]
  SoV_names <- names(df)[SoV_index[3:length(SoV_index)]]
  
  
  num_df <- df[,-factor_index]
  factor_df <- df[,factor_index]
  c_factor <- length(factor_index)
  rm(factor_index)
  
  #----finding absorbances with NA's (NOT to be included into final df)----
  
  absorbance_index <- which(grepl(x = names(num_df), pattern = "absorbance_"))
  absorbances <- unlist(strsplit(names(num_df)[absorbance_index], "absorbance_"))
  absorbances <- absorbances[-which(absorbances=="")]
  
  absorbance<- c()
  for(i in 1:length(absorbance_index)){
    if(length(which(is.na(num_df[,paste0("absorbance_",absorbances[i])])))!=0){
      absorbance[i] <- absorbances[i]
    }
  }
  rm(i)
  
  absorbance <-absorbance[!is.na(absorbance)]
  
  absorbance_index <- c() 
  
  for(i in 1:length(absorbance)){
    absorbance[i] <- paste0("absorbance_",absorbance[i])
    absorbance_index[i] <- which(names(num_df)==absorbance[i])
  }
  rm(i);rm(absorbance);rm(absorbances)
  c_absorbance <- length(absorbance_index)
  
  if(is.empty(absorbance_index)){
    absorbance_index <- NA
  }
  
  #----Comparing Rel_clo with SPAD_605, if identical (NOT to be included into final df)----
  
  if(identical(df$Relative.Chlorophyll,df$SPAD_650)){
    clo_index <- as.numeric(which(names(num_df)=="SPAD_650"))
    clo_rem <- 1
  } else {
    clo_index <- NA
    clo_rem <- 0
  }
  
  #----if NA's in spatial coordinates, removing the columns (NOT to be included into final df)----
  if(any(names(num_df) == "Latitude")){
    if(any(is.na(num_df$Latitude))){
      location_index <- c(which(names(num_df)=="Latitude"), which(names(num_df)=="Longitude"))
      c_lon_lat_a <- 2
      c_lon_lat <- 2
    } else {
      location_index <- NULL
      c_lon_lat <- 0
      c_lon_lat_a <- 0
    }
  } else {
    location_index <- NULL
    c_lon_lat <- 2
    c_lon_lat_a <- 0
  }
  
  f_indices <- unique(c(absorbance_index, clo_index,location_index))
  f_indices <- f_indices[!is.na(f_indices)]   
  
  if(!is.empty(f_indices)){
    factor_df <- cbind(factor_df,num_df[,f_indices])
    num_df <- num_df[,-f_indices]
  }
  
  rm(absorbance_index,clo_index, f_indices)
  
  #-----Finding and removing rows with NA's-----
  cat("Finding and removing rows with NA's\n")
  row_NA <- list()
  
  
  for(i in 1:ncol(num_df)){
    row_NA[[i]] <- num_df[,i] %>% 
      is.na() %>% 
      which() 
  }
  rm(i)
  
  row_NA <- unique(unlist(row_NA))
  r_rem <- length(row_NA)
  #-----Discarding rows with PAR Phi2, PhiNPQ, PhiNO and Rel_clo out of range-----
  cat("Discarding rows with PAR Phi2, PhiNPQ, PhiNO and Rel_clo out of range\n")
  
  par_out <- which(num_df$Light.Intensity..PAR. > 2500 | num_df$Light.Intensity..PAR. < 1)  
  phi2_out <- which(num_df$Phi2 > 0.85 | num_df$Phi2 < 0.05)
  phinpq_out <- which(num_df$PhiNPQ > 0.85 | num_df$PhiNPQ < 0)
  phino_out <- which(num_df$PhiNO > 0.5 | num_df$PhiNO < 0)
  
  if(clo_rem==1){
    rel_clo_out <- which(num_df$Relative.Chlorophyll > 75 | num_df$Relative.Chlorophyll < 0)
  } else if(clo_rem==0) {
    rel_clo_out <- which(num_df$SPAD_650 > 75 | num_df$SPAD_650 < 0)
  }
  psac <- which(num_df$PS1.Active.Centers > 50 | num_df$PS1.Active.Centers < -50)
  psopc <- which(num_df$PS1.Open.Centers > 50 | num_df$PS1.Open.Centers < -50)
  psorc <- which(num_df$PS1.Over.Reduced.Centers > 50 | num_df$PS1.Over.Reduced.Centers < -50)
  psoxc <- which(num_df$PS1.Oxidized.Centers > 50 | num_df$PS1.Oxidized.Centers < -50)
  npqt <- which(num_df$NPQt > 30 | num_df$NPQt < -30)
  vh <- which(num_df$vH. > 10 | num_df$vH. < -10)
  ecs_m <- which(num_df$ECSt.mAU > 5 | num_df$ECSt.mAU < -5)
  gh <- which(num_df$gH. > 1200 | num_df$gH. < -1200)
  ecs_t <- which(num_df$ECS_tau > 100 | num_df$ECS_tau < -100)
  
  rows_out <- unique(c(issues_index, row_NA))
  out_param <- unique(c(par_out, phi2_out, phinpq_out, phino_out, rel_clo_out, psac, psopc, psorc, psoxc, npqt, vh, ecs_m, gh, ecs_t))
  final_removals <- unique(c(rows_out, out_param))
  r_out_param <- length(out_param)
  
  #-----Making final df-----
  cat("Making final df\n")
  if(!is.empty(final_removals)){
    num_df <- num_df[-final_removals,]
    factor_df <- factor_df[-final_removals,]
    rm_data <- df[final_removals,]
  } else {
    rm_data <- NA
  }
  

  #-----Making summary table-----
  cat("Making summary table\n")
  summ_tab <- as.data.frame(matrix(NA, 1))
  names(summ_tab) <- "file_name"
  summ_tab$file_name <- data_name ###### index with list.files()
  summ_tab$total_obs <- r_full
  summ_tab$total_variables <- c_full
  summ_tab$non_numeric_variables_removed <- c_factor
  summ_tab$absorbance_variables_removed_by_NA <- c_absorbance
  summ_tab$LAT_LON_removed <- c_lon_lat
  summ_tab$SPAD_650_removed <- clo_rem
  summ_tab$obs_removed_by_Issues <- r_issue
  summ_tab$obs_removed_by_NA <- r_rem
  summ_tab$obs_removed_by_param_out <- r_out_param
  summ_tab$Total_obs_removed <- length(final_removals)
  summ_tab$perc_removed_obs <- round(((length(final_removals))*100/r_full), digits = 2)
  summ_tab$perc_removed_var <- round(((c_lon_lat_a+clo_rem+c_factor+c_absorbance)*100/c_full), digits = 1)
  
  
  #----Dropping empty levels in factors----
  num_df <- droplevels(num_df)
  factor_df <- droplevels(factor_df)
  if(is.data.frame(rm_data)){
  rm_data <- droplevels(rm_data)
  }
  #----output list----
  if(summ_tab$perc_removed_obs!=0){
    rm_table <- rm_data[,SoV_index[-2]] %>% table %>% as.data.frame
    if(any(rm_table$Freq==0)){
      rm_table <- rm_table[-which(rm_table[,"Freq"]==0),]
    }
    rm_table$date <- ymd(rm_table$date)
    rm_table %<>% arrange(Freq)
    
  } else {
    rm_table <- "No observations were removed"
  }
  summ_tab %<>% t
  
  out <- list(num_df, factor_df, summ_tab, rm_data, rm_table, SoV_names)
  names(out) <- c("numeric_dataset", "non_numeric_dataset", "summary", "removed_observations", "removed_freq", "Sources_ov_Variation")
  
  
  
  cat("Done!!!\n")
  return(out)
}

