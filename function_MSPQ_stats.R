MSPQ_stats <- function(out){ #### function arguments
    
if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
library(pacman)
p_load(kmlShape, lmerTest, lme4)
    
df_a <- out[["numeric_dataset"]]

df_a$Leaf.Temperature.Differential %<>% abs

SoV <- c("date", "time", out$Sources_ov_Variation)

arr_SoV <- SoV[-c(grep("date", SoV),grep("time", SoV),grep("Gen|gen", SoV))]

means <- df_a %>%              #Outputs 
  as_tibble %>%                #Outputs 
  group_by(.dots = SoV) %>%    #Outputs 
  summarise_all("mean") %>%    #Outputs 
  ungroup                      #Outputs 
                               #Outputs 
st_d <- df_a %>%               #Outputs 
  as_tibble %>%                #Outputs 
  group_by(.dots = SoV) %>%    #Outputs 
  summarise_all("sd") %>%      #Outputs 
  ungroup                      #Outputs 


# ANOVAs ------------------------------------------------------------------
resp <- c("Leaf.Temperature.Differential", "Relative.Chlorophyll", "phi_index","Phi2", "NPQt", "LEF")

formulas <- list()
# for(i in 1:length(resp)){
#   formulas[[i]] <- as.formula(paste0(resp[i],"~ 1+",paste0("as.factor(",SoV[3:length(SoV)],collapse = ")*") ,")*Light.Intensity..PAR.*Ambient.Humidity+as.factor(date)+as.factor(time)"))
#   names(formulas)[[i]] <- resp[i]
# }

for(i in 1:length(resp)){
if(any(grepl("Treat|Trat|treat|trat", SoV))){
  
  formulas[[i]] <- as.formula(paste0(resp[i],"~", paste0("as.factor(", c(arr_SoV),")", collapse = "+"), 
       "+as.factor(time)+(1|",SoV[grep("Gen|gen", SoV)], ")",
       "+(1|",paste0(SoV[grep("Gen|gen", SoV)]),":",SoV[grep("Treat|Trat|treat|trat", SoV)],")",
       "+(", paste0("date|", SoV[grep("Gen|gen", SoV)]), ")"))
  
} else if (!is.empty(arr_SoV)){
  
  formulas[[i]] <- as.formula(paste0(resp[i],"~", paste0("as.factor(", c(arr_SoV),")", collapse = "+"), 
         "+as.factor(time)+(1|",SoV[grep("Gen|gen", SoV)],")",
         "+(", paste0("date|", SoV[grep("Gen|gen", SoV)]), ")"))
  
} else {
  formulas[[i]] <- as.formula(paste0(resp[i],"~",  
                                     "as.factor(time)+(1|",SoV[grep("Gen|gen", SoV)],")",
                                     "+(", paste0("date|", SoV[grep("Gen|gen", SoV)]), ")"))
  }
}
rm(i)


ANOVAs <- list() #Outputs

for(i in 1:length(resp)){
  ANOVAs[[i]] <- lmer(formulas[[i]], data = means, control = lmerControl(optimizer ="Nelder_Mead")) %>% summary
  names(ANOVAs)[[i]] <- resp[i]
}

rm(formulas, resp, i)

# Ranks -------------------------------------------------------------------

dates_rank <- lapply(as.character(unique(means$date)), function(x)filter(means, date == x))

dates_rank <- lapply(1:length(dates_rank), function(x){
  as.data.frame(dates_rank[[x]][,-which(names(dates_rank[[x]])=="time")])
  })

for(i in 1:length(dates_rank)){
dates_rank[[i]] %<>% 
  group_by(.dots = SoV[3:length(SoV)]) %>% 
  summarise_all("mean") %>% 
  ungroup %>% 
  as.data.frame
}
rm(i)

x <- c("Leaf.Temperature.Differential", "Relative.Chlorophyll", "phi_index","Phi2", "LEF", "Thickness")
ranks <- list()

if(is.empty(arr_SoV)){
  
  for(i in 1:length(dates_rank)){
    per_date <- list()
    
    for(j in 1:length(x)){
      per_date[[j]] <- dates_rank[[i]][order(-dates_rank[[i]][,x[j]]),] %>% 
        droplevels %>% 
        as_tibble %>% 
        dplyr::select(contains(c(SoV[grep("Gen|gen", SoV)]))) %>% 
        dplyr::mutate(x = 1:dplyr::n()) %>% 
        arrange_(.dots = SoV[grep("Gen|gen", SoV)]) %>% 
        as.data.frame
      names(per_date)[[j]] <- x[j]
    }
  
  
  scores <- paste0(x,"_score")
  ranks[[i]] <- join_all(per_date, by= SoV[3:length(SoV)], type='left')
  names(ranks[[i]])[grep("x", names(ranks[[i]]))] <- scores
  ranks[[i]]$final_score <- rowSums(ranks[[i]][,grep("_score", names(ranks[[i]]))])
  ranks[[i]]$date <- dates_rank[[i]]$date[1]
  }
  
} else{
  
for(i in 1:length(dates_rank)){
  per_date <- list()
  
  for(j in 1:length(x)){
    per_date[[j]] <- dates_rank[[i]][order(dates_rank[[i]][,c(arr_SoV)], -dates_rank[[i]][,x[j]]),] %>% 
      droplevels %>% 
      as_tibble %>% 
      dplyr::select(contains(c(SoV[grep("Gen|gen", SoV)], arr_SoV))) %>% 
      group_by(.dots = arr_SoV) %>%
      dplyr::mutate(x = 1:dplyr::n()) %>% 
      arrange_(.dots = SoV[grep("Gen|gen", SoV)]) %>% 
      as.data.frame
      names(per_date)[[j]] <- x[j]
  }

  
  scores <- paste0(x,"_score")
  ranks[[i]] <- join_all(per_date, by= SoV[3:length(SoV)], type='left')
  names(ranks[[i]])[grep("x", names(ranks[[i]]))] <- scores
  ranks[[i]]$final_score <- rowSums(ranks[[i]][,grep("_score", names(ranks[[i]]))])
  ranks[[i]]$date <- dates_rank[[i]]$date[1]
}  
  }
rm(per_date)


# Rank Selection ----------------------------------------------------------

gen_le <- means[,grep("Gen|gen", SoV)] %>% as.data.frame() %>% .[,1] %>% unique %>% length

gen_le_5 <- round(gen_le*0.05)

if(gen_le_5==0){
  gen_le_5 <- 1
}

if(is.empty(arr_SoV)){
  
  selection <- list()
  
  for(s in 1:length(ranks)){
    good <- ranks[[s]] %>%
      arrange(final_score) %>% 
      .[,grep("Gen|gen", names(ranks[[s]]))] %>% head(.,n=gen_le_5) %>% as.character
    
    bad <- ranks[[s]] %>%
      arrange(final_score) %>% 
      .[,grep("Gen|gen", names(ranks[[s]]))] %>% tail(.,n=gen_le_5) %>% as.character
    
    selection[[s]] <- data.frame(date = ranks[[s]]$date[1:gen_le_5],
                                 select_good = good,
                                 select_bad = bad)
    rm(good,bad)
  }
} else if(length(arr_SoV)==1){
  if(length(unique(ranks[[1]][,arr_SoV]))==1){
    
    selection <- list()
    
    for(s in 1:length(ranks)){
    good <- ranks[[s]] %>%
    arrange(final_score) %>% 
    .[,grep("Gen|gen", names(ranks[[s]]))] %>% head(.,n=gen_le_5) %>% as.character
    
    bad <- ranks[[s]] %>%
      arrange(final_score) %>% 
      .[,grep("Gen|gen", names(ranks[[s]]))] %>% tail(.,n=gen_le_5) %>% as.character
    
    selection[[s]] <- data.frame(date = ranks[[s]]$date[1:gen_le_5], SoV = arr_SoV, 
                                 SoV_ID = unique(ranks[[s]][,arr_SoV]),
                                 select_good = good,
                                 select_bad = bad)
    rm(good,bad)
    }
  } else {
    
    selection <- list()
    
    
    for(m in 1:length(ranks)){
      a <- list()
    sel <- lapply(unique(ranks[[m]][,arr_SoV]), function(x) subset(ranks[[m]], ranks[[m]][,arr_SoV]==x))
    
    for(i in 1:length(sel)){
  
        good <- sel[[i]] %>%
        arrange(final_score) %>% 
        .[,grep("Gen|gen", names(sel[[i]]))] %>% head(.,n=gen_le_5) %>% as.character
      
      bad <- sel[[i]] %>%
        arrange(final_score) %>% 
        .[,grep("Gen|gen", names(sel[[i]]))] %>% tail(.,n=gen_le_5) %>% as.character
      
      a[[i]] <- data.frame(select_good = good,
                              select_bad = bad,
                              date = unique(sel[[i]]$date),
                              SoV = arr_SoV,
                              SoV_ID = unique(sel[[i]][,arr_SoV]))
      rm(good,bad)
    }
    selection[[m]] <- do.call(rbind, a)
    }
    rm(a,sel)
  }
} else{

  #### mas de una fuente de variacion diferente a genotipo??  
}


date_len <- length(dates_rank)

ranks <- do.call(rbind, ranks) #Outputs

ranks$date %<>% ymd

selection <- do.call(rbind, selection) #Outputs

sel_good <- selection$select_good %>% table %>% as.data.frame() %>% arrange(desc(Freq))
sel_bad <- selection$select_bad %>% table %>% as.data.frame() %>% arrange(desc(Freq))

names(sel_good)[which(names(sel_good)==".")] <- SoV[grep("Gen|gen", SoV)]
names(sel_bad)[which(names(sel_bad)==".")] <- SoV[grep("Gen|gen", SoV)]

sel_bad %<>% 
  mutate(Freq_rel = Freq/date_len)

sel_good %<>% 
  mutate(Freq_rel = Freq/date_len)

ranking <- list(number_of_samplings = date_len, 
                scores = ranks,
                selection_per_date = selection,
                Freq_good = sel_good,
                Freq_bad = sel_bad)      ###### Outputs

output <- list(means = means,
               std = st_d,
               ANOVA = ANOVAs,
               ranking = ranking)

return(output)
}

# cor_means <- cor(means[-c(1:4)])
# corrplot(cor_means, method = "number", type = "upper", tl.cex = 0.5, number.cex = 0.5)
# 
# times <- means$date %>% unique %>% length
#   
# AM <- means %>% filter(time == "AM") %>% arrange(Genotipo)
# 
# ######IDENTIFICAR Y COMPLETAR n PARTIENDO POR GENOTIPOS Y VALUE NA
# AM %>% 
#   group_by(Genotipo) %>% 
#   summarise(n = n()) %>% arrange(n) %>% View
# 
# 
# 
# 
# 
# AM_clds <- AM[,c("Genotipo", "Leaf.Temperature.Differential")] %>%
#   mutate(Genotipo = as.factor(Genotipo)) %>% as.data.frame
# 
# 
# gens <- lapply(unique(AM_clds$Genotipo), function(x)filter(AM_clds,Genotipo == x))
# 
# for(i in 1:length(gens)){
# gens[[i]]$times <- 0:(nrow(gens[[i]])-1)/10
# gens[[i]] %<>% dplyr::select(Genotipo, times, Leaf.Temperature.Differential)
# }
# 
# AM_clds <- do.call(rbind,gens)
# rm(gens)
# 
# set.seed(1)
# 
# nored <- AM_clds %>% as.data.frame %>% 
#   cldsLong %>% kmlShape(., toPlot = "none")
# 
# ###plotting K-means
# nored@trajMeans %>% 
#   ggplot(.,aes(x = times, y = traj, colour = factor(iCenters)))+
#   geom_line(size = 1)+
#   ylim(c(-10,2))
# 
# 
# 
# 
# df_a %>% names
# 
# mod <- df_a %>% 
#   as_tibble %>% 
#   dplyr::select(date,time, Genotipo,R,G,B,Phi2,PhiNO,PhiNPQ,phi_index, Thickness, LEF, Relative.Chlorophyll,
#                 Leaf.Temperature.Differential) %>% 
#   mutate(Genotipo = as.factor(Genotipo)) %>%
#   #group_by(date,time,Genotipo,Tratamiento) %>% 
#   #summarise_all("mean") %>% 
#   #ungroup %>% 
#   mutate(date = as.factor(date)) %>% 
#   aov(Phi2 ~ time + Genotipo + Error(date/Genotipo), data = .)
# 
# summary(mod)
# 
# dates <- list()
# 
# 
# fechas <- df_a$date %>% unique %>% as.character %>% sort
# for(i in 1:length(fechas)){
#   dates[[i]] <- filter(df_a, date == fechas[i])
#   names(dates)[[i]] <- fechas[i]
# }
# rm(i)
# 
# dates[[1]] %>% str
# 
# ?cldsLong
# 
# 
# g <- function(x)dnorm(0:100,runif(1,25,75),10)*rnorm(1,5,1)
# dn <- data.frame(id=rep(1:200,each=101),
#                  times=rep((0:100)/10,times=20),
#                  traj=as.numeric(sapply(1:200,g)))
# 
