# ## find best model for each country
# setwd("~/Downloads/pfpp-heapstates/pf_summarizer_age_results_sigma05_cauchy001_nonshift")
# ll <- list.files()
# setwd("~/Downloads/pfpp-heapstates/pf_summarizer_age_results_sigma05_cauchy001_shift")
# jj <- list.files()
# k <- c(ll,jj)
over5 <- read.csv('~/Desktop/work/over5-constant.csv', header=T) # reads in over 5 cfr
under5 <- read.csv('~/Desktop/work/under5-constant.csv', header=T) # reads in under 5 cfr

df_u5 <- matrix(NA,nrow=length(isos), ncol=2)
for (ii in 1:length(isos)){
  df_u5[ii,1] <- as.numeric(under5[which(under5$X==isos[ii]),2][1])
}
df_u5 <- as.data.frame(df_u5)
colnames(df_u5) <- c("cfr", "country") 
df_u5$country <- isos

# add in two countries that do not exist in cfr
df_u5[which(df_u5$country=="MNE"),1] <- 0.005

# create df of cfr values for over 5
df_o5 <- matrix(NA,nrow=length(isos), ncol=2)
for (ii in 1:length(isos)){
  df_o5[ii,1] <- as.numeric(over5[which(over5$X==isos[ii]),2][1])
}
df_o5 <- as.data.frame(df_o5)
colnames(df_o5) <- c("cfr", "country")
df_o5$country <- isos
# add in two countries that do not exist in cfr
df_o5[which(df_o5$country=="MNE"),1] <- 0.005/2
o5 <- df_o5 %>% add_row(cfr=0.005/2, country= "KAZ")

# create one cfr df with both under 5 and over 5 values 
cfr <- rbind(df_u5,df_o5)
cfr$u5 <- c(rep(TRUE,dim(cfr)[1]/2), rep(FALSE,dim(cfr)[1]/2))

setwd("~/Desktop/case_outputs_nLL_orig")
files <- list.files()
isos <- str_sub(files, 1, 3)

setwd("~/Downloads/pfpp-heapstates/pf_summarizer_result_fixedV2_all")
k <- list.files()

# normal = 0 / mcv2 = 1 / sia = 2 / mcv2sia = 3
iso <- c()
for (ii in 1:length(k)){
  c <- substr(k[ii],26,28)
  m <- substr(k[ii], 29,41)
  re <- c("_0.5_0", "0.01", "_0", "0_01", ".01", "_","[.]", "5")
  m <- str_remove_all(m, paste(re, collapse = "|"))
  temp <- cbind(c,m)
  iso <- rbind(iso,temp) 
}
iso <- data.frame(iso)
iso <- tibble(iso)

colnames(iso) <- c("iso", "model")
iso <- iso %>% mutate(modelnum = 
                case_when(model == "normal" | model== "normalshift" ~ 0, 
                model == "mcv2"| model == "mcv2shift" ~ 1,
                model == "sia" | model == "siashift" ~ 2,
                model == "mcv2sia"| model == "mcv2siashift" ~ 3))

setwd("~/Downloads/pfpp-heapstates/pf_summarizer_result_fixedV2_all")
df <- c()
for (ii in 1:dim(iso)[1]){
  temp <- k[ii]
  load(temp)
  a <- data.frame(temp, iso$iso[ii], iso$model[ii], iso$modelnum[ii], min(result$par_lik$nll_orig))
  df <- rbind(df,a)
}
df <- tibble(df)
colnames(df) <- c("file", "iso", "model", "modelnum", "nLL")

best_nll <- df %>% group_by(iso) %>% slice(which.min(nLL))

nll_shift <- best_nll %>% filter(model=="siashift" | model=="mcv2siashift" | model=="normalshift" | model==
                     "mcv2shift") %>% mutate(name = paste0(iso,"shift"), sianame = paste0("sia", iso,"shift"))

nll_non <- best_nll %>% filter(model=="sia" | model=="mcv2sia" | model=="normal" | model==
                      "mcv2") %>% mutate(name = iso, sianame =paste0("sia", iso))


best_nll <- full_join(nll_shift, nll_non)

ds <- c()
for (jj in 1:dim(best_nll)[1]){
  tmp_data <- load(file=best_nll$file[jj])
  tmp_data2 <- cbind(result$par_lik, rep(best_nll$iso[jj],100), rep(best_nll$model[jj],100), rep(best_nll$modelnum[jj],100))
  ds <- rbind(ds,tmp_data2)
}
ds <- tibble(ds)
colnames(ds)[7:9] <- c("iso", "model", "modelnum")
ds$number <- rep(1:100,length(unique(iso$iso)))
ds$shift <- ifelse(grepl("shift", ds$model), "shift", NA)

shell_names <- best_nll %>% select(iso, model, modelnum, name, sianame)
shell_names <- shell_names[order(shell_names$iso),]

shell_names <- shell_names %>% mutate(iest_name = 
                                        case_when(model== "normalshift" ~ "normal_shift", 
                                                  model == "mcv2shift" ~ "mcv2_shift",
                                                  model == "siashift" ~ "sia_shift",
                                                  model == "mcv2siashift" ~ "mcv2sia_shift",
                                                  model== "normal" ~ "normal",
                                                  model=="sia" ~ "sia",
                                                  model=="mcv2" ~ "mcv2",
                                                  model=="mcv2sia" ~ "mcv2sia"))
#### CHOOSE COUNTRY TO USE ####
# isoIdx <- unique(iso$iso)[ii]
# gg <- ds %>% filter(iso==isoIdx)
#### LOOP THROUGH ALL COUNTRIES ####
# for (ii in 1:length(unique(iso$iso))){
# for (ii in 21:193){
for (ii in 6:193){
  isoIdx <- unique(shell_names$iso)[ii]
  gg <- ds %>% filter(iso==isoIdx)

for (kk in 1:dim(gg)[1]){
  y <- gg[kk,]
  sig <- 0.100042
  h <- rbind(y$nll_orig, sig, y$beta0, y$beta1, y$pobs, y$pobshigh)
  setwd("~/Downloads/pfpp-heapstates/IEST_files_generated")
  write.table(h, paste0(y$iso,y$model, y$number,".txt"), row.names = F, col.names = F)
}

#######
# change Data_for_fitting2019 to Data_for_fitting2020
n <- c()
for (hh in 1:100){
   u <- paste0("./pfpp 5 ", paste0(shell_names$iso[ii],shell_names$model[ii],hh)," Data_for_fitting2020/", shell_names$name[ii],".txt", " Data_for_fitting2020/", shell_names$sianame[ii], ".txt ",shell_names$modelnum[ii], " 0.5 0.1 IEST_files_generated/",paste0(shell_names$iso[ii],shell_names$model[ii],hh,".txt"))
  n <- rbind(n,u)
}

setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates")
write.table(n, "shell_script_generated.sh", col.names = F, row.names = F, quote=F)

#######
setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates")
system("sh shell_script_generated.sh")

########
age <- rep(seq(1,100), 1000)
# SUM UP CASES FOR ONE COUNTRY
nYears <- 40
nSim <- 100
case_array <- array(NA, dim = c(nYear, 3, nSim))
for (qq in 1:100){
  # file <- paste(gg$iso[qq], gg$model[qq], gg$number[qq], "_i_particles_fxsave.txt", sep="")
  file <- paste(gg$iso[qq], gg$model[qq], gg$number[qq], "_iap_particles_fxsave.txt", sep="")
  best <- read.table(file=file, sep="", header=F, row.names = NULL)
  
  iap <- t(best)
  country_cfr <- cfr %>% filter(country %in% gg$iso[1])
  cfr_matrix <- matrix(rep(rep(country_cfr$cfr,c(5,95)),1000),nrow=100000,ncol=40,byrow=FALSE)
  iap_death <- iap*cfr_matrix
  
  iap_death <- data.frame(iap_death); iap_death <- tibble(iap_death)
  iap_death <- iap_death %>% add_column(age)
  
  by_mean_death <- iap_death %>% group_by(age) %>% summarise_all(list(~mean(.)))
  by_lb_death <- iap_death %>% group_by(age) %>% summarise_all(list(~quantile(.,probs = 0.025)))
  by_ub_death <- iap_death %>% group_by(age) %>% summarise_all(list(~quantile(.,probs = 0.975)))
  
  by_null_death <- by_mean_death %>% mutate(age=NULL); iap_death <- as.data.frame(t(by_null_death));
  iap_death <- iap_death %>% add_column(iso=rep(gg$iso[1],40)); iap_death <- iap_death[,c(101, 1:100)]
  by_null_lb_death <- by_lb_death %>% mutate(age=NULL); iap_lb_death <- as.data.frame(t(by_null_lb_death));
  iap_lb_death <- iap_lb_death %>% add_column(iso=rep(gg$iso[1],40)); iap_lb_death <- iap_lb_death[,c(101, 1:100)]
  by_null_ub_death <- by_ub_death %>% mutate(age=NULL); iap_ub_death <- as.data.frame(t(by_null_ub_death));
  iap_ub_death <- iap_ub_death %>% add_column(iso=rep(gg$iso[1],40)); iap_ub_death <- iap_ub_death[,c(101, 1:100)]
  
  iap_ages_death <- iap_death[,2:101]
  iap_ages_death <- matrix(unlist(t(iap_ages_death)), byrow=TRUE, (40*100), 1) #n countries x 40 yrs x 100 ages
  iap_ages_death <- as.data.frame(iap_ages_death)
  iap_death_a <- rep(seq(1981, 2020), each=100) #each= number age classes
  iap_ages_death <- iap_ages_death %>% add_column(year=rep(iap_death_a,1)) %>% add_column(age=rep(seq(0,99),1*40))
  iap_ages_death <- as.data.frame(iap_ages_death)
  iap_ages_death <- iap_ages_death %>% add_column(country=rep(gg$iso[1],each=4000)) #each= number age classes x 101
  
  iap_ages_death_lb <- iap_lb_death[,2:101]
  iap_ages_death_lb <- matrix(unlist(t(iap_ages_death_lb)), byrow=TRUE, (40*100), 1) #n countries x 40 yrs x 100 ages
  iap_ages_death_lb <- as.data.frame(iap_ages_death_lb)
  iap_death_a_lb <- rep(seq(1981, 2020), each=100) #each= number age classes
  iap_ages_death_lb <- iap_ages_death_lb %>% add_column(year=rep(iap_death_a,1)) %>% add_column(age=rep(seq(0,99),1*40))
  iap_ages_death_lb <- as.data.frame(iap_ages_death_lb)
  iap_ages_death_lb <- iap_ages_death_lb %>% add_column(country=rep(gg$iso[1],each=4000)) #each= number age classes x 101
  
  iap_ages_death_ub <- iap_ub_death[,2:101]
  iap_ages_death_ub <- matrix(unlist(t(iap_ages_death_ub)), byrow=TRUE, (40*100), 1) #n countries x 40 yrs x 100 ages
  iap_ages_death_ub <- as.data.frame(iap_ages_death_ub)
  iap_death_a_ub <- rep(seq(1981, 2020), each=100) #each= number age classes
  iap_ages_death_ub <- iap_ages_death_ub %>% add_column(year=rep(iap_death_a_ub,1)) %>% add_column(age=rep(seq(0,99),1*40))
  iap_ages_death_ub <- as.data.frame(iap_ages_death_ub)
  iap_ages_death_ub <- iap_ages_death_ub %>% add_column(country=rep(gg$iso[1],each=4000)) #each= number age classes x 101
  
  deaths_lb <- iap_ages_death_lb %>% group_by(year) %>% summarise(deaths_lb = sum(V1, na.rm=T))
  deaths_mean <- iap_ages_death %>% group_by(year) %>% summarise(deaths_mean = sum(V1, na.rm=T))
  deaths_ub <- iap_ages_death_ub %>% group_by(year) %>% summarise(deaths_ub = sum(V1, na.rm=T))
  
  
  case_array[,1,qq] <- as.matrix(deaths_lb$deaths_lb)
  case_array[,2,qq] <- as.matrix(deaths_mean$deaths_mean)
  case_array[,3,qq] <- as.matrix(deaths_ub$deaths_ub)
  
  # best <- as.matrix(best)
  # case_array[,1,qq] <- apply(best,1, quantile, probs = 0.025, na.rm = T)
  # case_array[,2,qq] <- apply(best, 1, quantile, probs=0.5, na.rm=T)
  # case_array[,3,qq] <- apply(best,1, quantile, probs = 0.975, na.rm = T)
  # case_array[,4,qq] <- apply(best,1, mean, na.rm = T)
  
}

lb <- apply(case_array[,1,],1, mean, na.rm = T)
# med <- apply(case_array[,2,],1, mean, na.rm = T)
mean <- apply(case_array[,2,],1, mean, na.rm = T)
ub <- apply(case_array[,3,],1, mean, na.rm = T)

quantiles <- cbind(lb, mean, ub)
quantiles <- as.data.frame(quantiles)
case_quantiles <- tibble(quantiles)

setwd("~/Desktop/case_outputs_nLL_orig_deaths")
#write.csv(case_quantiles, paste0(gg$iso[1],"_case_totals.csv"))
write.csv(case_quantiles, paste0(gg$iso[1],"_death_totals.csv"))

#### REMOVE UNNECESSARY FILES FOR NEXT COUNTRY ####
setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates")
remove <- paste0("rm *", paste0(gg$iso[1],gg$model), "*")
system(remove)

setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates/IEST_files_generated")
remove <- paste0("rm *", paste0(gg$iso[1],gg$model), "*")
system(remove)
}

# ################
# # CREATE SCRIPT TO RUN IAP FILES FOR BEST MODELS
# shell_names <- shell_names %>% mutate(iest_name =
#                         case_when(model== "normalshift" ~ "normal_shift",
#                                   model == "mcv2shift" ~ "mcv2_shift",
#                                   model == "siashift" ~ "sia_shift",
#                                   model == "mcv2siashift" ~ "mcv2sia_shift",
#                                   model== "normal" ~ "normal",
#                                   model=="sia" ~ "sia",
#                                   model=="mcv2" ~ "mcv2",
#                                   model=="mcv2sia" ~ "mcv2sia"))
# 
# n <- c()
# for (hh in 1:length(shell_names$iso)){
#   u <- paste0("./pfpp 5 ", paste0(shell_names$iso[hh],shell_names$model[hh])," Data_for_fitting2020/", shell_names$name[hh],".txt", " Data_for_fitting2020/",shell_names$sianame[hh],".txt ",shell_names$modelnum[hh]," 0.5 0.1 burden_2020_all_format_fix/Iest_",paste0(shell_names$iso[hh],shell_names$iest_name[hh],".txt"))
#   n <- rbind(n,u)
# }
# 
# setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates")
# write.table(n, "shell_script_runs.sh", col.names = F, row.names = F, quote=F)
# 
