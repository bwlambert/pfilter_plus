# ## find best model for each country
# setwd("~/Downloads/pfpp-heapstates/pf_summarizer_age_results_sigma05_cauchy001_nonshift")
# ll <- list.files()
# setwd("~/Downloads/pfpp-heapstates/pf_summarizer_age_results_sigma05_cauchy001_shift")
# jj <- list.files()
# k <- c(ll,jj)

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
  a <- data.frame(temp, iso$iso[ii], iso$model[ii], iso$modelnum[ii], min(result$par_lik$nll_new))
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
for (ii in 1:length(unique(iso$iso))){
# for (ii in 21:193){
## ii = 76 error
## start at ii = 130
## for (ii in 76){
  isoIdx <- unique(shell_names$iso)[ii]
  gg <- ds %>% filter(iso==isoIdx)

for (kk in 1:dim(gg)[1]){
  y <- gg[kk,]
  sig <- 0.100042
  h <- rbind(y$nll_new, sig, y$beta0, y$beta1, y$pobs, y$pobshigh)
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
# SUM UP CASES FOR ONE COUNTRY
nYears <- 41
nSim <- 100
case_array <- array(NA, dim = c(nYears, 4, nSim))
for (qq in 1:100){
  file <- paste(gg$iso[qq], gg$model[qq], gg$number[qq], "_i_particles_fxsave.txt", sep="")
  best <- read.table(file=file, sep="", header=F, row.names = NULL)
  
  best <- as.matrix(best)
  case_array[,1,qq] <- apply(best,1, quantile, probs = 0.025, na.rm = T)
  case_array[,2,qq] <- apply(best, 1, quantile, probs=0.5, na.rm=T)
  case_array[,3,qq] <- apply(best,1, quantile, probs = 0.975, na.rm = T)
  case_array[,4,qq] <- apply(best,1, mean, na.rm = T)
  
}

lb <- apply(case_array[,1,],1, mean, na.rm = T)
med <- apply(case_array[,2,],1, mean, na.rm = T)
ub <- apply(case_array[,3,],1, mean, na.rm = T)
mean <- apply(case_array[,4,],1, mean, na.rm = T)

quantiles <- cbind(lb, med, ub, mean)
quantiles <- as.data.frame(quantiles)
case_quantiles <- tibble(quantiles)

setwd("~/Downloads/pfpp-heapstates/case_outputs_nLL_new")
write.csv(case_quantiles, paste0(gg$iso[1],"_case_totals.csv"))

#### REMOVE UNNECESSARY FILES FOR NEXT COUNTRY ####
setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates")
remove <- paste0("rm *", paste0(gg$iso[1],gg$model), "*")
system(remove)

setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates/IEST_files_generated")
remove <- paste0("rm *", paste0(gg$iso[1],gg$model), "*")
system(remove)
}

################
# CREATE SCRIPT TO RUN IAP FILES FOR BEST MODELS
shell_names <- shell_names %>% mutate(iest_name =
                        case_when(model== "normalshift" ~ "normal_shift",
                                  model == "mcv2shift" ~ "mcv2_shift",
                                  model == "siashift" ~ "sia_shift",
                                  model == "mcv2siashift" ~ "mcv2sia_shift",
                                  model== "normal" ~ "normal",
                                  model=="sia" ~ "sia",
                                  model=="mcv2" ~ "mcv2",
                                  model=="mcv2sia" ~ "mcv2sia"))

n <- c()
for (hh in 1:length(shell_names$iso)){
  u <- paste0("./pfpp 5 ", paste0(shell_names$iso[hh],shell_names$model[hh])," Data_for_fitting2020/", shell_names$name[hh],".txt", " Data_for_fitting2020/",shell_names$sianame[hh],".txt ",shell_names$modelnum[hh]," 0.5 0.01 burden_2020_all_format_fix/Iest_",paste0(shell_names$iso[hh],shell_names$iest_name[hh],".txt"))
  n <- rbind(n,u)
}

setwd("/Users/sarahhauryski/Downloads/pfpp-heapstates")
write.table(n, "shell_script_runs.sh", col.names = F, row.names = F, quote=F)

