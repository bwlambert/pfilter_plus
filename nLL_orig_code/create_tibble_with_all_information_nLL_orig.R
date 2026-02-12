########################################
# create_tibble_with_all_information
# this code creates a tibble with all of the information necessary to calculate any metric necessary
# inputs: pf_summarizer_results_sigma05_cauchy001_all_countries (.dat files), over5-constant.csv, under5-constant.csv,
# output(s): 1) all_country_model: tibble containing the columns for years, age, mean_cases, lb_cases, ub_cases, country, model, nLL, nLL_new, under5, cfr, mean_deaths, lb_deaths, ub_deaths and best model (as determined by nLL), 2) best_sums and norm_sums: the sum of deaths using the numbers from the best model and the sum of deaths using the normal model, 4) weighted_country_model: tibble that accounts for model weights
########################################
setwd("~/Downloads/pfpp-heapstates/pf_summarizer_result_fixedV2_all") # load all results from all countries
files <- list.files() # all files
nfiles <- length(files) # number of files
#########################################
### SETS UP DF ALL
# set up first one df for first file
ii <- 1
load(files[ii])
df <- data.frame(result$out_data)
df$country <- str_sub(files[ii], 26, 28)
df$model <- str_sub(files[ii], 29, -14)
df$nll_orig <- mean(result$par_lik$nll_orig)
df$nll_new <- mean(result$par_lik$nll_new)
df$year <- rep(1981:2020,each=100)
rm("result")

# run loop for remaining files
for(ii in 2:nfiles){
  load(files[ii])
  dft <- data.frame(result$out_data)
  dft$country <- str_sub(files[ii], 26, 28)
  dft$model <- str_sub(files[ii], 29, -14)
  dft$nll_orig <- mean(result$par_lik$nll_orig)
  dft$nll_new <- mean(result$par_lik$nll_new)
  dft$year <- rep(1981:2020,each=100)
  df <- rbind(df,dft)
  rm("result")
}

# make tibble for df
df_all <- as_tibble(df)

# make df_all for shift and df_all for nonshift; run once for shift and once for non_shift
# df_all_nonshift <- df_all
# df_all_shift <- df_all
# df_all <- rbind(df_all_shift, df_all_nonshift) # create df will both shift and non shift

# df_all <- as_tibble(df_all)
df_all <- df_all %>% mutate(u5=age<5) # creates new column in df_all to set u5 as t/f

### CFR VALUES 
over5 <- read.csv('~/Desktop/work/over5-constant.csv', header=T) # reads in over 5 cfr
under5 <- read.csv('~/Desktop/work/under5-constant.csv', header=T) # reads in under 5 cfr
iso <- unique(df_all$country) # lists out unique iso in df_all

# create df of cfr values for under 5
df_u5 <- matrix(NA,nrow=length(iso), ncol=2)
for (ii in 1:length(iso)){
  df_u5[ii,1] <- as.numeric(under5[which(under5$X==iso[ii]),2][1])
}

df_u5 <- as.data.frame(df_u5)
colnames(df_u5) <- c("cfr", "country") 
df_u5$country <- iso

# add in two countries that do not exist in cfr
df_u5[which(df_u5$country=="MNE"),1] <- 0.005

# df_u5 <- tbl_df(df_u5)
# df_u5 <- df_u5 %>% add_row(cfr=0.005/2, country= "KAZ")


# create df of cfr values for over 5
df_o5 <- matrix(NA,nrow=length(iso), ncol=2)
for (ii in 1:length(iso)){
  df_o5[ii,1] <- as.numeric(over5[which(over5$X==iso[ii]),2][1])
}

df_o5 <- as.data.frame(df_o5)
colnames(df_o5) <- c("cfr", "country")
df_o5$country <- iso
# add in two countries that do not exist in cfr
df_o5[which(df_o5$country=="MNE"),1] <- 0.005/2

# df_o5 <- tbl_df(df_o5)
# df_o5 <- df_o5 %>% add_row(cfr=0.005/2, country= "KAZ")

# create one cfr df with both under 5 and over 5 values 
cfr <- rbind(df_u5,df_o5)
cfr$u5 <- c(rep(TRUE,dim(cfr)[1]/2), rep(FALSE,dim(cfr)[1]/2))

# JOIN df_all and cfr based on country and age (u5)
df_cfr <- full_join(df_all,cfr) # matches cfr to country and age group

# creates columns for deaths, ub and lb
df_cfr <- df_cfr %>%
    mutate(mean_deaths=(mean_cases*cfr), 
    lb_deaths=(lb_cases*cfr),
    ub_deaths=(ub_cases*cfr))

### finds best model
# df that determines if mcv2 and sia are present in a country
mcv2sia <- read.csv(file="~/Desktop/work/is_mcv2sia.csv")
mcv2sia <- as_tibble(mcv2sia)

# filter out models models that aren't correct for a country where there are no mcv2 or sia
# i.e. if a country has no mcv2, then mcv2/mcv2_shift/mcv2sia/mcv2sia_shift models are not possible
# if there is no sia, then no shift models will work
df_cfr <- df_cfr %>% 
  filter(!(!is.na(match(country,mcv2sia$country[which(!mcv2sia$is.mcv2)])) & (model == "mcv2" | model == "mcv2_shift" | model == "mcv2sia" | model == "mcv2sia_shift"))) %>%
  filter(!(!is.na(match(country,mcv2sia$country[which(!mcv2sia$is.sia)])) & (model == "sia" | model == "sia_shift" | model == "mcv2sia" | model == "mcv2sia_shift" | model == "normal_shift" | model == "mcv2_shift")))

# find the best model based on the lowest nLL
best_mod <- df_cfr %>% group_by(country, model) %>%  # group by country and model
    summarise(nll_orig = min(nll_orig)) %>%     # get the corresponding nll
    group_by(country) %>%                               # then group result by country
    summarize(nll_orig = min(nll_orig))         # and choose the model with the lowest nll

best_mod$best <- "best"                               # label model with "best"

# join with original datset so "best" labels appears on the model with lowest nll_orig
df_cfr_best <- full_join(df_cfr,best_mod,by=c("country","nll_orig")) 

## FIX WHO (NON-MODEL) COUNTRIES
# remove WHO countries and then put them back into the df
who_iso <- read.csv(file="~/Desktop/work/input-ignore.indices.csv")
who_reported <- read.csv(file="~/Desktop/work/input-cases.csv",row.names = 1)

non_who_countries <- df_cfr_best %>%        # subset out all models that are non-who method countries
  filter(is.na(match(country,who_iso$iso))) # filter to all models and non-who method countries

# add in WHO countries
# could just calculate these directly 
reported_cases <- as.vector(t(who_reported[which(!is.na(match(rownames(who_reported),who_iso$iso))),-1]))  # identify countries that use WHO method 2 calculation

# creates a data frame of all the models and cases and combines with non-who method countries
all_country_model <- rbind(non_who_countries, 
            data.frame(year = rep(1981:2020,length(who_iso$iso)), 
                       age = NA,
                       mean_cases = reported_cases/.2,    # assume 20% reporting 
                       lb_cases = reported_cases/.4,     # assume lower bound is 5% reporting
                       ub_cases = reported_cases/.05,      # assume uppper bound is 40% reporting
                       country = rep(who_iso$iso,each=40),
                       model = "who",
                       nll_orig = NA,
                       nll_new = NA,
                       u5 = NA,
                       cfr = NA,
                       mean_deaths = NA,
                       lb_deaths = NA,
                       ub_deaths = NA,
                       best = "who")
)

## FIND SUMS
# filter out only models designated as best and select who countries
# NA's in who countries -- remove NA's
best_models <- all_country_model %>% filter(best=="best" | best=="who")
# sum best model for each country for each year
best_sums <- best_models %>% group_by(year) %>%  # group best model by year
  summarise(sum_deaths = sum(mean_deaths,na.rm=T), 
            sum_lb_deaths = sum(lb_deaths, na.rm = T), 
            sum_ub_deaths = sum(ub_deaths, na.rm = T))

norm_models <- all_country_model %>% filter(model=="normal" | model=="who")
# sum normal model for each country for each year
norm_sums <- norm_models %>% group_by(year) %>% 
  summarise(sum_deaths = sum(mean_deaths, na.rm = T), 
            sum_lb_deaths = sum(lb_deaths, na.rm = T), 
            sum_ub_deaths = sum(ub_deaths, na.rm = T))

# normal model - deaths summed up -- sum across all contries for each year
# best model - deaths summed up -- sum across all countries for each year
# library(RColorBrewer)
# col=brewer.pal(10, "Spectral")

#############
#### weighted model
# logsummexpnorm function
logsumexpnorm <- function(X){
  c <- max(X)
  (exp(X - c)) / sum(exp(X - c))
}

# create df for country 1
isoIdx <- unique(non_who_countries$country)
df_country <- all_country_model %>% filter(country == isoIdx[1]) # subset to country == iso
df_summary <- all_country_model %>% filter(country == isoIdx[1]) %>% # subset to country == iso
  group_by(model) %>%                   # group by model
  summarise(nll = min(nll_orig)) %>%    # get the nll per model
  mutate(nllweight = logsumexpnorm(-nll)) # reweight likelihoods across models to sum to 1
df_nllweight <- inner_join(df_country,df_summary,by="model")        # join aa and model weights

for (ii in 2:length(isoIdx)){
  df_country <- all_country_model %>% filter(country == isoIdx[ii]) # subset to country == iso
  df_summary <- all_country_model %>% filter(country == isoIdx[ii]) %>% # subset to country == iso
    group_by(model) %>%                   # group by model
    summarise(nll = min(nll_orig)) %>%    # get the nll per model
    mutate(nllweight = logsumexpnorm(-nll)) # reweight likelihoods across models to sum to 1
  df_nllweight_new <- inner_join(df_country,df_summary,by="model")        # join df_country and model weights
  
  df_nllweight <- df_nllweight %>% add_row(df_nllweight_new)
}

weighted_country_model <- rbind(df_nllweight, 
                           data.frame(year = rep(1981:2019,length(who_iso$iso)), age = NA,
                                      mean_cases = reported_cases/.2,    # assume 20% reporting 
                                      lb_cases = reported_cases/.05,     # assume lower bound is 5% reporting
                                      ub_cases = reported_cases/.4,      # assume uppper bound is 40% reporting
                                      country = rep(who_iso$iso,each=39),
                                      model = "who",
                                      nll_orig = NA,
                                      nll_new = NA,
                                      u5 = NA,
                                      cfr = NA,
                                      mean_deaths = NA,
                                      lb_deaths = NA,
                                      ub_deaths = NA,
                                      best = "who",
                                      nll=NA, 
                                      nllweight= NA))


