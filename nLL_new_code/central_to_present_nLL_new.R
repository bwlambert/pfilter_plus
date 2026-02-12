library(dplyr)
library(tidyr)
library(tibble)
##### NOTES 8/2/21
## LEAVE DEATHS AND DALYS COLUMNS EMPTY

##### This script creates the central burden estimates files for VIMC/WHO. It's written to loop over
##### several scenarios. You can just comment those (lines 13, 15:17, 188) if you're doing them one at a time. 

##### Make sure you have all of the files needed for calculating deaths, dalys, and cohort size, 
##### and that those are reading in from the right place. (lines 54, 66/85 or 105, 117, 142:143, 146)

##### Also be mindful of which version of the deaths calculation is being used!! 

# scen <- c('scen1', 'scen2', 'scen3', 'scen4', 'scen5', 'scen6', 'scen7', 'scen8', 'scen9', 'scen10')


## create iap files for just best models

setwd("~/Downloads/pfpp-heapstates")
isos <- list.files(pattern="*iap_particles_fxsave.txt")
isos <- substr(isos, 1, 3)

template <- read.csv("~/Downloads/pfpp-heapstates/demographic data/central-burden-template.201910gavi-3.Measles_PSU-Ferrari_standard.csv")
template_countries <- unique(template$country)


# rename files if they are not in the correct format
file.rename(list.files(pattern="*iap_particles_fxsave.txt"), paste0(isos, '_iap_particles_fxsave.txt', sep=""))
isos <- intersect(isos, template_countries)


# for (s in isos) {
  # print(s)
  # age <- rep(seq(1,100), 1000)
  # grid <- data.frame("disease"=NA, "year"=NA, "age"=NA, "country"=NA, "country_name"=NA,
  #                    "cohort_size"=NA, "deaths"=NA, "cases"=NA, "dalys"=NA)
  # 
  # for (i in isos) {
  #   for (s in isos) {
  #     print(s)
  #     age <- rep(seq(1,100), 1000)
  #     grid <- data.frame("disease"=NA, "year"=NA, "age"=NA, "country"=NA, "country_name"=NA,
  #                        "cohort_size"=NA, "deaths"=NA, "cases"=NA, "dalys"=NA)
  #     
      for (i in isos) {
        age <- rep(seq(1,100), 1000)
        grid <- data.frame("disease"=NA, "year"=NA, "age"=NA, "country"=NA, "country_name"=NA,
                           "cohort_size"=NA, "deaths"=NA, "cases"=NA, "dalys"=NA)
        
        setwd("~/Downloads/pfpp-heapstates")
        print(i)
        iap <- read.delim(paste(i,"_iap_particles_fxsave.txt",sep=""), sep="", header=FALSE)
        # iap1 <- iap
        # iap <- iap[,1:100]
        iap <- t(iap)
        iap <- data.frame(iap)
        iap <- tibble(iap)
        iap <- iap %>% add_column(age)
        by_age <- iap %>% group_by(age)
        by_mean <- by_age %>% summarise_all(list(~mean(.)))
        by_null <- by_mean %>% mutate(age=NULL)
        iap <- as.data.frame(t(by_null))
        iap <- iap %>% add_column(iso=rep(i,40)) # add iso for each column... only 40 columns, one for each year
        iap <- iap[,c(101, 1:100)]
        colnames(iap) <- c('iso', 'X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15','X16',
                           'X17','X18','X19','X20','X21','X22','X23','X24','X25','X26','X27','X28','X29','X30','X31', 
                           'X32','X33','X34','X35','X36','X37','X38','X39','X40','X41','X42','X43','X44','X45','X46', 
                           'X47','X48','X49','X50','X51','X52','X53','X54','X55','X56','X57','X58','X59','X60','X61',
                           'X62','X63','X64','X65','X66','X67','X68','X69','X70','X71','X72','X73','X74','X75','X76',
                           'X77','X78','X79','X80','X81','X82','X83','X84','X85','X86','X87','X88','X89','X90','X91',
                           'X92','X93','X94','X95','X96','X97','X98','X99', 'X100')
        
        ## ADD TO TEMPLATE
        ages <- iap[,2:101]
        ages <- matrix(unlist(t(ages)), byrow=TRUE, (40*100), 1) #n countries x 39 yrs x 100 ages
        ages <- as.data.frame(ages)
        a <- rep(seq(1981, 2020), each=100) #each= number age classes
        ages <- ages %>% add_column(year=rep(a,1)) %>% add_column(age=rep(seq(0,99),1*40))
        ages <- as.data.frame(ages)
        ages <- ages %>% add_column(country=rep(i,each=4000)) #each= number age classes x 101
        
        template <- read.csv("~/Downloads/pfpp-heapstates/demographic data/central-burden-template.201910gavi-3.Measles_PSU-Ferrari_standard.csv")
        template <- template %>% filter(country %in% i)
        template <- unique(template$country_name)
        ages$country_name <- template
        ages$disease <- 'Measles'
        colnames(ages)[1] <- 'cases'
        ages$deaths <- NA
        ages$dalys <- NA
        ages$cohort_size <- NA
        ages <- ages[,c('disease', 'year', 'age', 'country', 'country_name', 'cohort_size', 'deaths', 'cases', 'dalys')]
        
        #CFR / DEATHS
        # under5 <- read.csv('../../output docs/input-cfr1-constant.csv')
        #
        # if (i == 'SSD') {
        #   under5 <- under5 %>% filter(X == 'SDN')
        # }
        #
        # if (i != 'SSD') {
        # under5 <- under5 %>% filter(X %in% i)
        # }
        #
        # under5[,40:122] <- under5[,'X2018']
        # under5[,2:21] <- NULL
        # colnames(under5)[2:102] <- seq(2000,2100)
        # under5 <- matrix(unlist(t(under5[,-c(1)])), byrow=TRUE, 101, 1)
        # under5 <- as.data.frame(under5)
        # under5 <- under5 %>% add_column(country=rep(i,each=101)) %>% add_column(year=rep(seq(2000,2100),1))
        # under5 <- under5[rep(seq_len(nrow(under5)), each=5),]
        # under5 <- under5 %>% add_column(age=rep(seq(0,4),101))
        #
        # over5 <- read.csv('../../output docs/input-cfr3-constant.csv')
        #
        # over5 <- over5 %>% filter(X %in% i)
        # over5[,40:122] <- over5[,'X2018']
        # over5[,2:21] <- NULL
        # colnames(over5)[2:102] <- seq(2000,2100)
        # over5 <- matrix(unlist(t(over5[,-c(1)])), byrow=TRUE, 101, 1)
        # over5 <- as.data.frame(over5)
        # over5 <- over5 %>% add_column(country=rep(i,each=101)) %>% add_column(year=rep(seq(2000,2100),1))
        # over5 <- over5[rep(seq_len(nrow(over5)), each=95),] #each= number age classes - 5
        # over5 <- over5 %>% add_column(age=rep(seq(5,99),101))
        #
        # colnames(over5)[1] <- 'cfr'
        # colnames(under5)[1] <- 'cfr'

        # t2 <- left_join(ages, over5, by=c('year','age','country')) %>% left_join(., under5, by=c('year','age','country'))
        # t2$cfr <- rowSums(t2[,c("cfr.x", "cfr.y")], na.rm=TRUE)
        # t2[,c('cfr.x','cfr.y')] <- NULL
        # t2$deaths <- t2$cases * t2$cfr
        #
        # cfr <- read.csv('~/Downloads/pfpp-heapstates/demographic data/cfr_bestcase_campaign_central_Measles-PSU-Ferrari.csv')
        # iso <- filter(cfr, country==i)
        # # iso <- iso[-c(6)]
        # iso <- filter(iso, year <= 2020)
        # y80 <- filter(iso, year==2000)
        # y80 <- y80[rep(seq_len(nrow(y80)),20),]
        # y80$year <- rep(seq(1980,1999), each=100)
        # cfr <- rbind(iso,y80)
        # cfr <- cfr[c(2:4,12)]
        #
        # # t2 <- left_join(cfr,ages,  by=c('year', 'age', 'country))
        # # t2 <- full_join(ages, iso)
        # t2 <- left_join(cfr,ages,  by=c('year', 'age', 'country'))
        # t2$deaths <- t2$cases * t2$cfr
        #
        #
        # ## DALYS
        # le <- read.csv("~/Downloads/pfpp-heapstates/demographic data/201910gavi-4_dds-201910_2_life_ex_both.csv")
        # le <- le %>% filter(between(year, 1980,2019)) %>% filter(age_to <= 99) %>% filter(country_code %in% i)
        # le0 <- le %>% filter(age_from == 0) %>% slice(rep(1:n(), each=5))
        # le0$year <- 1980:2020
        #
        # le <- le %>% filter(age_from > 0)
        # le <- le %>% arrange(age_from, year) %>% slice(rep(1:n(), each = 5))
        # le$year <- rep(1980:2020, 20)
        # le <- le %>% arrange(year, age_from) %>% slice(rep(1:n(), each=5))
        # le$age_from <- rep(1:100, 41)
        # le <- le %>% filter(age_from < 100)
        # le <- bind_rows(le, le0) %>% arrange(age_from, year)
        # #
        # yll <- (le$value* t2$deaths)
        # yld <- t2$cases * 1/52 * (0.051  + 0.133)/2
        # t2$dalys <- yll + yld

        ## COHORT SIZE
        births <- read.csv("~/Downloads/pfpp-heapstates/demographic data/201910gavi-2_dds-201910_births_both.csv")
        cdr <- read.csv("~/Downloads/pfpp-heapstates/demographic data/201910gavi-3_dds-201910_cdr_both.csv")
        cohort_size <- numeric()
        
        pop <- read.csv("~/Downloads/pfpp-heapstates/demographic data/201910gavi-2_dds-201910_2_tot_pop_both.csv")
        pop <- pop %>% filter(between(year,1980,2020))
        pop <- pop %>% filter(country_code==i) %>% select(year, value)
        
        births.coh <- births[which(births$country_code == i),]
        births.coh <- births.coh[which(births.coh$year>=1950 & births.coh$year<=2020),"value"]
        births.coh <- c(births.coh,births.coh[length(births.coh)])  # add one more year to get 2100
        births.coh <- c(rep(births.coh[1],70),births.coh)  # add one more year to get 2100
        cdr.coh <- cdr[which(cdr$country_code == i),]
        cdr.coh <- cdr.coh[which(cdr.coh$year>=1950 & cdr.coh$year<=2020), "value"]
        cdr.coh <- c(rep(cdr.coh[length(cdr.coh)],70),cdr.coh)  # add one more year to get 2100
        years <- 1980:2020
        srate <- 1-cdr.coh
        
        for(k in 1:length(years)){
          cohort_size <-
            c(cohort_size, births.coh[(100+k):(1+k)] * c(1,cumprod(srate[(99+k):(1+k)])))
        }
        cohort_size <- data.frame(cohort_size)
        cohort_size$year <- rep(years, keach=100)
        chs <- c()
        
        for (y in 1980:2020) {
          newy <- filter(pop, year==y)
          newy <- newy[,2]
          newc <- cohort_size %>% filter(year==y) %>% select(cohort_size)
          su <- sum(newc)
          su <- (newy/su)
          newc <- newc * su
          chs <- c(chs, newc)
        }
        cohort_size <- matrix(unlist(chs), byrow=T, nrow=4100)
        # t2 <- cbind(t2, cohort_size)
        cohort_size <- as.data.frame(cohort_size)
        cohort_size <- cohort_size %>% add_column(year = rep(1980:2020,100))
        
        # t2 <- t2[c(5,1:3,6,11,8:10)]
        # t2 <- filter(t2, year >= 1981)
        # t2 <- t2[order(t2$year),]
        
        cohort_size <- cohort_size %>% filter(year>=1981)
        cohort_size <- as_tibble(cohort_size)
        cohort_size <- cohort_size %>% arrange(year)
       
        set <- cbind(ages, cohort_size)
        set <- set[c(1,2,3,4,5,10,7,8,9)]
        
        names(set)[6] <- "cohort_size"
        grid <- rbind(grid, set)
        grid <- grid[-c(1),]
        setwd("~/Desktop/central_burden_portnoy_nLL_new")
        write.csv(grid, paste0(i,'central_burden_portnoy.csv'), row.names=F)
      }
  #  }
 # }
# }
  