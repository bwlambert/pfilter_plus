setwd("~/Desktop/case_outputs_nLL_orig")
files <- list.files()
isos <- str_sub(files, 1, 3)
nYear <- 40
nStats <- 3

array1 <- array(NA, dim=c(length(isos),nYear,nStats))
rownames(array1) <- c(isos)
colnames(array1) <- 1981:2020
for(ii in isos){
  tt <- read.csv(paste0(ii, "_case_totals.csv"))
  tt <- tt[-1,]
  array1[ii,,1] <- as.numeric(tt$lb)
  array1[ii,,2] <- as.numeric(tt$mean)
  array1[ii,,3] <- as.numeric(tt$ub)
}

lb <- array1[,,1]
mean <- array1[,,2]
ub <- array1[,,3]

setwd("~/Desktop")
write.csv(lb, "lb.csv")
write.csv(ub, "ub.csv")
write.csv(mean, "mean.csv")
