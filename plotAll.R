
cgridplot <- function(iso){
#library(fields)
#library(scales)

print(iso)
# regular
#params <- read.table(paste("~/Dropbox/GAVI/Code_May2017/Batch1/sim_parameters.txt", sep = ""),header = TRUE)[country,]
                                        #dir.create(fig.here <- paste0("./SIMcountry_",country, "_figures/"), showWarnings = FALSE)
fig.here <- getwd()#dirname(gridsearchinfo) #"." #paste0("~/Dropbox/GAVI/Code_May2017/Batch", batch,"/Data",country)
                                        #load(paste("~/Dropbox/GAVI/Code_January2017/", country, "-example-4par.Rdat", sep = ""))
print(fig.here)
s <- 0.5

palette(colorRampPalette(c("yellow", "green", "blue"))(n = 100))

##################################################################################
#Likelihood Figures
    #load(paste(fig.here,"/GridSearchInfo_",s,".Rdat", sep = ""))

#    load(gridsearchinfo)
                                        #for(r in c(1:2)){
irun1 <- read.csv(paste0(iso,"_run1.txt"),sep=" ",header=TRUE)
irun2 <- read.csv(paste0(iso,"_run2.txt"),sep=" ",header=TRUE)

inll1 <- as.vector(irun1[5])
irun1$nLL <- NULL

#irun2$nLL <- NULL
#colnames(irun2)[4] <- "nLL"
inll2 <- as.vector(irun2[5])
irun2$nLL <- NULL
#irun1 <- c(data.matrix(irun1))
#irun2 <-  c(data.matrix(irun2))
run1 <- c()
run2 <- c()
run1$params  <- data.matrix(irun1)
run2$params <- data.matrix(irun2)
run1$nLL <- inll1
run2$nLL <- inll2
 
print(run2$nLL)
print(run2$params)

# collect run 1 and run 2 togeter (run2 is smaller grid range)
 cparams <- rbind(run1[[1]],run2[[1]])
params <- cparams
 cnLL <- c(run1[[2]],run2[[2]])

run.all <- list(params = cparams, nLL = cnLL)
#rm.obs <- which(run.all$params[,2]< 0)
run <- run.all
run$nLL <- as.numeric(unlist(cnLL))
#run$params <- run$params[-rm.obs,]
#run$nLL <- run$nLL[-rm.obs]
best <- run$params[which.min(run$nLL),]
#print(paste0(c," ",best)#        run = eval(parse(text=paste0("run", r)))
        run$nLL <- as.numeric(run$nLL)
        col.val = ceiling((run$nLL-min(run$nLL, na.rm = T))/(max(run$nLL, na.rm = T)-min(run$nLL, na.rm = T))*100)
        col.val[which(col.val == 0)] = 1
### FIgure 1
        print(length(run$nLL))
        print(length(run[[1]][,1]))
        print(run[[1]][,1])
png(paste0("run_",iso,"_likelihood_s.png"), width = 1000, height = 720, units = 'px')
        #pdf(paste0(fig.here, "/likelihood_s", s, ".pdf"))
        par(mfrow = c(2,2))
        plot(run[[1]][,1],run$nLL, col = scales::alpha(col.val), pch = 19, xlab = "Beta0", ylab  = "Neg Log-Likelihood")
abline(v = params[1])
abline(v = best[1], col = "red")
                                        #        lines(supsmu(run[[1]][,1], run$LL), col = "red")
legend("topright", bty = "n", col = c("black", "red"), lwd = 2, c("Sim values", "Best Estimate"))

        plot(run[[1]][,2],run$nLL, col = scales::alpha(col.val), pch = 19, xlab = "Beta1", ylab  = "Neg Log-Likelihood")
abline(v = params[2])
abline(v = best[2], col = "red")

                                        #       lines(supsmu(run[[1]][,2], run$LL), col = "red")

        plot(run[[1]][,3],run$nLL, col = scales::alpha(col.val), pch = 19, xlab = "P(reporting)", ylab  = "Neg Log-Likelihood")
abline(v = params[3])
abline(v = best[3], col = "red")

                                        #      lines(supsmu(run[[1]][,3], run$LL), col = "red")

        plot(run[[1]][,1], run[[1]][,2],col = scales::alpha(col.val,.5), pch = 19,
             xlab = "", ylab = "Beta1",xaxt = "n")
        points(best[1], best[2], pch = "*", col = "red", cex = 2)
        abline(v = params[1])
abline(h = params[2])
abline(v = best[1],h = best[2], col = "red")

        axis(3)
        mtext(side = 3, "Beta0", line =2)
        text(best[1], best[2], pos = 3, col = "red",
             paste0("Beta0 = ", round(best[1],2), ", Beta1 = ", round(best[2],2)), cex = .8)
        text(best[1], best[2]*1.1, pos = 3,col = "red",
             paste0("P(Reporting) ", round(best[3],2)), cex = .8)
      #  text(best[1], best[2]+3, pos = 3,col = "red",
           #  paste0("P(High Rep) = ", round(best[4],2)), cex = .8)
        fields::image.plot( legend.only=TRUE,col = palette(), zlim= range(run$nLL, na.rm = T), horizontal = TRUE)
dev.off()


### FIGURE 2
       # pdf(paste0(fig.here, "/likelihoods2_", s, ".pdf"))
png(paste0("run_",c,"_likelihood2_s.png"), width = 1000, height = 720, units = 'px')
        par(mfrow = c(2,2))
        plot(run[[1]][,1],run[[1]][,3], col = scales::alpha(col.val), pch = 19, xlab = "Beta0", ylab  = "Prob(Reporting)")
abline(h = params[3], v = params[1])
abline(h = best[3], v = best[1], col = "red")
legend("topright", bty = "n", col = c("black", "red"), lwd = 2, c("Sim values", "Best Estimate"))


        plot(run[[1]][,2],run[[1]][,3], col = scales::alpha(col.val), pch = 19, xlab = "Beta1", ylab  = "Prob(Reporting)")
abline(h = params[3], v = params[2])
abline(h = best[3], v = best[2], col = "red")

        plot(-run[[1]][,1]/run[[1]][,2], run[[1]][,3],col = scales::alpha(col.val,.5), pch = 19,
             ylab = "P(reporting)", xlab = "-beta0/beta1", xlim = c(0,.25))
abline(h = params[3], v = -params[1]/params[2])
abline(h = best[3], v = -best[1]/best[2], col = "red")

        fields::image.plot( legend.only=TRUE,col = palette(), zlim= range(run$nLL, na.rm = T), horizontal = FALSE)
        plot(run[[1]][,1], run[[1]][,2],col = scales::alpha(col.val,.5), pch = 19,
             xlab = "Beta0", ylab = "Beta1")
abline(h = params[2], v = params[1])
abline(v = best[1],h = best[2], col = "red")

dev.off()

}

args = commandArgs(trailingOnly=TRUE)
fc = args[1]
c = fc
print(fc)
x <- paste0(fc,"_i_particles_fxsave.txt")
    if(file.exists(x)){



png(paste0("run_",fc,".png"), width = 1000, height = 720, units = 'px')

x <- read.csv(paste0(fc,"_i_particles_fxsave.txt"),sep=" ",header=FALSE)
mx <- as.matrix(x)
plot(mx[,1],main=paste0(c),ylab="I",xlab="Year",ylim=c(0,7000000))
#plot(mx[,1])
for(i in 1:1000){ lines(mx[,i],col="blue")}

#x <- read.table(paste0("~/code/PFBurdenEstimator/inst/extdata/",c,".txt"),header=TRUE)
x <- read.table(paste0("Data_for_fitting/",c,".txt"),header=TRUE)
lines(x$C,col="red",lwd=3)

sa <- x$SIA
print(x)
for(k in 1:36){

    if(sa[k]>0){
        abline(v=k, col="orange", lwd=3, lty=2)
    }
    if(is.na(x$high.years[k])){
        abline(v=k, col="purple", lwd=3, lty=2)
    }
}
dev.off()

if(file.exists(paste0(fc,"_run1.txt"))){
cgridplot(fc)
}
}
