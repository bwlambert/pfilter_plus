# like grid_plots_FAKE.R but for grid within grid
rm(list = ls())
library(fields)

                                        #country = "NGA"
#setwd("~/Dropbox/GAVI/Code_May2017/")
setwd(".")
batch = 1
country = 2
# regular
#params <- read.table(paste("~/Dropbox/GAVI/Code_May2017/Batch1/sim_parameters.txt", sep = ""),header = TRUE)[country,]
                                        #dir.create(fig.here <- paste0("./SIMcountry_",country, "_figures/"), showWarnings = FALSE)
fig.here <- "." #paste0("~/Dropbox/GAVI/Code_May2017/Batch", batch,"/Data",country)
        
                                        #load(paste("~/Dropbox/GAVI/Code_January2017/", country, "-example-4par.Rdat", sep = ""))
s = 0.4

palette(colorRampPalette(c("yellow", "green", "blue"))(n = 100))


library(scales)

##################################################################################
#Likelihood Figures
#for(s in c(.25, .5, .75)){
    load(paste(fig.here,"/GridSearchInfo_",s,".Rdat", sep = ""))

                                        #for(r in c(1:2)){

# collect run 1 and run 2 togeter (run2 is smaller grid range)
 cparams <- rbind(run1[[1]],run2[[1]])
 cnLL <- c(run1[[2]],run2[[2]])
run.all <- list(params = cparams, nLL = cnLL)
#rm.obs <- which(run.all$params[,2]< 0)
run <- run.all
#run$params <- run$params[-rm.obs,]
#run$nLL <- run$nLL[-rm.obs]
#best <- run$params[which.min(run$nLL),]
#        run = eval(parse(text=paste0("run", r)))
        #run$nLL <- as.numeric(run$nLL)
        col.val = ceiling((run$nLL-min(run$nLL, na.rm = T))/(max(run$nLL, na.rm = T)-min(run$nLL, na.rm = T))*100)
        col.val[which(col.val == 0)] = 1
        
### FIgure 1
        pdf(paste0(fig.here, "/likelihood_s", s, ".pdf"))
        par(mfrow = c(2,2))
        plot(run[[1]][,1],run$nLL, col = alpha(col.val), pch = 19, xlab = "Beta0", ylab  = "Neg Log-Likelihood")
abline(v = params[1])
abline(v = best[1], col = "red")
                                        #        lines(supsmu(run[[1]][,1], run$LL), col = "red")
legend("topright", bty = "n", col = c("black", "red"), lwd = 2, c("Sim values", "Best Estimate"))
        
        plot(run[[1]][,2],run$nLL, col = alpha(col.val), pch = 19, xlab = "Beta1", ylab  = "Neg Log-Likelihood")
abline(v = params[2])
abline(v = best[2], col = "red")

                                        #       lines(supsmu(run[[1]][,2], run$LL), col = "red")

        plot(run[[1]][,3],run$nLL, col = alpha(col.val), pch = 19, xlab = "P(reporting)", ylab  = "Neg Log-Likelihood")
abline(v = params[3])
abline(v = best[3], col = "red")

                                        #      lines(supsmu(run[[1]][,3], run$LL), col = "red")
        
        plot(run[[1]][,1], run[[1]][,2],col = alpha(col.val,.5), pch = 19,
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
        image.plot( legend.only=TRUE,col = palette(), zlim= range(run$nLL, na.rm = T), horizontal = TRUE)
dev.off()


### FIGURE 2
        pdf(paste0(fig.here, "/likelihoods2_", s, ".pdf"))
        par(mfrow = c(2,2))
        plot(run[[1]][,1],run[[1]][,3], col = alpha(col.val), pch = 19, xlab = "Beta0", ylab  = "Prob(Reporting)")
abline(h = params[3], v = params[1])
abline(h = best[3], v = best[1], col = "red")
legend("topright", bty = "n", col = c("black", "red"), lwd = 2, c("Sim values", "Best Estimate"))


        plot(run[[1]][,2],run[[1]][,3], col = alpha(col.val), pch = 19, xlab = "Beta1", ylab  = "Prob(Reporting)")
abline(h = params[3], v = params[2])
abline(h = best[3], v = best[2], col = "red")

        plot(-run[[1]][,1]/run[[1]][,2], run[[1]][,3],col = alpha(col.val,.5), pch = 19,
             ylab = "P(reporting)", xlab = "-beta0/beta1", xlim = c(0,.25))
abline(h = params[3], v = -params[1]/params[2])
abline(h = best[3], v = -best[1]/best[2], col = "red")
        
        image.plot( legend.only=TRUE,col = palette(), zlim= range(run$nLL, na.rm = T), horizontal = FALSE)
        plot(run[[1]][,1], run[[1]][,2],col = alpha(col.val,.5), pch = 19,
             xlab = "Beta0", ylab = "Beta1")
abline(h = params[2], v = params[1])
abline(v = best[1],h = best[2], col = "red")

dev.off()

########################################################################
##################################################################################
## clouds
#r = 2
    df <- read.table(paste0(fig.here, "/Data.txt"), header = TRUE)[1:33,]

pdf(paste0(fig.here, "/clouds_", s, ".pdf"))
#par(mfrow = c(3,1))
                                        #for(s in c(.5, .75)){

#    load(paste(fig.here,"/GridSearchInfo_",s,".Rdat", sep = ""))
                                        #    load(paste("KE_workflow2_",country, "_s",s,".Rdat", sep = ""))
#load(paste("KE_workflow3_",country, "_s",s,"_fixed.Rdat", sep = ""))
#    cov = sum(coverage, na.rm = TRUE)/32
#cparams <- rbind(run1[[1]],run2[[1]])
#    cnLL <- c(run1[[2]],run2[[2]])
#    run <- list(params = cparams, nLL = cnLL)
                                        #   best.params <- run$params[which.min(run$nLL),]
best.params <- best

    prob.rep <- rep(best.params[3], dim(df)[1])  #assume baseline prob of reprting
    #prob.rep[which(is.na(df$high.years))] = best.params[4] # replace high years with inflated prob of reporting
    I = Iest$Isum
    C <- df$C/prob.rep
    Iest95 <- apply(Iest$Isum,1,quantile,probs=c(.025,.975),na.rm=T)
    my.cov <- (C > Iest95[1,] & C < Iest95[2,])
    cov = sum(my.cov, na.rm = TRUE)/(length(my.cov) - sum(is.na(my.cov)))
    yrange = range(I, C, na.rm = TRUE)


    matplot(I,col = alpha("grey", .2), pch = 19, ylim = yrange, main = paste0("S = ", s))
    points(C, col = "red", pch = 19)
    points(Iest95[1,], pch = "-", col = "black")
    points(Iest95[2,], pch = "-", col = "black")
    legend("topright", c("Prediction Cloud", "Cases/P(Reporting)", paste0("Coverage = ", round(cov,3)*100, "%")), col = c("grey", "red", "white"), pch = 19)
dev.off()

