library(abind)
##########################################################################################
# re-weighting function to handle underflow
# can replace smp.wt line in fx.4par
logsumexpnorm <- function(X){
  c <- max(X)
  (exp(X - c)) / sum(exp(X - c))
}

##########################################################################################

# Read in needed arguments:
  #samplesize = number of draws from parameters
  #pars = matrix of parameters from gridserach
  #nll = vector of nll from pf run - length equals dim(pars)[1]
  #CIcover = CI coverage
  #df = input dataframe
  #sia.info = input sia file
  #S0 = input initial susceptibles
  #npart = number of particles for particle filter
  #n.ages = number of age classes
  #country = country
  #sigma = noise in the attack rate function
  #cauchy.weight = mixture weight on cauchy in the observation step

# Rscript pf_summarize_age_driver.R 10 AFGnormal_runs_showing_penalty_ab1bdcf3a_220714.txt 0.95 AFGnormal data_for_fitting2020/AFG.txt data_for_fitting2020/siaAFG.txt 0 ab1bdcf3a 220714 220714 0.5 0.001

#Rscript pf_summarize_age_particles.R 100 #{c}#{m}_runs_showing_penalty_#{sha}_#{dte}.txt 0.95 #{c}#{m} Data_for_fitting2022/#{c}.txt Data_for_fitting2022/sia#{c}.txt #{mn} #{sha} #{dte} #{today} #{sig} #{cau} Iest_#{c}#{m}_#{sha}_#{dte}.txt
# Rscript pf_summarize_age_particles.R 10 TZAnormal_runs_showing_penalty_abe912ff6_230726.txt 0.95 TZAnormal Data_for_fitting2022/TZA.txt Data_for_fitting2022/siaTZA.txt 0 abe912ff6 230726 231103 0.5 0.01 Iest_TZAnormal_abe912ff6_230726.txt
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    samplesize <- args[1] # number of draws from parameters
    runmatrix <- args[2] # full run matrix, parameters and nLLs
    CIcover <- args[3] # CI coverage
    country.model <- args[4] # country
    df <- args[5] # input dataframe
    sia.info <- args[6] # input sia file
    model.number <- args[7] # input sia file
    sha <- args[8] # input sia file
    shanv <- sha # input sia file
    #shanv <- args[9] # input sia file
    tstamp <- args[9] # input sia file
    today <- args[10] # 
    npart <- 1000 #  number of particles for particle filter
    n.ages <- 100 # number of age classes
    n.yrs <- 42 # number of years
    sigma <- args[12] # noise in the attack rate function
    cauchy.weight <- args[13] # mixture weight on cauchy in the observation step

    print(runmatrix)
    run.df <- read.csv(runmatrix,header=TRUE,sep=' ')
    nll <- run.df$nLL # vector of nll from pf run - length equals dim(pars)[1]
    keep <- c("b0", "b1","pr","prhigh")
    pars <- data.matrix(run.df[keep])
    n <- length(nll)

    llweight <- logsumexpnorm(-nll)

    ind <- sample(1:n, samplesize, replace=TRUE, prob=llweight)

    samplepar <- pars[ind,]
    aa <- matrix()
    nll_new <- rep(NA,dim(samplepar)[1])
    print(samplepar)

    for(ii in 1:dim(samplepar)[1]){
        best.trans <- c(samplepar[ii,1],samplepar[ii,2],samplepar[ii,3], samplepar[ii,4])
        print(best.trans)
        iestcontent = c(0.0,0.0,best.trans)

        iestfn = paste(country.model,"_",sha,"_",tstamp,"_",ii,".txt",sep="")
        zz <- file(iestfn, "wb")
        writeBin( paste(iestcontent, collapse="\n"), zz ) 
        close(zz)

        country.model.iter = paste(country.model,"_",ii,sep="")

        ccommand <- paste("./pfpp 5",
                          country.model.iter,
                          df, 
                          sia.info,
                          model.number,
                          sigma,
                          cauchy.weight,
                          iestfn,
                          sep = " "
        )

        print(ccommand)

        # run pfpp with appropriate arguments:
        res <- system(ccommand)
        print(res)
    #

    }

    # read nll_new
    for(ii in 1:dim(samplepar)[1]){
        con=file(paste("Iest_single_",country.model,"_",ii,"_",shanv,"_",today,".txt",sep=""),open="r")
        linn=readLines(con,1)
        # close the connection
        close(con)
        nll_new[ii] <- as.numeric(linn) 
   }
    print(dim(aa))
    print("nll_new") 
    print(nll_new) 

    # read iap and produce dat. 
    for(ii in 1:dim(samplepar)[1]){
        ab <- as.matrix(read.table(paste(country.model,"_",ii,"_iap_particles_fxsave_",shanv,"_",today,".txt",sep=""), header = FALSE, sep = " "))
         zxs = array(ab,,dim=c(1000, n.yrs, 100))
         tzxs <- aperm(zxs,c(3,2,1))
	      ifelse(ii==1, aa <- tzxs, aa <- abind(aa, tzxs, along=3))
     }

    mean_traj <- apply(aa,1:2,mean,na.rm=T)
    print(mean_traj)
    CIcover <- as.numeric(CIcover)
    print(CIcover)
    lim <- (1 - CIcover)/2
    CI_traj <- apply(aa,1:2,quantile,probs=c(lim,1-lim),na.rm=T)
    #return(list(mean = mean_traj,CI = CI_traj, pars = samplepar))  

    out_data <- data.frame(
			   cbind(
				 rep(1:n.yrs,each=n.ages), # column of year for reference
				 rep(1:n.ages,times=n.yrs), # column of age clases for reference
                         as.vector(mean_traj),     # column of means
                         as.vector(CI_traj[1,,]),  # column of lower bounds
                         as.vector(CI_traj[2,,]))) # colum of upper bounds
                         # order of years is rep(1980:2019,n.ages)

    par_lik <- data.frame(cbind(samplepar,nll[ind],nll_new))
    names(out_data) <- c("year","age","mean_cases","lb_cases","ub_cases")
    names(par_lik) <- c("beta0","beta1","pobs","pobshigh","nll_orig","nll_new")

    result <- list(out_data=out_data,par_lik = par_lik)

    save(result,file=paste0("pf_summarizer_result_iage_",country.model,"_",sha,"_",shanv,"_",tstamp,"_",today,"_",sigma,"_",cauchy.weight,".dat"))
    long_file=paste0("long_pf_summarizer_result_iage_",country.model,"_",sha,"_",shanv,"_",tstamp,"_",today,"_",sigma,"_",cauchy.weight,".csv")
    write.csv(result$out_data, long_file, row.names = FALSE) 

    aa <- matrix()
    # read iap and produce dat. 
    for(ii in 1:dim(samplepar)[1]){
        ab <- as.matrix(read.table(paste(country.model,"_",ii,"_sap_particles_fxsave_",shanv,"_",today,".txt",sep=""), header = FALSE, sep = " "))
         zxs = array(ab,,dim=c(1000, n.yrs, 100))
         tzxs <- aperm(zxs,c(3,2,1))
	      ifelse(ii==1, aa <- tzxs, aa <- abind(aa, tzxs, along=3))
     }

    mean_traj <- apply(aa,1:2,mean,na.rm=T)
    print(mean_traj)
    CIcover <- as.numeric(CIcover)
    print(CIcover)
    lim <- (1 - CIcover)/2
    CI_traj <- apply(aa,1:2,quantile,probs=c(lim,1-lim),na.rm=T)
    #return(list(mean = mean_traj,CI = CI_traj, pars = samplepar))  

    out_data <- data.frame(
			   cbind(
				 rep(1:n.yrs,each=n.ages), # column of year for reference
				 rep(1:n.ages,times=n.yrs), # column of age clases for reference
                         as.vector(mean_traj),     # column of means
                         as.vector(CI_traj[1,,]),  # column of lower bounds
                         as.vector(CI_traj[2,,]))) # colum of upper bounds
                         # order of years is rep(1980:2019,n.ages)

    par_lik <- data.frame(cbind(samplepar,nll[ind],nll_new))
    names(out_data) <- c("year","age","mean_susceptibles","lb_susceptibles","ub_susceptibles")
    names(par_lik) <- c("beta0","beta1","pobs","pobshigh","nll_orig","nll_new")

    result <- list(out_data=out_data,par_lik = par_lik)

    save(result,file=paste0("pf_summarizer_result_sage_",country.model,"_",sha,"_",shanv,"_",tstamp,"_",today,"_",sigma,"_",cauchy.weight,".dat"))
    long_file=paste0("long_pf_summarizer_result_sage_",country.model,"_",sha,"_",shanv,"_",tstamp,"_",today,"_",sigma,"_",cauchy.weight,".csv")
    write.csv(result$out_data, long_file, row.names = FALSE) 
}

# Sleep for a random number of seconds between 0 and 600
#Sys.sleep(runif(1, min = 0, max = 180))

main()


