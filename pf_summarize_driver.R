
##########################################################################################
# re-weighting function to handle underflow
# can replace smp.wt line in fx.4par
logsumexpnorm <- function(X){
  c <- max(X)
  (exp(X - c)) / sum(exp(X - c))
}

##########################################################################################
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


#  Rscript pf_summarize_driver.R 10 subset_sigma05_cauchy001/SOMnormal_runs_showing_penalty.txt 0.95 SOM data_for_fitting2019/SOM.txt data_for_fitting2019/siaSOM.txt 0.5 0.001


main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    samplesize <- args[1] # number of draws from parameters
    runmatrix <- args[2] # full run matrix, parameters and nLLs
    CIcover <- args[3] # CI coverage
    country.model <- args[4] # country
    df <- args[5] # input dataframe
    sia.info <- args[6] # input sia file
    model.number <- args[7] # input sia file
    npart <- 1000 #  number of particles for particle filter
    n.ages <- 100 # number of age classes
    sigma <- args[8] # noise in the attack rate function
    cauchy.weight <- args[9] # mixture weight on cauchy in the observation step
    iestfile <- args[10] # mixture weight on cauchy in the observation step

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
        #best.trans <- c(log(-samplepar[ii,1]),logit(samplepar[ii,2]/250),logit(samplepar[ii,3]), logit(samplepar[ii,4]))	# note that transformations aren't critical for grid search, so we should switch to all transformed or all natural scale
        best.trans <- c(samplepar[ii,1],samplepar[ii,2],samplepar[ii,3], samplepar[ii,4])
        print(best.trans)
        iestcontent = c(0.0,0.0,best.trans)

        iestfn = paste(country.model,"_",ii,".txt",sep="")
        zz <- file(iestfn, "wb")
        writeBin( paste(iestcontent, collapse="\n"), zz ) 
        close(zz)

        country.model.iter = paste(country.model,"_",ii,sep="")

        #    Iest.states<-fx.4par(best.trans, df = df,sia.info=sia.info,S0=S0, npart = npart, n.ages=n.ages,age.dist = NA, directory=state.loc, b=NA, yr.adj = df2$high.years,save=TRUE,country=country, sigma=s, cauchy.weight = .01)
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

        con=file(paste("Iest_single_",country.model,"_",ii,".txt",sep=""),open="r")
        linn=readLines(con,1)

        nll_new[ii] <- as.numeric(linn) 

        ab <- as.matrix(read.table(paste(country.model,"_",ii,"_i_particles_fxsave.txt",sep=""), header = FALSE, sep = " ", skip = 1))
        #
           #aa <- Iest.states$states[,2,,], # I states by age class
           #aa <- abind(aa, Iest.states$Isum,along = 3)) # append next run
        ifelse(ii==1,
               aa <- ab,
               aa <- cbind(aa, ab))
    }

    print("nll_new") 
    print(nll_new) 

    save(aa,file=paste0("aa_",country.model,"_",sigma,"_",cauchy.weight,".dat"))

    mean_traj <- apply(aa,1:2,mean,na.rm=T)
    CIcover <- as.numeric(CIcover)
    print(CIcover)
    lim <- (1 - CIcover)/2
    #CI_traj <- apply(aa,1,quantile,probs=c(lim,1-lim),na.rm=T)
    CI_traj <- apply(aa,1:2,quantile,probs=c(lim,1-lim),na.rm=T)
    #return(list(mean = mean_traj,CI = CI_traj, pars = samplepar))  

    #out_data <- data.frame(cbind(mean_traj,t(CI_traj)))
    out_data <- data.frame(cbind(rep(1:n.ages,each=samplesize), # column of age clases for reference
                         as.vector(mean_traj),     # column of means
                         as.vector(CI_traj[1,,]),  # column of lower bounds
                         as.vector(CI_traj[2,,]))) # colum of upper bounds
                         # order of years is rep(1980:2019,n.ages)

    par_lik <- data.frame(cbind(samplepar,nll[ind],nll_new))
    names(out_data) <- c("age","mean_cases","lb_cases","ub_cases")
    names(par_lik) <- c("beta0","beta1","pobs","pobshigh","nll_orig","nll_new")

    #return(list(out_data=out_data,par_lik = par_lik))
    result <- list(out_data=out_data,par_lik = par_lik)

    print(out_data)
    save(result,file=paste0("pf_summarizer_result_",country.model,"_",sigma,"_",cauchy.weight,".dat"))

}

main()


