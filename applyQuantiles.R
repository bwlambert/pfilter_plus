
args = commandArgs(trailingOnly=TRUE)
n = 37
countrytag = args[1]
nsim = 1000

Ipaths <- NULL

rs <- sample(1:1000, 100)

	print(countrytag)
	x <-read.csv(paste0(countrytag,"_i_particles_fxboot",1,".txt"), sep=' ', row.names=NULL)
	row.names(x) <- 1:37
	Ipaths <- x[,rs]

	for(b in 0:999){
			x <-read.csv(paste0(countrytag,"_i_particles_fxboot",b,".txt"), sep=' ', row.names=NULL)
			row.names(x) <- 1:37
			Ipaths <- cbind(Ipaths,x[,rs])
	}

	Ilb <- apply(Ipaths,1, quantile, probs = 0.025, na.rm = TRUE)
	Iub <- apply(Ipaths,1, quantile, probs = 0.975, na.rm = TRUE)
	Imed <-apply(Ipaths,1, quantile, probs = 0.5, na.rm = TRUE)
	Imean <-apply(Ipaths,1, mean, na.rm = TRUE)

	write.table(Ilb, file = paste0(countrytag, "_smooth_Ilb.csv"), sep = ",")
	write.table(Iub, file = paste0(countrytag, "_smooth_Iub.csv"), sep = ",")
	write.table(Imean, file = paste0(countrytag, "_smooth_Imean.csv"), sep = ",")
	write.table(Imed, file = paste0(countrytag, "_smooth_Imed.csv"), sep = ",")
