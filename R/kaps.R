##########################################################################
####
#### Multiway-splits adaptive partitioning with distribution free approach 
####
####		Soo-Heang EO and HyungJun CHO 
####
####		Version 1.0.0
####
####		30 Oct, 2013
####
##########################################################################

kaps <- function(formula, data, K = 2:5, mindat, type = c("perm", "Boot", "random","test", "NULL"), ...){

	# For debugging
	# K = 2:3; mindat = 15; CV = TRUE; minors = kaps.control(V = 5) # add CV option in kaps.control
	##########################
	##### pre-processing step
	#options(warn = -1)
	if(missing(mindat)) mindat = floor(nrow(data) * 0.05)
	if(!is.Formula(formula)) {
		formula <- as.Formula(formula)
		#class(formula) <- "Formula"
	}
	# for minor options used in kaps algorithm
	# minors = kaps.control()
	minors = kaps.control(...)

	if(any(K == 1)) stop("the number of subgroups (K) have to be greater than subjects.")

	n <- nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n

    if(n == 0L) stop("0 (non-NA) cases.")
	if(length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	if(length(K) > 10) stop("the maximum number of subgroups (K) is too large.")

	call <- match.call()

	######################################################
	###	Finding optimal subgroups
	type <- match.arg(type)
	if(length(K) == 1 & type == "NULL"){
		result = kaps.fit(formula = formula, data = data, K = K, mindat = mindat, minors = minors)
		test.stat2 = adj.test(result)
		result@call <- call
		result@Z <- test.stat2[1] #overall statistic 
		result@X <- test.stat2[2] #optimal pair p-value
		result@WH <- test.stat2[3] #cube-root transformation
		result@t <- test.stat2[4] #WH approximation statistic
		return(result)
	} else if(type == "perm"){
		# Fit KAPS with simple permutation test
		cat("Now, finding optimal number of subgroups (K) by simple permetation. \n")
		# have to build naive scripts for the selection of subgroups
		test.stat = kaps.CV(formula, data, K, V = minors@V, minors)
		# choose maximal pairwise permutation p-value
		# NEED TO MODIFY
	} else if(type == "test"){
		cat("Now, finding optimal number of subgroups (K) by test estimates. \n")
		#lapply(K, kaps.fit, formula = formula, data = data, mindat = mindat, minors = minors)
		test.stat = kaps.test(formula, data, K, minors)
		# choose maximal pairwise permutation p-value
		# NEED TO MODIFY
	} else if(type %in% c("Boot", "random")){
		cat("Now, finding optimal number of subgroups (K) by Double permuting kaps. \n")
		fit.boot = kaps.boot(formula = formula, data = data, K = K, B = minors@N.boot, minors = minors, type = type)
		test.stat = fit.boot$stat
		# choose maximal pairwise permutation p-value
		index = 1:nrow(test.stat)
		index.pair = test.stat[,3] <= 0.05
		#index.pair <- ifelse(is.na(index.pair), FALSE, index.pair)

		if(all(index.pair == FALSE)){
			fit = lapply(K, kaps.fit, formula = formula, data = data, mindat= mindat, minors = minors)
			test.stat2 = sapply(fit, adj.test)
			test.stat2 <- as.matrix(test.stat2)
			
			result = new("kaps")
			result@index <- as.integer(0)
			result@groups <- 1
			attr(result@groups,"names") <- paste("K=",1,sep="")
			result@Z <- as.vector(test.stat2[1, ]) #overall test statistic 
			result@X <- as.vector(test.stat2[2, ]) # test statistic for the worst pair
			result@WH <- as.vector(test.stat2[3, ]) #cube-root transformation
			result@t <- as.vector(test.stat2[4, ]) #WH approximation statistic
			result@pvalue <- test.stat[1,c(1,3)] # overall and worst pair p-values for the selected 
			result@elbow <- test.stat
			result@over.stat.sample <- fit.boot$boot.over.stat
			result@pair.stat.sample <- fit.boot$boot.pair.stat
			result@call <- call
			result@Options <- minors
			return(result)
		} else{
			index <- index[sum(index.pair)]
			# choose minimal overall permutation p-value
			#index = index[test.stat[index.pair,1] == max(test.stat[index.pair,1], na.rm = TRUE)]
			#index <- index[sum(index.pair)]
			#if(length(index) >= 2) index <- index[length(index)]
		}
	} else if (type == "NULL"){
		index = 1
	}

	# select optimal K
	K.sel = K[index]

	######################################################
	# Obtain log-rank statistics at K
	### parallel computing in order to find optimal k subgroups
	cat("Now, selecting a set of cut-off points...\n")
	#fit = kaps.fit(formula = formula, data = data, K = K.sel, mindat = mindat, minors = minors)
	fit = lapply(K, kaps.fit, formula = formula, data = data, mindat= mindat, minors = minors)
	test.stat2 = sapply(fit, adj.test)
	test.stat2 <- as.matrix(test.stat2)

	result <- fit[[index]]
	result@index <- as.integer(index)
	result@groups <- K
	attr(result@groups,"names") <- paste("K=",K,sep="")
	result@Z <- as.vector(test.stat2[1, ]) #overall statistic 
	result@X <- as.vector(test.stat2[2, ]) #optimal pair p-value
	result@WH <- as.vector(test.stat2[3, ]) #cube-root transformation
	result@t <- as.vector(test.stat2[4, ]) #WH approximation statistic
	result@results <- fit
	result@pvalue <- test.stat[index,c(1,3)] # overall and worst pair p-values for the selected group
	if(type != "NULL") {
		result@elbow <- test.stat
		result@over.stat.sample <- fit.boot$boot.over.stat
		result@pair.stat.sample <- fit.boot$boot.pair.stat
	}
	result@call <- call
	result@Options <- minors
	return(result)
}


############################################################
# Fit KAPS with Permutation test by cross validation approach
############################################################
kaps.CV <- function(formula, data, K = 2:5, V = 5, minors = kaps.control()){

	#options(warn = -1)
	#cat("Plz, wait some minutes.... \n")
	if(any(K == 1)) stop("the number of subgroups (K) is greater than 1.")

	n = nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n

    if(n == 0L) stop("0 (non-NA) cases.")
	if(length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	#if(length(K) > 10) stop("the maximum number of subgroups (K) is too large.")

	call = match.call()

	# kaps by V-fold cross-validation
	if(V < 2) stop("V must be an integer which is greater than and equal to 2.")
	
	fold.no = rep(1:V, nrow(data) %/% V) # generate new var.: fold number
	fold.no = c(fold.no, (1:V)[nrow(data) %% V])
	fold.over.stat = fold.pair.stat = matrix(0, nrow = V, ncol = length(K))
	result = matrix(0, nrow = length(K), ncol = 4)
	colnames(result) <- c("over_mean", "over_se", "pair_mean", "pair_se")
	rownames(result) <- paste("K=",K,sep="")
	## CHECK ME: reduce computing time by parallel computing
	if(minors@ncl == 1){
		for(i in 1:V){
			cat("V =",i,"\n")	
			learning = data[fold.no != i,, drop = FALSE]
			test = data[fold.no == i,, drop = FALSE]
			mindat = floor(nrow(learning) * 0.05)
		
			fit = lapply(K, kaps.fit, formula = formula, data = learning, mindat = mindat, minors = minors)
			test.stat = sapply(fit, kaps.perm, newdata = test)			
			rownames(test.stat) =  c( "perm_overall_pval", "perm_min_pair_pval")
			colnames(test.stat) = paste("K=",K,sep="")
			print(round(test.stat,3))

			fold.over.stat[i,] <- test.stat[1,]
			fold.pair.stat[i,] <- test.stat[2,]
			rm(fit, train, test, test.stat)
		}
	}
	else{
		ncl <- makeCluster(minors@ncl)
		Export.func <- c("as.Formula", "model.part", "survdiff", "Surv")
		clusterExport(ncl, Export.func)
		for(i in 1:V){	
			cat("V =",i,"\n")	
			train = data[fold.no != i,, drop = FALSE]
			test = data[fold.no == i,, drop = FALSE]
			fit = parLapply(ncl, K, kaps.fit, formula = formula, data = learning, mindat = mindat, minors = minors)
			test.stat = parSapply(ncl, fit, kaps.perm, newdata = test)			
			rownames(test.stat) =  c( "perm_overall_pval", "perm_min_pair_pval")
			colnames(test.stat) = paste("K=",K,sep="")
			print(round(test.stat,3))

			fold.over.stat[i,] <- test.stat[1,]
			fold.pair.stat[i,] <- test.stat[2,]
			rm(fit, train, test, test.stat)
		}
		stopCluster(ncl)
	}
	### output
	result[,1] <- apply(fold.over.stat, 2, mean, na.rm = TRUE)
	result[,2] <- apply(fold.over.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	result[,3] <- apply(fold.pair.stat, 2, mean, na.rm = TRUE)
	result[,4] <- apply(fold.pair.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	#result <- as.data.frame(result)
	return(result)
}

############################################################
# Fit KAPS with test estimates approach 
############################################################
kaps.test <- function(formula, data, K = 2:5, minors = kaps.control()){

	#options(warn = -1)
	
	if(any(K == 1)) stop("the number of subgroups (K) is greater than 1.")

	n = nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n

    if(n == 0L) stop("0 (non-NA) cases.")
	if(length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	#if(length(K) > 10) stop("the maximum number of subgroups (K) is too large.")

	#call = match.call()

	# kaps by test estimates approach
	fold.over.stat = fold.pair.stat = matrix(0, nrow = 1, ncol = length(K))
	result = matrix(0, nrow = length(K), ncol = 4)
	colnames(result) <- c("over_pval", "over_stat", "pair_pval", "pair_stat")
	rownames(result) <- paste("K=",K,sep="")

	## CHECK ME: reduce computing time by parallel computing
	index = sample(1:n, floor(n*0.7))
	learning = data[index,, drop = FALSE]
	test = data[-index,, drop = FALSE]
	mindat = floor(nrow(learning) * 0.05)
	
	fit = lapply(K, kaps.fit, formula = formula, data = learning, mindat = mindat, minors = minors)
	test.stat = sapply(fit, kaps.perm, newdata = test)			
	rownames(test.stat) =  c( "perm_overall_pval", "perm_min_pair_pval")
	colnames(test.stat) = paste("K=",K,sep="")
	print(round(test.stat,3))

	fold.over.stat[1,] <- test.stat[1,]
	fold.pair.stat[1,] <- test.stat[2,]
	### output
	result[,1] <- fold.over.stat
	#result[,2] <- apply(fold.over.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	result[,3] <- fold.pair.stat
	#result[,4] <- apply(fold.pair.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	#result <- as.data.frame(result)
	return(result)
}

# END @ OCT 19, 2013