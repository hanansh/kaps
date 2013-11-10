###################################################################
#### KAPS for change point detection
#### 	Soo-Heng EO and HyungJun CHO
#### 	Korea University, Seoul.
#### 	Last updated Oct 6, 2013
###################################################################
kaps.fit <- function(formula, data, K, mindat, minors) {
	#######################
	## pre-processing step
	n <- nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n
    if(n == 0L) stop("0 (non-NA) cases.")

	# original formula
    #mf <- match.call() # save a copy of the call
	#formula <- as.Formula(formula)
	#mf$formula <- f
	#mf <- eval(mf, parent.frame())
	#vars <- all.vars(formula)
	
	if(!is.data.frame(data)) data <- as.data.frame(data)
	## Set the object of result
    result <- new("kaps")
	result@formula <- formula
	
	##### Model Fitting
	X <- model.part(formula, data = data, rhs = 1, drop = FALSE)
	f <- update(formula, . ~ K)

	## treat pre-determined split points
	if(is.null(attr(minors@pre.pt,"names"))) { 
		pt.set <- lapply(X, function(x) sort(unique(x)) )
	} else{
		pre.name <- colnames(X) != attr(minors@pre.pt, "names")
		pt.set <- lapply(X[,pre.name, drop = FALSE], function(x) sort(unique(x)) )
		pt.set <- c(pt.set, minors@pre.pt)
		X <- X[,attr(pt.set, "names"), drop = FALSE]
	}

	## treat pre-determined ranges
	scope <- minors@scope
	if(!is.null(attr(scope, "names"))) {
		scope.vars <- which(names(pt.set) == attr(scope, "names"))
		pt.set.name <- names(pt.set)
	}
	
	result@X <- 0
	pt.set <- lapply(pt.set, function(x,upper, lower) x[x <= upper & x >= lower], upper = minors@upper.limit, lower = minors@lower.limit)
	v <- K-1
	####################################
	# NEED FORTRAN ROUTINE UPDATE
	# 1. group.sel() goes to FORTRAN 90 routine
	# 2. gorup.sel() with loop goes to FORTRAN 90 routine 

	for(i in 1: length(pt.set)){
		if(!is.null(attr(scope, "names"))){
			if(i %in% scope.vars) {
				rngs <- scope[[pt.set.name[i]]]  
				pt.set[[i]] <- pt.set[[i]][ pt.set[[i]] >= rngs[1] & pt.set[[i]] <= rngs[2] ]
			}
		}

		ord <- order(X[,i])
		data <- data[ord,]
		X <- X[ord, , drop = FALSE]
		test.where <- group.sel(X[,i], pt.set[[i]], K, mindat, data, f, minors)
		x.test <- test.where$test
		
		##### Result
		if(x.test[1,test.where$index] >= result@X) {
			index <- test.where$index
			result@X <- x.test[1,index] # pairwise test statistic X
			#result@WH <- (x.test[1,index] / v)^(1/3)
			#result@t <- (result@WH - (1- (2 / (9*v)))) / sqrt(2 / (9 * v))
			result@pvalue <-  x.test[2,index] # pairwise p-value
			result@pair <- x.test[3,index] # pair for selected
			result@groupID <- test.where$where
			result@split.pt <- sapply(unique(result@groupID), function(x,y,where) max(y[where == x]), y = X[,i], where = result@groupID) 
			result@split.pt <- result@split.pt[-length(result@split.pt)]
			result@split.var <- colnames(X[,i,drop =FALSE])
			result@data <- data
		}
	}
	####################################
	result@groups <- K
	attr(result@groups,"names") <- paste("G=",K,sep="")
	result@mindat <- mindat
	result@Options <- minors
	return(result)
}

############################################################
# Fit KAPS with Cross-Vaildation
############################################################
CV.kaps.fit <- function(K, V, formula, data, mindat, minors){
	## Finding optimal K by cross-validation
	mf <- match.call() 
	formula <- as.Formula(formula)
	## kaps by V-fold Cross-Validation
	if(V < 2) stop("V must be an integer which is greater than and equal to 2.")
	fold.no <- rep(1:V, nrow(data) %/% V) # generate new var.: fold number
	fold.no <- c(fold.no, (1:V)[nrow(data) %% V])
	fold.chi <- matrix(NA, nrow = V, ncol = length(K))
	fold.adj.chi <- matrix(NA, nrow = V, ncol = length(K))
	fold.WH <- matrix(NA, nrow = V, ncol = length(K))
	fold.t <- matrix(NA, nrow = V, ncol = length(K))
	result <- matrix(NA, nrow = 4, ncol = length(K))
	## CHECK ME: modify by parallel computing
	if(minors@ncl == 1){
		for(i in 1:V){	
			train <- data[fold.no != i,, drop = FALSE]
			test <- data[fold.no == i,, drop = FALSE]
			skaps <- lapply(K, kaps.fit, formula = formula, data = train, mindat= mindat, minors = minors)
			tmp <- sapply(skaps, predict, newdata = test, type = "kaps")
			fold.chi[i,] <- tmp[1,]
			fold.adj.chi[i,] <- tmp[2,]
			fold.WH[i,] <- tmp[3,]
			fold.t[i,] <- tmp[4,]
			rm(skaps, train, test, tmp)
		}
	}
	else{
		ncl <- makeCluster(minors@ncl)
		Export.func <- c("as.Formula", "model.part", "survdiff", "Surv")
		clusterExport(ncl, Export.func)
		for(i in 1:V){	
			train <- data[fold.no != i,, drop = FALSE]
			test <- data[fold.no == i,, drop = FALSE]
			skaps <- parLapply(ncl, K, kaps.fit, formula = formula, data = train, mindat= mindat, minors = minors)
			tmp <- parSapply(ncl, skaps, predict, newdata = test, type = "kaps")
			fold.chi[i,] <- tmp[1,]
			fold.adj.chi[i,] <- tmp[2,]
			fold.WH[i,] <- tmp[3,]
			fold.t[i,] <- tmp[4,]
			rm(skaps, train, test, tmp)
		}
		stopCluster(ncl)
	}
	### output
	result[1,] <- apply(fold.chi,2, mean, na.rm = TRUE)
	result[2,] <- apply(fold.adj.chi,2, mean, na.rm = TRUE)
	result[3,] <- apply(fold.WH,2, mean, na.rm = TRUE)
	result[4,] <- apply(fold.t,2, mean, na.rm = TRUE)
	return(result)
}


############################################################
# Fit KAPS with Permutation test by test estimates approach
# This algorithm select best subgroups by permutation test
# in which the data is randomly divided into 
# learning and validation set
############################################################
kaps.perm <- function(fit, newdata){
	#if(!inherits(fit, "kaps")) stop("This function requires the object of kaps class as a main object.")
	if(missing(newdata)){
		dat.perm = fit@data
		dat.perm$Group <- factor(fit@groupID)
	} else{
		if(!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
		dat.perm <-	predict(fit, newdata)
		dat.perm$Group <- factor(dat.perm$Group)
	}
	res.perm = matrix(1, nrow = 1, ncol = 2)
		
	if(length(levels(dat.perm$Group)) == 1){
		return(res.perm)
	} else{
		Ku = fit@groups
		formula = fit@formula
		minors = fit@Options
		
		formula1 = update(formula, ~ Group)
		formula1 <- as.formula(formula1)

		# Set storages
		perm.over.pval = perm.pair.stat = perm.pair.pval = c()
		res.perm = matrix(0, nrow = 1, ncol = 2)
		#res.perm <- as.data.frame(res.perm)
		#colnames(res.perm) <- c("overall.stat", "overall.pval", "pair.stat", "pair.pval")
		#colnames(res.perm) <- c( "perm_overall_pval", "perm_min_pair_pval")
		######################################################################
		# GOD by permutation test with asymptotic distribution
		# require(coin)
		# cat("Now, permutation kaps is working. Please, wait a minute.^^\n")
		# 1. overall permuatation test statistic
		tmp <- surv_test(formula = formula1, data = dat.perm, distribution = approximate(B = minors@N.perm))
		# Naive approach
		#tmp <- clinfun::permlogrank(formula1, dat.perm)
		#perm.over.stat <- statistic(tmp)
		perm.over.pval <- pvalue(tmp)
		#names(perm.over.stat) <- "overall_statistic"
		#surv_test(formula = Surv(Y,cens) ~ Group, data = dat.perm, distribution = exact())
		
		# 2. pairwise permutation test statistic
		pair = combn(Ku, 2)
		#pair = combat(Ku, 2)
		if(minors@shortcut) pair <- pair[, !duplicated(pair[1,]), drop = FALSE]

		for(i in 1:ncol(pair)){
			dat.perm.tmp = dat.perm[dat.perm$Group %in% pair[,i],]
			#if(length(levels(dat.perm.tmp$Group)) <= 1) {
			#	perm.pair.stat <- c(perm.pair.stat, 100)
			#	perm.pair.pval <- c(perm.pair.pval, 0)
			#} else if(nrow(dat.perm.tmp) == 0){
			#	perm.pair.stat <- c(perm.pair.stat, 100)
			#	perm.pair.pval <- c(perm.pair.pval, 0)
			#} else if(all(table(dat.perm.tmp$Group) >= fit@mindat)){
			#	perm.pair.stat <- c(perm.pair.stat, 100)
			#	perm.pair.pval <- c(perm.pair.pval, 0)
			#}	else {
			#	tmp <- surv_test(formula = formula1, data = dat.perm.tmp, distribution = approximate(B = minors@N.perm))
			#	perm.pair.stat <- c(perm.pair.stat, statistic(tmp))
			#	perm.pair.pval <- c(perm.pair.pval, pvalue(tmp))
			#}
			tmp <- try(surv_test(formula = formula1, data = dat.perm.tmp, distribution = approximate(B = minors@N.perm)), silent = TRUE)
			if(class(tmp) != "try-error"){
				perm.pair.stat <- c(perm.pair.stat, statistic(tmp))
				perm.pair.pval <- c(perm.pair.pval, pvalue(tmp))
			} 
		}
		#names(perm.pair.stat) <- apply(pair, 2, function(x) paste(x, collapse = "&"))

		##########
		# results
		min.ind = which.min(perm.pair.stat)

		#res.perm[1,1] <- perm.over.stat
		res.perm[1,1] <- perm.over.pval
		#res.perm[1,3] <- perm.pair.stat[min.ind]
		res.perm[1,2] <- ifelse(length(min.ind) == 0, 1.00, perm.pair.pval[min.ind])
		return(res.perm)
	}
}
# END @ Oct 20, 2013
