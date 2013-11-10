kapsNews <- function(){
	file.locate <- file.path(system.file(package = "kaps"), "NEWS")
	file.show(file.locate)
}

## Count minimum sample size
count.mindat <- function(formula, data, part = 10){
## output: the minimum number of samle size.

	formula <- as.Formula(formula)
	X <- model.part(formula, data, rhs = 1)
	subgr <- apply(X, 2, function(x) {
		res <- table(x) / length(x)
		res <- res[order(res, decreasing = FALSE)]
		return(res)
		}
	    )
	res <- apply(subgr, 2, function(x, data, part){
		dat.prop <- x[length(x) - part]
		res <- floor(nrow(data) * dat.prop)
		return(res)}, 
		data =  data, part = part
		)
	names(res) <- NULL
	return(res)
}

## Calculate all possible subset
# from gRbase package
combnat <- function(x,m){
  if (length(x)==1 && is.numeric(x))
    x <- seq(x)
  if (length(x) < m)
    stop("Error in combnat: n < m\n")
  NCAND = as.integer(length(x))
  NSEL  = as.integer(m)
  NSET <- as.integer(choose(NCAND,NSEL))
  ANS  <- as.integer(rep.int(0, NSET*NSEL))
  res <- .C("combnC", NSEL, NCAND, NSET, ANS, DUP=FALSE
            ,PACKAGE="kaps" )[[4]]
  res <- x[res]
  dim(res) <- c(NSEL, NSET)  
  return(res)
}

## Calculate Chi-square and Adj. chi-square test statistics among subgroups.
adj.test <- function(x){
	data <- x@data
	data$subgroups <- x@groupID
	#data$group <- gr.sel[,index]
	f <- update(x@formula, . ~ subgroups)		
	test.stat <- survdiff(f, data = data)$chisq
	v <- length(x@split.pt)
	x.mean <- 1 - (2 / (9 * v))
	x.std <- sqrt(2 / (9 * v))
	WH <- (test.stat / v)^(1/3)
	t <-(WH-x.mean)/(x.std)
	#WH <- ( (7/9) + sqrt(nu)*((test.stat/ nu )^(1/3) - 1 + (2 / (9 * nu))))^3
	#WH <- max(0, WH)
	### ouptut structure
	# over.stat = overall test statistic for selected candidate
	# pair.stat = pairwise test statistic for selected candidate
	# cube.stat = cube-root transformation
	# WH.stat = t statistic
	return(c(over.stat = test.stat, pair.stat = x@X, cube.stat = WH, WH.stat = t))
}

## Median Survival Time, refered to the print.survfit() in the survival package
survmed <- function(surv,time, tol= 1.0e-5) {
	keep <- (!is.na(surv) & surv <(.5 + tol))
	if (!any(keep)) NA
	else {
		time <- time[keep]
		surv <- surv[keep]
		if (abs(surv[1]-.5) <tol  && any(surv< surv[1])) 
			(time[1] + time[min(which(surv<surv[1]))])/2
		else time[1]
	}
}

## summary functions for adaptive staging algorithms
surv.yrs <- function(pt, surv, time){
	if(any(time == pt)) return(surv[time == pt])
	else {
		mod <- min(abs(time - pt))
		if(mod > 7) {
			return(0)
		}
		if(any(time == (pt + mod))) return(surv[time ==(pt + mod)])
		else return(surv[time ==(pt - mod)])
	}
}

## summary functions for kaps class
setGeneric("summary")
setMethod("summary","kaps", function(object,K){
	if(!missing(K)) object <- object@results[[which(object@groups == K)]]

	data <- object@data
	f <- update(object@formula, . ~ 1)
	surv.root <- survfit(f, data = data)
	rootS <- summary(surv.root)
	
	data$Group <- object@groupID 
	f <- update(object@formula, . ~ Group)
	surv.all <- survfit(f, data = data)
	objS <- summary(surv.all)
	 
	subgr.surv <- list()
	subgr.time <- list()
	level.objS <- levels(objS$strata)
	for(i in 1:length(level.objS)) {
		subgr.surv[[i]] <- objS$surv[objS$strata == level.objS[i]] 
		subgr.time[[i]] <- objS$time[objS$strata == level.objS[i]] 
	}
	gr.med <- objS$table[,"median"]
	root.med <-  rootS$table["median"]
	root.1 <- surv.yrs(pt = 12, surv = rootS$surv, time = rootS$time)
	root.3 <- surv.yrs(pt = 36, surv = rootS$surv, time = rootS$time)
	root.5 <- surv.yrs(pt = 60, surv = rootS$surv, time = rootS$time)
	root <- round(c(nrow(data), root.med, root.1, root.3, root.5), 3)
	names(root) <- NULL
	gr.1yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 12))
	gr.3yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 36))
	gr.5yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 60))
	obs <- table(object@groupID)
	obs <- as.vector(obs)
	names(obs) <- names(gr.med)
	res <- data.frame( 
		N = obs,
		Med = gr.med,
		yrs.1 = round(gr.1yrs,3),
		yrs.3 = round(gr.3yrs,3),
		yrs.5 = round(gr.5yrs,3)
	)
	rownames(res) <- names(obs)
	res <- res[order(res$Med, decreasing = TRUE),]
	res <- rbind(All = root, res)
	return(res)
	}
)

####################################################################
## print furnctions for an kaps S4 class
setGeneric("print")
setMethod("print","kaps", function(x,K){
	call <- match.call()
    if(!missing(K)) {
      index <- which(x@groups == K)
      Z <- x@Z[index]
      object <- x@results[[index]]
      object@results[[1]] <- x
      object@Z <- Z
      object@Options@fold <- FALSE
	  object@call <- call
    }
	else object <- x
	
	pts <- round(object@split.pt,2)
	cat("Selecting a set of cut-off points:\n")
	spt <- lapply(object@results, function(x) x@split.pt)
	spt <- sapply(spt, function(x) paste(x, collapse = ", "))
	Z.stat <- zapsmall(object@Z, digits = 3) # overall test statistic
	Z.stat.df <- (object@groups)-1
	Z.stat.p <- round(1 - pchisq(q = Z.stat, df= Z.stat.df),4) #overall p-value
	Signif <- symnum(Z.stat.p, corr = FALSE, na = FALSE, 
			  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1), 
			  symbols = c("***", "**", "*", "."," "))

	X.stat <- zapsmall(object@X,digits = 3) # optimal pair statistic
	test.stat <- data.frame(
		X.stat,
		df2 = 1,
		pvalue = round(object@pvalue,4),
		Z.stat,
		df1 = Z.stat.df,
		pvalue1 = Z.stat.p,
		pts = spt,
		sig = format(Signif))
	dimnames(test.stat) <- list(paste("K=",object@groups,sep=""), c("X", "df", "Pr(>|X|)", "Xk", "df", "Pr(>|k|)",  "cut-off points", ""))
	print(test.stat)
	
	cat("\nP-values of pairwise comparisons\n")
	data <- object@data
	data$groupID <- object@groupID
	pair <- combn(unique(object@groupID),2)
	f <- update(object@formula, . ~ groupID)
	res <- c()
	for(i in 1:ncol(pair)){
		data.tmp <- data[data$groupID %in% pair[,i],]
		tmp <- survdiff(f, data = data.tmp, rho = object@Options@rho)
		res[i] <- 1 - pchisq(tmp$chisq, 1)
	}
	# Adjsut P-values for Multiple Comparisons 
	res <- p.adjust(res, method = object@Options@p.adjust.methods)
	pair.mat <- matrix(NA, ncol = length(object@split.pt), nrow = length(pts))
	pair.mat[!upper.tri(pair.mat)] <- res
	pair.mat <- round(pair.mat, 4)
	## adjust split points
	tmps <- range(data[,object@split.var])
	pts <- c(tmps[1], pts, tmps[2])
	n.pts <- length(pts)
	colnames(pair.mat) <- paste(pts[-(n.pts:(n.pts-1))],"<",object@split.var,"<=", sep = "")
	colnames(pair.mat) <- paste(colnames(pair.mat),pts[-c(1,length(pts))], sep = "")
	colnames(pair.mat)[1] <- paste(pts[1],"<=",object@split.var,"<=",pts[2], sep = "")
	rownames(pair.mat) <- paste(pts[-c(1,length(pts))],"<",object@split.var,"<=", sep = "")
	rownames(pair.mat) <- paste(rownames(pair.mat), pts[-(1:2)],sep="")
	pair.mat <- ifelse(pair.mat < 1.0e-5, "<.0000", pair.mat)
	pair.mat <- ifelse(!is.na(pair.mat), pair.mat, "-")
	pair.mat <- as.data.frame(pair.mat)
	print(pair.mat)
	#cat("\nP-value adjustment method:" ,object@Options@p.adjust.methods, "\n")
	invisible(object)
	}
)

setGeneric("show")
setMethod("show", "kaps", function(object){
	cat("Call:\n")
	print(object@call)
	cat("\n")
	cat("	K-Adaptive Partitioning for Survival Data\n\n")
	if (object@Options@fold & length(object@groups) >= 2) cat("Samples=", nrow(object@data), "				Optimal K=",object@groups[object@index], "\n")
	else cat("Samples=", nrow(object@data), "\n")
	
	pts <- round(object@split.pt,2)
	
	cat("\n")
	cat("\nSelecting a set of cut-off points:\n")
	spt <- lapply(object@results, function(x) x@split.pt)
	spt <- sapply(spt, function(x) paste(x, collapse = ", "))
	Z.stat <- zapsmall(object@Z, digits = 3) # overall test statistic
	Z.stat.df <- (object@groups)-1
	Z.stat.p <- round(1 - pchisq(q = Z.stat, df= Z.stat.df),4) #overall p-value
	Signif <- symnum(Z.stat.p, corr = FALSE, na = FALSE, 
			  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1), 
			  symbols = c("***", "**", "*", "."," "))

	X.stat <- zapsmall(object@X,digits = 3) # optimal pair statistic
	test.stat <- data.frame(
		X.stat,
		df2 = 1,
		pvalue = round(object@pvalue,4),
		Z.stat,
		df1 = Z.stat.df,
		pvalue1 = Z.stat.p,
		pts = spt,
		sig = format(Signif))
	dimnames(test.stat) <- list(attr(object@groups,"names"), c("X", "df", "Pr(>|X|)", "Xk", "df", "Pr(>|k|)",  "cut-off points", ""))
	print(test.stat)
	
	if(object@Options@fold & length(object@groups) >= 2){
		cat("\nFinding an optimal K with cross-validation:\n")
		cv.pvalue <- round(1 - pnorm(q = object@elbow[4,]), 4)
		cv.xvalue <- round(1 - pchisq(q = object@elbow[2,],df=1), 4)
		cv.zvalue <- round(1 - pchisq(q = object@elbow[1,],df=Z.stat.df), 4)
		Signif <- symnum( cv.pvalue, corr = FALSE, na = FALSE, 
		  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
		  symbols = c("***", "**", "*", ".", " "))
		cv.stat <- data.frame(
			X = object@elbow[2,],
			df1 = 1,
			Z = object@elbow[1,],
			df2 = Z.stat.df,
			t = object@elbow[4,])
		cv.stat<- zapsmall(cv.stat,digits = 3)
		cv.stat <- cbind(
			cv.stat[,1:2], 
			xvalue = cv.xvalue,
			cv.stat[,3:4],
			zvalue = cv.zvalue,
			cv.stat[,5],
			pvalue = cv.pvalue,
			sig = format(Signif))
		dimnames(cv.stat) <- list(attr(object@groups,"names"), c("X","df","Pr(>|X|)","Xk", "df", "Pr(>|k|)",  "W","Pr(>|W|)", ""))
		print(cv.stat)
		#cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
	}
	cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")

	cat("\nP-values of pairwise comparisons\n")
	data <- object@data
	data$groupID <- object@groupID
	pair <- combn(unique(object@groupID),2)
	f <- update(object@formula, . ~ groupID)
	res <- c()
	for(i in 1:ncol(pair)){
		data.tmp <- data[data$groupID %in% pair[,i],]
		tmp <- survdiff(f, data = data.tmp, rho = object@Options@rho)
		res[i] <- 1 - pchisq(tmp$chisq, 1)
	}
	# Adjsut P-values for Multiple Comparisons 
	res <- p.adjust(res, method = object@Options@p.adjust.methods)
	pair.mat <- matrix(NA, ncol = length(object@split.pt), nrow = length(pts))
	pair.mat[!upper.tri(pair.mat)] <- res
	pair.mat <- round(pair.mat, 4)
	## adjust split points
	tmps <- range(data[,object@split.var])
	pts <- c(tmps[1], pts, tmps[2])
	n.pts <- length(pts)
	colnames(pair.mat) <- paste(pts[-(n.pts:(n.pts-1))],"<",object@split.var,"<=", sep = "")
	colnames(pair.mat) <- paste(colnames(pair.mat),pts[-c(1,length(pts))], sep = "")
	colnames(pair.mat)[1] <- paste(pts[1],"<=",object@split.var,"<=",pts[2], sep = "")
	rownames(pair.mat) <- paste(pts[-c(1,length(pts))],"<",object@split.var,"<=", sep = "")
	rownames(pair.mat) <- paste(rownames(pair.mat), pts[-(1:2)],sep="")
	pair.mat <- ifelse(pair.mat < 1.0e-5, "<.0000", pair.mat)
	pair.mat <- ifelse(!is.na(pair.mat), pair.mat, "-")
	pair.mat <- as.data.frame(pair.mat)
	print(pair.mat)
	#cat("\nP-value adjustment method:" ,object@Options@p.adjust.methods, "\n")
	invisible(object)
	}
)

####################################################################
### panel function for Kaplan-Meier curves in nodes
## A standard tree plot with KM curves at the terminal nodes
## Soo-Heang Eo
setGeneric("plot")
setMethod("plot","kaps",function(x, K, ...){
		if(!missing(K) | !x@Options@fold) {
			if (!missing(K)) x <- x@results[[which(x@groups == K)]]
			else {
				par(mfrow = c(1,2))

				fvars <- all.vars(x@formula)
				plot(x@data[,fvars[3]], x@data[,fvars[1]], pch = c(1,3)[x@data[,fvars[2]]+1] , 
				axes=FALSE, col= ifelse(x@data[,fvars[2]] == 1, "red2","blue"), xlab="", ylab="")
				mtext(fvars[3], side=1, line=3, cex=1.2)
				mtext("Survival months", side=2, line=2, cex=1.2)
				axis(1)
				axis(2, at= c(12,24,48,72,96,120, 144, 168, 192, 216, 240), labels=c(12,24,48,72,96,120, 144, 168, 192, 216, 240))
				legend("topright",c("Event","Censored"), col=c("blue","red2"), cex=1, pch=c(1,3), bty = "n")
				r <- locfit.censor(x = x@data[,fvars[3]], y = x@data[,fvars[1]], cens= (1- x@data[,fvars[2]])) 
				lines(r, lwd=2, lty=1)
			}
			km.curve(x)
			legend("topright", legend = paste("G", 1:(length(x@split.pt)+1),sep = ""), 
				bty = "n", col = unique(x@groupID), lty = unique(x@groupID))
			return(invisible(x))
		}
		else K <- x@groups[x@index]
		
		if( (length(x@groups) == 1) ){
			km.curve(x, main = paste("# of groups = ",K))
			legend("topright", legend = paste("G", 1:(length(x@split.pt)+1),sep = ""), 
				bty = "n", col = unique(x@groupID), lty = unique(x@groupID))
			return(invisible(x))
		}
	
		fvars <- all.vars(x@formula)
		par(mfrow = c(2,2))
		
		plot(x@data[,fvars[3]], x@data[,fvars[1]], pch = c(1,3)[x@data[,fvars[2]]+1] , 
			axes=FALSE, col= ifelse(x@data[,fvars[2]] == 1, "red2","blue"), xlab="", ylab="")
		mtext(fvars[3], side=1, line=3, cex=1.2)
		mtext("Survival months", side=2, line=2, cex=1.2)
		axis(1)
		axis(2, at= c(12,24,48,72,96,120, 144, 168, 192, 216, 240), labels=c(12,24,48,72,96,120, 144, 168, 192, 216, 240))
		legend("topright",c("Event","Censored"), col=c("blue","red2"), cex=1, pch=c(1,3), bty = "n")
		r <- locfit.censor(x = x@data[,fvars[3]], y = x@data[,fvars[1]], cens= (1-x@data[,fvars[2]])) 
		lines(r, lwd=2, lty=1)
		
		km.curve(x)
		legend("topright", legend = paste("G", 1:(length(x@split.pt)+1),sep = ""), 
			bty = "n", col = unique(x@groupID), lty = unique(x@groupID))

		plot(x@elbow[2,], type = "b", col = "blue", xlab = "", ylab = "",axes = FALSE, ylim = c(0,max(x@elbow[2,])),...)
		mtext("K", side=1, line=3, cex=1.2)
		mtext(expression(X), side=2, line=2, cex=1.2)
		abline(h = 3.84, col = "gray", lty = 2)
		axis(side = 1, at = 1:ncol(x@elbow), labels =  x@groups)
		axis(side = 2)

		plot(x@elbow[4,], type = "b", col = colors()[630], xlab = "", ylab = "", ylim = c(0,ceiling(max(x@elbow[4,]))), axes = FALSE,...)
		mtext("K", side=1, line=3, cex=1.2)
		mtext(expression(W), side=2, line=2, cex=1.2)
		abline(h = 1.64, col = "gray", lty = 2)
		axis(side = 1, at = 1:ncol(x@elbow), labels =  x@groups)
		axis(side = 2)
	}
)

## plot Kaplan-Meire survival curves for termninal nodes
km.curve <- function(object, 
	x.lab = c(0,24,48,72,96,120, 144, 168, 192, 216, 240), lwd = 1.5, ...){
	if(class(object) == "kaps"){
		object@data$groupID <- object@groupID 
		id.n <- length(unique(object@groupID))
	}
	else if(class(object) == "Tree"){
		object@data$groupID <- object@where 
		id.n <- length(unique(object@where))
	}
	f <- update(object@formula, . ~ groupID)
	surv.all <- survfit(f, data = object@data)
	plot(surv.all, col= 1:id.n, lty = 1:id.n, axes=FALSE, cex=1, lwd= lwd, ...)
	mtext("Survival months", side=1, line=3, cex=1.2)
	mtext("Survival probability", side=2, line=3, cex=1.2)
	axis(2)
	axis(1, at= x.lab, labels=x.lab)
}
###################################3
## update functions for kaps class
setGeneric("update")
setMethod("update","kaps", function(object, type, sel, alpha = 0.05){

	# Data manipulation
	call <- match.call()
	K = object@groups

	# sample statistic
	Boot.over.stat = object@over.stat.sample
	Boot.pair.stat = object@pair.stat.sample
	over.pair.stat = matrix(0, nrow = nrow(object@elbow), ncol = 4)

	if(sel == "mean"){
		# mean
		over.pair.stat[,1] <- apply(Boot.over.stat, 2, mean, na.rm = TRUE)
		over.pair.stat[,3] <- apply(Boot.pair.stat, 2, mean, na.rm = TRUE)
		#result[,2] <- apply(Boot.over.stat, 2, sd, na.rm = TRUE) / sqrt(B)
		#result[,4] <- apply(Boot.pair.stat, 2, sd, na.rm = TRUE) / sqrt(B)
	} else if(sel == "median"){
		# median
		over.pair.stat[,1] <- apply(Boot.over.stat, 2, median, na.rm = TRUE)
		over.pair.stat[,3] <- apply(Boot.pair.stat, 2, median, na.rm = TRUE)
		#result[,2] <- apply(Boot.over.stat, 2, sd, na.rm = TRUE) / sqrt(B)
		#result[,4] <- apply(Boot.pair.stat, 2, sd, na.rm = TRUE) / sqrt(B)
	} else if(sel == "trim"){
		# censoring-related trimmed mean
		CR = table(data[,all.vars(formula)[2]])[1]
		over.pair.stat[,1] <- apply(Boot.over.stat, 2, mean, na.rm = TRUE, trim = CR)
		over.pair.stat[,3] <- apply(Boot.pair.stat, 2, mean, na.rm = TRUE, trim = CR)
		#result[,2] <- apply(Boot.over.stat, 2, sd, na.rm = TRUE) / sqrt(B)
		#result[,4] <- apply(Boot.pair.stat, 2, sd, na.rm = TRUE) / sqrt(B)
	} else if(sel == "test"){
		# test of proportion
		# this algorithm is used for testing the null that the proportion are the same
		# H0: p >= 0.05
		# H1: p < 0.05
		over.sum = pair.sum = c()
		for(i in 1:ncol(Boot.over.stat)){
			over.sum[i] <- sum(Boot.over.stat[,i] > alpha)
			test.tmp <- prop.test(over.sum[i], n = B, alternative = "less")
			over.pair.stat[i,1] <- test.tmp$p.value
			over.pair.stat[i,2] <- test.tmp$statistic

			pair.sum[i] <- sum(Boot.pair.stat[,i] > alpha)
			test.tmp <- prop.test(pair.sum[i], n = B, alternative = "less")
			over.pair.stat[i,3] <- test.tmp$p.value
			over.pair.stat[i,4] <- test.tmp$statistic
		}
	}

	# Selection of optimal number of subgroups
	index = 1:nrow(over.pair.stat)
	index.pair = over.pair.stat[,3]
	if(all(index.pair == FALSE)){
		fit = lapply(K = index, kaps.fit, formula = object@formula, data = object@data, mindat= object@mindat, minors = object@Options)
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
		return(result)
	} else{
		index <- index[sum(index.pair)]
	}

	######################################################
	# Obtain log-rank statistics at K
	### parallel computing in order to find optimal k subgroups
	#cat("Now, selecting a set of cut-off points...\n")
	#fit = kaps.fit(formula = formula, data = data, K = K.sel, mindat = mindat, minors = minors)
	#fit = lapply(K, kaps.fit, formula = formula, data = data, mindat= mindat, minors = minors)
	#test.stat2 = sapply(fit, adj.test)
	#test.stat2 <- as.matrix(test.stat2)

	result <- object@results[[index]]
	result@index <- as.integer(index)
	result@groups <- K
	attr(result@groups,"names") <- paste("K=",K,sep="")
	result@Z <- object@Z
	result@X <- object@X
	result@WH <- object@WH
	result@t <- object@t
	result@results <- object@results
	result@pvalue <- c(object@Z[index], object@WH[index])
	if(type != "NULL") {
		result@elbow <- object@elbow
		result@over.stat.sample <- Boot.over.stat
		result@pair.stat.sample <- Boot.pair.stat
	}
	result@call <- call
	result@Options <- minors
	return(result)
	}
)
