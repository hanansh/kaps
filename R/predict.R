### predict functions for kaps class
setGeneric("predict")
setMethod("predict","kaps", function(object, newdata, type = c("predict", "kaps")){
	## predict the ID for terminal subgroups 
	type <- match.arg(type)
	if(missing(newdata)) {
		where <- object@groupID
		newdata <- object@data
	}
	else {
		where <- pred.kaps(object@split.pt, object@formula, newdata)
	}

	if( type == "predict"){
		#tmps <- range(object@data[,object@split.var])
		#match(, )
		#pts <- c(tmps[1], round(object@split.pt,2), tmps[2])
		result <- data.frame(newdata, Group = where)
		#rownames(result) <- "The estimated group="
		#result$Range <- paste("(", "<X<=",")", sep = "")
		#colnames(result) <- c("Newdata", "Group", "Range")
		#colnames(result) <- c("Group", colnames(newdata))
		return(result)
	}
	else if(type == "kaps"){
		newdata$subgroups <- where
		f <- update(object@formula, . ~  subgroups)
		test.chi <- survdiff(f, data = newdata)$chisq
		test.adj.chi <- pairwise.logrank.test(where, data = newdata, formula = object@formula, rho = object@Options@rho, adj = object@Options@p.adjust.methods, splits = object@Options@splits, shortcut = object@Options@shortcut)[1]
		v <- length(object@split.pt)
		x.mean <- 1 - 2 / (9 * v)
		x.std <- sqrt(2 / (9 * v))
		WH <- (test.chi / v)^(1/3)
		t <-abs((WH - x.mean)/(x.std))
		return(pred.stat = c(test.chi, test.adj.chi, WH, t))
	}
	}
)

pred.kaps <- function(split.pt, f, newdata){
#### find ID number for new data set
##Input
# split points by train data
# new covariates data with the type of data.frame
## output
# ID number (the observations that assins each terminal group)
## FIXME by more efficient way in the near future
	X <- model.part(f, data = newdata, rhs = 1, drop = FALSE)
	nc <- length(split.pt)
	gClass <- matrix(NA, ncol = nc, nrow = nrow(newdata))
	gClass <- sapply(split.pt, function(x,y) y > x, y = X)
	if(is.vector(gClass)) gClass <- t(gClass)
	where <- apply(gClass, 1, sum)
	where <- where + 1
	return(where)
}
