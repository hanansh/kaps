kaps.control <- function(pre.pt = list(), scope = list(), 
	rho = 0, ncl = 1, lower.limit = 0, upper.limit = 100,
	shortcut = TRUE, N.perm = 1999, N.boot = 200, alpha = 0.05, 
	splits = c("logrank", "exact"),
	sel = c("mean","median","trim","test"),
	p.adjust.methods = c("none", "holm", "hochberg", "hommel", "bonferroni",  "BH", "BY", "fdr")){
	p.adjust.methods <- match.arg(p.adjust.methods)
	splits <- match.arg(splits)
	sel <- match.arg(sel)
	res <- new("kapsOptions")
	#res@method <- method
	res@lower.limit <- lower.limit
	res@upper.limit <- upper.limit
	res@N.perm <- N.perm
	res@N.boot <- N.boot
	res@alpha <- alpha
	res@pre.pt <- pre.pt
	res@scope <- scope
	res@rho <- rho
	res@shortcut <- shortcut
	res@splits <- splits
	res@boot.sel <- sel
	res@ncl <- as.integer(ncl)
	res@p.adjust.methods <- p.adjust.methods
	return(res)
}

