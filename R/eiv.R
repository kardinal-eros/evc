#	clean string
eivClean <- function (x) {
	x[ c(which(x == "x"), which(x == "-"), which(x == "?"),  which(x == "n.a.")) ] <- NA
	r <- as.numeric(x)
	return(r)
}

#	calculate median
eivThreshold <- function (obj, plot, summary = TRUE) {
	xi <- obj[ plot, ]
	xx <- data.frame(t(as.numeric(xi)))
	names(xx) <- "cov"
	xx$cov.scaled <- xx$cov / sum(xx$cov) * 100
	xx$abbr <- decode(xx, xi)$abbr
	
	zz <- taxonomy(taxonomy(xi))
	zz <- zz[, c(grep("abbr", names(zz)), grep("taxon", names(zz)), grep("eiv", names(zz))) ]
	
	r <- cbind(xx, zz[ match(xx$abbr, zz$abbr), ])
	
	#	apply treshold
	r$threshold <- r$eiv.n <= 4
	r$threshold[ c(which(r$eiv.n == "x"), which(r$eiv.n == "-"), which(r$eiv.n == "?")) ] <- NA
		
	if (summary) {
		r <- data.frame(
			plot = plot,
			cov.threshold = sum(r[ which(r$threshold == TRUE), ]$cov),
			cov.scaled.threshold = round(sum(r[ which(r$threshold == TRUE), ]$cov.scaled), 2),
			cov.total = sum(r$cov),
			eiv.n.median = round(moeller(eivClean(r$eiv.n)), 2),
			eiv.f.median = round(moeller(eivClean(r$eiv.f)), 2),
			eiv.r.median = round(moeller(eivClean(r$eiv.r)), 2),
			eiv.t.median = round(moeller(eivClean(r$eiv.t)), 2)						
		)
	}
	return(r)
}

#	plot eiv
#plot.eiv <- function (obj, notch = TRUE) {
#	op <- par(mfrow = c(3,2)); on.exit(par(op))
#	hist(obj$richness[ obj$type == "pasture" ], main = "pasture", xlab = "richness", xlim = c(0, 80))
#	
#	if (any(names(obj) != "richness")) {
#		message("calculate richness")
#		obj$richness <- richness(obj, "sample")
#	}
#	scatter.smooth(obj$cov.scaled.threshold, obj$richness,
#		xlab = "N cover treshold (scaled)", ylab = "richness")
#	abline(v = 50, lty = "dotted")
#	abline(v = 75, lty = "dotted")
#	
#	plot(obj$eiv.n.median ~ obj$cov.scaled.threshold,
#		ylab = "N median", xlab = "N cover treshold (scaled)")
#	fit <- lm(obj$eiv.n.median ~ obj$cov.scaled.threshold)
#	abline(fit, col = 2)
#	abline(h = 4, lty = "dotted")
#	abline(v = 50, lty = "dotted")
#	abline(v = 75, lty = "dotted")
#
#	boxplot(obj$richness ~ obj$threshold.scaled, notch = notch,
#		xlab = "Threshold > 50", ylab = "richness")	
#	
#	boxplot(obj$richness ~ obj$threshold2.scaled, notch = notch,
#		xlab = "Threshold > 75", ylab = "richness")
#}