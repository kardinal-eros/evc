evcLookup <- function () {
	file <- system.file("extdata", "EuroVegChecklistClasses.csv",
                    package = "evc")
	r <- read.csv(file, colClasses = "character")
	return(r)	
}

evcCode <- function (obj) {
	r <- strsplit(taxonomy(obj)$evc.class.code, "+", fixed = TRUE)
	names(r) <- taxonomy(obj)$abbr
	return(r)
}

#	function to deparse evc class code string and cast to species by class matrix	
evcClass <- function (obj, plot) {
	labels <- evcLookup()$evc.class.code
	
	#	subset plot
	xi <- as(obj[ plot, ], "Vegsoup")
	decostand(xi) <- NULL
	
	#	decompose class membership string  
	ri <- evcCode(xi)

	#	matrix of class membership	
	r <- matrix(nrow = length(ri), ncol = length(labels))
	dimnames(r) <- list(names(ri), labels)
	r <- as.data.frame(r)
	r1 <- r
	r1[ ] <- FALSE
	
	for (i in 1:length(ri)) {
		r1[ i, unique(ri[[ i ]]) ] <- TRUE
	}
	
	#	matrix of cover value
	ri <- t(as.numeric(layers(xi, collapse = "ol")))
	#	this should hold true
	#	rownames(ri) <- decode(ri, xi)$abbr; match(rownames(ri), rownames(r1))	
	r2 <- r1 + ri
	r2[ !r1 ] <- 0
		
	#	integer matrix of taxon in classes
	r3 <- r1 * as.integer(1)
	n <- rowSums(r3)
	r3 <- r3 * n
	
	#	weighted cover
	r4 <- r2 / r3
	r4[ t(apply(r4, 1, is.nan)) ] <- 0

	r <- list(weighted.cover = r4, cover = r2, membership = r1)
	
	return(r)
}	

#	function to cast class memberships to matrix
evcMatrix <- function (obj, select, restrict, weighted = TRUE) {
		
	r <- sapply(rownames(obj), function (y) { 
		evcClass(obj, y)
		}, simplify = FALSE)

	if (!weighted) {
		r1 <- sapply(r, "[[", 3, simplify = FALSE) # membership
		r1 <- t(sapply(r1, colSums))
		r1 <- as.data.frame(r1)
		r1 <- r1[ restrict ]
		r1 <- r1 [ select ]
		r <- r1 / rowSums(r1)
	} else {
		r2 <- sapply(r, "[[", 2, simplify = FALSE) # weighted.cover
		r2 <- t(sapply(r2, colSums))
		r2 <- as.data.frame(r2)
		r2 <- r2[ restrict ]
		r2 <- r2 [ select ]		
		r <- r2 / rowSums(r2)
	}
	
	return(r)		
}

#	function to transform evc classes to species and build Vegsoup object
evc2vegsoup <- function (obj, restrict, select) {
	labels <- evcLookup()$evc.class.code
		
	X <- evcMatrix(obj, weighted = TRUE,
			restrict = restrict,
			select = select)
	X <- t(X)
	X <- data.frame(abbr = rownames(X), layer = "ol", X)
	
	X <- stackSpecies(x = X)[, 1:4]
			
	Z <- labels[ match(unique(X$abbr), labels), ]
	names(Z) <- c("abbr", "taxon")
	
	Y <- sites(obj)
	Y$plot <- rownames(obj)
	Y$latitude <- coordinates(obj)[, 2]
	Y$longitude <- coordinates(obj)[, 1]
	Y <- stackSites(Y)
	
	XZ <- SpeciesTaxonomy(X, Z)
	
	r <- Vegsoup(XZ, Y, coverscale = "percentage")
}