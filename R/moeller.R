#	calculate median according to Möller 1992
moeller <- function (x) {
	x <- x[ !is.na(x) ]
	x_M <- median(x)
	U <- x_M - 0.5
	b <- 1
	n <- length(x[ !is.na(x) ])
	B_u <- length(x[ x < x_M ])
	n_x <- length(x[ x == x_M ])	
	r <- U + b * ( (n/2-B_u) / ifelse(n_x != 0, n_x, 1))
	if (is.na(r)) r <- NA # too many NA species values to perform median calculation
	return(r)
}

#	function to round values preserving a sum of 100%
.round.preserve.sum <- function(x, digits = 0) {
	u <- 10^digits
	x <- x * u
	y <- floor(x)
	i <- tail(order(x - y), round(sum(x)) - sum(y))
	y[ i ] <- y[ i ] + 1
	r <- y / u
	return(r)
}