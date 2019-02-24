SPE.lims <- function (my.asca, alpha)
{
limits <- vector("list", length(my.asca))
names(limits) <- names(my.asca)
	for (i in 1:length(my.asca)) {
		assign ("model", my.asca[[i]])
		if (!is.null(model$SPE)) {
			m <- mean(model$SPE)
			v <- var(model$SPE)
			g <- v/(2*m)
			h <- 2*m*m/v
			limits[[i]] <- g*qchisq(1-alpha, df=h)
		}
	}
limits
}

