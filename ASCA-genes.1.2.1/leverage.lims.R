leverage.lims <- function(data = data, R = 100, FUN, Designa = Designa, Designb = Designb, Designc = NULL, Fac = c(1,2,2,2), type = 2, alpha = 0.01, showvar=FALSE, showscree=FALSE)
{
## Compute ASCA model for data
	Model <- FUN(X = t(data), Designa = Designa,Designb = Designb,Designc = Designc, Fac = Fac,type = type)
	n <- ncol(data)
	lim <- Selection <- vector(mode = "list", length = length(Model)-1)
	names(lim) <- names(Selection) <- names(Model)[1:length(lim)]

## Calculate the reference distribution of leverages
	for (i in 1:R) {
		place <- sample(1:n)
		permu <- data[,place]
		A <- FUN (X = t(permu), Designa = Designa, Designb = Designb, Designc = Designc, Fac = Fac, type = type, showvar=showvar, showscree=showscree)
		for (j in 1:length(lim)) {
			lim[[j]] <- cbind(lim[[j]],A[[j]]$leverage)
		}
	}
		rm(A)
		gc()
## Find alpha quantile value for reference distribution
	QC <- apply(sapply(lim,function(x) {apply(x,1,quantile,probs=1 - alpha, na.rm=T)}),2,quantile,probs=1 - alpha, na.rm=T)

# Gene selection
	for (h in 1:length(lim)) {
		Selection[[h]] <- rownames(data)[Model[[h]]$leverage > QC[[h]]]
	}
	output <- list(lim, QC, Selection)
	names(output) <- c("NullDistribution", "Cutoff", "Selection")
	output
}