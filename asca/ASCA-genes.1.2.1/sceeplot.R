screeplot <- function (asca) {
num <- length(asca)
r <- round(sqrt(num))
c <- num-r
n <- r*c
layout(matrix(c(1:num), 1, num, byrow = TRUE))

   for ( i in 1: length (asca)) {

      plot(asca[[i]]$var.exp[,1], type="l", main=paste("Screeplot", names(asca)[[i]]), xlab="Component", ylab="Expalined variability")

   }

}
