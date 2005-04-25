\name{imageplot.bma}
\alias{imageplot.bma}
\title{Images of models used in Bayesian model averaging }
\description{
  Creates an image of the models selected using \code{\link{bicreg}}, \code{\link{bic.glm}} or \code{\link{bic.surv}}.
}
\usage{
imageplot.bma(bicreg.out,color="default",...) 
}
\arguments{
  \item{bicreg.out}{An object of type 'bicreg', 'bic.glm' or 'bic.surv'}
  \item{color}{The color of the plot. The value "default" uses the current default R color scheme for \code{image}. 
The value "blackandwhite" produces a black and white image.}
  \item{...}{Other parameters to be passed to the \code{image} and \code{axis} functions.}
}
\details{ 
Creates an image of the models selected using \code{\link{bicreg}}, \code{\link{bic.glm}} or \code{\link{bic.surv}}.

  }
\author{Adrian E. Raftery \email{raftery@AT@stat.washington.edu} and Hana Sevcikova}
\references{Clyde, M. (1999) Bayesian Model Averaging and Model Search Strategies (with discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.
}
\seealso{ \code{\link{bicreg}}, \code{\link{bic.glm}}, \code{\link{bic.surv}} }
\examples{
# logistic regression using bic.glm
library("MASS")
data(birthwt)
y<- birthwt$lo
x<- data.frame(birthwt[,-1])
x$race<- as.factor(x$race)
x$ht<- (x$ht>=1)+0
x<- x[,-9]
x$smoke <- as.factor(x$smoke)
x$ptl<- as.factor(x$ptl)
x$ht <- as.factor(x$ht)
x$ui <- as.factor(x$ui)

glm.out1<- bic.glm(x, y, strict = TRUE, OR = 20, glm.family="binomial", factor.type=TRUE)
imageplot.bma(glm.out1)



# logistic regression using glib
library("MASS")
data(birthwt)
y<- birthwt$lo
x<- data.frame(birthwt[,-1])
x$race<- as.factor(x$race)
x$ht<- (x$ht>=1)+0
x<- x[,-9]
x$smoke <- as.factor(x$smoke)
x$ptl<- as.factor(x$ptl)
x$ht <- as.factor(x$ht)
x$ui <- as.factor(x$ui)

glib.birthwt<- glib(x,y, error="binomial", link = "logit")
glm.birthwt<- as.bic.glm(glib.birthwt)
imageplot.bma(glm.birthwt)

}
\keyword{regression}
\keyword{models}

