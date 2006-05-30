\name{summary.iBMA}
\alias{summary.iBMA.glm}
\alias{summary.iBMA.bicreg}
\alias{summary.iBMA.surv}
\alias{print.iBMA.glm}
\alias{print.iBMA.bicreg}
\alias{print.iBMA.surv}
\alias{summary.iBMA.intermediate.glm}
\alias{summary.iBMA.intermediate.bicreg}
\alias{summary.iBMA.intermediate.surv}
\alias{print.iBMA.intermediate.glm}
\alias{print.iBMA.intermediate.bicreg}
\alias{print.iBMA.intermediate.surv}\
\title{Summaries of iterated Bayesian model averaging objects}
\description{
     \code{summary} and \code{print} methods for iterated Bayesian model averaging objects.
}
\usage{
\method{summary}{iBMA.glm}(object, ...)
\method{summary}{iBMA.bicreg}(object, ...)
\method{summary}{iBMA.surv}(object, ...)
\method{print}{iBMA.glm}(x, ...)
\method{print}{iBMA.bicreg}(x, ...)
\method{print}{iBMA.surv}(x, ...)
\method{summary}{iBMA.intermediate.glm}(object, ...)
\method{summary}{iBMA.intermediate.bicreg}(object, ...)
\method{summary}{iBMA.intermediate.surv}(object, ...)
\method{print}{iBMA.intermediate.glm}(x, ...)
\method{print}{iBMA.intermediate.bicreg}(x, ...)
\method{print}{iBMA.intermediate.surv}(x, ...)
}

\arguments{
  \item{object}{object of type iBMA.glm, iBMA.bicreg or iBMA.surv, iBMA.intermediate.glm, iBMA.intermediate.bicreg or iBMA.intermediate.surv. }
  \item{x}{object of type iBMA.glm, iBMA.bicreg or iBMA.surv, iBMA.intermediate.glm, iBMA.intermediate.bicreg or iBMA.intermediate.surv. }
  \item{...}{other parameters to be passed to \code{print.bic.lmg}, \code{print.bicreg} or \code{print.bic.surv} }
}
\details{
  These methods provide concise and summarized information about the variables that have been examined up to the last iteration. 
  If the result is a final result then the methods also display the results of calling print or summary on the Bayesian model average object for the final set of variables.
}
\author{Ian Painter \email{ian.painter@AT@gmail.com}}
\note{The summary function does not create a summary object (unlike \code{summary.lm} or \code{summary.glm}), instead it directly prints the summary. Note that no calculations are done to create the summary.}
\examples{

}
\keyword{print}

