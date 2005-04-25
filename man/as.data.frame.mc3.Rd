\name{as.data.frame.mc3}
\alias{as.data.frame.mc3}
\alias{[.mc3}
\title{ Type conversion and indexing for mc3 object }
\description{
  Allows for converting an \code{mc3} object to a \code{.data.frame} and indexing on \code{mc3} objects
}
\usage{
as.data.frame.mc3(x, ...)

x[i,...]
x[i,j, ...]


}
\arguments{
  \item{x}{object of type \code{mc3}}
  \item{i,j}{elements to extract from an \code{mc3} object}
  \item{...}{other parameters to be passed to \code{[.data.frame}}
}
\details{
 \code{as.data.frame} converts an object of class "mc3" to a dataframe. 

 Indexing on an object of class "mc3" is achieved by first converting the object to a dataframe.
}
\value{
  an object of class .data.frame
}
\references{ }
\author{Ian Painter}
\note{ }
\seealso{ \code{\link{MC3.REG}}, \code{\link{[.data.frame}}}
\examples{

}
\keyword{ methods }

