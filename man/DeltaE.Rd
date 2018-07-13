\name{DeltaE}
\alias{DeltaE}

\title{Calculate the Color Difference between Two Colors}

\description{Calculate the 1976 Standard CIE Color Difference between two Colors.
}

\usage{DeltaE( Lab0, Lab1 ) }

\arguments{
\item{Lab0}{a numeric Nx3 matrix with Lab values in the rows, or a vector that can be converted to such a matrix, by row.
\code{Lab0} can also be a numeric 3-vector with a single Lab, and it is then replicated
to match the size of \code{Lab1}.
}

\item{Lab1}{a numeric Nx3 matrix with Lab values in the rows, or a vector that can be converted to such a matrix, by row.
\code{Lab1} can also be a numeric 3-vector with a single Lab, and it is then replicated
to match the size of \code{Lab0}.
}

}

\value{
A numeric N-vector with the pairwise differences, 
i.e. between row i of \code{Lab0} and  row i of \code{Lab1}.
The difference is simply the Euclidean distance between the two points in Lab.
}



\section{Future Work}{
The current DeltaE was approved by CIE in 1976.
Since then there have been many improved formulas.
}



\examples{
DeltaE( c(50,0,0),  c(51,2,2,  52,10,11,  46,-13,16)  )
## [1]  3 15 21
}
\keyword{datasets}

