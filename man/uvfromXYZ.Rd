\name{uvfromXYZ}
\alias{uvfromXYZ}


\title{Convert from XYZ to Uniform Chromaticity Spaces}

\description{Convert from XYZ to Uniform Chromaticity Spaces}

\usage{
uvfromXYZ( XYZ, version=1976 ) 
}

\arguments{
\item{XYZ}{ a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.}

\item{version}{the year of the version - valid versions are 1976 and 1960.
The default \code{version=1976} has largely superseded \code{version=1960},
but the latter is still used for computing Correlated Color Temperature; 
see \code{\link[colorSpec]{computeCCT}()}.
The 1976 version is used in the function \code{\link{LuvfromXYZ}()}.
The 1960 coordinates are usually denoted by uv,
and the 1976 coordinates by u'v'.
The only difference is that \eqn{v' = (3/2) v}.
}
}

\value{
a numeric Nx2 matrix with u'v' (or uv) coordinates in the rows.
If \eqn{X + 15Y + 3Z \le 0}, uv are set to \code{NA}.

The rownames are copied from input to output.
If the input is invalid, the function returns \code{NULL}.
}



\references{
Wikipedia.
CIE 1931 color space.
\url{https://en.wikipedia.org/wiki/CIE_1931_color_space}

Wikipedia.
CIE 1960 color space.
\url{https://en.wikipedia.org/wiki/CIE_1960_color_space}
}


\seealso{
\code{\link{LuvfromXYZ}()},
\code{\link{standardXYZ}()}
}

\examples{
# locate some standard illuminants on the 1976 UCS diagram
uvfromXYZ( standardXYZ( c('C','D50','D65','E')  ) )

##             u'        v'
##  C   0.2008921 0.4608838
##  D50 0.2091601 0.4880734
##  D65 0.1978398 0.4683363
##  E   0.2105263 0.4736842
}

\keyword{XYZ}

