\name{toXYZ}
\alias{XYZfromxyY}
\alias{XYZfromLab}
\alias{XYZfromLuv}

\title{Convert other Color Spaces to XYZ}

\description{Convert other Color Spaces to XYZ}

\usage{
XYZfromxyY( xyY ) 
XYZfromLab( Lab, white ) 
XYZfromLuv( Luv, white ) 
}

\arguments{
\item{xyY}{ a numeric Nx3 matrix with CIE xyY coordinates in the rows, or a vector that can be converted to such a matrix, by row.}

\item{Lab}{ a numeric Nx3 matrix with CIE Lab coordinates in the rows, or a vector that can be converted to such a matrix, by row.}

\item{Luv}{ a numeric Nx3 matrix with CIE Luv coordinates in the rows, or a vector that can be converted to such a matrix, by row.}

\item{white}{a numeric 3-vector giving the XYZ of reference white.
\code{white} can also be a character string with the name of a standard illuminant,
which is passed to \code{\link{standardXYZ}()} to get the XYZ.
}
}

\value{
a numeric Nx3 matrix with XYZ coordinates in the rows.
The rownames are copied from input to output.
\cr
For \code{XYZfromxyY()} an input xyY with Y=0 is treated as a special case (black).
In this case x and y are ignored, and all XYZ are set to 0.
}



\references{
Wikipedia.
CIE 1931 color space.
\url{https://en.wikipedia.org/wiki/CIE_1931_color_space}

}


\seealso{
\code{\link{standardXYZ}()}
}

\examples{
XYZfromxyY(c(0.310897, 0.306510, 74.613450))
##              X        Y        Z
##  [1,] 75.68137 74.61345 93.13427


XYZfromLab( c(50,2,-3), 'D50' )
##               X         Y         Z
##  [1,] 0.1813684 0.1841865 0.1643335
}

\keyword{XYZ}

