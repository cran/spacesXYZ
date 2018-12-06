\name{fromXYZ}
\alias{xyYfromXYZ}
\alias{LabfromXYZ}
\alias{LuvfromXYZ}

\title{Convert from XYZ to other Color Spaces}

\description{Convert from XYZ to other Color Spaces}

\usage{
xyYfromXYZ( XYZ ) 
LabfromXYZ( XYZ, white ) 
LuvfromXYZ( XYZ, white ) 
}

\arguments{
\item{XYZ}{ a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.}

\item{white}{a numeric 3-vector giving the XYZ of reference white; all 3 numbers must be positive.
\code{white} can also be a character string with the name of a standard illuminant,
which is passed to \code{\link{standardXYZ}()} to get the XYZ.
}
}

\value{
\item{xyYfromXYZ()}{ returns a numeric Nx3 matrix with CIE xyY coordinates in the rows.
If the sum \code{X+Y+Z==0}, \code{xy} are set to \code{NA}.
}

\item{LabfromXYZ()}{ returns a numeric Nx3 matrix with CIE Lab coordinates in the rows}

\item{LuvfromXYZ()}{ returns a numeric Nx3 matrix with CIE Luv coordinates in the rows}

In all cases, the rownames are copied from input to output.
If the input is invalid, the functions return \code{NULL}.
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
D65 = standardXYZ( 'D65' )

xyYfromXYZ( D65 )
##              x         y Y
##  D65 0.3127266 0.3290231 1     # probably not familiar

round( xyYfromXYZ(D65), 4 )
##           x     y Y
##  D65 0.3127 0.329 1        # probably more familiar


LabfromXYZ( 0.18*D65, D65 )   # 18\% gray card
##             L a b
##  D65 49.49611 0 0          # exactly neutral, and L is about 50


D50 = standardXYZ( 'D50' )

LabfromXYZ( D50, D65 )
##        L        a        b
##  D50 100 2.399554 17.65321   # D50 is far from neutral (yellowish) in D65 viewing environment
}

\keyword{XYZ}

