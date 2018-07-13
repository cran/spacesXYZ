\name{toPolar}
\alias{LCHabfromLab}
\alias{LCHuvfromLuv}

\title{Convert CIE Lab and Luv to Polar Form}

\description{Convert the Rectangular Form of CIE Lab and Luv to Polar Form}

\usage{
LCHabfromLab( Lab ) 
LCHuvfromLuv( Luv ) 
}

\arguments{
\item{Lab}{ a numeric Nx3 matrix with CIE Lab coordinates in the rows, or a vector that can be converted to such a matrix, by row.}

\item{Luv}{ a numeric Nx3 matrix with CIE Luv coordinates in the rows, or a vector that can be converted to such a matrix, by row.}

}

\value{
a numeric Nx3 matrix with LCH coordinates in the rows.
The lightness L is simply copied from input to output.
The chroma C corresponds to radius in polar coordinates,
and the hue H corresponds to theta.
H is in degrees, not radians.
The rownames are copied from input to output.
}



\references{
Wikipedia.
CIE 1931 color space.
\url{https://en.wikipedia.org/wiki/CIE_1931_color_space}

}


\seealso{
\code{\link{LabfromLCHab}()},
\code{\link{LuvfromLCHuv}()}
}

\examples{
LCHabfromLab( c(50,0,0) )   #  a neutral gray
##        L Cab Hab
##  [1,] 50   0   0         # Hue is undefined, but set to 0

LCHabfromLab( c(50,0,20) ) 
##        L Cab Hab
##  [1,] 50  20  90         # 90 degrees, on yellow axis


LCHabfromLab( c(50,0,-20) ) 
##        L Cab Hab
##  [1,] 50  20 270         # 270 degrees, on blue axis

}

\keyword{XYZ}

