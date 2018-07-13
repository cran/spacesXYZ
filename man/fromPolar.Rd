\name{fromPolar}
\alias{LabfromLCHab}
\alias{LuvfromLCHuv}


\title{Convert CIE Lab and Luv from Polar Form to Rectangular Form}

\description{Convert the Polar Form of CIE Lab and Luv to Rectangular Form}

\usage{
LabfromLCHab( LCHab ) 
LuvfromLCHuv( LCHuv ) 
}

\arguments{
\item{LCHab}{ a numeric Nx3 matrix with CIE LCHab coordinates in the rows, or a vector that can be converted to such a matrix, by row.  The hue angle H must be in degrees.}

\item{LCHuv}{ a numeric Nx3 matrix with CIE LCHuv coordinates in the rows, or a vector that can be converted to such a matrix, by row.  The hue angle H must be in degrees.}
}

\value{
\item{LabfromLCHab()}{ returns a numeric Nx3 matrix with CIE Lab coordinates in the rows.
}

\item{LuvfromLCHuv()}{ returns a numeric Nx3 matrix with CIE Luv coordinates in the rows.
}

In both cases, the rownames are copied from input to output.
If the input is invalid, the functions return \code{NULL}.
}



\references{
Wikipedia.
CIE 1931 color space.
\url{https://en.wikipedia.org/wiki/CIE_1931_color_space}

}


\seealso{
\code{\link{LCHabfromLab}()},
\code{\link{LCHuvfromLuv}()}
}

\examples{
LabfromLCHab( c(50,10,45) ) 
##        L        a        b
##  [1,] 50 7.071068 7.071068     #  on line with slope 1
}

\keyword{XYZ}

