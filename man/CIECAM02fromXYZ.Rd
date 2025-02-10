\name{CIECAM02fromXYZ}
\alias{CIECAM02fromXYZ}

\title{Compute CIECAM02 Quantities from XYZ}

\description{CIECAM02 is a widely used Color Appearance Model.
This function computes many of the quantities in this model from XYZ.
}

\usage{CIECAM02fromXYZ( XYZ, XYZ_w, L_A=100, Y_b=20, surround='Average', discount=TRUE ) }

\arguments{
\item{XYZ}{a numeric Nx3 matrix with XYZ values in the rows.
It can also be a vector that can be converted to such a matrix, by row.}

\item{XYZ_w}{XYZ of the adapting white point. It should be scaled so that \eqn{Y_w}=100.
All values must be positive.
\code{XYZ_w} can also be a character string with the name of a standard illuminant,
which is passed to \code{\link{standardXYZ}()} to get the XYZ.}

\item{L_A}{absolute luminance of the adapting field.}

\item{Y_b}{luminous factor of background.}

\item{surround}{surround conditions, valid values are \code{'Average'}, \code{'Dim'}, \code{'Dark'}.
Matching is partial and caseless.}

\item{discount}{whether or not observers are discounting the illuminant (color constancy is in effect).}
}


\value{
a \code{data.frame} with N rows and these columns:

\item{XYZ}{the input Nx3 matrix}

\item{h}{hue angle, in degrees}

\item{J}{correlate of Lightness}

\item{Q}{correlate of Brightness}

\item{C}{correlate of Chroma}

\item{M}{correlate of Colorfulness}

\item{Jp}{Lightness \eqn{J'} in CAM02-UCS space}

\item{abp}{rectangular coordinates \eqn{(a',b')}  in CAM02-UCS space}

The \code{rownames} are set to those of \code{XYZ}.
The other function arguments (\code{XYZ_w}, ...) are attached to the \code{data.frame} as attributes.
}

\references{
Fairchild, Mark D.
CIECAM02. in
\bold{Color Appearance Models} (3rd ed.).
pp. 287-302.
John Wiley & Sons.
2013.
}

\examples{
CIECAM02fromXYZ( c(19.01,20,21.78), c(95.05,100,108.88), L_A=318.31 )
 
# XYZ.1 XYZ.2 XYZ.3        h        J        Q         C         M       Jp       abp.1       abp.2
# 19.01 20.00 21.78 219.0484 41.73109 195.3713 0.1047078 0.1088422 54.90433 -0.03666471 -0.02974185
}


\seealso{
\code{\link{DeltaE_CAM02}()}
}

\author{Glenn Davis}
\keyword{ISCC-NBS}

