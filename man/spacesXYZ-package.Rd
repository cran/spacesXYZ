% File spacesXYZ-package.Rd
\docType{package}
\name{spacesXYZ-package}
\alias{spacesXYZ-package}
\title{CIE XYZ and some of Its Derived Color Spaces}



\description{
This package covers the basic CIE 1931 space, and derived spaces CIE xyY, Lab, and Luv.
The equations are taken from Bruce Lindbloom's \cite{CIE Color Calculator}.
Color areas that are *not* covered are:
\itemize{
\item spectral color data
\item device color spaces, e.g. RGB and CMYK
\item color order systems, e.g. Munsell, DIN, NCS, Ostwald, ...
}

}


\section{The API}{
The API is small.  There are functions to
\itemize{
\item convert between CIE XYZ and other CIE spaces
\item create and perform some standard chromatic adaptation transforms (CATs)
\item compute the CIE 1976 color difference \eqn{\Delta}E
\item retrieve XYZ of some standard illuminants
}
}


\section{Other Packages}{
Package \pkg{colorscience} is a superset of this one.
\cr
Package \pkg{colorspace} has similar functionality, and is much faster because it in compiled C.
But the whitepoint (e.g. used when converting XYZ to Lab) is a global variable
and cannot be passed as a parameter.
\cr
Package \pkg{grDevices} also has similar functionality
(in the function \code{\link[grDevices]{convertColor}()}),
but is missing chromaticities xy, uv, and u'v'.
}


\author{
Glenn Davis <gdavis@gluonics.com>
}

\references{
Lindbloom, Bruce.
\bold{CIE Color Calculator}.
\url{http://brucelindbloom.com/index.html?ColorCalculator.html}

Lindbloom, Bruce.
\bold{Color Difference Calculator}.
\url{http://brucelindbloom.com/index.html?ColorDifferenceCalc.html}
}



\keyword{package}

