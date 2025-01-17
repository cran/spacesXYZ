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
\item compute 3 standard color difference \eqn{\Delta}E metrics
\item retrieve XYZ and xy of some standard illuminants
}
}


\section{Other Packages}{
Package \pkg{colorscience} is a superset of this one.
\cr
Packages \pkg{colorspace} and  \pkg{farver} have similar functionality,
and are faster because they are in compiled C/C++.
\cr
Package \pkg{grDevices} also has similar functionality
(in the function \code{\link[grDevices]{convertColor}()}),
but is missing chromaticities xy, uv, and u'v'.
}

\section{Logging}{

Logging is performed using the package  \pkg{logger}.
This is a powerful package that allows a separate configuration
for logging from within \pkg{spacesXYZ}, and that is what I have done.
When \pkg{spacesXYZ} is loaded, the logging threshold is changed from \code{INFO} to \code{WARN}.
To change it back again, one can execute:

\code{logger::log_threshold( logger::INFO, namespace="spacesXYZ" )     # preferred}  \cr
or  \cr
\code{library( logger )} \cr
\code{log_threshold( INFO, namespace="spacesXYZ" )    # not preferred}

The layout callback function is customized;
it adds the name of the calling function to the message.
To change it back again, one can execute:  \cr
\code{log_layout( layout_simple, namespace="spacesXYZ" )} \cr
or to install ones own layout function, one can execute:  \cr
\code{log_layout( <your function>, namespace="spacesXYZ" )}.

The appender callback functions is also customized;
it comes to an immediate stop if the log event level is \code{FATAL} or \code{ERROR}.
To continue execution after such a log event, one can execute:  \cr
\code{log_appender( appender_console, namespace="spacesXYZ" )}

The formatter callback function is initialized to be \code{formatter_sprintf()};
this should not be changed.
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

