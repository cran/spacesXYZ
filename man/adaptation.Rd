\name{adaptation}
\alias{CAT}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chromatic Adaptation Transforms (CATs)
}
\description{
Construct transforms from a source viewing enviroment with a given illuminant,
to a target viewing environment with a different illuminant.
Some standard linear von-Kries-based CAT methods are available.
}
\usage{
CAT( source.XYZ, target.XYZ, method="Bradford" )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{source.XYZ}{
the XYZ of the illuminant in the source viewing environment.
\code{source.XYZ} can also be a string with the name of a standard illuminant
as in the function \code{standardXYZ()}.
}

\item{target.XYZ}{
the XYZ of the illuminant in the target viewing environment.
\code{target.XYZ} can also be a string with the name of a standard illuminant
as in the function \code{standardXYZ()}.
}

\item{method}{
the method used for the chromatic adaptation. Available methods are:
\code{"Bradford"}, \code{"VonKries"}, \code{"MCAT02"}, \code{"Bianco+Schettini"},
and \code{"scaling"};
see \bold{References}.
Partial matching is enabled, and matching is case-insensitive.
\cr
\code{method} can also be a 3x3 matrix, which is the cone response matrix
used to construct a von-Kries-based CAT.
The matrix must be invertible and map both \code{source.XYZ} and \code{target.XYZ}
to the positive octant.
}

}


\value{
\code{CAT()} returns an object with S3 class \bold{CAT}, which can be passed to 
\code{\link{adaptXYZ}()}, 
\code{\link{adaptxyY}()}, 
\code{\link{adaptLab}()}, or 
\code{\link{adaptLuv}()}.

An object with S3 class \bold{CAT} is a list with the following items:
\describe{
\item{method}{full name of the adaptation method, as in \bold{Arguments}.
If argument \code{method} is a 3x3 matrix, then this \code{method} is \code{NA}.}
\item{Ma}{3x3 \emph{cone response matrix} \eqn{M_A} for the method, as defined in \cite{Lindbloom}}
\item{source.XYZ}{XYZ of the illuminant in the source viewing environment}
\item{source.xyY}{xyY of the illuminant in the source viewing environment}

\item{target.XYZ}{XYZ of the illuminant in the target viewing environment}
\item{target.xyY}{xyY of the illuminant in the target viewing environment}

\item{M}{3x3 matrix defining the CAT.  
The matrix is written on the left and the source XYZ is
written as a column vector on the right.
This matrix depends continuously on source.XYZ and target.XYZ,
and when these are equal, M is the identity.
Therefore, when source.XYZ and target.XYZ are close, M is close to the identity.
Compare with \cite{Lindbloom}.}
}
}


\note{
Chromatic adaptation can be viewed as an Aristotelian Analogy of Proportions.
For more about this, see the vignette
\href{../doc/adaptation.html}{\bold{Chromatic Adaptation}}.
}


\references{


Bianco, Simone and Raimondo Schettini.
\bold{Two new von Kries based chromatic adaptation transforms found by numerical optimization}.
Color Research & Application.
v. 35. i. 3.
Jan 2010.


Hunt, R. W. G.
\bold{The Reproduction of Colour}.  6th Edition.
John Wiley & Sons.
2004.

International Color Consortium.
ICC.1:2001-04.
File Format for Color Profiles.
2001.

Lindbloom, Bruce.
Chromatic Adaptation.
\url{http://brucelindbloom.com/Eqn_ChromAdapt.html}

Pascale, Danny.
A Review of RGB Color Spaces ...from xyY to R'G'B'.
\url{https://babelcolor.com/index_htm_files/A\%20review\%20of\%20RGB\%20color\%20spaces.pdf}
2003.

Wikipedia.
CIECAM02.
\url{https://en.wikipedia.org/wiki/CIECAM02}
}

\seealso{
\code{\link{standardXYZ}()},
\code{\link{adaptXYZ}()}, 
\code{\link{adaptxyY}()}, 
\code{\link{adaptLab}()}, 
\code{\link{adaptLuv}()}
}

\examples{
D65toC = CAT( 'D65', 'C' )
D65toC

##  $method
##  [1] "Bradford"
##  
##  $Ma
##          X       Y       Z
##  L  0.8951  0.2664 -0.1614
##  M -0.7502  1.7135  0.0367
##  S  0.0389 -0.0685  1.0296
##  
##  $source.XYZ
##            X Y       Z
##  D65 0.95047 1 1.08883
##  
##  $source.xyY
##              x         y Y
##  D65 0.3127266 0.3290231 1
##  
##  $target.XYZ
##          X Y       Z
##  C 0.98074 1 1.18232
##  
##  $target.xyY
##            x         y Y
##  C 0.3100605 0.3161496 1
##  
##  $M
##              X            Y           Z
##  X 1.009778519  0.007041913 0.012797129
##  Y 0.012311347  0.984709398 0.003296232
##  Z 0.003828375 -0.007233061 1.089163878
##  
##  attr(,"class")
##  [1] "CAT"  "list"


adaptXYZ( D65toC, c(1,1,0.5) )
##              X         Y         Z
##  [1,] 1.023219 0.9986689 0.5411773
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

