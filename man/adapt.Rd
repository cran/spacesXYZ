\name{adapt}
\alias{adaptxyY.CAT}
\alias{adaptxyY}
\alias{adaptXYZ.CAT}
\alias{adaptXYZ}
\alias{adaptLab.CAT}
\alias{adaptLab}
\alias{adaptLuv.CAT}
\alias{adaptLuv}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chromatic Adaptation Functions
}
\description{
Adapt XYZ, xyY, Lab, or Luv from a source viewing enviroment with a given illuminant,
to a target viewing environment with a different illuminant.
}
\usage{
\S3method{adaptXYZ}{CAT}( x, XYZ.src )

\S3method{adaptxyY}{CAT}( x, xyY.src )
\S3method{adaptLab}{CAT}( x, Lab.src )
\S3method{adaptLuv}{CAT}( x, Luv.src )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
a \bold{CAT} object as returned from \code{\link{CAT}()}
}

\item{XYZ.src}{
an Nx3 matrix, or a vector that can be converted to such a matrix, by rows.
Each row has an XYZ in the source viewing environment.}

\item{xyY.src}{
an Nx3 matrix, or a vector that can be converted to such a matrix, by rows.
Each row has an xyY in the source viewing environment.}

\item{Lab.src}{
an Nx3 matrix, or a vector that can be converted to such a matrix, by rows.
Each row has an Lab in the source viewing environment.}

\item{Luv.src}{
an Nx3 matrix, or a vector that can be converted to such a matrix, by rows.
Each row has an Luv in the source viewing environment.}
}



\details{
\code{adaptXYZ()} is the most fundamental of the group;
it simply multiplies each of the input XYZs by \code{x$M}.

\code{adaptxyY()} converts xyY.src to XYZ,
calls \code{adaptXYZ()},
and then converts back to xyY.tgt.
And it does an addition check:
if the xy of xyY.src is equal to the xy of \code{x$source.xyY},
then the xy of the returned xyY.tgt is set to be the xy of \code{x$target.xyY}.

\code{adaptLab()} and \code{adaptLuv()} work in a similar way.
When \code{Lab.src} is transformed to XYZ, the whitepoint is set to \code{x$source.XYZ}.
And when the adapted XYZ is transformed to adapted Lab, the whitepoint is set to \code{target.XYZ}.
}

\value{
\code{adaptXYZ()} returns an Nx3 matrix with adapted XYZ.  
Each row has an XYZ in the target viewing environment.
\cr
\code{adaptxyY()} returns an Nx3 matrix with adapted xyY.
Each row has an xyY in the target viewing environment.

\code{adaptLab()} and \code{adaptLuv()} return adapted Lab and Luv respectively.
}




\references{
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

Wikipedia.
CIECAM02.
\url{https://en.wikipedia.org/wiki/CIECAM02}
}

\seealso{
\code{\link{CAT}()},
\code{\link{standardXYZ}()}
}

\examples{
# try the Bradford method
bCAT = CAT( 'D50', 'D65', method='bradford' )

adaptXYZ( bCAT, c(1,1,0.5) )
##               X         Y         Z
##  [1,] 0.9641191 0.9921559 0.6567701


adaptLab( bCAT, c(50,20,-10) )
##              L        a         b
##  [1,] 49.97396 20.84287 -10.19661      # as expected, there is a change

adaptLab( bCAT, c(40,0,0) )
##        L a b
##  [1,] 40 0 0   # but adaptLab() always preserves neutrals


adaptLuv( bCAT, c(40,0,0) )
##        L u v
##  [1,] 40 0 0   # and adaptLuv() also preserves neutrals



# try the scaling method - now XYZ are scaled independently
sCAT = CAT( 'D50', 'D65', method='scaling' )

adaptLab( sCAT, c(50,20,-10) )
##        L  a   b
##  [1,] 50 20 -10    with sCAT, adaptLab() is now the identity for *all* colors


adaptLuv( sCAT, c(50,-20,10) )
##        L         u        v
##  [1,] 50 -18.32244 11.29946    but adaptLuv() is NOT the identity for all colors
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

