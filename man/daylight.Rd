\name{daylight}
\alias{daylightLocus}

\title{Compute Points on the Daylight Locus, in multiple Chromaticity Spaces}

\description{Compute points on the daylight locus, in multiple chromaticity spaces}

\usage{
daylightLocus( temperature, space=1931 )
}

\arguments{

\item{temperature}{an M-vector of temperatures (in K) at which to compute points on the 
daylight locus.
The valid temperatures range is 4000K to 25000K;
outside this range the output is set to \code{NA}.
}

\item{space}{the year of the output chromaticity space desired - valid values are
1931, 1976 and 1960.
The 1931 is the original and denoted by xy;  the others are derived from it.
The 1960 coordinates are usually denoted by uv,
and the 1976 coordinates by u'v'.
The only difference is that \eqn{v' = (3/2) v}.
}
}

\value{
a numeric Mx2 matrix with xy, uv, or u'v' coordinates in the rows.
The colnames of the output are set appropriately.

The names of the temperature are copied to the rownames of the output,
unless these names are \code{NULL} when the temperatures followed by 'K' are used.

If the input is invalid, the function returns \code{NULL}.
}



\references{
Wyszecki, Günther and W. S. Stiles.
\bold{Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.}
John Wiley & Sons, 1982.
pp. 145-146.
}


\seealso{
\code{\link{uvfromxy}()}
}

\examples{
# find some points on the daylight locus, and then their CCT
temp = seq( 2000, 10000, by=1000 )
xy = daylightLocus( temp )
cbind( xy, CCT=CCTfromxy(xy) )
##                x         y      CCT
## D2000         NA        NA       NA
## D3000         NA        NA       NA
## D4000  0.3823436 0.3837663 4000.096
## D5000  0.3457410 0.3586662 4999.749
## D6000  0.3216915 0.3377984 5998.015
## D7000  0.3053570 0.3216459 6997.858
## D8000  0.2937719 0.3092195 7997.599
## D9000  0.2852645 0.2995816 8999.301
## D10000 0.2787996 0.2919672 9991.920
}

\keyword{XYZ}

