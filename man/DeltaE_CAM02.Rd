\name{DeltaE_CAM02}
\alias{DeltaE_CAM02}

\title{Calculate the Color Difference in CAM02-UCS between Two Colors}

\description{Calculate Standard CIE Color Differences between two colors
in uniform color space CAM02-UCS
}

\usage{DeltaE_CAM02( CAM_1, CAM_2 ) }

\arguments{
\item{CAM_1}{a \code{data.frame} returned from \code{\link{CIECAM02fromXYZ}()}, with 1 or N rows.}


\item{CAM_2}{a \code{data.frame} returned from \code{\link{CIECAM02fromXYZ}()}, with 1 or N rows.}

These 2 \code{data.frame}s must have an equal number of rows, *or*
one of them can have a single row which is then replicated to match the number
or rows in the other one.
The \eqn{\Delta}E is then computed between corresponding rows.
}

\value{
\code{DeltaE_CAM02()} returns a numeric vector of length N,
containing the computed \eqn{\Delta}E in uniform color space CAM02-UCS,
between corresponding rows of \code{CAM_1} and \code{CAM_2}.
}


\references{
Fairchild, Mark D.
CIECAM02. in
\bold{Color Appearance Models} (3rd ed.).
pp. 287-302.
John Wiley & Sons.
2013.
}

\seealso{
\code{\link{CIECAM02fromXYZ}()}
}


\examples{
CAM = CIECAM02fromXYZ( c(19.01,20,21.78, 25,20,15), c(95.05,100,108.88), L_A=100 )
CAM
## XYZ.1 XYZ.2 XYZ.3         h        J        Q        C          M       Jp      abp.1      abp.2
## 19.01 20.00 21.78 211.58993 41.57911 153.4480  0.80776  0.7624239 54.74946 -0.2796279 -0.1719605
## 25.00 20.00 15.00  22.98284 43.12872 156.2812 33.08790 31.2308204 56.31674  9.4290922  3.9990801

DeltaE_CAM02( CAM[1, ], CAM[2, ] )

## [1] 10.68238
}

