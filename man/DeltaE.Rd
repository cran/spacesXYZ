\name{DeltaE}
\alias{DeltaE}

\title{Calculate the Color Difference between Two Colors}

\description{Calculate Standard CIE Color Differences between two Colors 
}

\usage{DeltaE( Lab1, Lab2, metric=1976 ) }

\arguments{
\item{Lab1}{a numeric Nx3 matrix with Lab values in the rows, or a vector that can be converted to such a matrix, by row.
\code{Lab1} can also be a numeric 3-vector with a single Lab, and it is then replicated
to match the size of \code{Lab2}.
}

\item{Lab2}{a numeric Nx3 matrix with Lab values in the rows, or a vector that can be converted to such a matrix, by row.
\code{Lab2} can also be a numeric 3-vector with a single Lab, and it is then replicated
to match the size of \code{Lab1}.
}

\item{metric}{a vector of color metric specifiers.
Valid values are \code{'1976'}, \code{'1994'}, and \code{'2000'}, 
which refer to the year the metric was recommended by the CIE.
They can also be given as integers, as shown.
These metrics are often denoted 
\eqn{\Delta E_{1976}}, \eqn{\Delta E_{1994}}, and \eqn{\Delta E_{2000}}.
When more than one metric is given, a matrix is returned, see \bold{Value}.
}

}

\value{
\code{DeltaE()} returns a numeric vector of length N, or an NxM matrix.
It returns a matrix iff \code{length(metric)=} M\eqn{\ge} 2;
the column names are set to the metric names.
The elements of the output are the pairwise differences, 
i.e. between row i of \code{Lab1} and  row i of \code{Lab2}.
The names or rownames are set to the rownames of one of the input matrices.

For \code{metric=1976}
the distance is simply the Euclidean distance between the two points in Lab, see \cite{Hunt} p. 111.
\cr
For \code{metric=1994} the symmetric variant is used, see \cite{Hunt} p. 670.
There is an asymmetric variant which is not available in this package.
The weighting coefficients are for \bold{graphic arts} (not for \bold{textiles}).
\cr
For \code{metric=2000} the distance is insanely complicated, see \cite{Hunt} p. 671.
\cr
All these metrics are \emph{symmetric},
which means that swapping \code{Lab1} and \code{Lab2} does not change the result.
}


\references{
Hunt, R. W. G.
\bold{The Reproduction of Colour}.  6th Edition.
John Wiley & Sons.
2004.
}


\examples{
DeltaE( c(50,0,0),  c(51,2,2,  52,10,11,  46,-13,16)  )
## [1]  3 15 21

path = system.file( "extdata/ciede2000testdata.txt", package='spacesXYZ' )
df   = read.table( path, sep='\t', quote='', head=TRUE )
Lab1 = as.matrix( df[ , 1:3 ] )
Lab2 = as.matrix( df[ , 4:6 ] )
cbind( Lab1, Lab2, DeltaE( Lab1, Lab2, metric=c(1976,2000) ) )[ 1:10, ]

##       LAB_L_REF LAB_A_REF LAB_B_REF LAB_L_SAM LAB_A_SAM LAB_B_SAM DeltaE.1976 DeltaE.2000
##  [1,]   50.0000    2.6772  -79.7751   50.0000    0.0000  -82.7485   4.0010633   2.0424597
##  [2,]   50.0000    3.1571  -77.2803   50.0000    0.0000  -82.7485   6.3141501   2.8615102
##  [3,]   50.0000    2.8361  -74.0200   50.0000    0.0000  -82.7485   9.1776999   3.4411906
##  [4,]   50.0000   -1.3802  -84.2814   50.0000    0.0000  -82.7485   2.0627008   0.9999989
##  [5,]   50.0000   -1.1848  -84.8006   50.0000    0.0000  -82.7485   2.3695707   1.0000047
##  [6,]   50.0000   -0.9009  -85.5211   50.0000    0.0000  -82.7485   2.9152927   1.0000130
##  [7,]   50.0000    0.0000    0.0000   50.0000   -1.0000    2.0000   2.2360680   2.3668588
##  [8,]   50.0000   -1.0000    2.0000   50.0000    0.0000    0.0000   2.2360680   2.3668588
##  [9,]   50.0000    2.4900   -0.0010   50.0000   -2.4900    0.0009   4.9800004   7.1791720
## [10,]   50.0000    2.4900   -0.0010   50.0000   -2.4900    0.0010   4.9800004   7.1791626
}

