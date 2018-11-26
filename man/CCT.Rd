% File computeCCT.Rd
\encoding{UTF-8}

\name{Correlated Color Temperature}
\title{Compute Correlated Color Temperature (CCT) and Points on the Planckian Locus}
\alias{CCTfromXYZ}
\alias{CCTfromxy}
\alias{CCTfromuv}
\alias{planckLocus}



\description{
Compute the CCT, in K, of XYZ, xy, and uv, by multiple methods.
And compute points on the Planckian locus.

At this time only one version of the locus is available - a spline 
(using \code{\link[stats:splinefun]{stats::splinefun}()} with \code{method="fmm"})
through the 31 uv (1960) locus points in \cite{Robertson} and \cite{Wyszecki & Stiles}.
This spline does not appear in \cite{Robertson}, but I am sure he would approve of it.
It has \eqn{C^1} continuity and good agreement with the true locus computed
from Planck spectra and the CIE 1931 standard observer color matching functions.
The maximum RMS error in the uv-plane is about \eqn{7.3 \times 10^{-6}}
over the valid temperature range from 1667 to \eqn{\infty} K.
A similar piecewise-linear interpolating path has a maximum RMS error about 10 times larger.
The 31 uv values in the table are accurate to the given 5 decimal places
(but see \bold{Note}).
The locus is parameterized directly by reciprocal color temperature (\eqn{10^6/T}),
and therefore indirectly by \eqn{T}.
We call either of these the \emph{native pameterization} of the locus.
The lines that are perpendicular to the locus are the \emph{native isotherms}.

Two more families of isotherms are available.
The \cite{Robertson} isotherms are tabulated just like the points on the locus,
and linear interpolation is used for intermediate temperatures.
The \cite{McCamy} isotherms are defined by a single cubic rational function in xy,
and no interpolation is necessary. 
Each isotherm family induce a slightly different parameterization of the locus -
the temperature at a locus point is the temperature of the isotherm
passing through that point.
\cr
The Robertson parameterization is \emph{not} \eqn{C^1},
but the \emph{geometric continuity} class is \eqn{G^1}.
The McCamy parameterization is \eqn{C^1}.
\cr
For the Robertson and native parameterization the valid temperature range is 1667 to \eqn{\infty} K.
For the McCamy parameterization the valid temperature range is 1702 to 34539 K.
}

\usage{
CCTfromXYZ( XYZ, isotherms='robertson',  locus='robertson', strict=FALSE )
CCTfromxy( xy, isotherms='robertson',  locus='robertson', strict=FALSE )
CCTfromuv( uv, isotherms='robertson', locus='robertson', strict=FALSE )

planckLocus( temperature, locus='robertson', param='robertson', delta=0, space=1960 )
}

\arguments{

\item{XYZ}{a numeric Mx3 matrix with XYZ tristimulus values (CIE 1931) in the rows,
or a numeric vector that can be converted to such a matrix, by row.}

\item{xy}{a numeric Mx2 matrix with xy chromaticity values (CIE 1931) in the rows,
or a numeric vector that can be converted to such a matrix, by row.}

\item{uv}{a numeric Mx2 matrix with uv chromaticity values (CIE UCS 1960) in the rows,
or a numeric vector that can be converted to such a matrix, by row.}

\item{isotherms}{A character vector whose elements match one
of the available isotherm families:
\code{'robertson'}, \code{'mccamy'}, and \code{'native'}.
Matching is partial and case-insensitive.
When more than one family is given, a matrix is returned, see \bold{Value}.
When \code{isotherms='native'} the isotherms are defined implicitly
as lines perpendicular to the locus, see \bold{Details}.
}

\item{locus}{
The only locus available value at this time is the \code{'robertson'} locus
defined above.
It is present here now for possible future use.
}

\item{strict}{The CIE considers the CCT of a chromaticity \code{uv} to be meaningful only
if the distance from \code{uv} to the Planckian locus is less than or equal to 0.05 
[in CIE UCS 1960].
If \code{strict=FALSE}, then this condition is ignored.
Otherwise, the distance is computed along the corresponding isotherm,
and if it exceeds 0.05 the returned CCT is set to \code{NA}.}

\item{temperature}{a M-vector of temperatures (in K) at which to compute points on the 
Planckian locus, either for \code{uv}, \code{u'v'}, or \code{xy}; see \code{space}.
}

\item{param}{the desired parameterization of the locus.
It can be either \code{'native'}, or a parameterization induced by the
\code{'robertson'} or \code{'mccamy'} isotherms.
}

\item{delta}{a vector of offset distances in \code{uv} (CIE UCS 1960), 
along the corresponding isotherms, from the locus.
Positive offsets are above the locus, and negative are below.
\code{delta} can have length M or 1, where M is the length of \code{temperature}.
If \code{delta} has length 1, that value is replicated to length M.
}

\item{space}{the year of the chromaticity space to return.
Valid values are 1960 (the default \code{uv}),
1976 (\code{u'v'}), and 1931 (\code{xy}).
}
}



\details{
Each of the isotherm families correspond to a parameterization of the locus.
All this is designed so a round trip:
temperature \rarrow  uv   \rarrow  CCT
(with the same choice of isotherm/parameterization)
has neglible error.

When \code{isotherms='native'}
the tangent line at a point on the locus is computed using the \code{deriv=1} argument to
\code{\link[stats:splinefun]{stats::splinefun}()}
and the normal line - the isotherm at the point - is then easily computed from the tangent line.

When \code{isotherms='robertson'} or \code{isotherms='mccamy'} the locus curve
has no effect on the computed CCT.
The locus is only used when computing the distance from the given uv point to the locus
(along the isotherm),
and therefore only affects the decision whether the CCT is meaningful
when \code{strict=TRUE}.

For \code{planckLocus()} with \code{param='mccamy'} the valid temperature range
is [1702,34539] K.
The \code{CCTfrom*()} functions with \code{isotherms='mccamy'} may return CCTs
outside this range if \code{strict=FALSE}, 
so it may be wise to do a sanity check on the value.
But if \code{strict=TRUE} they will not because the corresponding isotherm
does not intersect the locus, and so the distance to the locus cannot be computed.
}

\value{
\code{CCTfromXYZ()}, \code{CCTfromxy()}, and \code{CCTfromuv()} 
return a numeric vector of length M, or an MxN matrix.
It returns a matrix iff \code{length(isotherms) = N} \eqn{\ge} 2,
and the column names are set to the isotherm family names.
The names or rownames are set to the rownames of the input.
In case of error, the element of the vector or matrix is set to \code{NA_real_}.
In case there is an error in the arguments, the functions return \code{NULL}.
In these functions, the locus is not used unless
\code{isotherms='native'} or \code{strict=TRUE}.

\code{planckLocus()} returns an Mx2 matrix with chromaticies in the rows.
The column names are set appropriately for the value of \code{space}.
The row names are set from \code{temperature}.
In case of a single error, both entries in the row are set to \code{NA_real_}.
In case there is an error in the arguments, the functions return \code{NULL}.
}

\note{
The lookup table on page 228 in \cite{Wyszecki & Stiles}
contains an error at 325 mired,
which was corrected by Bruce Lindbloom (see \bold{Source}).

}

\source{
\url{http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_T.html}

}

\references{
McCamy, C. S.
\emph{Correlated color temperature as an explicit function of chromaticity coordinates}.
Color Research & Application.
Volume 17.
Issue 2.
pages 142-144.
April 1992.

Robertson, A. R.
Computation of correlated color temperature and distribution temperature.
Journal of the Optical Society of America.
58. pp. 1528-1535 (1968).

Wyszecki, Günther and W. S. Stiles.
\bold{Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.}
John Wiley & Sons, 1982.
Table 1(3.11). pp. 227-228.
}

\seealso{
\code{\link[stats:splinefun]{stats::splinefun}()},
\code{\link[colorSpec:computeCCT]{colorSpec::computeCCT}()},
the vignette \bold{Correlated Color Temperature Isotherms}
}

\examples{
# do a round trip and then compare
temperature = c(5003,6504)
uv  = planckLocus( temperature, delta=0.05 )
CCTfromuv( uv ) - temperature  
##  2.772227e-05 5.094369e-05 

# find some points on the daylight locus, and then their CCT
temperature = seq( 2000, 10000, by=1000 )
xy = daylightLocus( temperature )
cbind( xy, CCT=CCTfromxy(xy,iso='mccamy') )
##                x         y      CCT
## 2000K         NA        NA       NA
## 3000K         NA        NA       NA
## 4000K  0.3823436 0.3837663 4005.717
## 5000K  0.3457410 0.3586662 4999.998
## 6000K  0.3216915 0.3377984 5999.437
## 7000K  0.3053570 0.3216459 6997.542
## 8000K  0.2937719 0.3092195 7985.318
## 9000K  0.2852645 0.2995816 8948.809
## 10000K 0.2787996 0.2919672 9881.115

# compare all 3 different isotherms
CCTfromxy( xy, isotherms=c('robertson','native','mccamy') )
##        Robertson native.Rob   McCamy
## 2000K         NA         NA       NA
## 3000K         NA         NA       NA
## 4000K   4000.096   4000.062 4005.717
## 5000K   4999.749   4999.608 4999.998
## 6000K   5998.015   5999.242 5999.437
## 7000K   6997.858   6998.258 6997.542
## 8000K   7997.599   7996.985 7985.318
## 9000K   8999.301   8993.811 8948.809
## 10000K  9991.920   9992.672 9881.115
}



