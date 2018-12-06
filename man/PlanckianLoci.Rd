\encoding{UTF-8}
\name{Planckian Loci}
\alias{RobertsonLocus}
\alias{PrecisionLocus}
\docType{data}

\title{Planckian Loci - stored as Lookup Tables}

\description{
\tabular{ll}{
\code{RobertsonLocus} \tab the table from Robertson, with 31 points from 0 to 600 mired \cr
\code{PrecisionLocus} \tab a precomputed table, with 65 points from 0 to 1000 mired \cr
}
}

\format{
Both objects are \code{data.frame}s with these columns
\tabular{ll}{
{\code{mired}}  \tab {the reciprocal temperature \eqn{10^6/T}} \cr
{\code{u}}      \tab {the u chromaticity, in 1960 CIE} \cr
{\code{v}}      \tab {the v chromaticity, in 1960 CIE} \cr
}
The \code{PrecisionLocus}  \code{data.frame}  has these additional columns:
\tabular{ll}{
{\code{up}}  \tab {the 1st derivative of \code{u} with respect to \code{mired}} \cr
{\code{vp}}  \tab {the 1st derivative of \code{v} with respect to \code{mired}} \cr
{\code{upp}} \tab {the 2nd derivative of \code{u} with respect to \code{mired}} \cr
{\code{vpp}} \tab {the 2nd derivative of \code{v} with respect to \code{mired}} \cr
}
}

\details{
For \code{RobertsonLocus}, the values are taken from
\cite{Wyszecki & Stiles}.
The lookup table on page 228 
contains an error at 325 mired,
which was corrected by Bruce Lindbloom (see \bold{Source}).

For \code{PrecisionLocus}, the chromaticity values \code{u} and \code{v}
are computed from first principles, 
from the famous equation for the Planckian radiator
(with \eqn{c_2 = 1.4388 \times 10^{-2}})
and from the tabulated CIE 1931 standard observer color matching functions,
by summing from 360 to 830nm.
Let \eqn{\beta} denote the reciprocal temperature \eqn{10^6/T}.
We think of \code{u} as a function \eqn{u(\beta)}.
The column \code{up} is \eqn{u'(\beta)}, and \code{upp} is \eqn{u''(\beta)}.
And similarly for \code{v}.
The derivatives are computed from first principles, by summing the derivatives
of the Planckian formula from 360 to 830nm.
This includes the limiting case \eqn{\beta=0}.

When this package is loaded (during \code{.onLoad()}),
cubic splines are computed from \code{RobertsonLocus},
using \code{\link[stats:splinefun]{stats::splinefun}()} with \code{method="fmm"}).
And quintic splines are computed from \code{PrecisionLocus}.
Both splines are \eqn{C^2} continuous.
}

\source{
\url{http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_T.html}

}

\references{
Robertson, A. R.
Computation of correlated color temperature and distribution temperature.
Journal of the Optical Society of America.
58. pp. 1528-1535. 1968.

Wyszecki, Günther and W. S. Stiles.
\bold{Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.}
John Wiley & Sons, 1982.
Table 1(3.11). pp. 227-228.

}

\seealso{
\code{\link{CCTfromuv}()},
\code{\link{planckLocus}()}
}

\examples{
RobertsonLocus[ 1:10, ]
##    mired       u       v
## 1      0 0.18006 0.26352
## 2     10 0.18066 0.26589
## 3     20 0.18133 0.26846
## 4     30 0.18208 0.27119
## 5     40 0.18293 0.27407
## 6     50 0.18388 0.27709
## 7     60 0.18494 0.28021
## 8     70 0.18611 0.28342
## 9     80 0.18740 0.28668
## 10    90 0.18880 0.28997

PrecisionLocus[ 1:10, ]
##    mired         u         v           up           vp          upp          vpp
## 1      0 0.1800644 0.2635212 5.540710e-05 0.0002276279 7.115677e-07 1.977793e-06
## 2     10 0.1806553 0.2658948 6.291429e-05 0.0002469232 7.900243e-07 1.873208e-06
## 3     20 0.1813253 0.2684554 7.120586e-05 0.0002649377 8.679532e-07 1.722425e-06
## 4     30 0.1820820 0.2711879 8.026143e-05 0.0002812384 9.423039e-07 1.531723e-06
## 5     40 0.1829329 0.2740733 9.002982e-05 0.0002954676 1.010028e-06 1.309700e-06
## 6     50 0.1838847 0.2770894 1.004307e-04 0.0003073613 1.068393e-06 1.066350e-06
## 7     60 0.1849432 0.2802122 1.113592e-04 0.0003167582 1.115240e-06 8.120582e-07
## 8     70 0.1861132 0.2834161 1.226923e-04 0.0003235990 1.149155e-06 5.566812e-07
## 9     80 0.1873980 0.2866757 1.342971e-04 0.0003279171 1.169532e-06 3.088345e-07
## 10    90 0.1887996 0.2899664 1.460383e-04 0.0003298241 1.176525e-06 7.543963e-08
}



\keyword{datasets}
