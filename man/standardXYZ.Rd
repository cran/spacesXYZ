% File standardXYZ.Rd
\encoding{UTF-8}

\name{standardXYZ}
\title{Query the Standardized XYZ and xy values of Standard Illuminants and Whitepoints}
\alias{standardXYZ}
\alias{standardxy}

\description{
In careful calcuations with standard illuminants and whitepoints,
it is often helpful to have the 'official' values of XYZ and xy,
i.e. with the right number of decimal places.

}

\usage{
standardXYZ( name )

standardxy( name )
}

\arguments{
\item{name}{a subvector of
\code{c('A','B','C','C.NBS','C.NTSC','C.JOSA','D50','D50.ICC',}\cr
\code{'D55','D60','D60.ACES','D65','D75','E','F2','F7','F11')}, 
which are the names of some standard illuminants and white points.
Matching is partial and case-insensitive.
}
}

\details{
All XYZ values are taken from the ASTM publication in \bold{References},
except \code{B} which is taken from \cite{Wyszecki & Stiles}
and \code{D50.ICC} which is taken from ICC publications.

xy values were taken from \cite{CIE}, \cite{BT.709}, and \cite{TB-2018-001}.
For D65 the values in \cite{CIE} and \cite{BT.709} disagree;
the former has 5 digits and the latter has 4.
We have selected the value in \cite{BT.709} (page 3) since is far more commonly used.
Three of the Illuminant C variants are rarely used and obsolete.
}

\section{Warning}{
The returned XYZs are normalized so that Y=1.
In other color domains, it is common to normalize so that Y=100;
in these cases be sure to multiply by 100.

Some illuminants have no XYZ available and some have no xy.
In these cases, the rows are filled with NAs.
}

\value{
\code{standardXYZ()} returns
an Mx3 matrix where M is the length of \code{name}.  
Each row filled with the official XYZ, but if the illuminant name is not recognized
or if there is no data, the row is all \code{NA}s.
The output XYZ is normalized so that \code{Y=1}.
The matrix \code{rownames} are set to the full illuminant names, and \code{colnames} to \code{c('X','Y','Z')}.

Similarly, \code{standardxy()} returns an Mx2 matrix 
with \code{colnames} set to \code{c('x','y')}.
}



\references{
ASTM E 308 - 01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
2001.

BT.709.
Parameter values for the HDTV standards for production and international programme exchange.
June 2015.

CIE 015:2004 - Colorimetry, 3rd edition.
International Commission on Illumination (CIE).
Vienna Austria. 
Technical Report.
2004.

Günther Wyszecki and W. S. Stiles.
Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.
John Wiley & Sons, 1982.
Table I(3.3.8) p. 769.

TB-2018-001.
\bold{Derivation of the ACES White Point CIE Chromaticity Coordinates}.
The Academy of Motion Picture Arts and Sciences.
Science and Technology Council.
Academy Color Encoding System (ACES) Project.
June 15, 2018.
}



\examples{
standardXYZ( c('a','d50','D50.ICC','D65') ) 
#                 X Y         Z
# A       1.0985000 1 0.3558500
# D50     0.9642200 1 0.8252100
# D50.ICC 0.9642029 1 0.8249054
# D65     0.9504700 1 1.0888300

standardxy( c('a','D65','D60','D60.ACES','E','F2') )
#                x       y
# A        0.44758 0.40745
# D65      0.31270 0.32900
# D60      0.32163 0.33774
# D60.ACES 0.32168 0.33767
# E        0.33333 0.33333
# F2            NA      NA
}


