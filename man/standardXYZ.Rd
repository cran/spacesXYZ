% File standardXYZ.Rd
\encoding{UTF-8}

\name{standardXYZ}
\title{Query the Standardized XYZ values for Standard Illuminants}
\alias{standardXYZ}


\description{
In careful calcuations with standard illuminants,
it is often helpful to have the 'official' values of XYZ,
i.e. with the right number of decimal places.

}

\usage{
standardXYZ( name )
}

\arguments{
\item{name}{a subvector of
\code{c('A','B','C','D50','D50.ICC','D55','D65','D75',}
\code{'E','F2','F7','F11')}, which are the names of some standard illuminants.
Matching is partial and case-insensitive.
}
}

\details{
All XYZ values are taken from the ASTM publication in \bold{References},
except \code{B} which is taken from \cite{Wyszecki & Stiles}
and \code{D50.ICC} which is taken from ICC publications.
}

\section{Warning}{
The returned XYZs are normalized so that Y=1.
In other color domains, it is common to normalize so that Y=100;
in these cases be sure to multiply by 100.
}

\value{
An Mx3 matrix where M is the length of \code{name}.  
Each row filled with the official XYZ, but if the illuminant name is not recognized the row is all \code{NA}s.
The output XYZ is normalized so that \code{Y=1}.
The matrix \code{rownames} are set to the full illuminant names, and \code{colnames} to \code{c('X','Y','Z')}.
}



\references{
ASTM E 308 - 01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
(2001).

Günther Wyszecki and W. S. Stiles.
Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.
John Wiley & Sons, 1982.
Table I(3.3.8) p. 769.
}



\examples{
standardXYZ( c('a','d50','D50.ICC','D65') ) 
#                 X Y         Z
# A       1.0985000 1 0.3558500
# D50     0.9642200 1 0.8252100
# D50.ICC 0.9642029 1 0.8249054
# D65     0.9504700 1 1.0888300
}


