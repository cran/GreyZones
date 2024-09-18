\name{detectGreyZones}
\alias{detectGreyZones}

\title{
Detection of Grey Zones in Two-way Inter-rater Agreement Tables
}
\description{
Detects the existence of grey zones in in two-way inter-rater agreement tables.
}
\usage{
detectGreyZones(table)
}

\arguments{
  \item{table}{A \code{matrix} object that has equal number of columns and rows showing the agreement table for two raters.}
}

\details{
See Demirhan and Yilmaz (2023) for the details of the approach implemented by this function. 

If a grey zone is detected in the agreement table, the function reports AC2 and BP with quadratic weights. 

If no grey zone is detected in the agreement table, the function reports kappa, AC2 and BP with linear and quadratic weights.
}

\value{

\item{delta}{A matrix that shows the relative magnitude of deviance from the perfect agreement to the level of agreement for each cell of the given agreement table. Corresponds to Eq. (6) of Demirhan and Yilmaz (2023).}

\item{Delta}{A value that shows the maximum of the values in \code{delta}. Corresponds to Eq. (7) of Demirhan and Yilmaz (2023).}

\item{tau_Delta}{The value of the threshold for the detection of grey zones. If \code{Delta} > \code{tau_Delta} then there is a grey zone exists in the agreement table. Corresponds to Eq. (10) of Demirhan and Yilmaz (2023).}

\item{result}{Verbal result of the comparison of \code{Delta} and \code{tau_Delta}.}

\item{agreement}{Weighted agreement coefficient with linear and quadratic weights. See details. }

}

\references{
Demirhan, H., Yilmaz, A.E., 2023, Detection of grey zones in inter-rater agreement studies, BMC Medical Research Methodology,  23, Article number: 3, DOI: 10.1186/s12874-022-01759-7.

Tran, D., Demirhan, H., Dolgun, A., 2021, Bayesian approaches to the weighted kappa-like inter-rater agreement measures, Statistical Methods in Medical Research,  30, 10, 2329-2351, DOI: 10.1177/09622802211037068.

Tran, D., Dolgun, A., Demirhan, H., 2021, The impact of grey zones on the accuracy of agreement measures for ordinal tables, BMC Medical Research Methodology, 21:70, DOI: 10.1186/s12874-021-01248-3.
}

\author{
Haydar Demirhan and Ayfer Ezi Yilmaz

Maintainer: Haydar Demirhan <haydar.demirhan@rmit.edu.au>
}

\examples{
table <- matrix(c(10, 2, 7, 1, 3, 3, 0, 8, 80), nrow = 3, ncol = 3, byrow = TRUE) 
detectGreyZones(table)
}
