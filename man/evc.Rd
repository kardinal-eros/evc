\encoding{UTF-8}

\name{EuroVeg Checklist (EVC)}

\alias{EuroVeg Checklist}
\alias{evc}

\alias{evcLookup}
\alias{evcCode}
\alias{evcClass}
\alias{evcMatrix}
\alias{evc2vegsoup}

\title{
EuroVegChecklist (EVC) and Diagnostic Species for Classes.
}

\description{
  This set of functions provides means to analyze diagnostic species of classes
  of the plant communities dominated by vascular plants (EVC1) as provided with
  Mucina et al. (2016) in Electronic Appendix S1 (EuroSpeciesList 1).
}

\usage{
evcLookup()
evcCode(obj)
evcClass(obj, plot)
evcMatrix(obj, select, restrict, weighted = TRUE)
evc2vegsoup(obj, restrict, select)
}

\arguments{
  \item{ obj }{
    A \code{Vegsoup*} object.
  }
  \item{ plot }{
    \code{character}. A row name identifying a particular sample in a
    \code{Vegsoup*} object.
  }
  \item{ select }{
    \code{character}. A vector of EVC class codes to be selected for output.
  }
  \item{ restrict }{
    \code{character}. A vector of EVC class codes to restrict the analysis.
  }
  \item{ weighted }{
    A \code{logical} object.
  }
}

\details{
  Given an object where species (taxa) are mapped to vegetation classes in the
  \code{taxonomy} in a single column named \code{"evc.class.code"} ...

  \code{evcCode} returns a list with class assignments.
  \code{evcClass} expands the \code{"evc.class.code"} string and casts a
  species by class matrix for a given plot.
  \code{evcMatrix} returns a matrix as in \code{evcClass} for the whole object.
  \code{evcLookup()} returns a look table of class codes as used in Mucina et
  al. (2016).
  \code{evc2vegsoup} aggregates species according to their class membership and
  uses these as instaed of the original species.
}

\value{
Nothing here yet.
}

\references{
Mucina, L. et. al (2016). Vegetation of Europe: hierarchical floristic classification system of vascular plant, bryophyte, lichen, and algal communities. Applied Vegetation Science, 19:3-264
}

\author{
Roland Kaiser
}

\seealso{
\code{\link{Ellenberg ecological indicator values}}
}

\keyword{ manip }

