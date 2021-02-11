\encoding{UTF-8}

\name{EuroVeg Checklist (EVC)}

\alias{EuroVeg Checklist}
\alias{evc}
\alias{evcCode}
\alias{evcClass}
\alias{evcMatrix}
\alias{evc2vegsoup}

\title{
EuroVeg Checklist (EVC) and Diagnostic Species for Classes.
}

\description{
This set of functions provides means to analyze diagnostic species of classes of the plant communities dominated by vascular plants (EVC1) as provided with Mucina et al. (2016) in Electronic Appendix S1 (EuroSpeciesList 1).
}

\usage{
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
    \code{character}. A rowname identifying a particular sample in a \code{Vegsoup*} object.
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
Nothing here yet.
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

