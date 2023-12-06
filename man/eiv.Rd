\encoding{UTF-8}

\name{Ecological indicator values (EIV)}

\alias{Ellenberg ecological indicator values}
\alias{eiv}
\alias{eivClean}
\alias{eivMedian}
\alias{eivThreshold}

\title{
Work with Ecological Indicator Values.
}

\description{
Summarize species cover within defined threshold in order to identify nutrient poor or wet habitats.
}

\usage{
eivClean(x)
eivMedian(obj, plot)
eivThreshold(obj, plot, threshold = 4, parameter = "eiv.n", greater = FALSE, summary = TRUE)
}

\arguments{
  \item{ x }{
     \code{character}. A vector of strings to be cleaned.
  }
  \item{ obj }{
    A \code{Vegsoup*} object.
  }
  \item{ plot }{
    \code{character}. A rowname identifying a particular sample in a \code{Vegsoup*} object.
  }
  \item{ threshold }{
    \code{integer}. Threshold value of ecological indicator value.
  }
  \item{ parameter }{
    \code{character}. A column name identifying ecological indicator value (e.g. eiv.n or eiv.f).
  }
  \item{ greater }{
    \code{logical}. TRUE equals \verb{>=} and FALSE \verb{<=}.
  }
  \item{ summary }{
    \code{logical}. Report summary.
  }
}

\details{
Nothing here yet.
}

\value{
\code{numeric} vector.
}

\references{
Möller, H (2013). Zur Verwendung des Medians bei Zeigerwertberechnungen nach ELLENBERG. Tuexenia, 12:25-28
}

\author{
Roland Kaiser
}

\seealso{
\code{\link{help}}
}

\keyword{ database }

